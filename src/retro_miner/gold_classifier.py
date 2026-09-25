"""Training labels and a gradient-boosted score for gold MEI calls.

Positives are 1000 Genomes overlaps in the high-rate prefix of one genome-wide
gold list per sample. That list is re-sorted with the same evidence priority as
a single-chromosome gold table. Known overlaps below the prefix are left out.
Negatives are a family-matched random sample of other gold calls below it.

Chromosome, position, catalog-overlap fields, and the composite rank scores are
not features.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

GOLD_REVIEW_NAME = "candidate_loci.mei.gold_review.tsv"
MEI_TABLE_NAME = "candidate_loci.mei.tsv"

# Composite scores that define table order. Keeping them would reprint the rank.
_RANK_SCORE_COLUMNS = (
    "read_support_heuristic_score",
    "insertion_model_score",
)

_COORDINATE_COLUMNS = (
    "chrom",
    "window_start",
    "window_end",
    "discovery_window_start",
    "discovery_window_end",
    "consensus_insertion_breakpoint_pos",
    "consensus_breakpoint_interval_start",
    "consensus_breakpoint_interval_end",
    "consensus_insertion_mei_5p_coord",
    "consensus_insertion_mei_3p_coord",
    "consensus_insertion_mei_5p_coord_full",
    "consensus_insertion_mei_3p_coord_full",
    "asm_insertion_mei_start",
    "asm_insertion_mei_end",
    "asm_non_mei_partner_chrom",
    "asm_non_mei_partner_pos",
)

NUMERIC_FEATURES = (
    "consensus_poly_at_min_bp",
    "consensus_tsd_len_estimate",
    "consensus_breakpoint_interval_width_bp",
    "consensus_insertion_mei_span",
    "local_bam_peak_depth_z",
    "coherence_score",
    "split_cluster_window_reads",
    "split_cluster_reads",
    "split_cluster_read_fraction",
    "split_cluster_binomial_z",
    "mei_consensus_overlap_reads",
    "flank_mei_left",
    "flank_mei_right",
    "flank_mei_balance",
    "flank_polya_left",
    "flank_polya_right",
    "poly_at_reads_max",
    "poly_at_max_run_max",
    "junk_flag_count",
    "mate_junk_flag_count",
    "flag_segdup",
    "flag_low_mappability",
    "flag_gap_region",
    "flag_encode_blacklist",
    "flag_outside_giab_highconf",
    "breakpoint_simple_repeat",
)

CATEGORICAL_FEATURES = (
    "mei_family",
    "insertion_event_class",
    "complex_sv_signature_label",
    "nested_in_same_MEI",
    "consensus_insertion_orientation",
    "consensus_breakpoint_confidence_tier",
    "breakpoint_l1_en_motif_type",
    "discordant_mei_majority",
    "two_sided_support",
    "poly_at_supported",
    "tsd_or_polyA_supported",
    "family_agreement",
    "strand_agreement",
    "complex_mei_event",
    "complex_ins_with_del",
)

_JUNK_FLAG_COLUMNS = (
    "flag_segdup",
    "flag_low_mappability",
    "flag_gap_region",
    "flag_encode_blacklist",
    "flag_outside_giab_highconf",
    "junk_flag_count",
    "mate_junk_flag_count",
    "nested_repeat_class",
)

_MERGE_KEYS = ("chrom", "window_start", "window_end")


def mei_family(value: object) -> str:
    """Collapse consensus family labels to Alu, L1, SVA, or other."""
    text = str(value or "").upper()
    if "ALU" in text:
        return "Alu"
    if "SVA" in text:
        return "SVA"
    if "LINE" in text or "L1" in text:
        return "L1"
    return "other"


def _bernoulli_ll(successes: np.ndarray, trials: np.ndarray) -> np.ndarray:
    """Bernoulli log-likelihood at the maximum-likelihood rate. 0*log(0) is 0."""
    successes = np.asarray(successes, dtype=float)
    trials = np.asarray(trials, dtype=float)
    ll = np.zeros(np.broadcast(successes, trials).shape, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.divide(successes, trials, out=np.zeros_like(successes, dtype=float), where=trials > 0)
        hit = successes > 0
        miss = (trials - successes) > 0
        ll[hit] += successes[hit] * np.log(rate[hit])
        ll[miss] += (trials[miss] - successes[miss]) * np.log(1.0 - rate[miss])
    return ll


def bernoulli_changepoint(known: Sequence[bool]) -> int:
    """Exclusive end of the high-overlap prefix on a ranked gold list.

    One cut per list. The cut is the split whose overlap rate above and below
    differ the most, and it is kept only when that two-rate model beats a
    single rate by a BIC penalty. A fixed run of non-overlaps is not used, so
    the cut does not depend on how many gold calls a chromosome happens to have.
    """
    labels = np.asarray(list(known), dtype=bool)
    n = int(labels.size)
    if n < 2 or not bool(labels.any()):
        return 0
    counts = np.cumsum(labels.astype(float))
    total = float(counts[-1])
    prefix_n = np.arange(1, n, dtype=float)
    prefix_k = counts[:-1]
    suffix_n = n - prefix_n
    suffix_k = total - prefix_k
    prefix_rate = prefix_k / prefix_n
    suffix_rate = suffix_k / suffix_n
    two_rate = _bernoulli_ll(prefix_k, prefix_n) + _bernoulli_ll(suffix_k, suffix_n)
    one_rate = float(_bernoulli_ll(np.array([total]), np.array([float(n)]))[0])
    improved = (prefix_rate > suffix_rate) & ((two_rate - one_rate) > np.log(n))
    if not bool(np.any(improved)):
        return 0
    two_rate = np.where(improved, two_rate, -np.inf)
    return int(np.argmax(two_rate)) + 1


def rank_gold_like_review(gold: pd.DataFrame) -> pd.DataFrame:
    """Sort gold rows with the single-chromosome review priority.

    Catalog membership is blanked for the sort. The review order otherwise
    exempts known overlaps from the low-complexity penalty, which would pull
    the label into the ranking.
    """
    from retro_miner.mei_support import _prioritize_mei_candidates

    work = gold.copy()
    if "known_mei_polymorphism" in work.columns:
        work["_catalog_known"] = _as_bool(work["known_mei_polymorphism"]).to_numpy()
    else:
        work["_catalog_known"] = False
    work["known_mei_polymorphism"] = False
    if "tsd_seq" not in work.columns and "consensus_tsd_seq" in work.columns:
        work["tsd_seq"] = work["consensus_tsd_seq"]
    ranked = _prioritize_mei_candidates(work, stage_first=False)
    ranked["known_mei_polymorphism"] = ranked["_catalog_known"].astype(bool)
    return ranked.drop(columns=["_catalog_known"]).reset_index(drop=True)


def _as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    text = series.fillna("").astype(str).str.strip().str.lower()
    return text.isin({"1", "true", "t", "yes"})


def _read_columns(path: Path, columns: Iterable[str]) -> pd.DataFrame:
    """Read a TSV, keeping the first copy of any duplicated header name."""
    wanted = set(columns)
    frame = pd.read_csv(path, sep="\t", usecols=lambda name: name in wanted, low_memory=False)
    return frame.loc[:, ~frame.columns.duplicated()].copy()


def _numeric_col(frame: pd.DataFrame, name: str) -> pd.Series:
    if name not in frame.columns:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    series = frame[name]
    if series.dtype == bool:
        return series.astype(float)
    text = series.fillna("").astype(str).str.strip().str.lower()
    mapped = text.map({"1": 1.0, "0": 0.0, "true": 1.0, "false": 0.0, "t": 1.0, "f": 0.0})
    numeric = pd.to_numeric(series, errors="coerce")
    return numeric.fillna(mapped)


def _side_max(frame: pd.DataFrame, left_name: str, right_name: str, suffix: str) -> None:
    frame[f"flank_{suffix}_left"] = pd.concat(
        [
            _numeric_col(frame, f"disease_{left_name}"),
            _numeric_col(frame, f"control_{left_name}"),
        ],
        axis=1,
    ).max(axis=1)
    frame[f"flank_{suffix}_right"] = pd.concat(
        [
            _numeric_col(frame, f"disease_{right_name}"),
            _numeric_col(frame, f"control_{right_name}"),
        ],
        axis=1,
    ).max(axis=1)


def _attach_junk_flags(gold: pd.DataFrame, mei_path: Path) -> pd.DataFrame:
    if not mei_path.exists():
        return gold
    flags = _read_columns(mei_path, list(_MERGE_KEYS) + list(_JUNK_FLAG_COLUMNS))
    if flags.empty or not set(_MERGE_KEYS).issubset(flags.columns):
        return gold
    flags = flags.drop_duplicates(_MERGE_KEYS, keep="first")
    keep_flags = [c for c in _JUNK_FLAG_COLUMNS if c in flags.columns]
    return gold.merge(flags[list(_MERGE_KEYS) + keep_flags], on=list(_MERGE_KEYS), how="left")


def load_gold_chromosome(path: Path) -> pd.DataFrame:
    """Gold rows from one review table, in the table's priority order."""
    needed = set(_MERGE_KEYS) | {
        "analysis_stage_tier",
        "known_mei_polymorphism",
        "consensus_mei_family",
        "consensus_poly_at_min_bp",
        "consensus_tsd_len_estimate",
        "consensus_breakpoint_interval_width_bp",
        "consensus_insertion_mei_span",
        "local_bam_peak_depth_z",
        "coherence_score",
        "split_cluster_window_reads",
        "split_cluster_reads",
        "split_cluster_binomial_z",
        "mei_consensus_overlap_reads",
        "disease_left_flank_mei_reads",
        "disease_right_flank_mei_reads",
        "control_left_flank_mei_reads",
        "control_right_flank_mei_reads",
        "disease_left_flank_polya_reads",
        "disease_right_flank_polya_reads",
        "control_left_flank_polya_reads",
        "control_right_flank_polya_reads",
        "disease_poly_at_reads",
        "control_poly_at_reads",
        "disease_poly_at_max_run",
        "control_poly_at_max_run",
        "insertion_event_class",
        "complex_sv_signature_label",
        "nested_in_same_MEI",
        "consensus_insertion_orientation",
        "consensus_breakpoint_confidence_tier",
        "breakpoint_l1_en_motif_type",
        "discordant_mei_majority",
        "two_sided_support",
        "poly_at_supported",
        "tsd_or_polyA_supported",
        "disease_family_agreement",
        "control_family_agreement",
        "disease_strand_agreement",
        "control_strand_agreement",
        "complex_mei_event",
        "complex_ins_with_del",
        "disease_supporting_reads",
        "control_supporting_reads",
        "consensus_tsd_seq",
        "coherence_score",
        "insertion_model_score",
        "event_clip_overlap_consistency",
        "poly_at_reads",
        "poly_at_max_run",
    }
    frame = _read_columns(path, needed)
    if "analysis_stage_tier" not in frame.columns:
        return pd.DataFrame()
    gold = frame.loc[frame["analysis_stage_tier"].astype(str).str.lower().eq("gold")].copy()
    gold["known_mei_polymorphism"] = _as_bool(gold.get("known_mei_polymorphism", False))
    mei_path = path.with_name(MEI_TABLE_NAME)
    gold = _attach_junk_flags(gold, mei_path)
    return gold.reset_index(drop=True)


def label_gold_table(gold: pd.DataFrame) -> pd.DataFrame:
    """Mark positives, excluded catalog outliers, and the negative pool.

    The frame must already be gold-only and in genome-wide review rank order.
    """
    out = gold.copy()
    known = out["known_mei_polymorphism"].fillna(False).astype(bool)
    end = bernoulli_changepoint(known.tolist())
    out["gold_rank"] = np.arange(1, len(out) + 1)
    out["in_top_cluster"] = out["gold_rank"] <= end
    role = np.full(len(out), "unlabeled", dtype=object)
    role[known.to_numpy() & out["in_top_cluster"].to_numpy()] = "positive"
    role[known.to_numpy() & ~out["in_top_cluster"].to_numpy()] = "excluded_outlier_known"
    role[~known.to_numpy() & ~out["in_top_cluster"].to_numpy()] = "negative_pool"
    out["train_role"] = role
    if "consensus_mei_family" in out.columns:
        out["mei_family"] = out["consensus_mei_family"].map(mei_family)
    elif "mei_family" not in out.columns:
        out["mei_family"] = "other"
    return out


def sample_matched_negatives(
    labeled: pd.DataFrame,
    *,
    seed: int = 13,
) -> pd.DataFrame:
    """Keep every positive and an equal number of same-family negative-pool rows.

    Matching is within sample and family. Catalog outliers stay out of both classes.
    """
    if "sample" not in labeled.columns:
        raise ValueError("labeled rows need a sample column")
    rng = np.random.default_rng(seed)
    parts: list[pd.DataFrame] = []
    for (_, _), sub in labeled.groupby(["sample", "mei_family"], sort=False):
        positives = sub.loc[sub["train_role"].eq("positive")]
        pool = sub.loc[sub["train_role"].eq("negative_pool")]
        n_pos = len(positives)
        if n_pos == 0:
            continue
        take = min(n_pos, len(pool))
        if take == 0:
            parts.append(positives)
            continue
        chosen_idx = rng.choice(pool.index.to_numpy(), size=take, replace=False)
        negatives = pool.loc[chosen_idx].copy()
        negatives["train_role"] = "negative"
        parts.append(positives)
        parts.append(negatives)
    if not parts:
        return labeled.iloc[0:0].copy()
    out = pd.concat(parts, ignore_index=True)
    out["train_label"] = out["train_role"].map({"positive": 1, "negative": 0}).astype(int)
    return out


def _max_pair(frame: pd.DataFrame, disease: str, control: str) -> pd.Series:
    return pd.concat([_numeric_col(frame, disease), _numeric_col(frame, control)], axis=1).max(axis=1)


def _agreement(frame: pd.DataFrame, disease: str, control: str) -> pd.Series:
    if disease not in frame.columns and control not in frame.columns:
        return pd.Series("", index=frame.index)
    parts = []
    for col in (disease, control):
        if col in frame.columns:
            parts.append(frame[col].fillna("").astype(str))
    stacked = pd.concat(parts, axis=1)
    numeric = stacked.apply(lambda col: pd.to_numeric(col, errors="coerce"))
    return numeric.max(axis=1).fillna(0).astype(int).astype(str)


def add_model_features(frame: pd.DataFrame) -> pd.DataFrame:
    """Columns the booster is allowed to see."""
    out = frame.copy()
    if "consensus_mei_family" in out.columns:
        out["mei_family"] = out["consensus_mei_family"].map(mei_family)
    elif "mei_family" not in out.columns:
        out["mei_family"] = "other"
    _side_max(out, "left_flank_mei_reads", "right_flank_mei_reads", "mei")
    _side_max(out, "left_flank_polya_reads", "right_flank_polya_reads", "polya")
    left = pd.to_numeric(out["flank_mei_left"], errors="coerce").fillna(0)
    right = pd.to_numeric(out["flank_mei_right"], errors="coerce").fillna(0)
    denom = pd.concat([left, right], axis=1).max(axis=1).replace(0, np.nan)
    out["flank_mei_balance"] = pd.concat([left, right], axis=1).min(axis=1) / denom
    window = pd.to_numeric(out.get("split_cluster_window_reads"), errors="coerce")
    clustered = pd.to_numeric(out.get("split_cluster_reads"), errors="coerce")
    out["split_cluster_read_fraction"] = clustered / window.replace(0, np.nan)
    out["poly_at_reads_max"] = _max_pair(out, "disease_poly_at_reads", "control_poly_at_reads")
    out["poly_at_max_run_max"] = _max_pair(out, "disease_poly_at_max_run", "control_poly_at_max_run")
    out["family_agreement"] = _agreement(out, "disease_family_agreement", "control_family_agreement")
    out["strand_agreement"] = _agreement(out, "disease_strand_agreement", "control_strand_agreement")
    for col in (
        "discordant_mei_majority",
        "two_sided_support",
        "poly_at_supported",
        "tsd_or_polyA_supported",
        "complex_mei_event",
        "complex_ins_with_del",
    ):
        if col in out.columns:
            out[col] = _numeric_col(out, col).fillna(0).astype(int).astype(str)
    repeat = out.get("nested_repeat_class", pd.Series("", index=out.index)).fillna("").astype(str)
    out["breakpoint_simple_repeat"] = repeat.str.contains("Simple_repeat|Low_complexity", case=False, regex=True).astype(int)
    for col in _JUNK_FLAG_COLUMNS:
        if col == "nested_repeat_class":
            continue
        if col not in out.columns:
            out[col] = np.nan
        elif col.startswith("flag_") or col.endswith("_count"):
            out[col] = _numeric_col(out, col)
    return out


def iter_gold_tables(sample_dir: Path) -> Iterable[Path]:
    """Per-chromosome gold-review tables under a genome output directory."""
    direct = sample_dir / GOLD_REVIEW_NAME
    if direct.exists():
        yield direct
        return
    for chrom_dir in sorted(p for p in sample_dir.iterdir() if p.is_dir()):
        path = chrom_dir / GOLD_REVIEW_NAME
        if path.exists():
            yield path


@dataclass(frozen=True)
class SampleGold:
    sample: str
    directory: Path
    exclude_chroms: frozenset[str] = frozenset()


def build_labeled_gold(samples: Sequence[SampleGold]) -> pd.DataFrame:
    """Load each sample's chromosomes, re-rank them as one table, and label."""
    parts: list[pd.DataFrame] = []
    for spec in samples:
        chrom_parts: list[pd.DataFrame] = []
        for path in iter_gold_tables(spec.directory):
            gold = load_gold_chromosome(path)
            if gold.empty:
                continue
            chrom = str(gold["chrom"].iloc[0]) if "chrom" in gold.columns and len(gold) else path.parent.name
            if chrom in spec.exclude_chroms or path.parent.name in spec.exclude_chroms:
                continue
            gold.insert(0, "source_table", str(path))
            chrom_parts.append(gold)
        if not chrom_parts:
            continue
        ranked = rank_gold_like_review(pd.concat(chrom_parts, ignore_index=True))
        labeled = label_gold_table(ranked)
        labeled.insert(0, "sample", spec.sample)
        parts.append(labeled)
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


def training_table_from_labels(labeled: pd.DataFrame, *, seed: int = 13) -> pd.DataFrame:
    sampled = sample_matched_negatives(labeled, seed=seed)
    return add_model_features(sampled)


def encode_feature_matrix(
    frame: pd.DataFrame,
    *,
    categorical_levels: dict[str, list[str]] | None = None,
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Numeric matrix plus integer codes for categorical columns.

    Levels are frozen on the first call and reused at score time. Unseen
    categories become missing.
    """
    featured = add_model_features(frame) if "flank_mei_left" not in frame.columns else frame
    levels = {} if categorical_levels is None else dict(categorical_levels)
    columns: dict[str, pd.Series] = {}
    for name in NUMERIC_FEATURES:
        if name not in featured.columns:
            columns[name] = pd.Series(np.nan, index=featured.index, dtype=float)
        else:
            columns[name] = pd.to_numeric(featured[name], errors="coerce")
    for name in CATEGORICAL_FEATURES:
        raw = featured[name].fillna("").astype(str) if name in featured.columns else pd.Series("", index=featured.index)
        if name not in levels:
            levels[name] = sorted({v for v in raw.unique() if v != ""})
        mapping = {level: i for i, level in enumerate(levels[name])}
        columns[name] = raw.map(mapping).astype("float")
    matrix = pd.DataFrame(columns, index=featured.index)
    return matrix, levels


def _categorical_mask() -> list[bool]:
    names = list(NUMERIC_FEATURES) + list(CATEGORICAL_FEATURES)
    categorical = set(CATEGORICAL_FEATURES)
    return [name in categorical for name in names]


def fit_gold_classifier(
    matrix: pd.DataFrame,
    labels: pd.Series,
    *,
    seed: int = 13,
):
    """Shallow regularized gradient-boosted trees for mixed feature types."""
    from sklearn.ensemble import HistGradientBoostingClassifier

    model = HistGradientBoostingClassifier(
        learning_rate=0.05,
        max_depth=3,
        max_iter=300,
        min_samples_leaf=20,
        l2_regularization=1.0,
        categorical_features=_categorical_mask(),
        random_state=seed,
    )
    model.fit(matrix.to_numpy(dtype=float), labels.to_numpy(dtype=int))
    return model


def leave_one_sample_out_metrics(
    training: pd.DataFrame,
    *,
    seed: int = 13,
) -> list[dict[str, object]]:
    """Fit on every sample but one and score the held-out sample."""
    from sklearn.metrics import average_precision_score, roc_auc_score

    rows: list[dict[str, object]] = []
    samples = list(training["sample"].astype(str).unique())
    for held in samples:
        train = training.loc[training["sample"].astype(str) != held]
        test = training.loc[training["sample"].astype(str) == held]
        if train["train_label"].nunique() < 2 or test["train_label"].nunique() < 2:
            rows.append({"held_out_sample": held, "n": int(len(test)), "roc_auc": None, "average_precision": None})
            continue
        matrix, levels = encode_feature_matrix(train)
        model = fit_gold_classifier(matrix, train["train_label"], seed=seed)
        scored, _ = encode_feature_matrix(test, categorical_levels=levels)
        proba = model.predict_proba(scored.to_numpy(dtype=float))[:, 1]
        y = test["train_label"].to_numpy(dtype=int)
        record: dict[str, object] = {
            "held_out_sample": held,
            "n": int(len(test)),
            "n_positive": int(y.sum()),
            "roc_auc": float(roc_auc_score(y, proba)),
            "average_precision": float(average_precision_score(y, proba)),
        }
        for family, sub_y, sub_p in _family_slices(test, proba):
            if len(set(sub_y)) < 2:
                continue
            record[f"roc_auc_{family}"] = float(roc_auc_score(sub_y, sub_p))
        rows.append(record)
    return rows


def _family_slices(test: pd.DataFrame, proba: np.ndarray):
    families = test["mei_family"].astype(str).to_numpy()
    y = test["train_label"].to_numpy(dtype=int)
    for family in ("Alu", "L1", "SVA"):
        mask = families == family
        yield family, y[mask], proba[mask]


def role_summary(labeled: pd.DataFrame) -> pd.DataFrame:
    """Counts of each training role by sample and family."""
    return (
        labeled.groupby(["sample", "mei_family", "train_role"], dropna=False)
        .size()
        .rename("n")
        .reset_index()
    )


def dump_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
