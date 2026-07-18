"""Family-consistent panel→full MEI coordinate projection.

Guards the regression where Alu tip intervals were unioned with L1 body
projections (or used alone when L1 fragments were unmapped), producing bogus
“twin-truncation near L1 5′” spans on an Alu-length axis.
"""

from __future__ import annotations

import pandas as pd
import pytest

from retro_miner.mei_support import (
    FragmentToFullMap,
    _full_axis_union_from_panel_fragments,
    _project_panel_coords_to_full,
)


def _entry(
    frag: str,
    full: str,
    *,
    frag_start: int,
    frag_end: int,
    full_start: int,
    full_end: int,
    family: str,
    frag_len: int = 0,
    full_len: int = 0,
) -> FragmentToFullMap:
    return FragmentToFullMap(
        fragment_name=frag,
        fragment_length=frag_len or frag_end,
        full_name=full,
        full_length=full_len or full_end,
        fragment_aln_start=frag_start,
        fragment_aln_end=frag_end,
        full_aln_start=full_start,
        full_aln_end=full_end,
        strand="+",
        family=family,
    )


@pytest.fixture
def mixed_frag_map() -> dict[str, FragmentToFullMap]:
    """Map with L1 and Alu axes (as in real mei_fragment_to_full_coords.tsv)."""
    entries = [
        _entry(
            "L1MC_orf2",
            "L1HS_full",
            frag_start=1,
            frag_end=1000,
            full_start=3000,
            full_end=3999,
            family="LINE1",
            full_len=6049,
        ),
        _entry(
            "AluSc5",
            "AluSc5_full",
            frag_start=1,
            frag_end=312,
            full_start=1,
            full_end=312,
            family="ALU",
            full_len=312,
        ),
        _entry(
            "L1MCa_5end",
            "L1MA1_full",
            frag_start=791,
            frag_end=2434,
            full_start=96,
            full_end=1691,
            family="LINE1",
            full_len=5544,
        ),
    ]
    out: dict[str, FragmentToFullMap] = {}
    for entry in entries:
        out[entry.fragment_name] = entry
    return out


class TestProjectPanelCoordsToFull:
    def test_linear_plus_strand(self, mixed_frag_map):
        got = _project_panel_coords_to_full(100, 200, "L1MC_orf2", mixed_frag_map)
        assert got == (3099, 3199, "L1HS_full")


class TestFullAxisUnionFamilySafe:
    def test_rejects_alu_when_call_is_line1(self, mixed_frag_map):
        """L1MC_orf2 body + Alu tip must not smash into 157–3971-style unions."""
        row = pd.Series(
            {
                "consensus_mei_family": "LINE1",
                "disease_L_mei_subfamily": "L1MC_orf2",
                "disease_L_mei_start": 120,
                "disease_L_mei_end": 971,
                "disease_R_mei_subfamily": "AluSc5",
                "disease_R_mei_start": 157,
                "disease_R_mei_end": 186,
            }
        )
        union = _full_axis_union_from_panel_fragments(row, mixed_frag_map)
        assert union is not None
        assert union == (3119, 3970)

    def test_unmapped_l1_does_not_fall_back_to_alu_axis(self, mixed_frag_map):
        """Empty L1 map must not project AluSc5 onto a fake L1 full span."""
        # Drop L1MCa mapping so only Alu can project.
        frag_map = {
            k: v for k, v in mixed_frag_map.items() if not k.startswith("L1MCa")
        }
        row = pd.Series(
            {
                "consensus_mei_family": "LINE1",
                "consensus_mei_subfamily": "L1MCa_5end",
                "disease_L_mei_subfamily": "L1MCa_5end",
                "disease_L_mei_start": 136,
                "disease_L_mei_end": 1009,
                "disease_R_mei_subfamily": "AluSc5",
                "disease_R_mei_start": 136,
                "disease_R_mei_end": 309,
            }
        )
        assert _full_axis_union_from_panel_fragments(row, frag_map) is None

    def test_mapped_l1_5end_uses_l1_axis_only(self, mixed_frag_map):
        row = pd.Series(
            {
                "consensus_mei_family": "LINE1",
                "disease_L_mei_subfamily": "L1MCa_5end",
                "disease_L_mei_start": 800,
                "disease_L_mei_end": 1200,
                "disease_R_mei_subfamily": "AluSc5",
                "disease_R_mei_start": 10,
                "disease_R_mei_end": 50,
            }
        )
        union = _full_axis_union_from_panel_fragments(row, mixed_frag_map)
        assert union is not None
        # 800–1200 on L1MCa_5end → offset from frag 791 → full 96+(800-791) .. 96+(1200-791)
        assert union == (105, 505)

    def test_alu_call_keeps_alu_axis(self, mixed_frag_map):
        row = pd.Series(
            {
                "consensus_mei_family": "ALU",
                "disease_L_mei_subfamily": "AluSc5",
                "disease_L_mei_start": 10,
                "disease_L_mei_end": 100,
                "disease_R_mei_subfamily": "L1MC_orf2",
                "disease_R_mei_start": 120,
                "disease_R_mei_end": 200,
            }
        )
        union = _full_axis_union_from_panel_fragments(row, mixed_frag_map)
        assert union == (10, 100)

    def test_prefers_consensus_subfamily_full_axis(self):
        """L1HS_5end + L1HS_3end must union on L1HS_full, not a sibling L1PA* axis."""
        frag_map = {
            "L1HS_5end": _entry(
                "L1HS_5end",
                "L1PA2_full",
                frag_start=1,
                frag_end=2136,
                full_start=1,
                full_end=2135,
                family="LINE1",
                full_len=6045,
            ),
            "L1HS_3end": _entry(
                "L1HS_3end",
                "L1HS_full",
                frag_start=1,
                frag_end=902,
                full_start=5130,
                full_end=6031,
                family="LINE1",
                full_len=6049,
            ),
        }
        row = pd.Series(
            {
                "consensus_mei_family": "LINE1",
                "consensus_mei_subfamily": "L1HS_5end#LINE/L1",
                "disease_L_mei_subfamily": "L1HS_5end",
                "disease_L_mei_start": 1,
                "disease_L_mei_end": 840,
                "disease_R_mei_subfamily": "L1HS_3end",
                "disease_R_mei_start": 483,
                "disease_R_mei_end": 889,
            }
        )
        # Only the L1HS_full axis matches the consensus subfamily preference.
        assert _full_axis_union_from_panel_fragments(row, frag_map) == (5612, 6018)
