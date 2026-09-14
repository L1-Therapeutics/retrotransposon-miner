# Scientific Architecture & Mathematical Foundation

`retrotransposon-miner` is a high-throughput genomic analysis pipeline designed to detect, refine, genotype, and classify Mobile Element Insertions (MEIs)—including LINE-1 (L1HS), Alu (AluYa5, AluYb8), and SVA elements—from paired-end high-throughput sequencing data (BAM/CRAM format).

---

## 1. Evidence Extraction & Streaming Architecture

### 1.1 One-Pass BAM Streaming Engine
Traditional structural variant discovery algorithms utilize multi-pass BAM parsing, iterating through alignment files multiple times to collect discordant read pairs ($DP$) and split-reads ($SR$). `retrotransposon-miner` implements a single-pass streaming architecture using `pysam` C-bindings to extract both evidence classes simultaneously.

- **Discordant Pairs ($DP$):** Paired reads where mate mapping distance $|TLEN| > \mu_{TLEN} + 3\sigma_{TLEN}$ or mates map to non-syntenic chromosomes.
- **Split-Reads ($SR$):** Alignments containing soft-clip CIGAR operators (`S` / `4`) with clip length $L_{\text{clip}} \ge 15\text{ bp}$.

### 1.2 Memory Complexity & Key Isolation
*Literature Anchor: Cameron et al. (2017) GRIDSS, Genome Research.*

To prevent $O(N^2)$ memory growth in single-pass mate tracking, unmatched mate pairs are stored in a bounded hash buffer indexed by composite tuple keys:

$$K_{\text{mate}} = (\text{QNAME}, \text{RG}_{\text{ID}})$$

This prevents query-name collisions across merged multiplexed BAM libraries while maintaining $O(N)$ space complexity:

$$\text{Memory Overhead} \le O(M_{\text{unmatched}}) \ll O(N_{\text{total}})$$

---

## 2. Target Site Duplication (TSD) & Poly(A) Tail Refinement

### 2.1 Shannon Entropy Poly(A) Tail Gating
*Literature Anchor: Gardner et al. (2017) MELT, Nucleic Acids Research.*

Target-Primed Reverse Transcription (TPRT) mediated integrations leave characteristic 3' poly(A) or 5' poly(T) homopolymer tracts. Spurious soft-clips from low-complexity reference genomic repeats are filtered by calculating Shannon sequence entropy over sliding windows of size $W = 10\text{ bp}$:

$$H(X) = -\sum_{i \in \{A,C,G,T\}} P(x_i) \log_2 P(x_i)$$

Where $P(x_i)$ represents the empirical nucleotide frequency within the window.

- **Low-Entropy Poly(A) Filter:** Sequences with $H(X) < 1.2$ and nucleotide fraction $\max(P(A), P(T)) \ge 0.75$ are flagged as canonical retrotransposition tails (`POLYA=1`).

### 2.2 Micro-Homology Aware TSD Boundary Resolution
TPRT generates 10–25 bp Target Site Duplications (TSDs) flanking the insertion locus. When micro-homology exists between the reference site and the non-reference MEI insertion, standard alignment engines misplace breakpoint coordinates by 2–10 bp.

The TSD refiner evaluates 5' (upstream) and 3' (downstream) soft-clips using dynamic programming matrix $M_{i,j}$ to locate the longest exact common substring:

$$M_{i,j} = \begin{cases} M_{i-1, j-1} + 1 & \text{if } S_1[i] = S_2[j] \\ 0 & \text{otherwise} \end{cases}$$

Candidate TSDs are validated under biological bounds $5\text{ bp} \le L_{\text{TSD}} \le 35\text{ bp}$.

---

## 3. Bayesian Diploid MEI Genotyping Engine

*Literature Anchor: Li (2011) SAMtools / Garrison & Marth (2012) FreeBayes.*

### 3.1 Diploid Likelihood Formulation
Let $k_{\text{alt}}$ be the count of non-reference supporting evidence reads (split-reads + discordant pairs) and $k_{\text{ref}}$ be the count of concordant reference reads spanning the breakpoint locus. Total depth is $n = k_{\text{alt}} + k_{\text{ref}}$.

We evaluate binomial log-likelihoods for candidate diploid genotypes $G \in \{0/0, 0/1, 1/1\}$:

$$P(k_{\text{alt}} \mid n, G) = \binom{n}{k_{\text{alt}}} p_G^{k_{\text{alt}}} (1 - p_G)^{n - k_{\text{alt}}}$$

Under expected alternate allele proportions:
- **Homozygous Reference ($G = 0/0$):** $p_{0/0} = \epsilon$ (sequencing/mapping error rate, default $\epsilon = 0.02$)
- **Heterozygous Insertion ($G = 0/1$):** $p_{0/1} = 0.50$
- **Homozygous Non-Reference ($G = 1/1$):** $p_{1/1} = 1 - \epsilon = 0.98$

### 3.2 Posterior Probability & Phred Genotype Quality (GQ)
Assuming uniform prior probabilities $P(G) = \frac{1}{3}$, normalized posteriors are computed via log-sum-exp:

$$P(G_i \mid k_{\text{alt}}, n) = \frac{P(k_{\text{alt}} \mid n, G_i)}{\sum_{j} P(k_{\text{alt}} \mid n, G_j)}$$

The Phred-scaled Genotype Quality ($GQ$) quantifies confidence in the maximum a posteriori (MAP) call $G_{\text{max}}$:

$$GQ = \min\left(99.0, -10 \cdot \log_{10}\left(1 - P(G_{\text{max}} \mid k_{\text{alt}}, n)\right)\right)$$

Variant Allele Frequency ($VAF$) is reported directly as:

$$VAF = \frac{k_{\text{alt}}}{k_{\text{alt}} + k_{\text{ref}}}$$

---

## 4. Subfamily Classification via Diagnostic k-mer LLR

Disease-causing somatic retrotransposition in humans is driven almost exclusively by active, young element subfamilies ($L1HS$, $AluYa5$, $AluYb8$, $SVA\text{-}F$), whereas the genome contains $> 10^6$ fixed, truncated ancestral elements ($L1PA2\text{--}L1PA16$, $AluS$, $AluJ$).

The voting engine computes Log-Likelihood Ratios ($LLR$) over diagnostic nucleotide substitution motifs (e.g., L1HS 3' UTR 6015 ACA/GAG polymorphism):

$$LLR = \ln \frac{\mathcal{L}(S_{\text{active}} \mid \text{Reads})}{\mathcal{L}(S_{\text{ancestral}} \mid \text{Reads})} = \sum_{k \in K_{\text{diag}}} \text{count}(k) \cdot \ln(2.0)$$

Calls achieving $LLR \ge 2.0$ are annotated with specific active subfamilies, distinguishing polymorphic insertions from legacy reference alignment noise.

---

## 5. De Bruijn Micro-Assembly of Junction Unitigs

*Literature Anchor: Cameron et al. (2017) GRIDSS, Genome Research.*

To confirm complex insertion junctions without long-read sequencing, soft-clipped read fragments at candidate loci are decomposed into a directed De Bruijn graph $G = (V, E)$ using $k$-mer size $k = 15$.

1. **Spectrum Extraction:** Extract $k$-mers $K_k(S)$ across all soft-clips at the locus.
2. **Error Pruning:** Remove nodes with $k$-mer coverage $C(k) < 2$ to eliminate PCR and sequencing noise.
3. **Unitig Traversal:** Traverse non-branching graph paths from source nodes ($\text{in-degree} = 0$) to construct contiguous breakpoint unitigs.

---

## 6. VCF v4.3 Structural Variant Serialization

Callsets are output according to standard VCF v4.3 structural variant specifications using symbolic ALT alleles (`<INS:MEI:L1HS>`, `<INS:MEI:ALU>`, `<INS:MEI:SVA>`).

### Standard Meta-Information Tags:
- `INFO/SVTYPE`: Structural variant type (`INS`).
- `INFO/MEI_TYPE`: Predicted subfamily call.
- `INFO/TSD`: Refined Target Site Duplication sequence.
- `INFO/TSDLEN`: TSD length in base pairs.
- `INFO/POLYA`: Flag indicating 3' poly(A) tail detection (`1` or `0`).
- `INFO/MEI_LLR`: Subfamily log-likelihood ratio.
- `FORMAT/GT`: Diploid genotype (`0/0`, `0/1`, `1/1`, `./.`).
- `FORMAT/GQ`: Phred-scaled Genotype Quality score.
- `FORMAT/VAF`: Alternate Variant Allele Frequency.
- `FORMAT/AD`: Reference and Alternate read depths ($k_{\text{ref}}, k_{\text{alt}}$).
