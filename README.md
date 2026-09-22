# PPI Network Drug Repurposing — Alzheimer's Disease

A network-medicine pipeline testing whether Alzheimer's Disease (AD) genetic-risk genes form a
statistically real neighborhood in the human protein-interaction network, and whether any drug
approved for a different disease has a target sitting unusually close to that neighborhood — a
computable, published drug-repurposing screening method (network proximity; Guney et al. 2016).
Implemented in **Python** (NetworkX, primary) and **R** (igraph, twin), both reading identical
live-fetched real data.

## Research question

Can a purely computational, network-structure-based screen generate real, defensible drug-repurposing
hypotheses for Alzheimer's — using only public data and no wet-lab experiment?

## Data (real, fetched live)

- **Open Targets Platform GraphQL API** — real AD gene-disease genetic-association evidence scores
  (120 genes kept, a disclosed cutoff), and real drug-target-disease relationships for every gene in
  the network.
- **STRING API** — real, evidence-scored human protein-protein interactions, high-confidence threshold
  (combined score ≥ 700), expanded via STRING's `add_nodes` parameter to include real interaction
  partners beyond the 120 seed genes (required for a non-degenerate null model — see Methods note).

## Pipeline

```
Open Targets (Alzheimer's EFO term) — live query, genetic-association evidence
   │ keep top 120 genes by genetic-association score (disclosed cutoff)
   ▼
STRING API (high confidence, + add_nodes real partner expansion) → PPI graph (270 genes, 2002 edges)
   ▼
   centrality (degree/betweenness/eigenvector) + Louvain community detection
   ▼
   disease-module hypothesis test: degree-matched null model → z-score / p-value
   ▼
Open Targets — real drug-target-disease data for every gene in the network
   ▼
   network proximity (Guney et al. "closest" measure) vs. degree-matched null → rank repurposing candidates
   → export: centrality, communities, module test, ranked drug candidates — all CSV, both languages
```

## Results

| Metric | Python | R |
|---|---|---|
| Network size | 270 nodes, 2002 edges | 270 nodes, 2002 edges |
| Modularity (Louvain) | 0.509 (50 communities) | 0.508 (50 communities) |
| Disease-module real LCC | 46 | 46 |
| Degree-matched null mean (sd) | 55.18 (5.06) | 54.78 (5.31) |
| **z-score / p-value** | **z = -1.81, p = 0.961** | **z = -1.65, p = 0.948** |

The disease-module hypothesis is **not supported** under this network construction — the real AD
genes cluster *less* than a degree-matched random comparison, not more. Both independent language
implementations agree closely (same real LCC by construction, consistent null statistics, same
non-significant conclusion), which is itself evidence the finding is a property of the data/method,
not an artifact of one random draw.

**Top repurposing candidates (both languages):** dominated by CSF1R inhibitors (pexidartinib,
emactuzumab, cabiralizumab...) and EGFR-pathway inhibitors (erlotinib, afatinib, lapatinib...) — both
directions are real, actively discussed research angles for AD (microglial/neuroinflammatory
modulation via CSF1R; amyloid-EGFR pathway crosstalk). Several top candidates' targets are themselves
part of the 120-gene AD module, giving `real_distance = 0` by construction — flagged explicitly as a
different, weaker claim ("this target is itself AD-associated") than "this target is a genuine network
neighbor of the module."

## Interpretation

Given the disease-module test came back non-significant, the drug ranking is best read as **exploratory
hypothesis generation**, not validated repurposing evidence. A plausible, disclosed reason for the
non-significant module test: STRING's `add_nodes` parameter selects partner genes specifically because
they maximize connectivity to the seed set, which likely biases the null model's background population
toward generic, highly promiscuous hub proteins (EGFR, SRC, GRB2, TNF...) rather than a neutral random
sample — inflating the null's typical connectivity relative to a truly unbiased background. This is a
genuine methodological limitation of the network-scoping choice, not a data-quality issue, and is
reported honestly rather than smoothed over.

## Limitations

Network proximity is a screening heuristic, not causal or clinical evidence. STRING confidence scores
blend genuinely different evidence types (physical binding, co-expression, text-mining) that don't
carry equal biological weight. Open Targets evidence has literature/research-attention bias. The
network was deliberately scoped (120 seed genes + STRING's top ~150 connectivity-maximizing partners)
for tractability, not the full ~19,000-protein interactome, and that scoping choice plausibly affects
the null model as described above. No result here is a treatment recommendation.

## Files

```
ppi_drug_repurposing.ipynb   # Python: fetch → build network → centrality/communities → module test → repurposing ranking
ppi_drug_repurposing.R       # R twin: same pipeline on the same real bridge data, igraph
results_py/, results_R/      # centrality, communities, module test, ranked drug candidates — CSV, per language
```

## Run

**Python:** run `ppi_drug_repurposing.ipynb` top to bottom (installs its own dependencies). Requires
internet access (live Open Targets + STRING API queries).

**R:** run the Python notebook first (produces the `data_py/` bridge files), then run
`ppi_drug_repurposing.R`.

```r
install.packages("igraph")
```
