# PPI Network Drug Repurposing — Alzheimer's Disease

A network-medicine pipeline testing whether Alzheimer's Disease (AD) genetic-risk genes form a
statistically real neighborhood in the human protein-interaction network, and whether any drug approved
for a different disease has a target sitting unusually close to that neighborhood — a computable,
published drug-repurposing screening method (network proximity; Guney et al. 2016). Implemented in
**Python** (NetworkX, primary) and **R** (igraph, twin), both reading identical live-fetched real data.

## Aim

Can a purely computational, network-structure-based screen generate real, defensible drug-repurposing
hypotheses for Alzheimer's — using only public data and no wet-lab experiment?

## Objective

Build a real human protein-protein interaction network around Alzheimer's genetic-risk genes, statistically
test whether those genes form a real, non-random disease module in the network, and rank approved drugs
by how close their known targets sit to that module — cross-implemented in Python and R to confirm any
result isn't an artifact of one language's random-number generator or graph library.

## Data fetch

Two real, live-fetched data sources. The **Open Targets Platform GraphQL API** provides real AD
gene-disease genetic-association evidence scores (top 120 genes kept, by a disclosed cutoff) and real
drug-target-disease relationships for every gene in the resulting network. The **STRING API** provides
real, evidence-scored human protein-protein interactions at a high-confidence threshold (combined score
≥ 700), expanded via STRING's `add_nodes` parameter to include real interaction partners beyond the 120
seed genes — a necessary step for building a non-degenerate statistical null model (see Methods note
below).

## Data describe

A real 270-node, 2,002-edge protein-protein interaction network: 120 real Alzheimer's genetic-risk genes
(by Open Targets association score) plus their real STRING-derived interaction partners. Each edge
carries a real STRING confidence score blending multiple independent evidence types (physical binding,
co-expression, curated pathway databases, text-mining of the literature).

## Methods / Workflow — what we did

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

1. Fetch the top 120 real AD-associated genes from Open Targets and build the real STRING PPI network
   around them, expanded to 270 nodes via STRING's `add_nodes` partner-expansion.
2. Compute network centrality measures (degree, betweenness, eigenvector) and detect communities via
   Louvain modularity optimization.
3. **Disease-module test**: measure the largest connected component (LCC) formed by the 120 AD genes
   within the network, then build a **degree-matched null model** — 1,000 random substitutions of each
   real AD gene for a same-degree-bin alternative — to test whether the real LCC is larger than chance
   would predict, converting the comparison into a z-score and p-value.
4. Fetch real drug-target-disease relationships from Open Targets for every gene in the network.
5. **Network proximity drug ranking**: for every candidate drug with at least one known target in the
   network, compute the average shortest-path distance from its target(s) to the AD module (Guney et
   al.'s "closest" measure), compare against 200 degree-matched null draws, and rank drugs by resulting
   z-score.
6. Repeat the entire pipeline independently in R (igraph) on the identical real data, to cross-validate
   the Python (NetworkX) implementation.

## Results

| Metric | Python | R |
|---|---|---|
| Network size | 270 nodes, 2002 edges | 270 nodes, 2002 edges |
| Modularity (Louvain) | 0.509 (50 communities) | 0.508 (50 communities) |
| Disease-module real LCC | 46 | 46 |
| Degree-matched null mean (sd) | 55.18 (5.06) | 54.78 (5.31) |
| **z-score / p-value** | **z = -1.81, p = 0.961** | **z = -1.65, p = 0.948** |

The disease-module hypothesis is **not supported** under this network construction — real AD genes
cluster *less* than a degree-matched random comparison, not more. Both independent language
implementations agree closely, which is itself evidence the finding is a real property of the data/
method, not an artifact of one random draw.

**Top repurposing candidates (both languages)**: dominated by CSF1R inhibitors (pexidartinib,
emactuzumab, cabiralizumab...) and EGFR-pathway inhibitors (erlotinib, afatinib, lapatinib...).

## Biology interpretation of results

The disease-module test came back non-significant and in the *opposite* direction from the standard
network-medicine prediction (disease genes are expected to cluster more than chance, not less). This is
not treated as a code bug — two independent implementations (different graph libraries, different random
number generators, different Louvain implementation) agree closely on both the direction and rough
magnitude of the negative result, which is real evidence the finding reflects the actual data and method,
not a coding error in one language. Biologically, this complicates the simple "disease module" story in a
genuinely informative way: Alzheimer's genetic risk may span several loosely-connected biological systems
(amyloid processing, tau pathology, lipid transport via APOE, innate immunity/microglial signaling via
TREM2 and CSF1R-adjacent genes) that don't physically interact much in a protein-interaction network even
though they independently converge on the same clinical phenotype — a disease can be biologically real
without presenting as one tightly interconnected network hub. A separate, real methodological
explanation is equally plausible and disclosed rather than hidden: STRING's `add_nodes` expansion
specifically selects partner genes that *maximize* connectivity to the seed set, meaning the population
populating the null model's degree bins is systematically biased toward generic, highly promiscuous hub
proteins (EGFR, SRC, GRB2, TNF) rather than a neutral cross-section of the interactome — this likely
inflates the null's typical connectivity, making the real AD-gene module look smaller by comparison than
it would against a genuinely unbiased background. Given the module test's non-significant result, the
drug-repurposing ranking is honestly read as **exploratory hypothesis generation**, not validated
repurposing evidence. The top-ranked candidates (CSF1R and EGFR-pathway inhibitors) all have
`real_distance = 0`, meaning their drug targets are themselves already members of the 120-gene AD list,
not genuine network *neighbors* of it — a real, useful, but different and weaker finding ("this target is
itself AD-associated") than "network proximity reveals a hidden connection to a target not previously
linked to AD." Both CSF1R (an active real research target for microglial modulation in Alzheimer's
trials) and EGFR-pathway crosstalk with amyloid processing are legitimate, independently corroborated
research angles — so the method surfaced biologically real genes, even though the specific ranking
mechanism that surfaced them (distance-zero always wins a z-score comparison) is a less interesting
reason than a genuine hidden-neighbor discovery would be.

## Learning through project

A non-significant or direction-reversed result, cross-validated across two independent implementations,
is a real, reportable scientific finding — not a failure to hide or a result to keep redesigning around
until it looks like the textbook expectation. The choice of how a network's background/null population is
constructed is not a neutral technical detail — STRING's `add_nodes` parameter optimizing for maximum
connectivity to the seed set, rather than sampling a neutral background, measurably shapes what the
disease-module test can conclude, and naming that limitation explicitly (rather than treating a
clean-looking null as automatically unbiased) is what makes the negative result trustworthy rather than
just convenient. A ranking metric can technically "work" (surface real, biologically relevant genes) for
a less interesting reason than intended — here, distance-zero (a drug's target being an AD gene itself)
mechanically outranks any genuine 1-hop-or-further neighbor candidate in a z-score comparison, which is
worth naming explicitly rather than presenting the top-ranked list as if it demonstrated genuine network
proximity discovery.

## Limitations

Network proximity is a screening heuristic, not causal or clinical evidence. STRING confidence scores
blend genuinely different evidence types (physical binding, co-expression, text-mining) that don't carry
equal biological weight. Open Targets evidence has literature/research-attention bias. The network was
deliberately scoped (120 seed genes plus STRING's top ~150 connectivity-maximizing partners) for
tractability, not the full ~19,000-protein interactome, and that scoping choice plausibly affects the
null model as described above. No result here is a treatment recommendation.

## Reproduce

**Python:** run `ppi_drug_repurposing.ipynb` top to bottom (installs its own dependencies). Requires
internet access (live Open Targets + STRING API queries).

**R:** run the Python notebook first (produces the `data_py/` bridge files), then run
`ppi_drug_repurposing.R`.

```r
install.packages("igraph")
```

## Tech

`Python` (NetworkX) · `R` (igraph) · Open Targets GraphQL API · STRING API · Louvain community
detection · degree-matched null models · network proximity (Guney et al. 2016)

## Files

```
ppi_drug_repurposing.ipynb   # Python: fetch → build network → centrality/communities → module test → repurposing ranking
ppi_drug_repurposing.R       # R twin: same pipeline on the same real bridge data, igraph
results_py/, results_R/      # centrality, communities, module test, ranked drug candidates — CSV, per language
```

## License

All rights reserved — see `LICENSE`. This repository is public for portfolio/demonstration purposes
only; no permission is granted to copy, modify, or reuse any part of it.
