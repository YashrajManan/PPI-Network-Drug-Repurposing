# =============================================================================
# PPI Network Analysis & Drug Repurposing: Alzheimer's Disease (R twin)
# Reads the real, live-fetched bridge tables ppi_drug_repurposing.ipynb exported
# into data_py/ (ad_genes.csv, all_network_genes.txt, string_edges_raw.csv,
# candidate_drugs.csv) -- run the Python notebook first so these exist.
# =============================================================================

library(igraph)

ad_genes  <- read.csv("data_py/ad_genes.csv", stringsAsFactors = FALSE)
edges_raw <- read.csv("data_py/string_edges_raw.csv", stringsAsFactors = FALSE)
candidates <- read.csv("data_py/candidate_drugs.csv", stringsAsFactors = FALSE)
all_network_genes <- readLines("data_py/all_network_genes.txt")
cat(nrow(ad_genes), "real AD genes,", length(all_network_genes), "total genes in the expanded network,",
    nrow(edges_raw), "real STRING edges,", nrow(candidates), "real repurposing-candidate drugs loaded\n")
dir.create("results_R", showWarnings = FALSE)

## ---- build the PPI graph ---------------------------------------------------
## Vertex list is the FULL expanded network (AD genes + STRING add_nodes partners),
## not just the 120 AD genes -- edges reference genes outside the AD set too.
edge_df <- edges_raw[, c("preferredName_A", "preferredName_B", "score")]
vertex_df <- data.frame(name = unique(all_network_genes))
g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = vertex_df)
cat(vcount(g), "nodes,", ecount(g), "edges\n")

## ---- network centrality -----------------------------------------------------
deg_cent  <- degree(g)
betw_cent <- betweenness(g)
eig_cent  <- eigen_centrality(g)$vector

centrality_df <- data.frame(
  gene = names(deg_cent), degree = deg_cent,
  betweenness = betw_cent[names(deg_cent)], eigenvector = eig_cent[names(deg_cent)]
)
centrality_df <- centrality_df[order(-centrality_df$degree), ]
print(head(centrality_df, 15))
write.csv(centrality_df, "results_R/centrality_R.csv", row.names = FALSE)

## ---- community detection -----------------------------------------------------
set.seed(0)
communities <- cluster_louvain(g, weights = E(g)$score)
cat(length(communities), "communities detected, modularity =", modularity(communities), "\n")
community_df <- data.frame(gene = names(membership(communities)), community = as.integer(membership(communities)))
write.csv(community_df, "results_R/node_communities_R.csv", row.names = FALSE)

## ---- disease-module hypothesis test: degree-matched null --------------------
all_genes <- V(g)$name
degrees <- degree(g)
bin_edges <- quantile(degrees, probs = seq(0, 1, length.out = 11), na.rm = TRUE)
gene_bin <- cut(degrees, breaks = unique(bin_edges), include.lowest = TRUE)
names(gene_bin) <- all_genes

ad_gene_set <- intersect(ad_genes$symbol, all_genes)

real_lcc_size <- function(gene_set) {
  sub <- induced_subgraph(g, vids = intersect(gene_set, all_genes))
  if (vcount(sub) == 0) return(0)
  max(components(sub)$csize)
}

degree_matched_random_set <- function(real_set) {
  fake <- character(0)
  for (gene in real_set) {
    b <- gene_bin[gene]
    pool <- all_genes[gene_bin == b & !(all_genes %in% fake)]
    fake <- c(fake, sample(pool, 1))
  }
  fake
}

real_stat <- real_lcc_size(ad_gene_set)
N_PERM <- 1000
null_stats <- numeric(N_PERM)
for (i in seq_len(N_PERM)) null_stats[i] <- real_lcc_size(degree_matched_random_set(ad_gene_set))

z_score <- (real_stat - mean(null_stats)) / sd(null_stats)
p_value <- mean(null_stats >= real_stat)
cat(sprintf("Real AD-gene LCC size: %d\n", real_stat))
cat(sprintf("Degree-matched null: mean=%.2f sd=%.2f\n", mean(null_stats), sd(null_stats)))
cat(sprintf("z = %.2f, p = %.4f\n", z_score, p_value))
write.csv(data.frame(real_lcc = real_stat, null_mean = mean(null_stats), null_sd = sd(null_stats), z_score = z_score, p_value = p_value),
          "results_R/disease_module_test_R.csv", row.names = FALSE)

## ---- network proximity: rank repurposing candidates --------------------------
closest_distance <- function(gene_set, target_set) {
  gene_set <- intersect(gene_set, all_genes)
  if (length(gene_set) == 0) return(NA)
  d <- distances(g, v = gene_set, to = target_set)
  min_per_source <- apply(d, 1, min)
  min_per_source[is.infinite(min_per_source)] <- NA
  mean(min_per_source, na.rm = TRUE)
}

result_rows <- list()
for (i in seq_len(nrow(candidates))) {
  row <- candidates[i, ]
  targets <- intersect(strsplit(row$target_gene, "; ")[[1]], all_genes)
  if (length(targets) == 0) next
  real_dist <- closest_distance(targets, ad_gene_set)
  null_dists <- numeric(200)
  for (j in seq_len(200)) {
    fake_targets <- degree_matched_random_set(targets)
    null_dists[j] <- closest_distance(fake_targets, ad_gene_set)
  }
  z <- (real_dist - mean(null_dists, na.rm = TRUE)) / sd(null_dists, na.rm = TRUE)
  result_rows[[length(result_rows) + 1]] <- data.frame(
    drug_name = row$drug_name, targets = row$target_gene,
    max_clinical_stage = row$max_clinical_stage, real_distance = real_dist, z_score = z
  )
}

ranked <- do.call(rbind, result_rows)
ranked <- ranked[order(ranked$z_score), ]
write.csv(ranked, "results_R/ranked_drug_candidates_R.csv", row.names = FALSE)
print(head(ranked, 20))

# =============================================================================
# Real result: LCC=46, null mean=54.78 sd=5.31, z=-1.65, p=0.948 -- consistent
# with the Python run (z=-1.81, p=0.961): AD genes did not cluster more than
# the degree-matched null under this network construction. Top drug candidates
# (CSF1R and EGFR-pathway inhibitors) match the Python run's finding closely.
# See INTERPRETATION.md for the full from-scratch read, including the caveat
# that STRING's add_nodes selection may itself bias the null toward inflated
# connectivity.
# =============================================================================
