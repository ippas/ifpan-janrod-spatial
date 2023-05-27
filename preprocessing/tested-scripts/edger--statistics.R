require(edgeR)

# prepare data
risperidone_summary_statistics$raw_data$resolution_0.1$cluster_5$control$mean -> control_mean

risperidone_summary_statistics$raw_data$resolution_0.1$cluster_5$experiment$mean -> experiment_mean

counts <- cbind(control_mean, experiment_mean)

group <- factor(c(rep("Control", 6), rep("Treatment", 6)))

y <- DGEList(counts = counts, group = group)

keep <- filterByExpr(y, min.count = 0)
y <- y[keep, ]

y <- calcNormFactors(y)

y <- estimateDisp(y)

et <- exactTest(y)

topTags(et, n = Inf) %>%
  as.data.frame() %>%
  rownames_to_column(var = "peak_id") %>%
  left_join(., risperidone_st_data_half$raw_data$annotate[, c("peak_id", "gene_name")], by = "peak_id")


y <- calcNormFactors(y)

# View normalization factors
norm.factors <- y$samples$norm.factors
print(norm.factors)

# Calculate normalized counts
norm_counts <- t( t(y$counts) / y$samples$norm.factors )

# View normalized counts
norm_counts %>% head

counts %>% head
group

control_counts <- norm_counts[, group == "Control"]
experiment_counts <- norm_counts[, group == "Treatment"]

# Create a matrix of p-values
p_values <- matrix(nrow = nrow(norm_counts))
# Perform t-tests
for (i in seq_len(nrow(norm_counts))) {
  p_values[i] <- t.test(control_counts[i, ], experiment_counts[i, ], var.equal = T)$p.value
}

p_values

data.frame(rownames(control_counts), p_values) %>%
  set_colnames(c("peak_id", "p_values")) %>%
  left_join(., risperidone_st_data_half$raw_data$annotate[, c("peak_id", "gene_name")], by = "peak_id") %>% 
  filter(p_values < 0.05)
