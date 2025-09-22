# ##############################################################################
# ---- investigation on example Arc ----
# ##############################################################################
res_all$results$cluster_1$anova$global %>% filter(gene_symbol == "Arc")

res_all$results$cluster_1$anova$global %>% filter(gene_symbol == "Arc")

res_all$results$cluster_1$anova$main_first %>% filter(gene_symbol == "Arc")

res_all$results$cluster_1$anova$main_second %>% filter(gene_symbol == "Arc")

res_all$results$cluster_1$posthoc$`SECOND_risperidone-saline|FIRST_wtwt` %>% filter(gene_symbol == "Arc")
res_all$results$cluster_1$posthoc$`SECOND_risperidone-saline|FIRST_wtdel` %>% filter(gene_symbol == "Arc")
res_all$results$cluster_1$posthoc$`FIRST_wtdel-wtwt|SECOND_risperidone` %>% filter(gene_symbol == "Arc")
res_all$results$cluster_1$posthoc$`FIRST_wtdel-wtwt|SECOND_saline`%>% filter(gene_symbol == "Arc")


res_edger_MarginalProportional$results$cluster_1$anova$global %>% filter(gene_symbol == "Arc")
res_edger_MarginalProportional$results$cluster_1$anova$main_first %>% filter(gene_symbol == "Arc")
res_edger_MarginalProportional$results$cluster_1$anova$main_second %>% filter(gene_symbol == "Arc")

res_edger_MarginalProportional$results$cluster_1$posthoc$`SECOND_risperidone-saline|FIRST_wtwt` %>% filter(gene_symbol == "Arc")
res_edger_MarginalProportional$results$cluster_1$posthoc$`SECOND_risperidone-saline|FIRST_wtdel` %>% filter(gene_symbol == "Arc")
res_edger_MarginalProportional$results$cluster_1$posthoc$`FIRST_wtdel-wtwt|SECOND_risperidone` %>% filter(gene_symbol == "Arc")
res_edger_MarginalProportional$results$cluster_1$posthoc$`FIRST_wtdel-wtwt|SECOND_saline`%>% filter(gene_symbol == "Arc")



res_edger_MarginalProportional$results$cluster_1$anova$global %>% filter(FDR < 0.05)
