res_all$results$cluster_0$anova$global %>% 
  ggplot(aes(x = ))


res_all$results$cluster_0$posthoc$`SECOND_risperidone-saline|FIRST_wtdel`%>% 
  set_colnames(c("name_id", "peak_id", "gene_symbol", "salRisDel_logFC", "salRisDel_logCPM", "salRisDel_F", "salRisDel_PValue", "salRisDel_FDR"))

res_all$results$cluster_0$posthoc$`SECOND_risperidone-saline|FIRST_wtwt` %>% 
  set_colnames(c("name_id", "peak_id", "gene_symbol", "salRisWt_logFC", "salRisWt_logCPM", "salRisWt_F", "salRisWt_PValue", "salRisWt_FDR")) %>% dim


df %>% 
  select(cluster, name_id) %>% 
  filter(cluster == "cluster_2")




library(dplyr)
library(ggplot2)
library(ggrepel)

res_all$results$cluster_0$posthoc$`SECOND_risperidone-saline|FIRST_wtdel` %>% 
  set_colnames(c("name_id", "peak_id", "gene_symbol", 
                 "salRisDel_logFC", "salRisDel_logCPM", "salRisDel_F", 
                 "salRisDel_PValue", "salRisDel_FDR")) %>% 
  left_join(
    res_all$results$cluster_0$posthoc$`SECOND_risperidone-saline|FIRST_wtwt` %>% 
      set_colnames(c("name_id", "peak_id", "gene_symbol", 
                     "salRisWt_logFC", "salRisWt_logCPM", "salRisWt_F", 
                     "salRisWt_PValue", "salRisWt_FDR")),
    by = c("name_id", "peak_id", "gene_symbol")
  ) %>% 
  # usuń punkty ze środka
  filter(!(salRisDel_logFC > -0.2 & salRisDel_logFC < 0.2 &
             salRisWt_logFC > -0.2 & salRisWt_logFC < 0.2)) %>% 
  ggplot(aes(x = salRisDel_logFC, y = salRisWt_logFC)) +
  
  # 1. wszystkie punkty szare
  geom_point(color = "grey70", size = 1) +
  
  # 2. czarne punkty (tylko df$name_id, ale bez bordowych)
  geom_point(
    data = . %>% filter(name_id %in% df$name_id &
                          !name_id %in% c("peak-247180-Slc9a8",
                                          "peak-127298-Lats2",
                                          "peak-171286-Smim11")),
    color = "black", size = 1.5
  ) +
  
  # 3. bordowe punkty (tylko trzy wybrane)
  geom_point(
    data = . %>% filter(name_id %in% c("peak-247180-Slc9a8",
                                       "peak-127298-Lats2",
                                       "peak-171286-Smim11")),
    color = "darkred", size = 2
  ) +
  
  # 4. etykiety czarne
  geom_text_repel(
    data = . %>% filter(name_id %in% df$name_id &
                          !name_id %in% c("peak-247180-Slc9a8",
                                          "peak-127298-Lats2",
                                          "peak-171286-Smim11")),
    aes(label = gene_symbol),
    size = 3, color = "black"
  ) +
  
  # 5. etykiety bordowe
  geom_text_repel(
    data = . %>% filter(name_id %in% c("peak-247180-Slc9a8",
                                       "peak-127298-Lats2",
                                       "peak-171286-Smim11")),
    aes(label = gene_symbol),
    size = 3, color = "darkred"
  ) +
  
  # linie pomocnicze
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", color = "grey50") +
  
  # tytuł i osie
  labs(
    title = "cluster_0",
    x = "3q29 deletion: Risperidone vs Saline (log2FC)",
    y = "Wt/Wt: Risperidone vs Saline (log2FC)"
  ) +
  
  # wymuszenie kwadratu
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  theme_classic() -> scatterplot_cluster0

res_all$results$cluster_2$posthoc$`SECOND_risperidone-saline|FIRST_wtdel` %>% 
  set_colnames(c("name_id", "peak_id", "gene_symbol", 
                 "salRisDel_logFC", "salRisDel_logCPM", "salRisDel_F", 
                 "salRisDel_PValue", "salRisDel_FDR")) %>% 
  left_join(
    res_all$results$cluster_2$posthoc$`SECOND_risperidone-saline|FIRST_wtwt` %>% 
      set_colnames(c("name_id", "peak_id", "gene_symbol", 
                     "salRisWt_logFC", "salRisWt_logCPM", "salRisWt_F", 
                     "salRisWt_PValue", "salRisWt_FDR")),
    by = c("name_id", "peak_id", "gene_symbol")
  ) %>% 
  filter(!(salRisDel_logFC > -0.2 & salRisDel_logFC < 0.2 &
             salRisWt_logFC > -0.2 & salRisWt_logFC < 0.2)) %>% 
  ggplot(aes(x = salRisDel_logFC, y = salRisWt_logFC)) +
  
  # 1. wszystkie punkty szare
  geom_point(color = "grey70", size = 1) +
  
  # 2. czarne punkty (df$name_id, ale bez bordowego)
  geom_point(
    data = . %>% filter(name_id %in% df$name_id &
                          name_id != "peak-25960-Dcaf6"),
    color = "black", size = 1.5
  ) +
  
  # 3. bordowy punkt (tylko Dcaf6)
  geom_point(
    data = . %>% filter(name_id == "peak-25960-Dcaf6"),
    color = "darkred", size = 2
  ) +
  
  # etykiety czarne
  geom_text_repel(
    data = . %>% filter(name_id %in% df$name_id &
                          name_id != "peak-25960-Dcaf6"),
    aes(label = gene_symbol), size = 3, color = "black"
  ) +
  
  # etykieta bordowa
  geom_text_repel(
    data = . %>% filter(name_id == "peak-25960-Dcaf6"),
    aes(label = gene_symbol), size = 3, color = "darkred"
  ) +
  
  # linie pomocnicze
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", color = "grey50") +
  
  # tytuł i osie
  labs(
    title = "cluster_2",
    x = "3q29 deletion: Risperidone vs Saline (log2FC)",
    y = "Wt/Wt: Risperidone vs Saline (log2FC)"
  ) +
  
  # kwadratowy układ współrzędnych
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  theme_classic()  -> scatterplot_cluster2



scatterplot_cluster0 + scatterplot_cluster2


df <- read_excel("data/risperidone-3q29/tmp/3q29-ris-spatial-edger-gene-list-pInteraction0.01-all.xlsx")


df %>% filter(cluster == "cluster_0")


df %>% 
  select(cluster, name_id)
