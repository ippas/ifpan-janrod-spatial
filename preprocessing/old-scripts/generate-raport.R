install.packages(c("rtf"), repos = "http://cran.r-project.org")
require(rtf)



output<-"/projects/spatial-raport.doc" # although this is RTF, we can use the
# .doc extension so it opens in MS Word
rtf<-RTF(output,width=8.5,height=11,font.size=10,omi=c(1,1,1,1))
# Other rtf commands here...
done(rtf) # writes and closes the file


for(cluster in c(0:24)){
  addHeader(rtf,title=paste("Cluster", cluster))
  spatial_interest_cluster(cluster = cluster,
                           spatial_data = spatial_transcriptomic_data,
                           resolution = 1,
                           samples = {meta_data %>% filter(mouse_genotype == "tif-mutant" & sample_ID != "S5295Nr2") %>% .[order(.$treatment),] %>%.[,1]},
                           size= 1, 
                           ncol = 3) -> tmp_plot
  
  addPlot(rtf,plot.fun=print,width=4,height=4,res=50, tmp_plot)
  
  addNewLine(rtf)
  significant_statistics(spatial_data = spatial_transcriptomic_data,
                         type_data = "range_normalize",
                         stat_test = "t.test",
                         resolution = 1,
                         cluster = cluster,
                         per = "sample",
                         pvalue_threshold = 0.05,
                         log2ratio_threshold = 0.5) -> tmp_tab
  
  addTable(rtf,tmp_tab,font.size=9,row.names=FALSE,NA.string="-")
  addNewLine(rtf)
  addText(rtf, paste("pemrmutate fdr for cluster", cluster, ": ", as.numeric(permutate_median_vector[cluster + 1])))
  addNewLine(rtf)
  
  
}

done(rtf)



output<-"/projects/spatial-raport-pvalue0.05-permutate1000.doc" # although this is RTF, we can use the
# .doc extension so it opens in MS Word
rtf<-RTF(output,width=8.5,height=11,font.size=10,omi=c(1,1,1,1))
# Other rtf commands here...
done(rtf) # writes and closes the file


for(cluster in c(0:24)){
  addHeader(rtf,title=paste("Cluster", cluster))
  spatial_interest_cluster(cluster = cluster,
                           spatial_data = spatial_transcriptomic_data,
                           resolution = 1,
                           samples = {meta_data %>% filter(mouse_genotype == "tif-mutant" & sample_ID != "S5295Nr2") %>% .[order(.$treatment),] %>%.[,1]},
                           size= 1, 
                           ncol = 3) -> tmp_plot
  
  addPlot(rtf,plot.fun=print,width=6,height=6,res=150, tmp_plot)
  
  addNewLine(rtf)
  significant_statistics(spatial_data = spatial_transcriptomic_data,
                         type_data = "range_normalize",
                         stat_test = "t.test",
                         resolution = 1,
                         cluster = cluster,
                         per = "sample",
                         pvalue_threshold = 0.05,
                         log2ratio_threshold = 0.5) -> tmp_tab
  
  addTable(rtf,tmp_tab,font.size=9,row.names=FALSE,NA.string="-")
  addNewLine(rtf)
  # addText(rtf, paste("pemrmutate fdr for cluster", cluster, ": ", 
  #                    as.numeric(spatial_transcriptomic_data$range_normalize$statistics$cluster_resolution_1$permuatate_100$pvalue_0.05[cluster + 1])))
  addText(rtf, paste("permutate fdr for cluster", cluster, ": ",
                     permutate_median_vector[cluster + 1]))
  addNewLine(rtf)
  addPageBreak(rtf)
  
  
}

done(rtf)



