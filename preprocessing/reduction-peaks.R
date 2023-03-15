#############################################################
## Stuff I have to do
## mateuszzieba97@gmail.com - Mar 2022
## Script written in 4.1.2 r version
## Responsible for the reduction of peaks
#############################################################

######################################
# install and require needed package #
######################################

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicAlignments")
BiocManager::install("rtracklayer")
BiocManager::install("Gviz")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("bamsignals")
install.packages("tidyr")
install.packages("foreach")
install.packages("doParallel")

require(biomaRt)
require(GenomicAlignments)
require(rtracklayer)
require(Gviz)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(BSgenome.Mmusculus.UCSC.mm10)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(magrittr)
require(dplyr)
require(stringr)
require(tidyr)
require(purrr)
require(foreach)
require(doParallel)
require(bamsignals)



####################
# prepare function #
####################

# function to visualize interest peak
plot_peak_bam <- function(peaks_data, peak_id, bam_file){
  peak_data <- peaks_data[peaks_data$peak_id == peak_id, ]
  chromosome = peak_data[, 1]
  start = peak_data[, 2]
  end = peak_data[, 3]
  
  side_title <- paste(peak_data[,4],
                      "\n",
                      peak_data[, 5],
                      peak_data[, 7],
                      peak_data[, 8],
                      peak_data[, 9],
                      peak_data[, 10],
                      peak_data[, 11],
                      peak_data[, 16],
                      end - start,
                      peak_data[, 18],
                      sep = " ")
  
  ideogram_plot <- IdeogramTrack(genome = "mm10", chromosome = chromosome)
  
  gtrack <- GenomeAxisTrack()
  
  coverage_plot <- DataTrack(range = bam_file, genome = "mm10",
                             type = "horizon", name = "horizon coverage", window = -2, chromosome = chromosome,
                             fill.mountain=c("lighgrey", "lightgrey"), col.mountain="grey")
  
  polygon_plot <- DataTrack(range = bam_file, genome = "mm10",
                            type = "histogram", window = -1, chromosome = chromosome,
                            name = "historgram coverage",
                            fill.mountain=c("lighgrey", "lightgrey"), col.mountain="grey")
  
  heatmap_plot <- DataTrack(range = bam_file, genome = "mm10",
                            type = "heatmap", window = -1, chromosome = chromosome,
                            fill.mountain=c("lighgrey", "lightgrey"))
  
  
  reads_plot <- AnnotationTrack(range = bam_file, genome = "mm10",
                                name = "Reads", chromosome = chromosome,
                                stacking = "dense", alpha = 0.5,
                                max.height = 0.5)
  
  
  plotTracks(list(ideogram_plot, gtrack, coverage_plot, polygon_plot, heatmap_plot, reads_plot),
             from = start,
             to = end, main = side_title, col.main="black")
  
}

# function to visualize all peaks for interest gene
plot_gene_bam <-function(data, gene, bam_file){
  data %>% 
    filter(gene_name == gene) %>%
    .[, 4] -> peak_vector
  
  print(peak_vector)
  
  list_plot <- list()
  
  for(peak_id in peak_vector){
    plot_peak_bam(peaks_data = data, peak_id = peak_id, bam_file = bam_file)
  }
}


# function to visualize range for all peaks describe to interest gene
plot_gene_long_bam <- function(data = info_peaks, gene, bam_file){
  # function to visualize all peaks describe to interest gene
  # create GRanges
  data[data$gene_name == gene, ] %>%
    .[, c(1, 2, 3, 6)] %>% 
    rename(chr = chr_peak,
           start = start_peak,
           end = end_peak,
           strand = strand_coverage) %>% makeGRangesFromDataFrame() -> gene_granges
  
  peaks_name <- data[data$gene_name == gene, ] %>% 
    mutate(peak_id = str_replace(peak_id, "merged-samples_peak_", "")) %>% 
    .[, 4]
  
  peak_data <- data[data$gene_name == gene, ]
  
  chromosome = peak_data[, 1] %>% unique()
  
  start = peak_data[, 2] %>% min()
  end = peak_data[, 3] %>% max()
  
  ideogram_plot <- IdeogramTrack(genome = "mm10", chromosome = chromosome)
  gene_atrack <- AnnotationTrack(gene_granges, 
                                 id = peaks_name,
                                 name = "peaks")
  gtrack <- GenomeAxisTrack()
  coverage_plot <- DataTrack(range = bam_file, genome = "mm10",
                             type = "horizon", name = "coverage", window = -2, chromosome = chromosome,
                             fill.mountain=c("lighgrey", "lightgrey"), col.mountain="grey",
                             cex.title = 1.6,  background.title = "brown")
  polygon_plot <- DataTrack(range = bam_file, genome = "mm10",
                            type = "polygon", window = -1, chromosome = chromosome,
                            fill.mountain=c("lighgrey", "lightgrey"), col.mountain="grey")
  
  heatmap_plot <- DataTrack(range = bam_file, genome = "mm10",
                            type = "heatmap", window = -1, chromosome = chromosome,
                            fill.mountain=c("lighgrey", "lightgrey"))
  
  plotTracks(list(ideogram_plot, gtrack, polygon_plot, heatmap_plot, gene_atrack),
             shape = "box", featureAnnotation = "id", 
             rotation = 90, rotation.item=90,
             cex.feature = 0.7,
             fontcolor.feature = "black",
             from = start,
             to = end, main = gene)
  
}


# function to calculate max value from region or reads for region
parameter_region_multicore <- function(data, number_core, bam_file, parameter){

  # calculate using multi core
  name_chunk <- paste(rep("part", number_core), seq(1:number_core), sep = "_")
  name_chunk
  
  number_row <- nrow(data)
  number_row_chunk <- ceiling(nrow(data) / number_core)
  
  # split df
  split_df <- split(data, 
                    f = rep(seq_len(ceiling(number_row / number_row_chunk)), 
                            each = number_row_chunk, 
                            length.out = number_row))
  
  names(split_df) <- name_chunk
  
  # extract information from bam file
  registerDoParallel(number_core)
  foreach(chunk=name_chunk)%dopar% {
    
    split_df[[chunk]] %>% 
      .[, c(1, 2, 3, 6)] %>% 
      rename(chr = chr_peak,
             start = start_peak,
             end = end_peak,
             strand = strand_coverage) %>% 
      makeGRangesFromDataFrame() %>% 
      {if (parameter == "amplitude") {
        bamCoverage(bam_file, .) %>%
          lapply(., max) %>% 
          unlist()
      } else if(parameter == "counts"){
        bamCount(bam_file, .) 
      }} 
  } -> peak_data
}


#############
# read data #
#############
bam_file <- "data/gene-annotation/merged-samples.bam"
file_info_peaks <- "data/gene-annotation/peaks-annotate-sort.bed"
file_info_peaks <- "data/ldopa/gene-annotation/peaks-annotate-filt-nanopore.bed" # for ldopa
file_info_peaks <- "data/tmp-antypsychotics//gene-annotation/peaks-annotate-filt-nanopore.bed" # for tmp-antypsychotics
bam_file <- "data/tmp-antypsychotics/gene-annotation/merged-samples.bam" # for tmp-antypsychotics



info_peaks <- read.table(file_info_peaks, 
                       header = FALSE,
                       sep = "\t",
                       col.names = c('chr_peak',
                                     'start_peak',
                                     'end_peak',
                                     'peak_id',
                                     'score_int(-10*log10pvalue)',
                                     'strand_coverage',
                                     'fold_change_peak_summit',
                                     '-log10pvalue_peak_summit',
                                     '-log10qvalue_peak_summit',
                                     'relative_summit_position_peak_start',
                                     'type_peak',
                                     'chr_gene',
                                     'start_gene',
                                     'end_gene',
                                     'gene_id',
                                     'gene_name',
                                     'strand_gene')) 


#############################
# reduction number of peaks #
#############################
info_peaks %>%
  # calculate distance peaks from gene
  # if peak is inside gene distance = 0
  mutate(distance = ifelse(strand_gene == "+",
                           {
                             ifelse(start_peak < end_gene & start_peak > start_gene,
                                    0,
                                    ifelse(
                                      start_peak > end_gene,
                                      (start_peak - end_gene),
                                      (end_peak - start_gene)
                                    ))
                           },
                           {
                             ifelse(end_peak < end_gene & end_peak > start_gene,
                                    0,
                                    ifelse(
                                      end_peak < start_gene,
                                      (start_gene - end_peak),
                                      (end_gene - start_peak)
                                    ))
                           })) %>%
  # remove peaks with negative distance to gene
  filter(distance >= 0) %>%
  group_by(peak_id) %>%
  nest() %>%
  # calculate min value of distance for peak to describe gene
  mutate(min_distance = map(data, ~ min(.$distance))) %>%
  unnest(min_distance, data) %>%
  as.data.frame() %>%
  # filter peaks which are nearest gene
  filter(distance == min_distance) %>%
  # calculate abs distance for peaks inner gene
  # if peak is not inner, set value 0
  mutate(inner_closest = ifelse(strand_gene == "+",
                                {
                                  ifelse(distance == 0,
                                         abs(start_peak - end_gene), 0)
                                },
                                {
                                  ifelse(distance == 0,
                                         abs(start_gene - end_peak), 0)
                                })) %>%
  group_by(peak_id) %>%
  nest() %>%
  # calculate min for inner_closest
  mutate(min_inner_closest = map(data, ~ min(.$inner_closest))) %>%
  unnest(data, min_inner_closest) %>%
  as.data.frame() %>%
  # filter gene which are closest inner peak
  filter(inner_closest == min_inner_closest) %>%
  # reorder columns
  .[, c(2:4, 1, 5:18)] -> info_peaks_closest



# calculate amplitude for peaks
parameter_region_multicore(data = info_peaks_closest,
                        number_core = 16,
                        bam_file = bam_file,
                        parameter = "amplitude") -> amplitude_peaks

# calculate number reads for peak
parameter_region_multicore(data = info_peaks_closest,
                        number_core = 16,
                        bam_file = bam_file,
                        parameter = "counts") -> peak_counts

# add column with amplitude and reads to dataframe
info_peaks_closest$peak_counts <- unlist(peak_counts)

info_peaks_closest$amplitude_peak <- unlist(amplitude_peaks)

# filtering peaks
info_peaks_closest %>% 
mutate(
  ratio_counts_amplitude = .$peak_counts/.$amplitude_peak
) %>%
filter(
  score_int..10.log10pvalue. > 350,
  # peak_counts > 1200,
  peak_counts > 800,
  X.log10pvalue_peak_summit > 35,
  amplitude_peak > 600,
  ratio_counts_amplitude > 1.4
) -> info_peaks_reduction

###################################
# collect row to one for + strand #
###################################

# find duplicate peaks
info_peaks_reduction %>% 
group_by(peak_id) %>% 
count() %>% filter(n >1) %>% as.data.frame() %>%
.[, 1] -> duplicate_peak

# vector of column to group_by for plus strand
column_names <- info_peaks_reduction %>% 
colnames() %>%
.[!. %in% c("start_gene", "end_gene", "gene_id", "gene_name")] 

# collect many row to one for plus strand, 
# for minus strand not detect duplicate peak_id
info_peaks_reduction %>% unique() %>% 
  group_by(across(all_of(column_names))) %>%
  summarise(start_gene = str_c(start_gene, collapse = ", "),
            end_gene = str_c(end_gene, collapse = ", "),
            gene_id = str_c(gene_id, collapse = ", "),
            gene_name = str_c(gene_name, collapse = ", ")) %>% 
  as.data.frame() %>% 
  .[, c(1:12, 18, 19, 20, 21, 13, 14:17)] %>% 
  # save to file,
  write.table(., 
              file='data/tmp-antypsychotics//gene-annotation/peaks-annotate-reduction.tsv', 
              quote=FALSE, 
              sep='\t', 
              col.names = TRUE, 
              row.names = FALSE)
  # write.table(., 
  #             file='data/ldopa/gene-annotation/peaks-annotate-reduction.tsv', 
  #             quote=FALSE, 
  #             sep='\t', 
  #             col.names = TRUE, 
  #             row.names = FALSE)


###################################################
# Visualization interest peaks - check filtration #
###################################################

plot_peak_bam(peaks_data = info_peaks_reduction,
              peak_id = "merged-samples_peak_172",
              bam_file = bam_file)


plot_gene_bam(data = info_peaks_reduction, 
              gene = "Junb", 
              bam_file = bam_file)


plot_gene_long_bam(data = info_peaks_reduction,
                   gene = "Ttr", 
                   bam_file = bam_file)

plot_gene_long_bam(data = info_peak_reduction, 
                   gene="Junb", 
                   bam_file = bam_file)



