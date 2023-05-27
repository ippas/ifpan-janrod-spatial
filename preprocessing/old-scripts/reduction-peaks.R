######################################
# install and load needed packages #
######################################
# install.packages(package, lib = Sys.getenv("R_LIBS_USER"))
# List of CRAN packages

cran_packages <-
  c("tidyr",
    "optparse",
    "foreach",
    "doParallel",
    "magrittr",
    "dplyr",
    "stringr",
    "purrr")

# List of Bioconductor packages
bioconductor_packages <- c(
  "S4Vectors",
  "GenomicAlignments",
  "rtracklayer",
  "Gviz",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "BSgenome.Mmusculus.UCSC.mm10",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "bamsignals"
)

# Install and load CRAN packages
for (package in cran_packages) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install and load Bioconductor packages
for (package in bioconductor_packages) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(package)
    library(package, character.only = TRUE)
  }
}


####################
# prepare function #
####################
# Function to visualize the region of interest in a BAM file
plot_peak_bam <- function(peaks_data, peak_id, bam_file) {
  # Extract data for the selected peak
  peak_data <- peaks_data[peaks_data$peak_id == peak_id,]
  
  # Define variables from peak_data
  chromosome <- peak_data[, 1]
  start <- peak_data[, 2]
  end <- peak_data[, 3]
  
  # Create a side title with relevant information
  side_title <- paste(
    peak_data[, 4],
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
    sep = " "
  )
  
  # Create ideogram track (visual representation of the chromosome)
  ideogram_plot <-
    IdeogramTrack(genome = "mm10", chromosome = chromosome)
  
  # Create genome axis track (coordinate axis for the genome)
  gtrack <- GenomeAxisTrack()
  
  # Create horizon coverage track (visual representation of read coverage)
  coverage_plot <- DataTrack(
    range = bam_file,
    genome = "mm10",
    type = "horizon",
    name = "horizon coverage",
    window = -2,
    chromosome = chromosome,
    fill.mountain = c("lighgrey", "lightgrey"),
    col.mountain = "grey"
  )
  
  # Create histogram coverage track (bar plot of read coverage)
  polygon_plot <- DataTrack(
    range = bam_file,
    genome = "mm10",
    type = "histogram",
    window = -1,
    chromosome = chromosome,
    name = "historgram coverage",
    fill.mountain = c("lighgrey", "lightgrey"),
    col.mountain = "grey"
  )
  
  # Create heatmap coverage track (visual representation of read coverage as a heatmap)
  heatmap_plot <- DataTrack(
    range = bam_file,
    genome = "mm10",
    type = "heatmap",
    window = -1,
    chromosome = chromosome,
    fill.mountain = c("lighgrey", "lightgrey")
  )
  
  # Create annotation track (visual representation of individual reads)
  reads_plot <- AnnotationTrack(
    range = bam_file,
    genome = "mm10",
    name = "Reads",
    chromosome = chromosome,
    stacking = "dense",
    alpha = 0.5,
    max.height = 0.5
  )
  
  # Plot tracks together with the main title and side title
  plotTracks(
    list(
      ideogram_plot,
      gtrack,
      coverage_plot,
      polygon_plot,
      heatmap_plot,
      reads_plot
    ),
    from = start,
    to = end,
    main = side_title,
    col.main = "black"
  )
}


# Function to visualize all peaks for an interested gene in a BAM file
plot_gene_bam <- function(data, gene, bam_file_path) {
  # Filter data for the selected gene and extract peak IDs
  peak_vector <- data %>%
    filter(gene_name == gene) %>%
    pull(4)
  
  # Print the peak IDs
  print(peak_vector)
  
  # Initialize an empty list to store plots
  list_plot <- list()
  
  # Loop through peak IDs and call plot_peak_bam function for each peak ID
  for (peak_id in peak_vector) {
    plot_peak_bam(peaks_data = data,
                  peak_id = peak_id,
                  bam_file = bam_file_path)
  }
}


# Function to visualize the range of all peaks associated with an interested gene in a BAM file
plot_gene_long_bam <-
  function(data = info_peaks, gene, bam_file_path) {
    # Create GRanges object for the selected gene
    gene_granges <- data %>%
      filter(gene_name == gene) %>%
      dplyr::select(chr_peak, start_peak, end_peak, strand_coverage) %>%
      dplyr::rename(
        chr = chr_peak,
        start = start_peak,
        end = end_peak,
        strand = strand_coverage
      ) %>%
      makeGRangesFromDataFrame()
    
    # Extract peak names
    peaks_name <- data %>%
      filter(gene_name == gene) %>%
      mutate(peak_id = str_replace(peak_id, "merged-samples_peak_", "")) %>%
      pull(4)
    
    # Define variables for plotting
    peak_data <- data %>% filter(gene_name == gene)
    chromosome <- unique(peak_data[, 1])
    start <- min(peak_data[, 2])
    end <- max(peak_data[, 3])
    
    # Create tracks for plotting
    ideogram_plot <-
      IdeogramTrack(genome = "mm10", chromosome = chromosome)
    gene_atrack <-
      AnnotationTrack(gene_granges, id = peaks_name, name = "peaks")
    gtrack <- GenomeAxisTrack()
    coverage_plot <- DataTrack(
      range = bam_file_path,
      genome = "mm10",
      type = "horizon",
      name = "coverage",
      window = -2,
      chromosome = chromosome,
      fill.mountain = c("lighgrey", "lightgrey"),
      col.mountain = "grey",
      cex.title = 1.6,
      background.title = "brown"
    )
    polygon_plot <- DataTrack(
      range = bam_file_path,
      genome = "mm10",
      type = "polygon",
      window = -1,
      chromosome = chromosome,
      fill.mountain = c("lighgrey", "lightgrey"),
      col.mountain = "grey"
    )
    heatmap_plot <- DataTrack(
      range = bam_file_path,
      genome = "mm10",
      type = "heatmap",
      window = -1,
      chromosome = chromosome,
      fill.mountain = c("lighgrey", "lightgrey")
    )
    
    # Plot tracks together with the main title
    plotTracks(
      list(
        ideogram_plot,
        gtrack,
        polygon_plot,
        heatmap_plot,
        gene_atrack
      ),
      shape = "box",
      featureAnnotation = "id",
      rotation = 90,
      rotation.item = 90,
      cex.feature = 0.7,
      fontcolor.feature = "black",
      from = start,
      to = end,
      main = gene
    )
  }

# Function to calculate max value from region or reads for region using multiple cores
parameter_region_multicore <-
  function(data,
           number_core,
           bam_file_path,
           parameter) {
    # Create chunk names for parallel processing
    name_chunk <-
      paste(rep("part", number_core), seq(1:number_core), sep = "_")
    
    # Calculate the number of rows and chunks
    number_row <- nrow(data)
    number_row_chunk <- ceiling(number_row / number_core)
    
    # Split the data into chunks
    split_df <- split(data,
                      f = rep(seq_len(ceiling(
                        number_row / number_row_chunk
                      )),
                      each = number_row_chunk,
                      length.out = number_row))
    names(split_df) <- name_chunk
    
    # Register parallel backend
    registerDoParallel(number_core)
    
    # Extract information from the BAM file using multiple cores
    peak_data <- foreach(chunk = name_chunk) %dopar% {
      # Create GRanges object for the current chunk
      granges <- split_df[[chunk]] %>%
        dplyr::select(chr_peak, start_peak, end_peak, strand_coverage) %>%
        dplyr::rename(
          chr = chr_peak,
          start = start_peak,
          end = end_peak,
          strand = strand_coverage
        ) %>%
        makeGRangesFromDataFrame()
      
      # Calculate parameter value based on the input parameter
      if (parameter == "amplitude") {
        bamCoverage(bam_file_path, granges) %>%
          lapply(max) %>%
          unlist()
      } else if (parameter == "counts") {
        bamCount(bam_file_path, granges)
      }
    }
    
    return(peak_data)
  }

#####################
# prepare variables #
#####################

# Define command-line options
option_list <- list(
  make_option(c("-f", "--file_info_peaks"), type = "character", default = NULL,
              help = "File containes information about peaks", metavar = "PEAK_FILE"),
  make_option(c("-b", "--bam_file"), type = "character", default = NULL,
              help = "Path to the BAM file", metavar = "BAM_FILE"),
  make_option(c("-o", "--output_tsv"), type = "character", default = "output.tsv",
              help = "Path to the output TSV file [default: %default]", metavar = "OUTPUT_FILE"),
  make_option(c("-c", "--number_cores"), type = "integer", default = 1,
              help = "Number of cores to use for analysis [default: %default]", metavar = "NUM_CORES"),
  make_option(c("-p", "--peak_counts"), type = "integer", default = 1200,
              help = "Minimum count of reads in BAM file for each peak [default: %default]", metavar = "PEAK_COUNTS"),
  make_option(c("-a", "--amplitude_peak"), type = "numeric", default = 600,
              help = "Value of amplitude for each peak in BAM file [default: %default]", metavar = "AMPLITUDE_PEAK"),
  make_option(c("-r", "--ratio_counts_amplitude"), type = "numeric", default = 1.4,
              help = "Ratio to calculate counts/amplitude for filtering peaks [default: %default]", metavar = "RATIO_COUNTS_AMPLITUDE")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Check if file_info_peaks argument is provided
if (is.null(args$file_info_peaks)) {
  cat("Please provide the path to the peak information file using the -f option.\n")
  q(status = 1)
}

# Access and use the file_info_peaks argument
file_info_peaks <- args$file_info_peaks

# Check if bam_file argument is provided
if (is.null(args$bam_file)) {
  cat("Please provide the path to the BAM file using the -b option.\n")
  q(status = 1)
}

# Access and use the bam_file argument
bam_file <- args$bam_file


# Access the arguments
output_tsv_path <- args$output_tsv
number_cores <- args$number_cores
peak_counts_filt <- args$peak_counts
amplitude_peak_filt <- args$amplitude_peak
ratio_counts_amplitude_filt <- args$ratio_counts_amplitude

# # example of variables
# bam_file <- "/ifpan-janrod-spatial/data/risperidone/gene-annotation/merged-samples.bam"
# file_info_peaks <-"/ifpan-janrod-spatial/data/risperidone/gene-annotation/peaks-annotate-filt-nanopore.bed"
# output_tsv_path <- '/ifpan-janrod-spatial/data/risperidone/gene-annotation/peaks-annotate-reduction3.tsv'
# number_cores <- 40
# peak_counts_filt <- 1200
# amplitude_peak_filt <- 600
# ratio_counts_amplitude_filt <- 1.4

############
# Analysis #
############

# Read the file containing information about peaks
print("Read the file containing information about peaks")

info_peaks <- read.table(
  file_info_peaks,
  header = FALSE,
  sep = "\t",
  col.names = c(
    'chr_peak',
    # Chromosome of the peak
    'start_peak',
    # Start position of the peak
    'end_peak',
    # End position of the peak
    'peak_id',
    # ID of the peak
    'score_int(-10*log10pvalue)',
    # Score based on -10 * log10(p-value)
    'strand_coverage',
    # Strand coverage
    'fold_change_peak_summit',
    # Fold change at the peak summit
    '-log10pvalue_peak_summit',
    # -log10(p-value) at the peak summit
    '-log10qvalue_peak_summit',
    # -log10(q-value) at the peak summit
    'relative_summit_position_peak_start',
    # Relative summit position from the peak start
    'type_peak',
    # Type of the peak
    'chr_gene',
    # Chromosome of the gene
    'start_gene',
    # Start position of the gene
    'end_gene',
    # End position of the gene
    'gene_id',
    # ID of the gene
    'gene_name',
    # Name of the gene
    'strand_gene'                      # Strand of the gene
  )
)

#############################
# reduction number of peaks #
#############################
print("Preselection peaks")
info_peaks %>%
  # Calculate the distance of peaks from genes
  # If the peak is inside the gene, set the distance to 0
  mutate(distance = ifelse(
    strand_gene == "+",
    {
      ifelse(
        start_peak < end_gene & start_peak > start_gene,
        0,
        ifelse(
          start_peak > end_gene,
          (start_peak - end_gene),
          (end_peak - start_gene)
        )
      )
    },
    {
      ifelse(
        end_peak < end_gene & end_peak > start_gene,
        0,
        ifelse(
          end_peak < start_gene,
          (start_gene - end_peak),
          (end_gene - start_peak)
        )
      )
    }
  )) %>%
  filter(distance >= 0) %>% # Remove peaks with negative distance
  group_by(peak_id) %>%
  nest() %>%
  # Calculate the minimum distance for the nearest peak to the gene
  mutate(min_distance = map(data, ~ min(.$distance))) %>%
  unnest(c(min_distance, data)) %>%
  as.data.frame() %>%
  filter(distance == min_distance) %>% # Filter for the nearest peak
  # Calculate the absolute distance for peaks inside the gene
  # If the peak is not inside, set the value to 0
  mutate(inner_closest = ifelse(
    strand_gene == "+",
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
  # Calculate the minimum value for inner_closest
  mutate(min_inner_closest = map(data, ~ min(.$inner_closest))) %>%
  unnest(c(data, min_inner_closest)) %>%
  as.data.frame() %>%
  # Filter for genes that are closest to the inner peak
  filter(inner_closest == min_inner_closest) %>% 
  # Reorder columns
  .[, c(2:4, 1, 5:18)]  -> info_peaks_closest


# Calculate amplitude for peaks
# This function calculates the amplitude of each peak in parallel using multiple cores
# The resulting values are stored in the variable "amplitude_peaks"
print("Extracting the amplitude for peaks")
parameter_region_multicore(
  data = info_peaks_closest,
  number_core = number_cores,
  bam_file = bam_file,
  parameter = "amplitude"
) -> amplitude_peaks

# Calculate number of reads for each peak
# This function calculates the number of reads for each peak in parallel using multiple cores
# The resulting values are stored in the variable "peak_counts"
print("Extracting the number of reads per peaks")
parameter_region_multicore(
  data = info_peaks_closest,
  number_core = number_cores,
  bam_file = bam_file,
  parameter = "counts"
) -> peak_counts

# Add columns with amplitude and read counts to the "info_peaks_closest" dataframe
# The "unlist" function is used to convert the list of values returned by the parallel functions into a vector
info_peaks_closest$peak_counts <- unlist(peak_counts)
info_peaks_closest$amplitude_peak <- unlist(amplitude_peaks)


# filter peaks based on certain criteria
info_peaks_closest %>%
  # calculate the ratio of peak counts to peak amplitude
  mutate(
    ratio_counts_amplitude = peak_counts/amplitude_peak
  ) %>%
  filter(
    score_int..10.log10pvalue. > 350, # filter based on score
    peak_counts > peak_counts_filt, # filter based on peak counts
    X.log10pvalue_peak_summit > 35, # filter based on summit height
    amplitude_peak > amplitude_peak_filt, # filter based on peak amplitude
    ratio_counts_amplitude > ratio_counts_amplitude_filt # filter based on peak shape
  ) -> info_peaks_reduction

###################################
# collect row to one for + strand #
###################################

# find duplicate peaks based on peak ID
info_peaks_reduction %>%
  group_by(peak_id) %>%
  dplyr::count() %>%
  filter(n > 1) %>%
  as.data.frame() %>%
  .[, 1] -> duplicate_peak

# vector of column names to group by for plus strand
column_names <- info_peaks_reduction %>%
  colnames() %>%
  # exclude non-relevant columns from grouping
  .[! . %in% c("start_gene", "end_gene", "gene_id", "gene_name")]

print("Selecting peaks")
# collect multiple rows into one for plus strand
info_peaks_reduction %>%
  unique() %>%
  group_by(across(all_of(column_names))) %>%
  summarise(
    start_gene = str_c(start_gene, collapse = ", "),
    end_gene = str_c(end_gene, collapse = ", "),
    gene_id = str_c(gene_id, collapse = ", "),
    gene_name = str_c(gene_name, collapse = ", ")
  ) %>%
  as.data.frame() %>% 
  # reorder columns to desired order
  .[, c(1:12, 18, 19, 20, 21, 13, 14:17)] %>%
  write.table(.,
              # file='/ifpan-janrod-spatial/data/risperidone/gene-annotation/peaks-annotate-reduction.tsv',
              file = output_tsv_path,
              quote=FALSE,
              sep='\t',
              col.names = TRUE,
              row.names = FALSE)



# ##################################################
# Visualization interest peaks - check filtration #
# ##################################################
# 
# plot_peak_bam(peaks_data = info_peaks_reduction,
#               peak_id = "risperidone_peak_173",
#               bam_file = bam_file)
# 
# 
# plot_gene_bam(data = info_peaks_reduction,
#               gene = "Fkbp5",
#               bam_file = bam_file)
# 
# 
# plot_gene_long_bam(data = info_peaks_reduction,
#                    gene = "Ttr",
#                    bam_file = bam_file)
# 
# plot_gene_long_bam(data = info_peak_reduction,
#                    gene="Junb",
#                    bam_file = bam_file)
