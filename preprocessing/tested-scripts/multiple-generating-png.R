c("clozapine-peak-172728",
  "clozapine-peak-214362",
  "clozapine-peak-94947",
  "clozapine-peak-101104",
  "clozapine-peak-103990",
  "clozapine-peak-172728",
  "clozapine-peak-197060",
  "clozapine-peak-71605",
  "clozapine-peak-70652"
) -> tmp_vpeak

c("Tcf4",
  "Sdc4",
  "Mgll",
  "Slc38a3",
  "Lars2",
  "Lars2",
  "Paqr8",
  "Eef1b2",
  "Tmcc2",
  "Gpr37l1",
  "Csrp1",
  "Glrx2",
  "Pbx1",
  "Fabp7",
  "Ramp3",
  "Rps27a",
  "Fnip1",
  "Zbtb4",
  "Aldoc",
  "Tmem121",
  "Shc3",
  "Fam107a",
  "Chd8",
  "Gjb6",
  "Rrm2b",
  "Fam91a1",
  "Rpl8",
  "Cacng2",
  "Rbfox1",
  "Naa50",
  "Lix1",
  "Dlgap1",
  "Celf4",
  "Camk2a",
  "Tcf4",
  "Rnaseh2c",
  "Zfp91",
  "Psat1",
  "Gnaq",
  "Mxi1",
  "Epb41l1",
  "Rps3a1",
  "Naxe",
  "Pip5k1a",
  "Gstm5",
  "S1pr1",
  "Ppp3ca",
  "Rps20",
  "Tesk1",
  "Nfib",
  "Rps6",
  "Rps8",
  "Ptp4a2",
  "Trnp1",
  "Camk2n1",
  "Fgfr3",
  "Ldb2",
  "Ube3b",
  "Pdgfa",
  "Ptprz1",
  "Mgll",
  "Rps19",
  "Igf1r",
  "Rab8a",
  "Zfpm1",
  "Ntm",
  "Hepacam",
  "Rplp1",
  "Lars2",
  "Lars2",
  "Igfbp2",
  "Tmcc2",
  "Prdx6",
  "S100b",
  "Nefh",
  "Wwc1",
  "Zbtb4",
  "Timm22",
  "Aldoc",
  "Aldoc",
  "Asic2",
  "Rasl10b",
  "Slc9a3r1",
  "Shc3",
  "Gjb6",
  "Cldn10",
  "Selenop",
  "Eif3h",
  "Chrac1",
  "Nrbp2",
  "Mroh1",
  "Thap7",
  "Hcfc1r1",
  "Metrn",
  "Srsf3",
  "Dpp9",
  "Cox7a2l",
  "Asrgl1",
  "Fads1",
  "Gnaq",
  "Mxi1",
  "Lhx2",
  "Epc2",
  "Arhgap1",
  "Zfp106",
  "Trp53inp2",
  "Sdc4",
  "Spata2",
  "Bcas1",
  "Mtx1",
  "Gstm5",
  "Ube2r2",
  "Tmod1",
  "Rps6",
  "Camk2n1",
  "Camkk2",
  "Gna12",
  "Sfxn5",
  "Arhgap35",
  "Zfp260",
  "Luzp2",
  "Apba2",
  "Hs3st2",
  "Fam57b",
  "Htra1",
  "Gnao1",
  "6430548M08Rik",
  "Endod1",
  "Ecsit",
  "Pccb",
  "Lars2",
  "Sacm1l",
  "Dynlt3",
  "Pin4",
  "Ndrg1",
  "Tcf4",
  "Ssbp3",
  "Camk2n1",
  "Arpc1b",
  "Parva",
  "Lars2",
  "Arhgap5",
  "Otub2",
  "Fubp1",
  "Maged2"
) -> clozapine_genes_vector



# Create a directory to store the png files
dir.create("plots")

clozapine_st_data_half$raw_data$annotate %>% 
  # select(c(gene_name, peak_id)) %>% 
  filter(gene_name %in% clozapine_genes_vector) %>% 
  select(c(peak_id, gene_name)) -> peaks_df
  
pb <- txtProgressBar(min = 0, max = nrow(peaks_df[1:10,]), style = 3)

# Loop through each peak_id in peak and generate plots
for(i in 1:nrow(peaks_df[1:10,])) {
  peak <- peaks_df[i, 1]
  gene <- peaks_df[i, 2]
  
  png_filename <- paste0("plots/plot_", gene, "-", peak, ".png")
  
  # Start the png device
  png(filename = png_filename, width = 1100, height = 750)
  
  plot <- spatial_feature_plot(spatial_data = clozapine_st_data_half,
                               type_data = "raw_data",
                               peak_id = peak,
                               samples = c(samples_saline, samples_clozapine),
                               normalization = T,
                               tif_image = T) +
    plot_layout(ncol = 4) + 
    plot_annotation(title = paste(gene, peak, sep = ": "))
  
  print(plot)  # Explicitly print the plot
  dev.off()  # Close the png device
  
  # Update progress bar
  setTxtProgressBar(pb, i)
}

# Close the progress bar
close(pb)



start_time <- Sys.time()
zip("plots.zip", files = "plots")
end_time <- Sys.time()

end_time - start_time



# Delete the directory with all the PNG files
unlink("plots", recursive = TRUE)


# Create a zip file
start_time <- Sys.time()
tar("plots.tar", files = "plots")
end_time <- Sys.time()

end_time - start_time





spatial_feature_plot(spatial_data = clozapine_st_data_half,
                     type_data = "raw_data",
                     peak_id = peak,
                     samples = c(samples_saline, samples_clozapine),
                     # min_percentile = input$min_percentile,
                     # max_percentile = input$max_percentile,
                     normalization = T,
                     # size = input$spot_size,
                     tif_image = T) +
  plot_layout(ncol = 4)

for(i in 1:nrow(peaks_df)) {
  peak <- peaks_df[i, 1]
  # gene <- peaks_df[i, 2]
  # print(gene, peak)
  print(peak)
}


# Create a directory to store the png files
dir.create("plots")

# Get a subset of the data that contains the genes of interest.
# The selected dataframe is then stored in 'peaks_df'.
clozapine_st_data_half$raw_data$annotate %>% 
  filter(gene_name %in% clozapine_genes_vector) %>% 
  select(c(peak_id, gene_name)) -> peaks_df

# Setting up a text progress bar to track the progress of the loop.
pb <- txtProgressBar(min = 0, max = nrow(peaks_df[1:10,]), style = 3)

# Iterate over each row in the 'peaks_df' dataframe.
# Generate a plot for each peak_id in peaks_df and save it as a png file.
for(i in seq_len(nrow(peaks_df[1:10,]))) {
  # Extract the peak and gene values from the current row
  peak <- peaks_df[i, "peak_id"]
  gene <- peaks_df[i, "gene_name"]
  
  # Generate the filename for the png file
  png_filename <- file.path("plots", paste0("plot_", gene, "-", peak, ".png"))
  
  # Initialize the png device with the generated filename
  png(filename = png_filename, width = 1100, height = 750)
  
  # Generate the spatial feature plot
  plot <- spatial_feature_plot(
    spatial_data = clozapine_st_data_half,
    type_data = "raw_data",
    peak_id = peak,
    samples = c(samples_saline, samples_clozapine),
    normalization = TRUE,
    tif_image = TRUE) +
    plot_layout(ncol = 4) + 
    plot_annotation(title = paste0(gene, ": ", peak))
  
  # Print the plot to the png file
  print(plot)
  
  # Close the png device
  dev.off()
  
  # Update the progress bar
  setTxtProgressBar(pb, i)
}

# Close the progress bar after the loop finishes
close(pb)

# Record the time before zipping files
start_time <- Sys.time()

# Zip all the generated png files into 'plots.zip'
zip("plots.zip", files = "plots")

# Record the time after zipping files
end_time <- Sys.time()

# Print out the time taken to zip the files
print(end_time - start_time)

# Delete the 'plots' directory and all its contents
unlink("plots", recursive = TRUE)

