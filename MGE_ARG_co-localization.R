# ARG - MGE Co-localization Analysis ####
library(readxl)
library(openxlsx)
library(dplyr)

# 1. Load data
df_arg <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = "ARG")
df_mge <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = "All_unique_ARG-carrying_MGE")

# 2. Define the linkage function with specific MGE column selection
link_MGE_ARG <- function(df_mge, df_arg) {
  
  # Step 1: Pre-select specific columns from the MGE table
  # This ensures we only carry relevant MGE metadata
  mge_subset <- df_mge %>%
    select(
      specific_contig, pident, bitscore, subject_sequence_length, 
      e_value, query_sequence_length, mobileOG_ID, gene_name, 
      best_hit_accession_ID, major_mobileOG_category, minor_mobileOG_category, 
      source_database, start, stop, strand, gc_content, coverage, midpoint_MGE
    )
  
  # Step 2: Prepare the ARG data for linkage (make sure not replicate column name)
  # We extract identity info and taxonomy to attach to the MGEs
  arg_subset <- df_arg %>%
    select(
      Contig, class, fa_name, midpoint_ARG, 
      Kingdom, Phylum, Class, Order, Family, Genera, Species, Genome
    )
  
  # Step 3: Join tables and calculate distance
  # Using many-to-many relationship to capture all ARG neighbors for each MGE
  linked_df <- mge_subset %>%
    inner_join(arg_subset, by = c("specific_contig" = "Contig"), relationship = "many-to-many") %>%
    # Calculate physical gap between midpoints
    mutate(gap_distance = abs(midpoint_MGE - midpoint_ARG)) %>%
    # Filter for proximity (within 10kb)
    filter(gap_distance <= 10000)
  
  return(linked_df)
}

# 3. Execute the linkage
df_final_linked <- link_MGE_ARG(df_mge, df_arg)

# 4. Export the results
write.xlsx(df_final_linked, 'C:/Users/CFL/Desktop/UK/PhD/output_file_MGE_ARG_links.xlsx')
