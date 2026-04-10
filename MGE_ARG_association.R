# ARG - MGE Linkage Analysis ####
library(readxl)
library(openxlsx)
library(dplyr)

# Load data
df1 <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = "ARG")
df_final <- read_excel('C:/Users/CFL/Desktop/UK/PhD/AMR/Archaea/Archaea_output.xlsx', sheet = "All_unique_ARG-carrying_MGE")

# Define the co-localization linkage function
link_MGE_ARG <- function(mge_df, arg_df) {
  
  # 1. Extract required information from the ARG dataframe
  # Note: utilizing the standardized 'midpoint_ARG'
  arg_subset <- arg_df %>%
    select(Contig, class, fa_name, midpoint_ARG, Kingdom, Phylum, Class, Order, Family, Genera, Species, Genome) %>%
    # Rename columns to prevent naming conflicts and ambiguity after joining
    rename(ARG_class = class, ARG_fa_name = fa_name)
  
  # 2. Generate all potential MGE-ARG combinations on the same contig
  # This creates a many-to-many relationship, capturing all possible pairs
  linked_df <- mge_df %>%
    inner_join(arg_subset, by = c("Specific_Contig" = "Contig"), relationship = "many-to-many")
  
  # 3. Calculate the physical distance and apply the 10 kb (10,000 bp) proximity threshold
  filtered_links <- linked_df %>%
    # Calculate the absolute distance between the midpoints of the MGE and the ARG
    mutate(gap_distance = abs(midpoint_MGE - midpoint_ARG)) %>%
    # Retain only the pairs that are co-localized within 10 kb
    filter(gap_distance <= 10000)
  
  return(filtered_links)
}

# Execute the linkage function
df_final_linked <- link_MGE_ARG(df_final, df1)

# Export the co-localization results
write.xlsx(df_final_linked, 'C:/Users/CFL/Desktop/UK/PhD/output_file_MGE_ARG_links.xlsx')
