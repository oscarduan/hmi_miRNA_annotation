setwd("C:/Users/duany/OneDrive/Desktop/Srivastava Lab/sRNA-seq/")
library(dplyr)
library(readxl)
library(tidyr)
library(Biostrings)
library(ggplot2)
library(patchwork)
library(stringr)
library(rtracklayer)
library(ggrepel)
library(purrr)
library(readr)
`%notin%` <- Negate(`%in%`)

mGDB.all <- read_excel("Isana.2023.11/mGDB.all.xlsx", col_names = F) %>% 
  select(c(1:4,7:10)) %>%
  mutate(species = substr(.[[1]], 1,3))
colnames(mGDB.all)[1:8] <- c("mGDB.miR","mB.miR","family","seed","chr","start","end","strand")

all.diff <- mGDB.all %>%
  filter(family %in% c("LET-7", "MIR-10")) %>%
  # filter(seed %in% c("GAGGUAG", "CCCUGAG")) %>%
  filter(chr != "None" | strand != "None") %>%
  group_by(species, chr, strand) %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end)) %>%
  arrange(species, chr, strand, start) %>%
  mutate(next_start = lead(start), # Calculate next start within the group
         diff = next_start - start - 1 ,
         min_diff_1 = ifelse(row_number() == n(), NA_real_, diff),
         # Calculate distance with specified conditions
         distance = case_when(
           row_number() == 1 ~ min_diff_1, # First row in the group
           row_number() == n() ~ lag(diff, default = NA_real_), # Last row in the group
           n() >= 3 ~ pmin(min_diff_1, lag(min_diff_1, default = NA_real_)),
           # Group has more than 3 members
           TRUE ~ min_diff_1 # Default case
         )) %>%
  select(-next_start, -diff, - min_diff_1)  # Optionally remove intermediate columns

# add assembly name
read.assembly <- function(species){
  gff_lines <- read_lines(paste("reference1/MirGeneDB/GFF/", species, ".gff", sep= ""))
  third_line <- gff_lines[3]
  assembly <- gsub("# genome-build-id:  ", "", third_line)
  assembly
}

all.diff$assembly <- lapply(all.diff$species, read.assembly)
# all.diff %>% filter(family %in% c("LET-7", "MIR-10")) %>% pull(species) %>% unique() %>% length()

cluster.diff00 <- all.diff %>% 
  filter(distance < 10000 & distance >= 45) %>%
  group_by(species, chr, strand) %>%
  mutate(cluster_span = if_else(n() > 1, max(end) - min(start), NA_real_)) %>%
  mutate(group_id = cur_group_id(),) %>%
  ungroup() 


## Adjust clusters---stretch----clusters into .1, .2 id
cluster.diff0 <- cluster.diff00 %>% filter(cluster_span > 10000) %>%
  group_by(species, chr, strand) %>%
  mutate(min_start = min(start),  # Calculate minimum start for each group
         group_id = ifelse(start - min_start < 10000, 
                           paste(group_id, ".1", sep = ""), 
                           paste(group_id, ".2", sep = ""))) %>%
  mutate(group_id = as.numeric(group_id)) %>%
  select(-min_start) %>%  # Remove the temporary min_start column
  
  ungroup() %>%
  bind_rows(cluster.diff00  %>% filter(cluster_span <= 10000)) %>%
  group_by(group_id) %>%  # fix the two species with 3 sub-clusters
  mutate(min_start = min(start),  
         group_id = ifelse(start - min_start < 10000, 
                           group_id, 
                           group_id + 0.1)) %>%
  mutate(group_id = as.numeric(group_id)) %>%
  ungroup %>%
  select(-min_start) %>%
  
  group_by(group_id) %>%
  mutate(family_count = n_distinct(family),
         seed_count = n_distinct(seed),
         member_count = n(),
         cluster_span = if_else(n() > 1, max(end) - min(start), NA_real_)) %>%
  filter(seed_count > 1) %>%
  ungroup()



### heterogenous in family
cluster.diff1 <- cluster.diff0 %>% filter(family_count == 1)
cluster.diff2 <- cluster.diff0 %>% filter(family_count == 2)
# 
# cluster.diff2 %>% pull(species) %>% unique() %>% length()

#Add the 10 non-cluster species
non.cluster.list <- all.diff %>% filter(species %notin% as.vector(cluster.diff2$species)) %>% pull(species) %>% unique()
non.cluster.df <- data.frame(species = non.cluster.list,
                             group_id = c(900:909),
                             cluster_span = 1, 
                             member_count = 1, 
                             min_cluster_span = 1)

plot1.df <- cluster.diff2 %>% 
  group_by(group_id, species) %>%
  summarize(cluster_span = max(cluster_span),
            member_count = max(member_count),
            chr = dplyr::first(chr),
            cluster_start = min(start),
            cluster_end = max(end),
            strand = dplyr::first(strand), 
            assembly =dplyr::first(assembly)) %>%
  ungroup() 

species.names <- read.csv("intron/cluster.all.adj.csv", row.names = NULL) %>%
  dplyr::select(species, species.full, species.description) %>%
  dplyr::mutate(
    species.full = dplyr::case_when(
      species == "Dno" ~ "Dasypus novemcinctus",
      species == "Gja" ~ "Gekko japonicus",
      species == "Mmu" ~ "Mus musculus",
      species == "Rno" ~ "Rattus norvegicus",
      species == "Sha" ~ "Sarcophilus harrisii",
      TRUE ~ species.full
    ),
    species.description = dplyr::case_when(
      species == "Dno" ~ "Nine-banded armadillo",
      species == "Gja" ~ "Schlegels Japanese gecko",
      species == "Mmu" ~ "House mouse",
      species == "Rno" ~ "Norway rat",
      species == "Sha" ~ "Tasmanian devil",
      TRUE ~ species.description
    )
  ) %>%
  dplyr::distinct()

plot2.df <- plot1.df %>%
  group_by(species) %>%
  mutate(min_cluster_span = min(cluster_span)) %>%
  ungroup() %>%
  mutate(species = factor(species, levels = unique(species))) %>% # Reorder species based on min_cluster_span
  bind_rows(non.cluster.df) %>%
  left_join(species.names, by = "species") %>%
  arrange(min_cluster_span) %>%
  mutate(species = factor(species, levels = unique(species)),
         species.full = factor(species.full, levels = unique(species.full)))  
# write.csv(plot2.df, "intron/cluster.all.csv", row.names = F)

manual_p1_rows <- tibble::tibble(
  species = c("Hmi", "Xbo", "Sro"),
  group_id = c(9011, 9012, 9013),
  cluster_span = c(984, 2216, 1),
  member_count = NA_real_,
  min_cluster_span = cluster_span,
  species.full = c(
    "Hofstenia miamia",
    "Xenoturbella bocki",
    "Symsagittifera roscoffensis"
  ),
  species.description = NA_character_
)

plot2.df_p1 <- dplyr::bind_rows(plot2.df, manual_p1_rows) %>%
  dplyr::arrange(min_cluster_span) %>%
  dplyr::mutate(
    species = factor(species, levels = unique(species)),
    species.full = factor(species.full, levels = unique(species.full))
  )

# plot2.df_p1 %>% pull(species) %>% unique() %>% length() #74

p1 <- ggplot(plot2.df_p1, aes(x = species.full, y = log10(cluster_span))) +
  geom_point(aes(size = log10(cluster_span) ), shape = 21, size = 2.75) +  # Assuming size by cluster_span makes more sense
  scale_size_continuous(range = c(1, 1.1)) +
  scale_y_continuous(limits = c(2,4), breaks = c(2,3,4)) +
  guides(size = "none") +
  theme_test() +
  labs(x = "Species", y = "Log10 Cluster Span (BP)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p1


# Fetch Fasta file for cluster ----
ncbi.list <- c("Asu","Cli","Cgi","Aae")
ensembl.list <- c("Bla", "Cte","Gga","Dme", "Dsi","Lgi", "Dan","Dya","Dmo", "Gmo","Dpu", "Dca", "Tni","Tca")
other.list <- c("Xbo")
error.list <- c()

manual.list <- c(ncbi.list, error.list, ensembl.list, other.list)


cluster.df1 <- plot2.df %>%
  mutate(search_id = paste(species, group_id, sep = "_")) %>%
  mutate(common = case_when(
    grepl("chr", chr, ignore.case = TRUE) ~ "common",
    grepl("LG", chr, ignore.case = TRUE) ~ "common",
    grepl("2L|2R", chr, ignore.case = TRUE) ~ "common",
    grepl("scaff", chr, ignore.case = TRUE) ~ "common",
    grepl("supercont", chr, ignore.case = TRUE) ~ "common",
    grepl("^[0-9]+$", chr) ~ "common",  # Checks if the column contains only numbers
    grepl("^[0-9]", chr) ~ "common",  # Checks if the column starts with a number
    TRUE ~ "non.common"  # Default case if none of the above conditions are met
  )) %>%
  mutate(common = ifelse(is.na(chr), "non_clustered", common))

non.common.cmd <- cluster.df1 %>% filter(common == "non.common") %>%
  filter(species %notin% error.list) %>%
  mutate(cmd= paste("efetch -db nucleotide -id ", chr, " -format fasta ",
                    "-seq_start ", cluster_start,
                    " -seq_stop ", cluster_end,
                    " | sed '1s/>.*/>", search_id, "/' > ", search_id, ".fasta",
                    sep = ""
  ))

#### Common chr names 
common.cmd <- cluster.df1 %>% filter(common == "common") %>%
  filter(species %notin% manual.list) %>%
  mutate( chr = case_when(
    assembly == "ci3" ~ paste("chr", chr, sep = ""), 
    species == "Tgu" ~ paste("chr", chr, sep = ""),
    species == "Aca" ~ paste("chr", chr, sep = ""),
    species == "Loc" ~ paste("chr", chr, sep = ""),
    species == "Mdo" ~ paste("chr", chr, sep = ""),
    TRUE ~ chr
  ) ) %>%
  mutate(assembly = str_replace_all(assembly, "_add", ""),
         assembly = str_replace_all(assembly, "_trace", "")) %>%
  mutate(assembly = str_sub(assembly, 1, 1) %>%
           str_to_lower() %>%
           paste0(str_sub(assembly, 2))) %>%
  mutate(cmd = paste("twoBitToFa http://hgdownload.cse.ucsc.edu/gbdb/", assembly,"/", assembly, ".2bit ",
                     search_id, ".fasta",
                     " -seq=", chr, 
                     " -start=", cluster_start,
                     " -end=", cluster_end, 
                     " && sed -i '1s/>.*/>", search_id, "/' ", search_id, ".fasta",
                     sep = ""))

#Export cmd
# writeLines(common.cmd$cmd, "intron/common.cmd.sh")
# writeLines(non.common.cmd$cmd, "intron/non.common.cmd.sh")

ensembl.df <- cluster.df1 %>% filter(species %in% ensembl.list)
ncbi.df <- cluster.df1 %>% filter(species %in% ncbi.list | species %in% error.list)

# write.csv(ensembl.df, "intron/ensembl.df.csv", row.names = F)
# write.csv(ncbi.df, "intron/ncbi.and.error.df.csv", row.names = F)


# Generate sequence for splicing site ----
dna_seqs <- readDNAStringSet("intron/let7.cluster.merge.fasta", format = "fasta")

# Create a dataframe
fasta_df <- data.frame(search_id = names(dna_seqs), sequence = as.character(dna_seqs), stringsAsFactors = FALSE)

seq.df <- cluster.df1 %>% 
  select(search_id, group_id, species, cluster_span, 
         species.full, species.description, 
         common, strand) %>%
  left_join(fasta_df, by = "search_id")


# Function to get reverse complement if strand is "-"
get_rev_comp <- function(strand, sequence) {
  if (is.na(sequence)) {
    return(NA)
  } else if (strand == "-") {
    return(as.character(reverseComplement(DNAString(sequence))))
  } else {
    return(sequence)
  }
}

seq.df <- seq.df %>%
  rowwise() %>% mutate(sequence = get_rev_comp(strand, sequence))


# Function to write each row to a FASTA format
write_fasta <- function(search_id, sequence, file_path) {
  # Open the file in append mode
  file_conn <- file(file_path, open = "a")
  
  # Write the header and sequence
  writeLines(paste0(">", search_id), file_conn)
  writeLines(sequence, file_conn)
  
  # Close the file connection
  close(file_conn)
}

seq.export.df <- seq.df %>% filter(!is.na(sequence))
apply(seq.export.df, 1, function(x) write_fasta(x["search_id"], x["sequence"], "intron/cluster.export.fasta"))

# Analyze ASSP ----
# Read the file
ASSP.high.lines <- readLines("intron/ASSP_6.8.txt")
ASSP.low.lines <- readLines("intron/ASSP_2.2_4.5.txt")

# Combine keywords into a regular expression
keywords <- "Sequence:|Constitutive donor|Constitutive acceptor"

# Find lines containing any of the keywords
indices <- grep(keywords, ASSP.high.lines)
indices1 <- grep(keywords, ASSP.low.lines)
# Subset lines based on indices
ASSP.high.filtered <- ASSP.high.lines[indices]
ASSP.low.filtered <- ASSP.low.lines[indices1]
# Optionally, write the filtered linindices1# Optionally, write the filtered lines to a new file
writeLines(ASSP.high.filtered, "intron/ASSP.high.fitered.txt")
writeLines(ASSP.low.filtered, "intron/ASSP.low.fitered.txt")

## Generate the df
# Load the data
data1 <- read_lines("intron/ASSP.high.fitered.txt")
data2 <- read_lines("intron/ASSP.low.fitered.txt")


# Initialize dataframe with an explicit structure matching the loop's output
library(stringr) # Required for str_remove and str_split

assp.high.df <- data.frame(position = integer(),
                           constitutive_type = character(),
                           site_sequence = character(),
                           score = numeric(),
                           cluster = character(),
                           stringsAsFactors = FALSE)

current_cluster <- NA

for (line in data1) {
  if (startsWith(trimws(line), "Sequence:")) {
    current_cluster <- str_remove(trimws(line), "Sequence: ")
  } else if (grepl("Constitutive", line)) {
    parts <- str_split(line, "\\|", simplify = TRUE)
    
    if (length(parts) >= 4) { 
      position <- as.integer(trimws(parts[1]))
      constitutive_type <- str_extract(trimws(parts[2]), "acceptor|donor")
      site_sequence <- trimws(parts[3])
      score <- as.numeric(trimws(parts[4]))
      
      # Bind the new row to the dataframe
      assp.high.df <- rbind(assp.high.df, 
                            data.frame(position = position, 
                                       constitutive_type = constitutive_type, 
                                       site_sequence = site_sequence, 
                                       score = score, 
                                       cluster = current_cluster, 
                                       stringsAsFactors = FALSE))
    }
  }
}

assp.high.intron <- assp.high.df %>%
  group_by(cluster) %>%
  summarise(
    min_donor_position = ifelse(any(constitutive_type == "donor"), 
                                as.character(min(position[constitutive_type == "donor"], na.rm = TRUE)), 
                                "NO Donor"),
    max_acceptor_position = ifelse(any(constitutive_type == "acceptor"), 
                                   as.character(max(position[constitutive_type == "acceptor"], na.rm = TRUE)), 
                                   "NO Acceptor")
  ) %>%
  mutate(contain_intron = ifelse(min_donor_position != "NO Donor" & 
                                   max_acceptor_position != "NO Acceptor" & 
                                   as.numeric(max_acceptor_position) > as.numeric(min_donor_position), TRUE, FALSE))

# combine
cluster.df2 <- cluster.df1 %>%
  left_join(assp.high.intron %>% dplyr::rename(search_id = cluster), by = "search_id") %>%
  mutate(contain_intron = ifelse(is.na(min_donor_position), FALSE, contain_intron))

mir.df0 <- cluster.diff0 %>% # No_Intron
  mutate(search_id = paste(species, group_id, sep = "_")) %>%
  left_join(assp.high.intron %>% dplyr::rename(search_id = cluster), by = "search_id") %>%
  mutate(contain_intron = ifelse(is.na(min_donor_position), FALSE, contain_intron)) %>%
  mutate(Intron_Type = ifelse(contain_intron == FALSE, "No_Intron", 
                              ifelse(member_count <3, "Shorten", "TBD")))


cluster.removal.tbd.df <- mir.df0 %>% filter(Intron_Type == "TBD") %>%
  group_by(search_id) %>%
  arrange(search_id, start) %>%
  summarise(middle.start = nth(start, 2),
            first.mGDB = nth(mGDB.miR,1),
            second.mGDB = nth(mGDB.miR,2),
            third.mGDB = nth(mGDB.miR,3),
            first.mB = nth(mB.miR,1),
            second.mB = nth(mB.miR,2),
            third.mB = nth(mB.miR,3)) 

cluster.removal <- cluster.removal.tbd.df %>%
  left_join(cluster.df2 %>% 
              filter(search_id %in% as.vector(cluster.removal.tbd.df$search_id)),
            by = "search_id" ) %>%
  mutate(middle.start = as.numeric(middle.start)) %>%
  mutate(middle.position = case_when(
    strand == "+"~ middle.start- cluster_start +1,
    strand == "-"~ cluster_end - middle.start +1
  )) %>%
  
  mutate(Intron_Type = ifelse(min_donor_position <= (middle.position + 21) &
                                max_acceptor_position >= (middle.position +21), "Removal", "Shorten"))

cluster.merge.1 <- cluster.removal %>% select(search_id, Intron_Type)
cluster.merge.2 <- mir.df0 %>% filter(Intron_Type != "TBD") %>% select(search_id, Intron_Type)

cluster.merge <- rbind(cluster.merge.1, cluster.merge.2)

cluster.df3 <- cluster.df2 %>% left_join(cluster.merge, by= "search_id") %>%
  mutate(Intron_Type = ifelse(is.na(Intron_Type), "No_Intron", Intron_Type))

## remove Xenoburbella

cluster.df4 <- cluster.df3 %>%
  filter(species != "Xbo")

get_opentree_phylogeny <- function(species_names, cache_path) {
  species_names <- sort(unique(as.character(species_names)))

  graft_known_unplaced_species <- function(phylo_tree, unplaced_species) {
    schmidtea_species <- intersect("Schmidtea mediterranea", unplaced_species)
    if (length(schmidtea_species) == 0) {
      attr(phylo_tree, "unplaced_species") <- unplaced_species
      return(phylo_tree)
    }

    protostomia_tip_candidates <- c(
      "Aedes aegypti",
      "Anopheles gambiae",
      "Apis mellifera",
      "Bombyx mori",
      "Caenorhabditis elegans",
      "Capitella teleta",
      "Crassostrea gigas",
      "Drosophila ananassae",
      "Drosophila melanogaster",
      "Drosophila mojavensis",
      "Drosophila pseudoobscura",
      "Drosophila simulans",
      "Drosophila yakuba",
      "Lottia gigantea",
      "Octopus bimaculoides",
      "Tribolium castaneum"
    )
    present_protostomes <- intersect(protostomia_tip_candidates, phylo_tree$tip.label)

    if (length(present_protostomes) >= 2) {
      graft_node <- ape::getMRCA(phylo_tree, present_protostomes)
    } else if (length(present_protostomes) == 1) {
      graft_node <- which(phylo_tree$tip.label == present_protostomes)
    } else {
      graft_node <- length(phylo_tree$tip.label) + 1
    }

    for (species_to_graft in schmidtea_species) {
      graft_tip <- ape::read.tree(text = paste0("('", species_to_graft, "');"))
      phylo_tree <- ape::bind.tree(
        phylo_tree,
        graft_tip,
        where = graft_node,
        position = 0
      )
      phylo_tree$tip.label <- stringr::str_remove_all(phylo_tree$tip.label, "^'|'$")
    }

    unplaced_species <- setdiff(unplaced_species, schmidtea_species)
    attr(phylo_tree, "unplaced_species") <- unplaced_species
    phylo_tree
  }

  if (file.exists(cache_path)) {
    cached_tree <- ape::read.tree(cache_path)
    if (setequal(cached_tree$tip.label, species_names)) {
      return(ape::ladderize(cached_tree, right = FALSE))
    }
  }

  tnrs_response <- httr::POST(
    "https://api.opentreeoflife.org/v3/tnrs/match_names",
    body = list(names = species_names),
    encode = "json"
  )
  if (httr::status_code(tnrs_response) != 200) {
    stop("OpenTree TNRS request failed with status: ", httr::status_code(tnrs_response))
  }

  tnrs <- httr::content(tnrs_response, as = "parsed", simplifyVector = FALSE)
  match.df <- purrr::map_dfr(
    tnrs$results,
    function(result) {
      if (length(result$matches) == 0) {
        return(data.frame(species.full = result$name, ott_id = NA_real_))
      }
      data.frame(
        species.full = result$name,
        ott_id = result$matches[[1]]$taxon$ott_id
      )
    }
  )

  unmatched_species <- match.df %>%
    dplyr::filter(is.na(ott_id)) %>%
    dplyr::pull(species.full)
  if (length(unmatched_species) > 0) {
    stop("OpenTree could not match species: ", paste(unmatched_species, collapse = ", "))
  }

  unplaced_species <- character()
  subtree_response <- NULL

  repeat {
    subtree_response <- httr::POST(
      "https://api.opentreeoflife.org/v3/tree_of_life/induced_subtree",
      body = list(ott_ids = as.list(match.df$ott_id)),
      encode = "json"
    )

    if (httr::status_code(subtree_response) == 200) {
      break
    }

    subtree_error <- httr::content(subtree_response, as = "text", encoding = "UTF-8")
    pruned_ott_ids <- unique(as.numeric(stringr::str_remove_all(
      stringr::str_extract_all(subtree_error, "ott[0-9]+")[[1]],
      "ott"
    )))

    if (length(pruned_ott_ids) == 0 || !any(pruned_ott_ids %in% match.df$ott_id)) {
      stop(
        "OpenTree induced subtree request failed with status: ",
        httr::status_code(subtree_response),
        "\n",
        subtree_error
      )
    }

    newly_unplaced <- match.df$species.full[match(pruned_ott_ids, match.df$ott_id)]
    unplaced_species <- unique(c(unplaced_species, newly_unplaced))
    match.df <- match.df %>%
      dplyr::filter(!ott_id %in% pruned_ott_ids)

    if (nrow(match.df) < 2) {
      stop("OpenTree returned fewer than two placeable species after removing pruned taxa.")
    }
  }

  subtree <- httr::content(subtree_response, as = "parsed", simplifyVector = FALSE)
  phylo_tree <- ape::read.tree(text = subtree$newick)
  phylo_tree <- ape::collapse.singles(phylo_tree)

  original_tip_labels <- phylo_tree$tip.label
  tip_ott <- as.numeric(stringr::str_match(original_tip_labels, "_ott([0-9]+)$")[, 2])
  mapped_tip_labels <- match.df$species.full[match(tip_ott, match.df$ott_id)]

  if (any(is.na(mapped_tip_labels))) {
    cleaned_tip_labels <- original_tip_labels %>%
      stringr::str_remove_all("^'|'$") %>%
      stringr::str_remove("_ott[0-9]+$") %>%
      stringr::str_remove("\\s+ott[0-9]+$") %>%
      stringr::str_remove("\\s+\\([^\\)]+\\)$") %>%
      stringr::str_replace_all("_", " ")
    mapped_tip_labels[is.na(mapped_tip_labels)] <- cleaned_tip_labels[is.na(mapped_tip_labels)]
  }

  if (any(!mapped_tip_labels %in% match.df$species.full)) {
    stop(
      "Could not map all OpenTree tips back to species names: ",
      paste(original_tip_labels[!mapped_tip_labels %in% match.df$species.full], collapse = ", ")
    )
  }

  phylo_tree$tip.label <- mapped_tip_labels

  if (length(unplaced_species) == 0) {
    ape::write.tree(phylo_tree, file = cache_path)
  }

  phylo_tree <- ape::ladderize(phylo_tree, right = FALSE)
  phylo_tree <- graft_known_unplaced_species(phylo_tree, unplaced_species)
  phylo_tree <- ape::ladderize(phylo_tree, right = FALSE)
  phylo_tree
}

make_phylo_tree_segments <- function(phylo_tree, tip.df) {
  n_tip <- length(phylo_tree$tip.label)
  n_node <- phylo_tree$Nnode
  all_nodes <- seq_len(n_tip + n_node)
  root_node <- n_tip + 1

  node_x <- rep(NA_real_, length(all_nodes))
  node_y <- rep(NA_real_, length(all_nodes))
  node_depth <- rep(NA_real_, length(all_nodes))

  node_x[seq_len(n_tip)] <- tip.df$species_x[match(phylo_tree$tip.label, tip.df$species.full)]
  node_y[seq_len(n_tip)] <- tip.df$tip_y[match(phylo_tree$tip.label, tip.df$species.full)]
  node_depth[root_node] <- 0

  for (i in seq_len(nrow(phylo_tree$edge))) {
    parent <- phylo_tree$edge[i, 1]
    child <- phylo_tree$edge[i, 2]
    node_depth[child] <- node_depth[parent] + 1
  }

  get_descendant_tips <- function(node) {
    children <- phylo_tree$edge[phylo_tree$edge[, 1] == node, 2]
    if (length(children) == 0) {
      return(node)
    }
    unlist(lapply(children, get_descendant_tips))
  }

  internal_nodes <- seq.int(n_tip + 1, n_tip + n_node)
  for (node in internal_nodes) {
    node_x[node] <- mean(node_x[get_descendant_tips(node)])
  }

  tree_y_root <- 0.88
  tree_y_internal_top <- 1.22
  max_internal_depth <- max(node_depth[internal_nodes], na.rm = TRUE)
  if (max_internal_depth == 0) {
    node_y[internal_nodes] <- tree_y_root
  } else {
    node_y[internal_nodes] <- tree_y_root +
      (node_depth[internal_nodes] / max_internal_depth) * (tree_y_internal_top - tree_y_root)
  }

  dplyr::bind_rows(
    data.frame(
      x = node_x[phylo_tree$edge[, 2]],
      xend = node_x[phylo_tree$edge[, 2]],
      y = node_y[phylo_tree$edge[, 1]],
      yend = node_y[phylo_tree$edge[, 2]]
    ),
    data.frame(
      x = node_x[phylo_tree$edge[, 1]],
      xend = node_x[phylo_tree$edge[, 2]],
      y = node_y[phylo_tree$edge[, 1]],
      yend = node_y[phylo_tree$edge[, 1]]
    )
  )
}

p2_excluded_species <- c("Nve")

manual_cluster.df4_short_rows <- tibble::tibble(
  species.full = c(
    "Hofstenia miamia",
    "Xenoturbella bocki",
    "Symsagittifera roscoffensis"
  ),
  species_x = c(71, 72, 73),
  cluster_span = c(
    984, 2216, 1
  ),
  Intron_Type = c("Shorten","Shorten","No_Intron")
)

current_species_names <- dplyr::bind_rows(
  cluster.df4 %>%
    dplyr::filter(!species %in% p2_excluded_species) %>%
    dplyr::select(species.full),
  manual_cluster.df4_short_rows %>%
    dplyr::select(species.full)
) %>%
  dplyr::filter(!is.na(species.full)) %>%
  dplyr::pull(species.full) %>%
  as.character() %>%
  unique() %>%
  sort()

clade <- function(...) list(...)
newick_tip <- function(label) paste0("'", label, "'")
clade_to_newick <- function(x) {
  if (is.character(x)) {
    return(newick_tip(x))
  }
  paste0("(", paste(vapply(x, clade_to_newick, character(1)), collapse = ","), ")")
}

drosophila_clade <- clade(
  clade(clade("Drosophila melanogaster", "Drosophila simulans"), "Drosophila yakuba"),
  clade("Drosophila ananassae", "Drosophila mojavensis")
)
hexapoda_clade <- clade(
  "Blattella germanica",
  clade("Tribolium castaneum", clade("Heliconius melpomene", clade("Aedes aegypti", drosophila_clade)))
)
arthropoda_clade <- clade(
  clade("Ixodes scapularis", clade("Centruroides sculpturatus", "Limulus polyphemus")),
  clade(clade("Daphnia pulex", "Daphnia magna"), hexapoda_clade)
)
ecdysozoa_clade <- clade(
  clade("Ascaris suum", clade("Caenorhabditis briggsae", "Caenorhabditis elegans")),
  arthropoda_clade
)

mollusca_clade <- clade(
  clade("Nautilus pompilius", clade("Euprymna scolopes", clade("Octopus bimaculoides", "Octopus vulgaris"))),
  clade("Lottia gigantea", "Crassostrea gigas")
)
spiralia_clade <- clade(
  clade("Schmidtea mediterranea", "Brachionus plicatilis"),
  clade("Lingula anatina", clade(clade("Capitella teleta", "Eisenia fetida"), mollusca_clade))
)
protostomia_clade <- clade(ecdysozoa_clade, spiralia_clade)

ambulacraria_clade <- clade(
  clade("Saccoglossus kowalevskii", "Ptychodera flava"),
  clade("Strongylocentrotus purpuratus", "Patiria miniata")
)
actinopterygii_clade <- clade(
  "Lepisosteus oculatus",
  clade("Danio rerio", clade("Monopterus albus", clade("Tetraodon nigroviridis", "Gadus morhua")))
)
mammalia_clade <- clade(
  "Ornithorhynchus anatinus",
  clade(
    clade("Monodelphis domestica", "Sarcophilus harrisii"),
    clade(
      clade("Dasypus novemcinctus", "Echinops telfairi"),
      clade(
        clade("Bos taurus", "Canis familiaris"),
        clade(clade("Homo sapiens", "Macaca mulatta"), clade(clade(clade("Mus musculus", "Rattus norvegicus"), "Cavia porcellus"), "Oryctolagus cuniculus"))
      )
    )
  )
)
sauropsida_clade <- clade(
  clade(clade("Gallus gallus", clade("Columba livia", "Taeniopygia guttata")), "Alligator mississippiensis"),
  clade("Chrysemys picta bellii", clade("Sphenodon punctatus", clade("Anolis carolinensis", clade("Python bivittatus", "Gekko japonicus"))))
)
tetrapoda_clade <- clade(
  clade(clade("Xenopus tropicalis", "Xenopus laevis"), "Microcaecilia unicolor"),
  clade(mammalia_clade, sauropsida_clade)
)
sarcopterygii_clade <- clade("Latimeria chalumnae", tetrapoda_clade)
osteichthyes_clade <- clade(actinopterygii_clade, sarcopterygii_clade)
gnathostomata_clade <- clade(clade("Callorhinchus milii", "Scyliorhinus torazame"), osteichthyes_clade)
vertebrata_clade <- clade(clade("Eptatretus burgeri", "Petromyzon marinus"), gnathostomata_clade)
chordata_clade <- clade(
  clade("Branchiostoma floridae", "Branchiostoma lanceolatum"),
  clade("Ciona intestinalis", vertebrata_clade)
)
deuterostomia_clade <- clade(ambulacraria_clade, chordata_clade)
xenacoelomorpha_clade <- clade(
  clade("Hofstenia miamia", "Symsagittifera roscoffensis"),
  "Xenoturbella bocki"
)

mirgenedb_species_clade <- clade(
  "Nematostella vectensis",
  clade(xenacoelomorpha_clade, clade(protostomia_clade, deuterostomia_clade))
)
mirgenedb_species_newick <- paste0(clade_to_newick(mirgenedb_species_clade), ";")

mirgenedb_species_tree <- ape::read.tree(text = mirgenedb_species_newick)
mirgenedb_species_tree$tip.label <- stringr::str_remove_all(mirgenedb_species_tree$tip.label, "^'|'$")
missing_from_mirgenedb_tree <- setdiff(current_species_names, mirgenedb_species_tree$tip.label)
if (length(missing_from_mirgenedb_tree) > 0) {
  warning(
    "Species missing from MirGeneDB-compatible tree and omitted from p2/species tree PDF: ",
    paste(missing_from_mirgenedb_tree, collapse = ", ")
  )
}

species_tree <- ape::drop.tip(
  mirgenedb_species_tree,
  setdiff(mirgenedb_species_tree$tip.label, current_species_names)
)
species_tree <- ape::ladderize(species_tree, right = FALSE)
unplaced_species <- character()

phylo_species_order <- if (length(current_species_names) >= 2) {
  species_tree$tip.label
} else {
  current_species_names
}

phylo.species.df <- data.frame(
  species.full = phylo_species_order,
  species_x = seq_along(phylo_species_order),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    species_label = species.full,
    species_label_y = 1.50 + pmax(0, 0.34 - 0.012 * nchar(species_label)),
    tip_y = species_label_y
  )

cluster.df4 <- cluster.df4 %>%
  dplyr::left_join(
    phylo.species.df %>% dplyr::select(species.full, species_x),
    by = "species.full"
  )

cluster.df4_short <- cluster.df4 %>%
  # dplyr::filter(!species %in% c("Sme", "Nve")) %>%
  dplyr::filter(!species %in% p2_excluded_species) %>%
  dplyr::select(species.full, species_x, cluster_span, Intron_Type)

manual_cluster.df4_short_rows <- manual_cluster.df4_short_rows %>%
  dplyr::select(species.full, cluster_span, Intron_Type) %>%
  dplyr::left_join(
    phylo.species.df %>% dplyr::select(species.full, species_x),
    by = "species.full"
  ) %>%
  dplyr::select(species.full, species_x, cluster_span, Intron_Type)

cluster.df4_short <- dplyr::bind_rows(
  cluster.df4_short,
  manual_cluster.df4_short_rows
)

# cluster.df4_short %>% filter(Intron_Type != "No_Intron") %>% pull(species.full) %>% unique() %>% length() #38
# cluster.df4_short %>% filter(Intron_Type == "Removal") %>% pull(species.full) %>% unique() %>% length() #10



phylo.tree.df <- if (length(current_species_names) >= 2) {
  make_phylo_tree_segments(species_tree, phylo.species.df)
} else {
  phylo.species.df %>%
    dplyr::transmute(x = species_x, xend = species_x, y = 0.88, yend = tip_y)
}

if (length(unplaced_species) > 0) {
  phylo.tree.df <- dplyr::bind_rows(
    phylo.tree.df,
    phylo.species.df %>%
      dplyr::filter(species.full %in% unplaced_species) %>%
      dplyr::transmute(x = species_x, xend = species_x, y = 0.88, yend = tip_y)
  )
}
##
p2 <- ggplot(cluster.df4_short, aes(x = species_x, y = log10(cluster_span))) +
  geom_segment(
    data = phylo.species.df,
    aes(x = species_x, xend = species_x, y = 0.78, yend = 4),
    inherit.aes = FALSE,
    linewidth = 0.18,
    color = "grey88"
  ) +
  geom_segment(
    data = phylo.tree.df,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    linewidth = 0.28,
    color = "grey35"
  ) +
  geom_text(
    data = phylo.species.df,
    aes(x = species_x, y = species_label_y, label = species_label),
    inherit.aes = FALSE,
    angle = 90,
    hjust = 0,
    vjust = 0.5,
    size = 4
  ) +
  geom_point(aes(size = log10(cluster_span), fill = Intron_Type ), shape = 21, size = 2.5) +  # Assuming size by cluster_span makes more sense
  scale_size_continuous(range = c(1, 1.1)) +
  scale_x_continuous(
    breaks = phylo.species.df$species_x,
    labels = NULL,
    expand = ggplot2::expansion(add = 0.6)
  ) +
  scale_y_continuous(limits = c(0.78, 4), breaks = c(2,3,4)) +
  coord_cartesian(ylim = c(0.78, 4), clip = "off") +
  guides(size = "none") +
  theme_test() +
  labs(x = "Species", y = "Log10 Cluster Span (BP)") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(5.5, 5.5, 18, 5.5)
  ) +
  scale_fill_manual(values = c("Removal" =  "red2", "Shorten" = "blue3", "No_Intron" = "grey75"))

p2

# p2 with the low-threshold ASSP input ----
assp.low.df <- data.frame(position = integer(),
                          constitutive_type = character(),
                          site_sequence = character(),
                          score = numeric(),
                          cluster = character(),
                          stringsAsFactors = FALSE)

current_cluster <- NA

for (line in data2) {
  if (startsWith(trimws(line), "Sequence:")) {
    current_cluster <- str_remove(trimws(line), "Sequence: ")
  } else if (grepl("Constitutive", line)) {
    parts <- str_split(line, "\\|", simplify = TRUE)

    if (length(parts) >= 4) {
      position <- as.integer(trimws(parts[1]))
      constitutive_type <- str_extract(trimws(parts[2]), "acceptor|donor")
      site_sequence <- trimws(parts[3])
      score <- as.numeric(trimws(parts[4]))

      assp.low.df <- rbind(assp.low.df,
                           data.frame(position = position,
                                      constitutive_type = constitutive_type,
                                      site_sequence = site_sequence,
                                      score = score,
                                      cluster = current_cluster,
                                      stringsAsFactors = FALSE))
    }
  }
}

assp.low.intron <- assp.low.df %>%
  group_by(cluster) %>%
  summarise(
    min_donor_position = ifelse(any(constitutive_type == "donor"),
                                as.character(min(position[constitutive_type == "donor"], na.rm = TRUE)),
                                "NO Donor"),
    max_acceptor_position = ifelse(any(constitutive_type == "acceptor"),
                                   as.character(max(position[constitutive_type == "acceptor"], na.rm = TRUE)),
                                   "NO Acceptor")
  ) %>%
  mutate(contain_intron = ifelse(min_donor_position != "NO Donor" &
                                   max_acceptor_position != "NO Acceptor" &
                                   as.numeric(max_acceptor_position) > as.numeric(min_donor_position), TRUE, FALSE))

cluster.df2_low <- cluster.df1 %>%
  left_join(assp.low.intron %>% dplyr::rename(search_id = cluster), by = "search_id") %>%
  mutate(contain_intron = ifelse(is.na(min_donor_position), FALSE, contain_intron))

mir.df0_low <- cluster.diff0 %>%
  mutate(search_id = paste(species, group_id, sep = "_")) %>%
  left_join(assp.low.intron %>% dplyr::rename(search_id = cluster), by = "search_id") %>%
  mutate(contain_intron = ifelse(is.na(min_donor_position), FALSE, contain_intron)) %>%
  mutate(Intron_Type = ifelse(contain_intron == FALSE, "No_Intron",
                              ifelse(member_count <3, "Shorten", "TBD")))

cluster.removal.tbd.df_low <- mir.df0_low %>% filter(Intron_Type == "TBD") %>%
  group_by(search_id) %>%
  arrange(search_id, start) %>%
  summarise(middle.start = nth(start, 2),
            first.mGDB = nth(mGDB.miR,1),
            second.mGDB = nth(mGDB.miR,2),
            third.mGDB = nth(mGDB.miR,3),
            first.mB = nth(mB.miR,1),
            second.mB = nth(mB.miR,2),
            third.mB = nth(mB.miR,3))

cluster.removal_low <- cluster.removal.tbd.df_low %>%
  left_join(cluster.df2_low %>%
              filter(search_id %in% as.vector(cluster.removal.tbd.df_low$search_id)),
            by = "search_id" ) %>%
  mutate(middle.start = as.numeric(middle.start)) %>%
  mutate(middle.position = case_when(
    strand == "+"~ middle.start- cluster_start +1,
    strand == "-"~ cluster_end - middle.start +1
  )) %>%
  mutate(Intron_Type = ifelse(min_donor_position <= (middle.position + 21) &
                                max_acceptor_position >= (middle.position +21), "Removal", "Shorten"))

cluster.merge.1_low <- cluster.removal_low %>% select(search_id, Intron_Type)
cluster.merge.2_low <- mir.df0_low %>% filter(Intron_Type != "TBD") %>% select(search_id, Intron_Type)

cluster.merge_low <- rbind(cluster.merge.1_low, cluster.merge.2_low)

cluster.df3_low <- cluster.df2_low %>% left_join(cluster.merge_low, by= "search_id") %>%
  mutate(Intron_Type = ifelse(is.na(Intron_Type), "No_Intron", Intron_Type))

cluster.df4_low <- cluster.df3_low %>%
  filter(species != "Xbo")

cluster.df4_low <- cluster.df4_low %>%
  dplyr::left_join(
    phylo.species.df %>% dplyr::select(species.full, species_x),
    by = "species.full"
  )

cluster.df4_short_low <- cluster.df4_low %>%
  dplyr::filter(!species %in% p2_excluded_species) %>%
  dplyr::select(species.full, species_x, cluster_span, Intron_Type)

manual_cluster.df4_short_rows_low <- manual_cluster.df4_short_rows %>%
  dplyr::select(species.full, species_x, cluster_span, Intron_Type)

cluster.df4_short_low <- dplyr::bind_rows(
  cluster.df4_short_low,
  manual_cluster.df4_short_rows_low
)


cluster.df4_short_low %>% filter(Intron_Type != "No_Intron") %>% pull(species.full) %>% unique() %>% length() #43
cluster.df4_short_low %>% filter(Intron_Type == "Removal") %>% pull(species.full) %>% unique() %>% length() #11



p2_low <- ggplot(cluster.df4_short_low, aes(x = species_x, y = log10(cluster_span))) +
  geom_segment(
    data = phylo.species.df,
    aes(x = species_x, xend = species_x, y = 0.78, yend = 4),
    inherit.aes = FALSE,
    linewidth = 0.18,
    color = "grey88"
  ) +
  geom_segment(
    data = phylo.tree.df,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    linewidth = 0.28,
    color = "grey35"
  ) +
  geom_text(
    data = phylo.species.df,
    aes(x = species_x, y = species_label_y, label = species_label),
    inherit.aes = FALSE,
    angle = 90,
    hjust = 0,
    vjust = 0.5,
    size = 4
  ) +
  geom_point(aes(size = log10(cluster_span), fill = Intron_Type ), shape = 21, size = 2.5) +
  scale_size_continuous(range = c(1, 1.1)) +
  scale_x_continuous(
    breaks = phylo.species.df$species_x,
    labels = NULL,
    expand = ggplot2::expansion(add = 0.6)
  ) +
  scale_y_continuous(limits = c(0.78, 4), breaks = c(2,3,4)) +
  coord_cartesian(ylim = c(0.78, 4), clip = "off") +
  guides(size = "none") +
  theme_test() +
  labs(x = "Species", y = "Log10 Cluster Span (BP)") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(5.5, 5.5, 18, 5.5)
  ) +
  scale_fill_manual(values = c("Removal" =  "red2", "Shorten" = "blue3", "No_Intron" = "grey75"))

p2_low

# Export plots ----
ggsave("intron/Export/splice.site.p1.tiff", p1, width = 455, height = 165, dpi = 300, units = "mm")
ggsave("intron/Export/splice.site.p1.pdf", p1, width = 455, height = 165, units = "mm")
ggsave("intron/Export/splice.site.tiff", p2, width = 455, height = 165, dpi = 300, units = "mm")
ggsave("intron/Export/splice.site.phylogeny.pdf", p2, width = 455, height = 165, units = "mm")
ggsave("intron/Export/splice.site.low.tiff", p2_low, width = 455, height = 165, dpi = 300, units = "mm")
ggsave("intron/Export/splice.site.phylogeny.low.pdf", p2_low, width = 455, height = 165, units = "mm")
