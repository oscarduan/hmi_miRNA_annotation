setwd("C:/Users/duany/OneDrive/Desktop/Srivastava Lab/sRNA-seq/")
library(dplyr)
library(readxl)
library(stringr)
library(Biostrings)
library(tibble)
library(gridExtra)
library(ggplot2)
library(DT)
library(GenomicRanges)
`%notin%` <- Negate(`%in%`)
source("C:/Users/duany/OneDrive/Desktop/Srivastava Lab/sRNA-seq/Isana.2023.11/scripts/B00.head.R")

# 1 two functions ----
mGDB.prep <- function(species1,file){
  input <- readRNAStringSet(file)
  df <- data.frame(mGDB.miR = names(input), sequence = as.character(input), row.names = NULL)
  df$species <- species1
  df$seed = substr(df$sequence, start = 2, stop = 8)
  df$nonseed = substr(df$sequence, start = 11, stop = length(df$sequence))
  
  match.ref <- read.csv("reference1/mGDB.mB.name.csv", row.names = NULL) %>% select(everything()[-6]) %>%
    dplyr::rename(mGDB.fam = family, species = sepecis)
  
  df1 <- df  %>% 
    mutate(mGDB.3p5p = substr(mGDB.miR, nchar(mGDB.miR)-1, nchar(mGDB.miR))) %>%
    mutate(mGDB.miR.full = mGDB.miR) %>%
    mutate(mGDB.miR = substr(mGDB.miR, 1, nchar(mGDB.miR)-3)) %>%# Remove the "-3p" or "-5p"
    left_join(match.ref %>% select(mGDB.miR, mB.miR, mGDB.fam), by = "mGDB.miR") %>%
    mutate(mGDB.fam = tolower(mGDB.fam)) %>%
    select(mB.miR, mGDB.miR, mGDB.fam, everything())
  df1
}

nonseed.blast <- function(bait.seed, bait.nonseed, ref.miR.df){
  homolog.df <- ref.miR.df %>% filter(seed == bait.seed)  # subset hsa miRs with same seed as the cel bait
  bait.str <- RNAString(bait.nonseed) # get ready to blast
  ref.nonseed.pool <- RNAStringSet(homolog.df$nonseed)  #generate the non-seed sequence pool to blast (ref)
  
  score <- sapply(ref.nonseed.pool, function(seq){     # blast
    alignment <- pairwiseAlignment(bait.str, seq)
    score(alignment)
  })
  
  index_most_similar <- which.max(score) # highest scored sequence
  most_similar_sequence_with_context <- homolog.df[index_most_similar,]
  return(most_similar_sequence_with_context$mGDB.miR)
}

# 2 load Hmia input ----
filter.df9 <- read.csv("PB_sRNA/Key_files/filter.df10.csv", header = T)

plot.df9 <- filter.df9 %>%
  select(-X3p.5p_ziv) %>%
  left_join(count.df, by = "PBid") %>%
  select(PBid, X3p.5p_ziv, X3p5p_mod, total_mature, total_star, Mature_ziv, Star_ziv, Hairpin_seq_trimmed_ziv, Chr, Start, End ) %>%
  mutate(total_reads = total_mature + total_star) 

# modify mature strand & Extract seed
revise.df <- plot.df9 %>%
  mutate(X3p5p_swap = ifelse(X3p.5p_ziv == X3p5p_mod, FALSE, TRUE)) %>%
  mutate(Mature_mod = ifelse(X3p5p_swap == TRUE, Star_ziv, Mature_ziv),
         Star_mod = ifelse(X3p5p_swap == TRUE, Mature_ziv, Star_ziv)) %>% 
  mutate(
    # seed6      = substr(Mature_mod, 2, 7),
    seed      = substr(Mature_mod, 2, 8),
    # seed6_alt  = substr(Star_mod,   2, 7),
    seed_alt  = substr(Star_mod,   2, 8)
  ) 

# working df
hmia.df <- revise.df %>% select(-X3p.5p_ziv, -X3p5p_swap, - Hairpin_seq_trimmed_ziv, -Mature_ziv, -Star_ziv) 

revise.df00 <- plot.df9 %>%
  select(PBid, X3p.5p_ziv, X3p5p_mod, total_mature, total_star, Mature_ziv, Star_ziv, Hairpin_seq_trimmed_ziv, Chr, Start, End ) %>%
  mutate(total_reads = total_mature + total_star) 

isomiR.revise_coord.df <- read.csv("PB_sRNA/File_Export/isomiR_custom/revise_coord.df.csv", header = T)

revise.df0a <- revise.df00 %>% 
  filter(!PBid %in% isomiR.revise_coord.df$PBid)

revise.df0b <- revise.df00 %>% 
  filter(PBid %in% isomiR.revise_coord.df$PBid) %>%
  left_join(isomiR.revise_coord.df %>% select(PBid, isomiR_seq), by = "PBid") %>%
  mutate(Mature_ziv = isomiR_seq) %>% 
  select(-isomiR_seq)

revise.df0 <- bind_rows(revise.df0a, revise.df0b)

# modify mature strand & Extract seed
revise.df <- revise.df0 %>%
  mutate(X3p5p_swap = ifelse(X3p.5p_ziv == X3p5p_mod, FALSE, TRUE)) %>%
  mutate(Mature_mod = ifelse(X3p5p_swap == TRUE, Star_ziv, Mature_ziv),
         Star_mod = ifelse(X3p5p_swap == TRUE, Mature_ziv, Star_ziv)) %>% 
  mutate(
    # seed6      = substr(Mature_mod, 2, 7),
    seed      = substr(Mature_mod, 2, 8),
    # seed6_alt  = substr(Star_mod,   2, 7),
    seed_alt  = substr(Star_mod,   2, 8)
  ) 

# working df
hmia.df <- revise.df %>% select(-X3p.5p_ziv, -X3p5p_swap, - Hairpin_seq_trimmed_ziv, -Mature_ziv, -Star_ziv) 


# 3 Hsa 7mer ----
hsa.mGDB.df <- mGDB.prep("Hsa", "reference1/MirGeneDB/mGDB_3/mature/hsa.fas") 

hsa.ref <- hsa.mGDB.df %>% 
  group_by(mGDB.miR.full) %>%  # Group by PBid
  dplyr::slice_max(order_by = mGDB.fam, n = 1,
                   with_ties = F) %>% 
  ungroup() %>% 
  as.data.frame() %>%
  group_by(mGDB.fam, seed) %>%
  mutate(seed_count = n()) %>%
  ungroup() 

seed28.hsa <- hmia.df %>% 
  # filter(PBid %notin% match.xacl$PBid) %>%
  left_join(hsa.ref, by = "seed") %>% 
  filter(is.na(mGDB.miR) == F ) %>%
  mutate(seed.match = "7mer")

unique28.hsa <- seed28.hsa %>%
  group_by(PBid) %>%  # Group by PBid
  dplyr::slice_max(order_by = total_reads, n = 1,
                   with_ties = F) %>%  # Select the row with the highest Score.after.pen in each group
  ungroup() %>% # Remove the grouping
  select(PBid, X3p5p_mod, Mature_mod, seed, nonseed, total_reads, seed.match) %>%
  mutate(human.homolog =  mapply(nonseed.blast, seed, nonseed, MoreArgs = list(ref.miR.df = hsa.ref))) %>%
  mutate(mGDB.miR = human.homolog) %>%
  left_join(hsa.ref %>% 
              dplyr::rename(ref.seed = seed) %>%
              select(mGDB.miR, mGDB.miR.full, mGDB.fam, mB.miR, mGDB.3p5p, ref.seed, species), by = "mGDB.miR") %>%
  dplyr::rename(hmi.seed = seed,
                naming.species = species) %>%
  filter(hmi.seed == ref.seed)

unique.hsa <- rbind(unique28.hsa)

match.hsa <- unique.hsa %>% 
  select(PBid, naming.species, mGDB.miR, mGDB.miR.full, mGDB.fam, mB.miR, seed.match) %>%
  mutate(
    mGDB.trim = sub("^(.{4})(.{3}-\\d+).*", "\\2", mGDB.miR.full),
    mGDB.trim = casefold(mGDB.trim, upper = FALSE),
    mB.trim   = sub("^(.{4})(.{3}-\\d+).*", "\\2", mB.miR),
    mGDB_match_miRBase = ifelse(
      !is.na(mGDB.trim) & !is.na(mB.trim) & mGDB.trim == mB.trim,
      TRUE,
      FALSE
    )) %>%
  mutate(
    mGDB.fam = if_else(
      !is.na(mGDB.fam), 
      mGDB.fam, 
      tolower(
        if_else(
          str_detect(mGDB.miR, regex("novel", ignore_case = TRUE)),
          str_remove(mGDB.miR, "(-P|-o|-a).*"),  # Removes only if -P, -o, or -a appears as a suffix
          str_sub(str_remove(mGDB.miR, "(-P|-o|-a|-v).*"), 5)  )))) %>%
  left_join(revise.df %>% select(PBid, total_reads, total_mature, total_star, X3p5p_mod, 
                                 Mature_mod, Star_mod, Hairpin_seq_trimmed_ziv,  Chr, Start, End), by = "PBid") 

nrow(match.hsa) # 38
match.hsa %>% pull(mGDB.fam) %>% unique() # 38 miR, 22 Families

# 3-1 Hsa: process mir-125 seperately ----
# this match.hsa.0 will be used for later
match.hsa.0 <- match.hsa %>% filter(mGDB.fam == "mir-10") %>%
  arrange(PBid) %>%
  mutate(mGDB.fam = c("mir-10","mir-125","mir-125"), 
         mGDB.trim = c("mir-10", "mir-125","mir-125")) 
 
match.hsa.1 <- match.hsa %>% filter(mGDB.fam != "mir-10") %>%
  arrange(desc(seed.match), mGDB.fam) 

# 4 Model 7mer ----
#prep the references
cel.mGDB.df <- mGDB.prep("Cel", "reference1/MirGeneDB/mGDB_3/mature/cel.fas") 
mmu.mGDB.df <- mGDB.prep("Mmu", "reference1/MirGeneDB/mGDB_3/mature/mmu.fas") 
dme.mGDB.df <- mGDB.prep("Dme", "reference1/MirGeneDB/mGDB_3/mature/dme.fas") 
dre.mGDB.df <- mGDB.prep("Dre", "reference1/MirGeneDB/mGDB_3/mature/dre.fas") 

model.mGDB.df <- rbind(cel.mGDB.df, mmu.mGDB.df, dme.mGDB.df, dre.mGDB.df)

model.ref <- model.mGDB.df %>% 
  group_by(mGDB.miR.full) %>%  # Group by ivl.pre.id
  dplyr::slice_max(order_by = mGDB.fam, n = 1,
                   with_ties = F) %>% 
  ungroup()  %>%
  as.data.frame() %>%
  group_by(mGDB.fam, seed) %>%
  mutate(seed_count = n()) %>%
  ungroup() 

seed28.model <- hmia.df %>% 
  filter(PBid %notin% match.hsa$PBid) %>% 
  # filter(PBid %notin% match.xacl$PBid) %>%
  left_join(model.ref, by = "seed") %>%
  filter(is.na(mGDB.miR) == F ) %>%
  mutate(seed.match = "7mer")

unique28.model <- seed28.model %>%
  group_by(PBid) %>%  # Group by ivl.pre.id
  dplyr::slice_max(order_by = total_reads, n = 1,
                   with_ties = F) %>%  # Select the row with the highest read in each group
  ungroup() %>% # Remove the grouping
  select(PBid, Mature_mod, seed, nonseed, total_reads, seed.match) %>%
  mutate(homolog =  mapply(nonseed.blast, seed, nonseed, MoreArgs = list(ref.miR.df = model.ref))) %>%
  mutate(mGDB.miR = homolog)  %>%
  left_join(model.ref %>% select(mGDB.miR.full, mGDB.miR, mB.miR, mGDB.3p5p, mGDB.fam, species), by = "mGDB.miR") %>%
  dplyr::rename(hmi.seed = seed,
                naming.species = species) %>% 
  mutate(ref.seed = hmi.seed) %>%
  group_by(PBid) %>%  # Group by ivl.pre.id
  dplyr::slice_max(order_by = total_reads, n = 1,
                   with_ties = F) %>%  # Select the row with the highest Score.after.pen in each group
  ungroup()

unique.model <- rbind(unique28.model)

match.model <- unique.model %>% 
  select(PBid, naming.species, mGDB.miR, mGDB.miR.full, mGDB.fam, mB.miR, seed.match) %>%
  mutate(
    mGDB.trim = sub("^(.{4})(.{3}-\\d+).*", "\\2", mGDB.miR.full),
    mGDB.trim = casefold(mGDB.trim, upper = FALSE),
    mB.trim   = sub("^(.{4})(.{3}-\\d+).*", "\\2", mB.miR),
    mGDB_match_miRBase = ifelse(
      !is.na(mGDB.trim) & !is.na(mB.trim) & mGDB.trim == mB.trim,
      TRUE,
      FALSE
    )) %>%
  mutate(
    mGDB.fam = if_else(
      !is.na(mGDB.fam), 
      mGDB.fam, 
      tolower(
        if_else(
          str_detect(mGDB.miR, regex("novel", ignore_case = TRUE)),
          str_remove(mGDB.miR, "(-P|-o|-a).*"),  # Removes only if -P, -o, or -a appears as a suffix
          str_sub(str_remove(mGDB.miR, "(-P|-o|-a|-v).*"), 5)  )))) %>%
  left_join(revise.df %>% select(PBid, total_reads, total_mature, total_star, X3p5p_mod, 
                                 Mature_mod, Star_mod, Hairpin_seq_trimmed_ziv,  Chr, Start, End), by = "PBid") 

match.model.1 <- match.model %>% 
  arrange(mGDB.fam) 

# nrow(match.model.1) #19
# match.model.1 %>% pull(mGDB.fam) %>% unique() # 19 miR, 9 Families


# 5 Other Bilateria 7mer ----
source("C:/Users/duany/OneDrive/Desktop/Srivastava Lab/sRNA-seq/PB_sRNA/script/Head_02_Ref_ziv.R")

all.mGDB.df0 <- mGDB.prep("All", "reference1/MirGeneDB/mGDB_3/mature/all.fas") %>% 
  mutate(species = substr(mGDB.miR, 1, 3)) %>%
  # filter(species %notin% non_bilateria) %>%
  mutate(
    mGDB.fam = if_else(
      !is.na(mGDB.fam), 
      mGDB.fam, 
      tolower(
        if_else(
          str_detect(mGDB.miR, regex("novel", ignore_case = TRUE)),
          str_remove(mGDB.miR, "(-P|-o|-a|-v).*"),  # Removes only if -P, -o, -v or -a appears as a suffix
          str_sub(str_remove(mGDB.miR, "(-P|-o|-a|-v).*"), 5)  ))))

mGDB_all0 <- read.csv("Isana.2023.11/mGDB_30/mGDB.all.miR.html.csv", header = T)
colnames(mGDB_all0) <- c("MirGeneDB_ID", "MiRBase_ID", "Family", "Seed", "Accession_5p", "Accession_3p",
                         "Chromosome", "Start", "End", "Strand", "Node_of_origin_locus",
                         "Node_of_origin_family", "Iso.miR", "NTU_3", "UG", "UGUG", "CNNC")
mGDB_all <- mGDB_all0 %>%  
  mutate(species = substr(MirGeneDB_ID, 1, 3)) %>%
  filter(species %notin% non_bilateria) 

other.mGDB.df <- all.mGDB.df0 %>%
  filter(species %notin% c("Hsa","Mmu","Cel","Dme","Dre", "Xbo", "Sro", "Hmi")) %>%
  filter(species %notin% non_bilateria)

# other.mGDB.df %>% pull(species) %>% unique() %>% length() #100 species

other.ref <- other.mGDB.df %>% 
  group_by(mGDB.miR.full) %>%  # Group by PBid
  dplyr::slice_max(order_by = mGDB.fam, n = 1,
                   with_ties = F) %>% 
  ungroup() %>% 
  as.data.frame() %>%
  group_by(mGDB.fam) %>%
  mutate(mGDB_seed_count = n()) %>%
  ungroup() %>%
  group_by(seed) %>%
  mutate(distinct_seed_count = n()) %>%
  ungroup() %>%
  group_by(mGDB.fam) %>%
  filter(if (n_distinct(mGDB.3p5p) > 1) mGDB_seed_count == max(mGDB_seed_count) else TRUE) %>%
  ungroup() %>%# Select the row with the highest Score.after.pen in each group
  filter(mGDB.miR != "Ofu-Mir-210" & mGDB.miR != "Pca-Mir-210" & mGDB.miR != "Gsa-Mir-124-o7")


seed28.other <- hmia.df %>% 
  filter(PBid %notin% match.hsa$PBid) %>%
  # filter(PBid %notin% match.xacl.1$PBid) %>%
  filter(PBid %notin% match.model.1$PBid) %>%
  left_join(other.ref, by = "seed") %>% 
  filter(is.na(mGDB.miR) ==F ) %>%
  mutate(seed.match = "7mer")

unique28.other <- seed28.other %>%
  group_by(PBid) %>%  # Group by PBid
  dplyr::slice_max(order_by = total_reads, n = 1,
                   with_ties = F) %>%  # Select the row with the highest read in each group
  ungroup() %>% # Remove the grouping
  select(PBid, Mature_mod, seed, nonseed, total_reads, seed.match) %>%
  mutate(homolog =  mapply(nonseed.blast, seed, nonseed, MoreArgs = list(ref.miR.df = other.ref))) %>%
  mutate(mGDB.miR = homolog)  %>%
  left_join(other.ref %>% select(mGDB.miR.full, mGDB.miR, mB.miR, mGDB.3p5p, mGDB.fam, species), by = "mGDB.miR") %>%
  dplyr::rename(hmi.seed = seed,
                naming.species = species) %>% 
  mutate(ref.seed = hmi.seed) %>%
  group_by(PBid) %>%  # Group by ivl.pre.id
  dplyr::slice_max(order_by = total_reads, n = 1,
                   with_ties = F) %>%  # Select the row with the highest Score.after.pen in each group
  ungroup()

unique.other <- rbind(unique28.other)

match.other <- unique.other %>% 
  select(PBid, naming.species, mGDB.miR, mGDB.miR.full, mGDB.fam, mB.miR, seed.match) %>%
  mutate(
    mGDB.trim = sub("^(.{4})(.{3}-\\d+).*", "\\2", mGDB.miR.full),
    mGDB.trim = casefold(mGDB.trim, upper = FALSE),
    mB.trim   = sub("^(.{4})(.{3}-\\d+).*", "\\2", mB.miR),
    mGDB_match_miRBase = ifelse(
      !is.na(mGDB.trim) & !is.na(mB.trim) & mGDB.trim == mB.trim,
      TRUE,
      FALSE
    )) %>%
  mutate(
    mGDB.fam = if_else(
      !is.na(mGDB.fam), 
      mGDB.fam, 
      tolower(
        if_else(
          str_detect(mGDB.miR, regex("novel", ignore_case = TRUE)),
          str_remove(mGDB.miR, "(-P|-o|-a).*"),  # Removes only if -P, -o, or -a appears as a suffix
          str_sub(str_remove(mGDB.miR, "(-P|-o|-a|-v).*"), 5)  )))) %>%
  left_join(revise.df %>% select(PBid, total_reads, total_mature, total_star, X3p5p_mod, 
                                 Mature_mod, Star_mod, Hairpin_seq_trimmed_ziv,  Chr, Start, End), by = "PBid") 

match.other.1 <- match.other %>% 
  arrange(desc(seed.match), mGDB.fam) 

# nrow(match.other.1) # 62
# match.other.1 %>% pull(mGDB.fam) %>% unique() # 62 miRs, 24 Families


# 6 Xenacoelamorpha 7mer ----
#prep the references
sro.mGDB.df <- mGDB.prep("Sro", "reference1/MirGeneDB/mGDB_3/mature/sro.fas") 
xbo.mGDB.df <- mGDB.prep("Xbo", "reference1/MirGeneDB/mGDB_3/mature/xbo.fas") 

xacl.mGDB.df <- rbind(sro.mGDB.df, xbo.mGDB.df)

xacl.ref <- xacl.mGDB.df %>% 
  group_by(mGDB.miR.full) %>%  # Group by PBid
  dplyr::slice_max(order_by = mGDB.fam, n = 1,
                   with_ties = F) %>% 
  ungroup()  %>%
  as.data.frame() %>%
  group_by(mGDB.fam, seed) %>%
  mutate(seed_count = n()) %>%
  ungroup() 


seed28.xacl <- hmia.df %>% 
  filter(PBid %notin% match.hsa$PBid) %>%
  filter(PBid %notin% match.model.1$PBid) %>%
  filter(PBid %notin% match.other.1$PBid) %>%
  left_join(xacl.ref, by = "seed") %>%
  filter(is.na(mGDB.miR) == F ) %>%
  mutate(seed.match = "7mer",
         mGDB.miR = as.character(mGDB.miR))

unique28.xacl <- seed28.xacl %>%
  group_by(PBid) %>%  # Group by PBid
  dplyr::slice_max(order_by = total_reads, n = 1,
                   with_ties = F) %>%  # Select the row with the highest read in each group
  ungroup() %>% # Remove the grouping
  select(PBid, Mature_mod, seed, nonseed, total_reads, seed.match) %>%
  mutate(homolog =  mapply(nonseed.blast, seed, nonseed, MoreArgs = list(ref.miR.df = xacl.ref))) %>%
  mutate(mGDB.miR = homolog)  %>%
  left_join(xacl.ref %>% select(mGDB.miR.full, mGDB.miR, mB.miR, mGDB.3p5p, mGDB.fam, species), by = "mGDB.miR") %>%
  dplyr::rename(hmi.seed = seed,
                naming.species = species) %>% 
  mutate(ref.seed = hmi.seed) %>%
  group_by(PBid) %>%  # Group by ivl.pre.id
  dplyr::slice_max(order_by = total_reads, n = 1,
                   with_ties = F) %>%  # Select the row with the highest Score.after.pen in each group
  ungroup()

unique.xacl <- rbind(unique28.xacl)

match.xacl <- unique.xacl %>% 
  select(PBid, naming.species, mGDB.miR, mGDB.miR.full, mGDB.fam, mB.miR, seed.match) %>%
  mutate(
    mGDB.trim = sub("^(.{4})(.{3}-\\d+).*", "\\2", mGDB.miR.full),
    mGDB.trim = casefold(mGDB.trim, upper = FALSE),
    mB.trim   = sub("^(.{4})(.{3}-\\d+).*", "\\2", mB.miR),
    mGDB_match_miRBase = ifelse(
      !is.na(mGDB.trim) & !is.na(mB.trim) & mGDB.trim == mB.trim,
      TRUE,
      FALSE
    )) %>%
  mutate(
    mGDB.fam = if_else(
      !is.na(mGDB.fam), 
      mGDB.fam, 
      tolower(
        if_else(
          str_detect(mGDB.miR, regex("novel", ignore_case = TRUE)),
          str_remove(mGDB.miR, "(-P|-o|-a).*"),  # Removes only if -P, -o, or -a appears as a suffix
          str_sub(str_remove(mGDB.miR, "(-P|-o|-a|-v).*"), 5)  )))) %>%
  left_join(revise.df %>% select(PBid, total_reads, total_mature, total_star, X3p5p_mod, 
                                 Mature_mod, Star_mod, Hairpin_seq_trimmed_ziv,  Chr, Start, End), by = "PBid") 


match.xacl.1 <- match.xacl %>% 
  arrange(desc(seed.match), mGDB.fam) 

# nrow(match.xacl.1) # 10
# match.xacl.1 %>% pull(mGDB.fam) %>% unique() # 10 miR, 5 Families




# 7 non-bilateria 7mer ----
mGDB_all <- mGDB_all0 %>%  
  mutate(species = substr(MirGeneDB_ID, 1, 3)) 

non_bilateria.mGDB.df <- all.mGDB.df0 %>%
  filter(species %notin% c("Hsa","Mmu","Cel","Dme","Dre", "Xbo", "Sro", "Hmi")) %>%
  filter(species %in% non_bilateria)

# non_bilateria.mGDB.df %>% pull(species) %>% unique()  # 6 species

non_bilateria.ref <- non_bilateria.mGDB.df %>% 
  group_by(mGDB.miR.full) %>%  # Group by PBid
  dplyr::slice_max(order_by = mGDB.fam, n = 1,
                   with_ties = F) %>% 
  ungroup() %>% 
  as.data.frame() %>%
  group_by(mGDB.fam) %>%
  mutate(mGDB_seed_count = n()) %>%
  ungroup() %>%
  group_by(seed) %>%
  mutate(distinct_seed_count = n()) %>%
  ungroup() %>%
  group_by(mGDB.fam) %>%
  filter(if (n_distinct(mGDB.3p5p) > 1) mGDB_seed_count == max(mGDB_seed_count) else TRUE) %>%
  ungroup() %>%# Select the row with the highest Score.after.pen in each group
  filter(mGDB.miR != "Ofu-Mir-210" & mGDB.miR != "Pca-Mir-210")


seed28.non_bilateria <- hmia.df %>% 
  filter(PBid %notin% match.hsa$PBid) %>%
  filter(PBid %notin% match.xacl.1$PBid) %>%
  filter(PBid %notin% match.other.1$PBid) %>%
  filter(PBid %notin% match.model.1$PBid) %>%
  left_join(non_bilateria.ref, by = "seed") %>% 
  filter(is.na(mGDB.miR) ==F ) %>%
  mutate(seed.match = "7mer")

unique28.non_bilateria <- seed28.non_bilateria %>%
  group_by(PBid) %>%  # Group by PBid
  dplyr::slice_max(order_by = total_reads, n = 1,
                   with_ties = F) %>%  # Select the row with the highest read in each group
  ungroup() %>% # Remove the grouping
  select(PBid, Mature_mod, seed, nonseed, total_reads, seed.match) %>%
  mutate(homolog =  mapply(nonseed.blast, seed, nonseed, MoreArgs = list(ref.miR.df = non_bilateria.ref))) %>%
  mutate(mGDB.miR = homolog)  %>%
  left_join(non_bilateria.ref %>% select(mGDB.miR.full, mGDB.miR, mB.miR, mGDB.3p5p, mGDB.fam, species), by = "mGDB.miR") %>%
  dplyr::rename(hmi.seed = seed,
                naming.species = species) %>% 
  mutate(ref.seed = hmi.seed) %>%
  group_by(PBid) %>%  # Group by ivl.pre.id
  dplyr::slice_max(order_by = total_reads, n = 1,
                   with_ties = F) %>%  # Select the row with the highest Score.after.pen in each group
  ungroup()

unique.non_bilateria <- rbind(unique28.non_bilateria)

match.non_bilateria <- unique.non_bilateria %>% 
  select(PBid, naming.species, mGDB.miR, mGDB.miR.full, mGDB.fam, mB.miR, seed.match) %>%
  mutate(
    mGDB.trim = sub("^(.{4})(.{3}-\\d+).*", "\\2", mGDB.miR.full),
    mGDB.trim = casefold(mGDB.trim, upper = FALSE),
    mB.trim   = sub("^(.{4})(.{3}-\\d+).*", "\\2", mB.miR),
    mGDB_match_miRBase = ifelse(
      !is.na(mGDB.trim) & !is.na(mB.trim) & mGDB.trim == mB.trim,
      TRUE,
      FALSE
    )) %>%
  mutate(
    mGDB.fam = if_else(
      !is.na(mGDB.fam), 
      mGDB.fam, 
      tolower(
        if_else(
          str_detect(mGDB.miR, regex("novel", ignore_case = TRUE)),
          str_remove(mGDB.miR, "(-P|-o|-a).*"),  # Removes only if -P, -o, or -a appears as a suffix
          str_sub(str_remove(mGDB.miR, "(-P|-o|-a|-v).*"), 5)  )))) %>%
  left_join(revise.df %>% select(PBid, total_reads, total_mature, total_star, X3p5p_mod, 
                                 Mature_mod, Star_mod, Hairpin_seq_trimmed_ziv,  Chr, Start, End), by = "PBid") 

match.non_bilateria.1 <- match.non_bilateria %>% 
  arrange(desc(seed.match), mGDB.fam) 

# nrow(match.non_bilateria.1) # 3
# match.non_bilateria.1 %>% pull(mGDB.fam) %>% unique() # 3 miRs, 2 families


# 8 Combine All 7mer seed ----
all.7mer <- rbind(match.hsa.1, match.model.1, match.xacl.1, match.other.1, match.non_bilateria.1) %>%
  mutate(mGDB.fam = ifelse(is.na(mGDB.fam), mGDB.trim, mGDB.fam))

no_name.df7 <- hmia.df %>%
  filter(PBid %notin% as.vector(all.7mer$PBid) &
           PBid %notin% as.vector(match.hsa.0$PBid))

# 9-0 Isolate un-assigned to check strand ---- 
no_name.df8_process <- no_name.df7 %>%
  dplyr::mutate(
    total_mature_new = if_else(total_mature < total_star, total_star, total_mature),
    total_star_new   = if_else(total_mature < total_star, total_mature, total_star)
  ) %>%
  dplyr::select(-total_mature, -total_star) %>%
  dplyr::rename(
    total_mature = total_mature_new,
    total_star   = total_star_new
  ) %>% 
  mutate(m_length = nchar(Mature_mod), 
         s_length = nchar(Star_mod)) %>%
  mutate(m_perBaseCov = total_mature / m_length,
         s_perBaseCov = total_star / s_length) %>%
  mutate(ms_ratio = abs(total_mature/total_star))
  
no_strand_preference.id <- no_name.df8_process %>%
  filter(ms_ratio <= 2) %>%
  filter(s_perBaseCov >= 1) %>%
  filter(m_perBaseCov >= 5) %>%
  pull(PBid)

hmia.df8 <- hmia.df %>% filter(PBid %in% no_strand_preference.id) %>% 
  mutate(seed = seed_alt)  %>%
  select(-seed_alt)

# 9-1 seed_alt Hsa ----
seed_alt.hsa <- hmia.df8 %>% 
  # dplyr::filter(PBid %notin% match.xacl.2$PBid) %>%
  dplyr::left_join(hsa.ref, by = "seed") %>% 
  dplyr::filter(!is.na(mGDB.miR)) %>%
  dplyr::mutate(seed.match = "7mer_alt")

if (nrow(seed_alt.hsa) == 0) {
  
  message("seed_alt.hsa is EMPTY — skipping unique_alt.hsa and match.hsa.2.")
  
  ## create empty match.hsa.2 with the columns you will rely on later
  match.hsa.2 <- tibble::tibble(
    PBid               = character(),
    naming.species     = character(),
    mGDB.miR           = character(),
    mGDB.miR.full      = character(),
    mGDB.fam           = character(),
    mB.miR             = character(),
    seed.match         = character(),
    mGDB.trim          = character(),
    mB.trim            = character(),
    mGDB_match_miRBase = logical(),
    total_reads        = numeric(),
    total_mature       = numeric(),
    total_star         = numeric(),
    X3p5p_mod          = character(),
    Mature_mod         = character(),
    Star_mod           = character(),
    Hairpin_seq_trimmed_ziv = character(),
    Chr                = character(),
    Start              = numeric(),
    End                = numeric()
  )
  
} else {
  unique_alt.hsa <- seed_alt.hsa %>%
    group_by(PBid) %>%  # Group by PBid
    dplyr::slice_max(order_by = total_reads, n = 1,
                     with_ties = F) %>%  # Select the row with the highest read in each group
    ungroup() %>% # Remove the grouping
    select(PBid, Mature_mod, seed, nonseed, total_reads, seed.match) %>%
    mutate(homolog =  mapply(nonseed.blast, seed, nonseed, MoreArgs = list(ref.miR.df = hsa.ref))) %>%
    mutate(mGDB.miR = homolog)  %>%
    left_join(hsa.ref %>% select(mGDB.miR.full, mGDB.miR, mB.miR, mGDB.3p5p, mGDB.fam, species), by = "mGDB.miR") %>%
    dplyr::rename(hmi.seed = seed,
                  naming.species = species) %>% 
    mutate(ref.seed = hmi.seed) %>%
    group_by(PBid) %>%  # Group by ivl.pre.id
    dplyr::slice_max(order_by = total_reads, n = 1,
                     with_ties = F) %>%  # Select the row with the highest Score.after.pen in each group
    ungroup()
  
  match.hsa.2 <- unique_alt.hsa %>% 
    select(PBid, naming.species, mGDB.miR, mGDB.miR.full, mGDB.fam, mB.miR, seed.match) %>%
    mutate(
      mGDB.trim = sub("^(.{4})(.{3}-\\d+).*", "\\2", mGDB.miR.full),
      mGDB.trim = casefold(mGDB.trim, upper = FALSE),
      mB.trim   = sub("^(.{4})(.{3}-\\d+).*", "\\2", mB.miR),
      mGDB_match_miRBase = ifelse(
        !is.na(mGDB.trim) & !is.na(mB.trim) & mGDB.trim == mB.trim,
        TRUE,
        FALSE
      )) %>%
    mutate(
      mGDB.fam = if_else(
        !is.na(mGDB.fam), 
        mGDB.fam, 
        tolower(
          if_else(
            str_detect(mGDB.miR, regex("novel", ignore_case = TRUE)),
            str_remove(mGDB.miR, "(-P|-o|-a).*"),  # Removes only if -P, -o, or -a appears as a suffix
            str_sub(str_remove(mGDB.miR, "(-P|-o|-a|-v).*"), 5)  )))) %>%
    left_join(revise.df %>% select(PBid, total_reads, total_mature, total_star, X3p5p_mod, 
                                   Mature_mod, Star_mod, Hairpin_seq_trimmed_ziv,  Chr, Start, End), by = "PBid") 
  
  
  match.hsa.2 <- match.hsa.2 %>% 
    arrange(desc(seed.match), mGDB.fam) 
  
}

# nrow(match.hsa.2) # 3
# match.hsa.2 %>% pull(mGDB.fam) %>% unique() # 3 miR, 2 Families


# 9.2 seed_alt Models ----
seed_alt.model <- hmia.df8 %>% 
  dplyr::filter(PBid %notin% match.hsa.2$PBid) %>%
  dplyr::left_join(model.ref, by = "seed") %>% 
  dplyr::filter(!is.na(mGDB.miR)) %>%
  dplyr::mutate(seed.match = "7mer_alt")

if (nrow(seed_alt.model) == 0) {
  
  message("seed_alt.model is EMPTY — skipping unique_alt.model and match.model.2.")
  
  ## create empty match.model.2 with the columns you will rely on later
  match.model.2 <- tibble::tibble(
    PBid               = character(),
    naming.species     = character(),
    mGDB.miR           = character(),
    mGDB.miR.full      = character(),
    mGDB.fam           = character(),
    mB.miR             = character(),
    seed.match         = character(),
    mGDB.trim          = character(),
    mB.trim            = character(),
    mGDB_match_miRBase = logical(),
    total_reads        = numeric(),
    total_mature       = numeric(),
    total_star         = numeric(),
    X3p5p_mod          = character(),
    Mature_mod         = character(),
    Star_mod           = character(),
    Hairpin_seq_trimmed_ziv = character(),
    Chr                = character(),
    Start              = numeric(),
    End                = numeric()
  )
  
} else {
  
  unique_alt.model <- seed_alt.model %>%
    group_by(PBid) %>%  # Group by PBid
    dplyr::slice_max(order_by = total_reads, n = 1,
                     with_ties = F) %>%  # Select the row with the highest read in each group
    ungroup() %>% # Remove the grouping
    select(PBid, Mature_mod, seed, nonseed, total_reads, seed.match) %>%
    mutate(homolog =  mapply(nonseed.blast, seed, nonseed, MoreArgs = list(ref.miR.df = model.ref))) %>%
    mutate(mGDB.miR = homolog)  %>%
    left_join(model.ref %>% select(mGDB.miR.full, mGDB.miR, mB.miR, mGDB.3p5p, mGDB.fam, species), by = "mGDB.miR") %>%
    dplyr::rename(hmi.seed = seed,
                  naming.species = species) %>% 
    mutate(ref.seed = hmi.seed) %>%
    group_by(PBid) %>%  # Group by ivl.pre.id
    dplyr::slice_max(order_by = total_reads, n = 1,
                     with_ties = F) %>%  # Select the row with the highest Score.after.pen in each group
    ungroup()
  
  match.model.2 <- unique_alt.model %>% 
    select(PBid, naming.species, mGDB.miR, mGDB.miR.full, mGDB.fam, mB.miR, seed.match) %>%
    mutate(
      mGDB.trim = sub("^(.{4})(.{3}-\\d+).*", "\\2", mGDB.miR.full),
      mGDB.trim = casefold(mGDB.trim, upper = FALSE),
      mB.trim   = sub("^(.{4})(.{3}-\\d+).*", "\\2", mB.miR),
      mGDB_match_miRBase = ifelse(
        !is.na(mGDB.trim) & !is.na(mB.trim) & mGDB.trim == mB.trim,
        TRUE,
        FALSE
      )) %>%
    mutate(
      mGDB.fam = if_else(
        !is.na(mGDB.fam), 
        mGDB.fam, 
        tolower(
          if_else(
            str_detect(mGDB.miR, regex("novel", ignore_case = TRUE)),
            str_remove(mGDB.miR, "(-P|-o|-a).*"),  # Removes only if -P, -o, or -a appears as a suffix
            str_sub(str_remove(mGDB.miR, "(-P|-o|-a|-v).*"), 5)  )))) %>%
    left_join(revise.df %>% select(PBid, total_reads, total_mature, total_star, X3p5p_mod, 
                                   Mature_mod, Star_mod, Hairpin_seq_trimmed_ziv,  Chr, Start, End), by = "PBid") 
  
  
  match.model.2 <- match.model.2 %>% 
    arrange(desc(seed.match), mGDB.fam) 
}

# nrow(match.model.2) #1
# match.model.2 %>% pull(mGDB.fam) %>% unique() # 1 miR, 1 Families

# 9.3 seed_alt Other ----
seed_alt.other <- hmia.df8 %>% 
  # dplyr::filter(PBid %notin% match.xacl.2$PBid) %>%
  dplyr::filter(PBid %notin% match.hsa.2$PBid) %>%
  dplyr::filter(PBid %notin% match.model.2$PBid) %>%
  dplyr::left_join(other.ref, by = "seed") %>% 
  dplyr::filter(!is.na(mGDB.miR)) %>%
  dplyr::mutate(seed.match = "7mer_alt")

if (nrow(seed_alt.other) == 0) {
  
  message("seed_alt.other is EMPTY — skipping unique_alt.other and match.other.2.")
  
  ## create empty match.other.2 with the columns you will rely on later
  match.other.2 <- tibble::tibble(
    PBid               = character(),
    naming.species     = character(),
    mGDB.miR           = character(),
    mGDB.miR.full      = character(),
    mGDB.fam           = character(),
    mB.miR             = character(),
    seed.match         = character(),
    mGDB.trim          = character(),
    mB.trim            = character(),
    mGDB_match_miRBase = logical(),
    total_reads        = numeric(),
    total_mature       = numeric(),
    total_star         = numeric(),
    X3p5p_mod          = character(),
    Mature_mod         = character(),
    Star_mod           = character(),
    Hairpin_seq_trimmed_ziv = character(),
    Chr                = character(),
    Start              = numeric(),
    End                = numeric()
  )
  
} else {
  
  unique_alt.other <- seed_alt.other %>%
    group_by(PBid) %>%  # Group by PBid
    dplyr::slice_max(order_by = total_reads, n = 1,
                     with_ties = F) %>%  # Select the row with the highest read in each group
    ungroup() %>% # Remove the grouping
    select(PBid, Mature_mod, seed, nonseed, total_reads, seed.match) %>%
    mutate(homolog =  mapply(nonseed.blast, seed, nonseed, MoreArgs = list(ref.miR.df = other.ref))) %>%
    mutate(mGDB.miR = homolog)  %>%
    left_join(other.ref %>% select(mGDB.miR.full, mGDB.miR, mB.miR, mGDB.3p5p, mGDB.fam, species), by = "mGDB.miR") %>%
    dplyr::rename(hmi.seed = seed,
                  naming.species = species) %>% 
    mutate(ref.seed = hmi.seed) %>%
    group_by(PBid) %>%  # Group by ivl.pre.id
    dplyr::slice_max(order_by = total_reads, n = 1,
                     with_ties = F) %>%  # Select the row with the highest Score.after.pen in each group
    ungroup()
  
  match.other.2 <- unique_alt.other %>% 
    select(PBid, naming.species, mGDB.miR, mGDB.miR.full, mGDB.fam, mB.miR, seed.match) %>%
    mutate(
      mGDB.trim = sub("^(.{4})(.{3}-\\d+).*", "\\2", mGDB.miR.full),
      mGDB.trim = casefold(mGDB.trim, upper = FALSE),
      mB.trim   = sub("^(.{4})(.{3}-\\d+).*", "\\2", mB.miR),
      mGDB_match_miRBase = ifelse(
        !is.na(mGDB.trim) & !is.na(mB.trim) & mGDB.trim == mB.trim,
        TRUE,
        FALSE
      )) %>%
    mutate(
      mGDB.fam = if_else(
        !is.na(mGDB.fam), 
        mGDB.fam, 
        tolower(
          if_else(
            str_detect(mGDB.miR, regex("novel", ignore_case = TRUE)),
            str_remove(mGDB.miR, "(-P|-o|-a).*"),  # Removes only if -P, -o, or -a appears as a suffix
            str_sub(str_remove(mGDB.miR, "(-P|-o|-a|-v).*"), 5)  )))) %>%
    left_join(revise.df %>% select(PBid, total_reads, total_mature, total_star, X3p5p_mod, 
                                   Mature_mod, Star_mod, Hairpin_seq_trimmed_ziv,  Chr, Start, End), by = "PBid") 
  
  
  match.other.2 <- match.other.2 %>% 
    arrange(desc(seed.match), mGDB.fam) 
}

# nrow(match.other.2) #21
# match.other.2 %>% pull(mGDB.fam) %>% unique() # 21 miR, 7 Families


# 9-4 seed_alt Xacl----
seed_alt.xacl <- hmia.df8 %>% 
  dplyr::filter(PBid %notin% match.hsa.2$PBid) %>%
  dplyr::filter(PBid %notin% match.model.2$PBid) %>%
  dplyr::filter(PBid %notin% match.other.2$PBid) %>%
  dplyr::left_join(xacl.ref, by = "seed") %>% 
  dplyr::filter(!is.na(mGDB.miR)) %>%
  dplyr::mutate(seed.match = "7mer_alt")

if (nrow(seed_alt.xacl) == 0) {
  
  message("seed_alt.xacl is EMPTY — skipping unique_alt.xacl and match.xacl.2.")
  
  ## create empty match.xacl.2 with the columns you will rely on later
  match.xacl.2 <- tibble::tibble(
    PBid               = character(),
    naming.species     = character(),
    mGDB.miR           = character(),
    mGDB.miR.full      = character(),
    mGDB.fam           = character(),
    mB.miR             = character(),
    seed.match         = character(),
    mGDB.trim          = character(),
    mB.trim            = character(),
    mGDB_match_miRBase = logical(),
    total_reads        = numeric(),
    total_mature       = numeric(),
    total_star         = numeric(),
    X3p5p_mod          = character(),
    Mature_mod         = character(),
    Star_mod           = character(),
    Hairpin_seq_trimmed_ziv = character(),
    Chr                = character(),
    Start              = numeric(),
    End                = numeric()
  )
  
} else {
  
  unique_alt.xacl <- seed_alt.xacl %>%
    group_by(PBid) %>%  # Group by PBid
    dplyr::slice_max(order_by = total_reads, n = 1,
                     with_ties = F) %>%  # Select the row with the highest read in each group
    ungroup() %>% # Remove the grouping
    select(PBid, Mature_mod, seed, nonseed, total_reads, seed.match) %>%
    mutate(homolog =  mapply(nonseed.blast, seed, nonseed, MoreArgs = list(ref.miR.df = xacl.ref))) %>%
    mutate(mGDB.miR = homolog)  %>%
    left_join(xacl.ref %>% select(mGDB.miR.full, mGDB.miR, mB.miR, mGDB.3p5p, mGDB.fam, species), by = "mGDB.miR") %>%
    dplyr::rename(hmi.seed = seed,
                  naming.species = species) %>% 
    mutate(ref.seed = hmi.seed) %>%
    group_by(PBid) %>%  # Group by ivl.pre.id
    dplyr::slice_max(order_by = total_reads, n = 1,
                     with_ties = F) %>%  # Select the row with the highest Score.after.pen in each group
    ungroup()
  
  match.xacl.2 <- unique_alt.xacl %>% 
    select(PBid, naming.species, mGDB.miR, mGDB.miR.full, mGDB.fam, mB.miR, seed.match) %>%
    mutate(
      mGDB.trim = sub("^(.{4})(.{3}-\\d+).*", "\\2", mGDB.miR.full),
      mGDB.trim = casefold(mGDB.trim, upper = FALSE),
      mB.trim   = sub("^(.{4})(.{3}-\\d+).*", "\\2", mB.miR),
      mGDB_match_miRBase = ifelse(
        !is.na(mGDB.trim) & !is.na(mB.trim) & mGDB.trim == mB.trim,
        TRUE,
        FALSE
      )) %>%
    mutate(
      mGDB.fam = if_else(
        !is.na(mGDB.fam), 
        mGDB.fam, 
        tolower(
          if_else(
            str_detect(mGDB.miR, regex("novel", ignore_case = TRUE)),
            str_remove(mGDB.miR, "(-P|-o|-a).*"),  # Removes only if -P, -o, or -a appears as a suffix
            str_sub(str_remove(mGDB.miR, "(-P|-o|-a|-v).*"), 5)  )))) %>%
    left_join(revise.df %>% select(PBid, total_reads, total_mature, total_star, X3p5p_mod, 
                                   Mature_mod, Star_mod, Hairpin_seq_trimmed_ziv,  Chr, Start, End), by = "PBid") 
  
  
  match.xacl.2 <- match.xacl.2 %>% 
    arrange(desc(seed.match), mGDB.fam) }
  
# nrow(match.xacl.2) # 0
# match.xacl.2 %>% pull(mGDB.fam) %>% unique() # 0 miR, 0 Families


# 9-5 seed_alt Xacl----
seed_alt.non_bilateria <- hmia.df8 %>% 
  dplyr::filter(PBid %notin% match.hsa.2$PBid) %>%
  dplyr::filter(PBid %notin% match.model.2$PBid) %>%
  dplyr::filter(PBid %notin% match.other.2$PBid) %>%
  dplyr::filter(PBid %notin% match.xacl.2$PBid) %>%
  dplyr::left_join(non_bilateria.ref, by = "seed") %>% 
  dplyr::filter(!is.na(mGDB.miR)) %>%
  dplyr::mutate(seed.match = "7mer_alt")

if (nrow(seed_alt.non_bilateria) == 0) {
  
  message("seed_alt.non_bilateria is EMPTY — skipping unique_alt.non_bilateria and match.non_bilateria.2.")
  
  ## create empty match.non_bilateria.2 with the columns you will rely on later
  match.non_bilateria.2 <- tibble::tibble(
    PBid               = character(),
    naming.species     = character(),
    mGDB.miR           = character(),
    mGDB.miR.full      = character(),
    mGDB.fam           = character(),
    mB.miR             = character(),
    seed.match         = character(),
    mGDB.trim          = character(),
    mB.trim            = character(),
    mGDB_match_miRBase = logical(),
    total_reads        = numeric(),
    total_mature       = numeric(),
    total_star         = numeric(),
    X3p5p_mod          = character(),
    Mature_mod         = character(),
    Star_mod           = character(),
    Hairpin_seq_trimmed_ziv = character(),
    Chr                = character(),
    Start              = numeric(),
    End                = numeric()
  )
  
} else {
  
  unique_alt.non_bilateria <- seed_alt.non_bilateria %>%
    group_by(PBid) %>%  # Group by PBid
    dplyr::slice_max(order_by = total_reads, n = 1,
                     with_ties = F) %>%  # Select the row with the highest read in each group
    ungroup() %>% # Remove the grouping
    select(PBid, Mature_mod, seed, nonseed, total_reads, seed.match) %>%
    mutate(homolog =  mapply(nonseed.blast, seed, nonseed, MoreArgs = list(ref.miR.df = non_bilateria.ref))) %>%
    mutate(mGDB.miR = homolog)  %>%
    left_join(non_bilateria.ref %>% select(mGDB.miR.full, mGDB.miR, mB.miR, mGDB.3p5p, mGDB.fam, species), by = "mGDB.miR") %>%
    dplyr::rename(hmi.seed = seed,
                  naming.species = species) %>% 
    mutate(ref.seed = hmi.seed) %>%
    group_by(PBid) %>%  # Group by ivl.pre.id
    dplyr::slice_max(order_by = total_reads, n = 1,
                     with_ties = F) %>%  # Select the row with the highest Score.after.pen in each group
    ungroup()
  
  match.non_bilateria.2 <- unique_alt.non_bilateria %>% 
    select(PBid, naming.species, mGDB.miR, mGDB.miR.full, mGDB.fam, mB.miR, seed.match) %>%
    mutate(
      mGDB.trim = sub("^(.{4})(.{3}-\\d+).*", "\\2", mGDB.miR.full),
      mGDB.trim = casefold(mGDB.trim, upper = FALSE),
      mB.trim   = sub("^(.{4})(.{3}-\\d+).*", "\\2", mB.miR),
      mGDB_match_miRBase = ifelse(
        !is.na(mGDB.trim) & !is.na(mB.trim) & mGDB.trim == mB.trim,
        TRUE,
        FALSE
      )) %>%
    mutate(
      mGDB.fam = if_else(
        !is.na(mGDB.fam), 
        mGDB.fam, 
        tolower(
          if_else(
            str_detect(mGDB.miR, regex("novel", ignore_case = TRUE)),
            str_remove(mGDB.miR, "(-P|-o|-a).*"),  # Removes only if -P, -o, or -a appears as a suffix
            str_sub(str_remove(mGDB.miR, "(-P|-o|-a|-v).*"), 5)  )))) %>%
    left_join(revise.df %>% select(PBid, total_reads, total_mature, total_star, X3p5p_mod, 
                                   Mature_mod, Star_mod, Hairpin_seq_trimmed_ziv,  Chr, Start, End), by = "PBid") 
  
  
  match.non_bilateria.2 <- match.non_bilateria.2 %>% 
    arrange(desc(seed.match), mGDB.fam) }

nrow(match.non_bilateria.2) # 0
match.non_bilateria.2 %>% pull(mGDB.fam) %>% unique() # 0 miR, 0 Families


# 10 Combine ref-named ----
all.7mer_alt <- rbind(match.hsa.2, match.model.2, match.xacl.2, match.other.2, match.non_bilateria.2) %>%
  mutate(mGDB.fam = ifelse(is.na(mGDB.fam), mGDB.trim, mGDB.fam))

# name.df0 <- rbind(all.7mer, all.7mer_alt, match.hsa.0)

# this is to manually add mir-32 to mir-92
name.df0 <- rbind(all.7mer, all.7mer_alt, match.hsa.0) %>%
  mutate(mGDB.fam = ifelse(PBid %in% c("PB_574", "PB_96"),
                           "mir-92", mGDB.fam))

name.df <- name.df0  %>%
  arrange(desc(seed.match), mGDB.fam, Mature_mod) %>%
  group_by(mGDB.fam) %>%
  mutate(fam.size = n()) %>%
  ungroup() %>%
  group_by(Mature_mod) %>%
  mutate(dup.num = n()) %>%
  ungroup() %>%
  group_by(Mature_mod) %>%
  mutate(suffix2 = as.character(row_number())) %>%
  ungroup() %>%
  mutate(suffix2 = ifelse(dup.num == 1, "", paste("-",suffix2, sep = ""))) 

# Convert the dataframe into a list of dataframes, each corresponding to a unique mGDB.fam
list_by_fam <- split(name.df, name.df$mGDB.fam)

list_by_fam_modified <- lapply(list_by_fam, function(df) {
  # Create a mapping from unique Mature values to letters
  unique_mature <- unique(df$Mature_mod)
  mature_to_letter <- setNames(letters[seq_along(unique_mature)], unique_mature)
  
  df %>%
    arrange(seed.match,desc(total_reads), Chr, Start) %>%
    group_by(Mature_mod) %>%
    mutate(suffix1 = mature_to_letter[Mature_mod]) %>%
    mutate(suffix1 = paste("-", suffix1, sep='')) %>%
    ungroup()
})

# View(list_by_fam_modified[["mir-34"]])

name.merge <- bind_rows(list_by_fam_modified) %>%
  mutate(suffix1 = ifelse(fam.size == 1 | dup.num == fam.size, "", suffix1)) %>%
  mutate(pre_miRNA = paste("hmi-", mGDB.fam, suffix1, suffix2, sep = "")) %>%
  select(-suffix1, -suffix2)

# This is the final ref_named
ref_named.df <- name.merge %>%
  mutate(mB.miR = ifelse(is.na(mB.miR), "None", mB.miR) ) %>%
  dplyr::rename(ref_mGDB = mGDB.miR.full,
                naming_species = naming.species,
                ref_miRBase_to_mGDB = mB.miR,
                seed_match_type = seed.match,
                family = mGDB.fam) %>%
  select(-mGDB.miR, -mGDB.trim, -mB.trim, -mGDB_match_miRBase, -ref_miRBase_to_mGDB)

# 11 Assign hmi-novel- ----
novel.df0 <- revise.df %>%
  filter(PBid %notin% as.vector(ref_named.df$PBid)) %>%
  mutate(seed6  = substr(Mature_mod, 2, 7),
         seed6_alt = substr(Star_mod, 2, 7)) %>%
  select(-X3p.5p_ziv, -X3p5p_swap, -Mature_ziv, -Star_ziv) %>%
  mutate(naming_species = "novel",
         ref_mGDB = "N/A")

novel_fam.df1 <- novel.df0 %>%
  group_by(seed6)  %>%
  summarize(family_reads = sum(total_reads),
            fam.size = n(),
            iso.num = n_distinct(Mature_mod),
            dup.num = fam.size - iso.num) %>%
  ungroup() %>%
  arrange(desc(family_reads))

novel.list <- novel_fam.df1 %>%
  mutate(suffix = paste("-", 1:nrow(novel_fam.df1), sep = "")) %>%
  mutate(family = paste("hmi-novel", suffix, sep = "")) %>%
  select(-family_reads, -suffix)

generate_suffix <- function(n) {
  base_letters <- letters
  if (n <= 26) {
    return(base_letters[1:n])
  } else {
    extra_letters <- unlist(lapply(1:ceiling(n / 26), function(x) paste0(letters[x], base_letters)))
    return(c(base_letters, extra_letters)[1:n])
  }
}

novel.df1 <- novel.df0 %>% 
  left_join(novel.list, by = "seed6") %>%
  dplyr::rename(prefix = family) %>%
  group_by(prefix) %>%
  arrange(desc(total_reads), .by_group = TRUE) %>%
  mutate(suffix1 = paste("-", generate_suffix(n()), sep = '')) %>%
  ungroup() %>%
  mutate(suffix1 = ifelse(iso.num == 1, "", suffix1)) %>%
  group_by(Mature_mod) %>%
  mutate(suffix2 = paste("-", row_number(), sep = ''),
         suffix1 = dplyr::first(suffix1)) %>%
  ungroup() %>%
  mutate(suffix2 = ifelse(dup.num == 0, "", suffix2)) %>%
  mutate(pre_miRNA = paste(prefix, suffix1, suffix2, sep = '')) 

# 11.1 Assign alt_novel ----
# no_strand_preference.id_str <- no_name.df8_process %>%
#   filter(ms_ratio <= 1.5) %>%
#   filter(s_perBaseCov >= 10) %>%
#   filter(m_perBaseCov >= 10) %>%
#   pull(PBid)

no_strand_preference.id_str <- no_name.df8_process %>%
  filter(ms_ratio <= 2) %>%
  filter(s_perBaseCov >= 1) %>%
  filter(m_perBaseCov >= 5) %>%
  pull(PBid)

novel_alt.list <- novel.df0 %>%
  filter(PBid %in% no_strand_preference.id_str) %>%
  pull(PBid)  # 47 miRs

novel_alt.df0 <-novel.df0 %>% 
  mutate(seed6 = ifelse(PBid %in% novel_alt.list, seed6_alt, seed6))

novel_fam_alt.df1 <- novel_alt.df0 %>%
  group_by(seed6)  %>%
  summarize(family_reads = sum(total_reads),
            fam.size = n(),
            iso.num = n_distinct(Mature_mod),
            dup.num = fam.size - iso.num) %>%
  ungroup() %>%
  arrange(desc(family_reads))

novel_alt.list <- novel_fam_alt.df1 %>%
  mutate(suffix = paste("-", 1:nrow(novel_fam_alt.df1), sep = "")) %>%
  mutate(family = paste("hmi-novel", suffix, sep = "")) %>%
  select(-family_reads, -suffix)

novel_alt.df1 <- novel_alt.df0 %>% 
  left_join(novel_alt.list, by = "seed6") %>%
  dplyr::rename(prefix = family) %>%
  group_by(prefix) %>%
  arrange(desc(total_reads), .by_group = TRUE) %>%
  mutate(suffix1 = paste("-", generate_suffix(n()), sep = '')) %>%
  ungroup() %>%
  mutate(suffix1 = ifelse(iso.num == 1, "", suffix1)) %>%
  group_by(Mature_mod) %>%
  mutate(suffix2 = paste("-", row_number(), sep = ''),
         suffix1 = dplyr::first(suffix1)) %>%
  ungroup() %>%
  mutate(suffix2 = ifelse(dup.num == 0, "", suffix2)) %>%
  mutate(pre_miRNA = paste(prefix, suffix1, suffix2, sep = '')) 

# 11.2 compare dominant and alt ----
compare.df1 <- novel.df1 %>% select(PBid, prefix) %>%
  mutate(prefix_rank = as.integer(stringr::str_extract(prefix, "(?<=novel-)\\d+")))

compare.df2 <- novel_alt.df1 %>% select(PBid, prefix) %>%
  mutate(prefix_rank = as.integer(stringr::str_extract(prefix, "(?<=novel-)\\d+"))) %>%
  dplyr::rename(prefix_alt = prefix, 
                prefix_rank_alt = prefix_rank)

compare_novel_strand <- compare.df1 %>% 
  left_join(compare.df2, by = "PBid") %>%
  mutate(alt_higher = ifelse(prefix_rank_alt < prefix_rank -1 , TRUE, FALSE))

novel_swap.list <- compare_novel_strand %>% 
  filter(alt_higher == T) %>%
  pull(PBid)

# 12 re-assign novel_rank----
novel.df.swap <- novel.df0 %>%
  mutate(seed6 = ifelse(PBid %in% novel_swap.list, seed6_alt, seed6)) %>%
  select(-seed6_alt)

novel_fam.df2 <- novel.df.swap %>%
  group_by(seed6)  %>%
  summarize(family_reads = sum(total_reads),
            fam.size = n(),
            iso.num = n_distinct(Mature_mod),
            dup.num = fam.size - iso.num) %>%
  ungroup() %>%
  arrange(desc(family_reads))

novel.list2 <- novel_fam.df2 %>%
  mutate(suffix = paste("-", 1:nrow(novel_fam.df2), sep = "")) %>%
  mutate(family = paste("hmi-novel", suffix, sep = "")) %>%
  select(-family_reads, -suffix)

generate_suffix <- function(n) {
  base_letters <- letters
  if (n <= 26) {
    return(base_letters[1:n])
  } else {
    extra_letters <- unlist(lapply(1:ceiling(n / 26), function(x) paste0(letters[x], base_letters)))
    return(c(base_letters, extra_letters)[1:n])
  }
}

novel.df2 <- novel.df.swap %>% 
  left_join(novel.list2, by = "seed6") %>%
  dplyr::rename(prefix = family) %>%
  group_by(prefix) %>%
  arrange(desc(total_reads), .by_group = TRUE) %>%
  mutate(suffix1 = paste("-", generate_suffix(n()), sep = '')) %>%
  ungroup() %>%
  mutate(suffix1 = ifelse(iso.num == 1, "", suffix1)) %>%
  group_by(Mature_mod) %>%
  mutate(suffix2 = paste("-", row_number(), sep = ''),
         suffix1 = dplyr::first(suffix1)) %>%
  ungroup() %>%
  mutate(suffix2 = ifelse(dup.num == 0, "", suffix2)) %>%
  mutate(pre_miRNA = paste(prefix, suffix1, suffix2, sep = '')) %>%
  mutate(seed_match_type = ifelse(PBid %in% novel_swap.list, "6mer_alt", "6mer"))

novel_named.df <- novel.df2 %>%
  mutate(
    Mature_mod_orig = Mature_mod,
    Star_mod_orig = Star_mod
  ) %>%
  mutate(
    Mature_mod = ifelse(PBid %in% novel_swap.list, Star_mod_orig, Mature_mod_orig),
    Star_mod = ifelse(PBid %in% novel_swap.list, Mature_mod_orig, Star_mod_orig)
  ) %>%
  select(-iso.num, -suffix1, -suffix2, -seed6, -seed_alt, -seed) %>%
  select(-Mature_mod_orig, -Star_mod_orig) %>%
  dplyr::rename(family = prefix)

# 13 Combine all and export ----
all_combine.df <- rbind(ref_named.df, novel_named.df) %>%
  left_join(filter.df9 %>% select(PBid, X3p.5p_ziv), by = "PBid") %>%
  dplyr::arrange(as.integer(sub("PB_", "", PBid)))

all_combine.short <- all_combine.df %>%
  select(PBid, `pre_miRNA`, family, seed_match_type, naming_species, total_reads)

iso_combine.df <- all_combine.df %>% 
  filter(PBid %in% isomiR.revise_coord.df$PBid)

iso_combine.short <- all_combine.short %>% 
  filter(PBid %in% isomiR.revise_coord.df$PBid)


write.csv(all_combine.df, "PB_sRNA/Key_files/name.full_iso_revise.csv", row.names = F)
write.csv(all_combine.short, "PB_sRNA/Key_files/name.short_iso_revise.csv", row.names = F)










