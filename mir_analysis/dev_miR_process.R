source("C:/Users/duany/OneDrive/Desktop/Srivastava Lab/sRNA-seq/Isana.2023.11/scripts/B00.head.R")

count.df_iso <- read.csv("PB_sRNA/STAR_recount3_all_iso.csv", header = T) %>%
  filter(count_method == "frac") %>%
  select(-count_method)

desired_order1 <- c("EC", "GA", "DI", "PDi", "PDii", "PL", "PH", "HL", "IST", "AMP", "SMA")
sample_order1 <- unlist(lapply(desired_order1, function(stage) {
  sort(unique(count.df_iso$Sample[grepl(paste0("^", stage, "[0-9]+$"), count.df_iso$Sample)]))
}))

build_dev_raw <- function(df, features = NULL) {
  if (!is.null(features)) {
    df <- df %>% filter(feature %in% features)
  }

  df %>%
    group_by(Geneid, Sample) %>%
    summarise(raw_reads = sum(dedup), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Sample, values_from = raw_reads, values_fill = 0) %>%
    as.data.frame() %>%
    {
      rownames(.) <- .$Geneid
      .[, c("Geneid", sample_order1), drop = FALSE]
    }
}

build_dev_rpm <- function(df, features = NULL) {
  if (!is.null(features)) {
    df <- df %>% filter(feature %in% features)
  }

  sample_rpm <- df %>%
    group_by(Sample) %>%
    mutate(rpm = dedup / sum(dedup) * 1e6) %>%
    ungroup() %>%
    mutate(dev_stage = sub("[0-9]+$", "", Sample))

  sample_rpm %>%
    group_by(Geneid, dev_stage) %>%
    summarise(mean_rpm = mean(rpm), .groups = "drop") %>%
    mutate(dev_stage = factor(dev_stage, levels = desired_order1)) %>%
    arrange(Geneid, dev_stage) %>%
    tidyr::pivot_wider(names_from = dev_stage, values_from = mean_rpm, values_fill = 0) %>%
    as.data.frame() %>%
    {
      rownames(.) <- .$Geneid
      .[, c("Geneid", desired_order1), drop = FALSE]
    }
}

dev.mature.raw <- build_dev_raw(count.df_iso, "mature")
dev.star.raw <- build_dev_raw(count.df_iso, "star")
dev.ms.raw <- build_dev_raw(count.df_iso, c("mature", "star"))

dev.mature.rpm <- build_dev_rpm(count.df_iso, "mature")
dev.star.rpm <- build_dev_rpm(count.df_iso, "star")
dev.ms.rpm <- build_dev_rpm(count.df_iso, c("mature", "star"))

output_xlsx <- "PB_sRNA/dev.count_ms.xlsx"

if (nzchar(output_xlsx)) {
  openxlsx::write.xlsx(
    list(
      mature_raw = dev.mature.raw,
      star_raw = dev.star.raw,
      ms_raw = dev.ms.raw,
      mature_rpm = dev.mature.rpm,
      star_rpm = dev.star.rpm,
      ms_rpm = dev.ms.rpm
    ),
    file = output_xlsx,
    rowNames = FALSE,
    overwrite = TRUE
  )
}

# by 3p 5p ----
count.df_iso1 <- count.df_iso %>%
  left_join(name.full %>%
              select(PBid, pre_miRNA, X3p5p_mod) %>%
              dplyr::rename(Geneid = PBid), 
            by = "Geneid")

m.df <- count.df_iso1 %>%
  filter(feature == "mature") %>%
  mutate(miRNA = paste0(pre_miRNA, "-", X3p5p_mod))

s.df <- count.df_iso1 %>%
  filter(feature == "star") %>%
  mutate(X3p5p_mod0 = ifelse(X3p5p_mod == "5p", "3p", "5p")) %>%
  mutate(miRNA = paste0(pre_miRNA, "-", X3p5p_mod0)) %>%
  select(-X3p5p_mod0)

ms.df <- rbind(m.df, s.df)

build_ms_stage_raw <- function(df) {
  df %>%
    mutate(dev_stage = sub("[0-9]+$", "", Sample)) %>%
    group_by(Geneid, miRNA, pre_miRNA, dev_stage) %>%
    summarise(raw_reads = sum(dedup), .groups = "drop") %>%
    mutate(dev_stage = factor(dev_stage, levels = desired_order1)) %>%
    arrange(Geneid, miRNA, pre_miRNA, dev_stage) %>%
    tidyr::pivot_wider(names_from = dev_stage, values_from = raw_reads, values_fill = 0) %>%
    dplyr::rename(PBid = Geneid) %>%
    select(PBid, miRNA, pre_miRNA, all_of(desired_order1))
}

build_ms_stage_rpm <- function(df) {
  df %>%
    group_by(Sample) %>%
    mutate(rpm = dedup / sum(dedup) * 1e6) %>%
    ungroup() %>%
    mutate(dev_stage = sub("[0-9]+$", "", Sample)) %>%
    group_by(Geneid, miRNA, pre_miRNA, dev_stage) %>%
    summarise(mean_rpm = mean(rpm), .groups = "drop") %>%
    mutate(dev_stage = factor(dev_stage, levels = desired_order1)) %>%
    arrange(Geneid, miRNA, pre_miRNA, dev_stage) %>%
    tidyr::pivot_wider(names_from = dev_stage, values_from = mean_rpm, values_fill = 0) %>%
    dplyr::rename(PBid = Geneid) %>%
    select(PBid, miRNA, pre_miRNA, all_of(desired_order1))
}

ms.stage.raw <- build_ms_stage_raw(ms.df)
ms.stage.rpm <- build_ms_stage_rpm(ms.df)

output_xlsx_ms_stage <- "PB_sRNA/dev.count_3p5p.xlsx"

if (nzchar(output_xlsx_ms_stage)) {
  openxlsx::write.xlsx(
    list(
      raw = ms.stage.raw,
      rpm = ms.stage.rpm
    ),
    file = output_xlsx_ms_stage,
    rowNames = FALSE,
    overwrite = TRUE
  )
}

