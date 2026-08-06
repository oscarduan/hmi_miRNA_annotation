source("C:/Users/duany/OneDrive/Desktop/Srivastava Lab/sRNA-seq/Isana.2023.11/scripts/B00.head.R")
library(dplyr)
library(readxl)
library(tidyr)
library(Biostrings)
library(ggplot2)
library(stringr)
library(patchwork)
library(ggrepel)
library(data.table)
`%notin%` <- Negate(`%in%`)
library(reshape2)
library(gplots)
library(DESeq2)
library(tibble)

# 1 process read file ----
count1.miR <- function(subread.file, name){
  
  all.df <- fread(file = subread.file, 
                  skip = 1, header = TRUE,  select =  c(1,7:10)) # change this to include gene_name and counts (last col)
  colnames(all.df)[c(1,2,5)] <- c("miRNA","PBid","raw")
  # all.df <- all.df %>% select(c(1,2))
  
  miR.df <- all.df
  
  miR.df$rpm <- format(sweep(miR.df[,5], 2, colSums(miR.df[,5])/1000000, '/'),
                       digit = 3, scientific = F)
  miR.df$name <- as.character(name)
  miR.df
}

df.miR.H00_1 <- count1.miR("240126 regeneration whole/PB/subread/A1_H00.miRNA.txt", "H00_1") %>% mutate(rpm = as.numeric(rpm))  
df.miR.H00_2 <- count1.miR("240126 regeneration whole/PB/subread/A2_H00.miRNA.txt", "H00_2") %>% mutate(rpm = as.numeric(rpm))  
df.miR.H00_3 <- count1.miR("240126 regeneration whole/PB/subread/A3_H00.miRNA.txt", "H00_3") %>% mutate(rpm = as.numeric(rpm))  

df.miR.T00_1 <- count1.miR("240126 regeneration whole/PB/subread/B1_T00.miRNA.txt", "T00_1") %>% mutate(rpm = as.numeric(rpm))  
df.miR.T00_2 <- count1.miR("240126 regeneration whole/PB/subread/B2_T00.miRNA.txt", "T00_2") %>% mutate(rpm = as.numeric(rpm))  
df.miR.T00_3 <- count1.miR("240126 regeneration whole/PB/subread/B3_T00.miRNA.txt", "T00_3") %>% mutate(rpm = as.numeric(rpm))  

miR.raw.indv.df <- df.miR.H00_1 %>% transmute(miRNA, PBid, pre_miRNA, strand, H00_1 = raw) %>%
  left_join(df.miR.H00_2 %>% transmute(miRNA, PBid, pre_miRNA, strand, H00_2 = raw),
            by = c("miRNA", "PBid", "pre_miRNA", "strand")) %>%
  left_join(df.miR.H00_3 %>% transmute(miRNA, PBid, pre_miRNA, strand, H00_3 = raw),
            by = c("miRNA", "PBid", "pre_miRNA", "strand")) %>%
  left_join(df.miR.T00_1 %>% transmute(miRNA, PBid, pre_miRNA, strand, T00_1 = raw),
            by = c("miRNA", "PBid", "pre_miRNA", "strand")) %>%
  left_join(df.miR.T00_2 %>% transmute(miRNA, PBid, pre_miRNA, strand, T00_2 = raw),
            by = c("miRNA", "PBid", "pre_miRNA", "strand")) %>%
  left_join(df.miR.T00_3 %>% transmute(miRNA, PBid, pre_miRNA, strand, T00_3 = raw),
            by = c("miRNA", "PBid", "pre_miRNA", "strand"))

miR.rpm.indv.df <- df.miR.H00_1 %>% transmute(miRNA, PBid, pre_miRNA, strand, H00_1 = rpm) %>%
  left_join(df.miR.H00_2 %>% transmute(miRNA, PBid, pre_miRNA, strand, H00_2 = rpm),
            by = c("miRNA", "PBid", "pre_miRNA", "strand")) %>%
  left_join(df.miR.H00_3 %>% transmute(miRNA, PBid, pre_miRNA, strand, H00_3 = rpm),
            by = c("miRNA", "PBid", "pre_miRNA", "strand")) %>%
  left_join(df.miR.T00_1 %>% transmute(miRNA, PBid, pre_miRNA, strand, T00_1 = rpm),
            by = c("miRNA", "PBid", "pre_miRNA", "strand")) %>%
  left_join(df.miR.T00_2 %>% transmute(miRNA, PBid, pre_miRNA, strand, T00_2 = rpm),
            by = c("miRNA", "PBid", "pre_miRNA", "strand")) %>%
  left_join(df.miR.T00_3 %>% transmute(miRNA, PBid, pre_miRNA, strand, T00_3 = rpm),
            by = c("miRNA", "PBid", "pre_miRNA", "strand"))


write.csv(miR.raw.indv.df, "PB_sRNA/Key_files/head_tail/miR_raw.indv.df", row.names = F)
write.csv(miR.rpm.indv.df, "PB_sRNA/Key_files/head_tail/miR_rpm.indv.df", row.names = F)

merge.ms <- function(count1.miR.df){
  count1.miR.df %>% 
    group_by(PBid) %>%
    summarise(
      pre_miRNA = first(.data$pre_miRNA),
      name      = first(.data$name),
      raw       = sum(.data$raw, na.rm = TRUE),
      rpm       = sum(.data$rpm, na.rm = TRUE),
      .groups   = "drop"
    )
}

df.ms.H00_1 <- merge.ms(df.miR.H00_1)
df.ms.H00_2 <- merge.ms(df.miR.H00_2)
df.ms.H00_3 <- merge.ms(df.miR.H00_3)

df.ms.T00_1 <- merge.ms(df.miR.T00_1)
df.ms.T00_2 <- merge.ms(df.miR.T00_2)
df.ms.T00_3 <- merge.ms(df.miR.T00_3)

ms.raw.indv.df <- df.ms.H00_1 %>%
  transmute(Geneid = pre_miRNA, H00_1 = raw) %>%
  left_join(df.ms.H00_2 %>% transmute(Geneid = pre_miRNA, H00_2 = raw), by = "Geneid") %>%
  left_join(df.ms.H00_3 %>% transmute(Geneid = pre_miRNA, H00_3 = raw), by = "Geneid") %>%
  left_join(df.ms.T00_1 %>% transmute(Geneid = pre_miRNA, T00_1 = raw), by = "Geneid") %>%
  left_join(df.ms.T00_2 %>% transmute(Geneid = pre_miRNA, T00_2 = raw), by = "Geneid") %>%
  left_join(df.ms.T00_3 %>% transmute(Geneid = pre_miRNA, T00_3 = raw), by = "Geneid") %>%
  mutate(across(-Geneid, ~ as.integer(replace_na(., 0))))

ms.rpm.indv.df <- df.ms.H00_1 %>%
  transmute(Geneid = pre_miRNA, H00_1 = rpm) %>%
  left_join(df.ms.H00_2 %>% transmute(Geneid = pre_miRNA, H00_2 = rpm), by = "Geneid") %>%
  left_join(df.ms.H00_3 %>% transmute(Geneid = pre_miRNA, H00_3 = rpm), by = "Geneid") %>%
  left_join(df.ms.T00_1 %>% transmute(Geneid = pre_miRNA, T00_1 = rpm), by = "Geneid") %>%
  left_join(df.ms.T00_2 %>% transmute(Geneid = pre_miRNA, T00_2 = rpm), by = "Geneid") %>%
  left_join(df.ms.T00_3 %>% transmute(Geneid = pre_miRNA, T00_3 = rpm), by = "Geneid")

raw.indv.df <- ms.raw.indv.df
rpm.indv.df <- ms.rpm.indv.df

write.csv(ms.raw.indv.df, "PB_sRNA/Key_files/head_tail/ms_raw.indv.df", row.names = FALSE)
write.csv(ms.rpm.indv.df, "PB_sRNA/Key_files/head_tail/ms_rpm.indv.df", row.names = FALSE)

# 2 Volc plot ----
raw.input <- raw.indv.df %>%
  dplyr::rename(pre_miRNA = Geneid) %>%
  mutate(pre_miRNA = gsub("hmi-","",pre_miRNA)) %>%
  mutate(across(-pre_miRNA, ~ as.integer(replace_na(., 0))))
# left_join(miR_pre.list %>% select(ivl.pre.id, pre_miRNA), by = "pre_miRNA") %>%
# group_by(pre_miRNA) %>%
# summarise(across(2:7, sum, .names = "sum_{col}")) %>%
# ungroup()


count_data <- as.matrix(raw.input[, 2:7])
rownames(count_data) <- raw.input$pre_miRNA
col_data <- data.frame(condition = rep(c("H00", "T00"), each = 3))
rownames(col_data) <- colnames(count_data)
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "T00", "H00"))
res_ordered <- res[order(res$padj, na.last = NA), ]
output.df <- as.data.frame(res_ordered)
output.df$pre_miRNA <- rownames(res_ordered)

# write.csv(output.df, "Isana.2023.11/mGDB_30/head_tail/AP.deseq2.csv", row.names = F)
hi_conf_pre.list <- name.short %>% 
  filter(PBid %in% hi_conf.list) %>% 
  mutate(pre_miRNA = gsub("hmi-","", pre_miRNA)) %>%
  pull(pre_miRNA)

volc1.df <- output.df %>%
  filter(pre_miRNA %in% as.vector(hi_conf_pre.list)) %>%
  mutate(sig = case_when(log2FoldChange > log2(1.5) & padj <= 0.05 ~ "Posterior",
                         log2FoldChange < -log2(1.5) & padj <= 0.05 ~ "Anterior",
                         TRUE ~ "Not Significant"))

volc1.sig.df <- volc1.df %>% 
  filter(sig != "Not Significant") %>%
  mutate(pre_miRNA = gsub("hmi-", "", pre_miRNA))

volc.p1 <- ggplot(volc1.df, aes(x = -log2FoldChange, y = -log10(padj))) +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = log2(1.5), linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log2(1.5), linetype = "dashed", color = "black") +
  geom_point(aes(color = sig), size = 2.5, shape = 16, alpha = 1) +
  geom_text_repel(data = volc1.sig.df, aes(label = pre_miRNA, color = sig),
                  size = 6.5,
                  max.overlaps = 50,
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0, "lines"),
                  show.legend = F)  +
  scale_x_continuous(limits = c(-3, 4), breaks = c(-3, -2, -1, -0.58, 0, 0.58, 1, 2, 3, 4, 5), oob = scales::squish) +
  scale_y_continuous(limits = c(0, 10), oob = scales::squish) +
  scale_color_manual(values = c(
    "Anterior" = "firebrick",   # Replace "sig1" with the actual value for sig
    "Posterior" = "steelblue3",    # Replace "sig2" with the actual value for sig
    "Not Significant" = "grey70" # Replace "non-sig" with the actual value for non-significant points
  )) +
  labs(x = "log2FC (Anterior v.s. Posterior)") + 
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14.5),
        legend.title = element_blank(), 
        legend.text = element_text(size = 16, face = "bold"))

volc.p1
ggsave(plot = volc.p1, filename = "PB_sRNA/Key_files/key_pic1/PDF_pic/volcano.pdf", width = 275, height = 150, units = "mm", dpi = 400)

################################ Table S4 ############################
export.S4A <- volc1.df

export.S4B <- miR.rpm.indv.df

export.S4C <- miR.raw.indv.df %>%
  mutate(across(where(is.numeric), round))

export.S4D <- ms.rpm.indv.df %>%
  mutate(across(where(is.numeric), round))

export.S4E <- ms.raw.indv.df %>%
  mutate(across(where(is.numeric), round))


library(openxlsx)
wb4 <- createWorkbook()

addWorksheet(wb4, "A.DESeq2_results")
writeData(wb4, sheet = "A.DESeq2_results", export.S4A)

addWorksheet(wb4, "B.RPM_3p5p_miR")
writeData(wb4, sheet = "B.RPM_3p5p_miR", export.S4B)

addWorksheet(wb4, "C.raw_3p5p_miR")
writeData(wb4, sheet = "C.raw_3p5p_miR", export.S4C)

addWorksheet(wb4, "D.RPM_3p+5p")
writeData(wb4, sheet = "D.RPM_3p+5p", export.S4D)

addWorksheet(wb4, "E.raw_3p+5p")
writeData(wb4, sheet = "E.raw_3p+5p", export.S4E)

saveWorkbook(wb4, "PB_sRNA/Key_files/Table_S4.xlsx", overwrite = TRUE)

###############################################################################################

# 3 MA plot ---
ma.plot.df <- rpm.indv.df %>%
  dplyr::rename(pre_miRNA = Geneid) %>%
  mutate(pre_miRNA = gsub("hmi-","",pre_miRNA)) %>%
  filter(pre_miRNA %in% hi_conf_pre.list) %>% 
  rowwise() %>% 
  mutate(mean_rpm = mean(c_across(2:7), na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(volc1.df %>% select(pre_miRNA, log2FoldChange, pvalue, padj, sig), by = "pre_miRNA") %>%
  filter(if_all(everything(), ~ !is.na(.)))

ma.sig.df <- ma.plot.df %>% 
  filter(sig != "Not Significant") %>%
  mutate(pre_miRNA = gsub("hmi-", "", pre_miRNA))

ma.p1 <- ggplot(ma.plot.df, aes(x = log10(mean_rpm), y = -log2FoldChange, color = sig)) +
  geom_point(size = 2.5, shape = 16, alpha = 0.79) +
  geom_hline(yintercept = -log2(1.5), linetype = "dashed", color = "grey30") +
  geom_hline(yintercept = log2(1.5), linetype = "dashed", color = "grey30") +
  theme_classic() +
  geom_text_repel(data = ma.sig.df, aes(label = pre_miRNA, color = sig),
                  size = 6,
                  max.overlaps = 50,
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0, "lines"),
                  show.legend = F)  +
  scale_y_continuous(limits = c(-2, 4), oob = scales::squish) +
  scale_x_continuous(limits = c(1, 5), oob = scales::squish) +
  scale_color_manual(values = c(
    "Anterior" = "firebrick",   # Replace "sig1" with the actual value for sig
    "Posterior" = "steelblue3",    # Replace "sig2" with the actual value for sig
    "Not Significant" = "grey70" # Replace "non-sig" with the actual value for non-significant points
  )) +
  labs(x = "Expression", y = "log2 FC") + 
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16),
        legend.title = element_blank(), 
        legend.text = element_text(size = 16, face = "bold"))

ma.p1  
ggsave(plot = ma.p1, filename = "PB_sRNA/Key_files/key_pic1/PDF_pic/MA.pdf", width = 295, height = 150, units = "mm", dpi = 400)


# let7 mir-125a
let7_125_head_tail_rpm.df <- ms.rpm.indv.df %>%
  dplyr::rename(pre_miRNA = Geneid) %>%
  dplyr::filter(pre_miRNA %in% c("hmi-let-7", "hmi-mir-125-a", "let-7", "mir-125-a")) %>%
  dplyr::mutate(
    pre_miRNA = dplyr::case_when(
      pre_miRNA == "let-7" ~ "hmi-let-7",
      pre_miRNA == "mir-125-a" ~ "hmi-mir-125-a",
      TRUE ~ pre_miRNA
    )
  ) %>%
  tidyr::pivot_longer(
    cols = -pre_miRNA,
    names_to = "sample",
    values_to = "rpm_ms"
  ) %>%
  dplyr::mutate(
    tissue = dplyr::case_when(
      stringr::str_detect(sample, "^H00") ~ "Head",
      stringr::str_detect(sample, "^T00") ~ "Tail",
      TRUE ~ NA_character_
    ),
    replicate = stringr::str_extract(sample, "[0-9]+$"),
    tissue = factor(tissue, levels = c("Head", "Tail")),
    pre_miRNA = factor(pre_miRNA, levels = c("hmi-let-7", "hmi-mir-125-a"))
  ) %>%
  dplyr::filter(!is.na(tissue))

let7_125_head_tail_mean.df <- let7_125_head_tail_rpm.df %>%
  dplyr::group_by(pre_miRNA, tissue) %>%
  dplyr::summarise(mean_rpm_ms = mean(rpm_ms, na.rm = TRUE), .groups = "drop")

p_let7_125_head_tail_rpm <- ggplot(
  let7_125_head_tail_rpm.df %>% dplyr::filter(rpm_ms > 0),
  aes(x = tissue, y = rpm_ms, color = pre_miRNA, group = pre_miRNA)
) +
  geom_line(
    data = let7_125_head_tail_mean.df %>% dplyr::filter(mean_rpm_ms > 0),
    aes(y = mean_rpm_ms),
    linewidth = 0.6
  ) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 2,
    alpha = 0.85
  ) +
  theme_classic() +
  labs(x = NULL, y = expression(log[10]("mature + star RPM")), color = "miRNA") +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("hmi-let-7" = "firebrick3", "hmi-mir-125-a" = "dodgerblue3")) +
  scale_y_log10()

p_let7_125_head_tail_rpm
ggsave(
  plot = p_let7_125_head_tail_rpm,
  filename = "PB_sRNA/Key_files/key_pic1/PDF_pic/let7_125_head_tail_ms_rpm.pdf",
  width = 85,
  height = 75,
  units = "mm",
  dpi = 300
)
