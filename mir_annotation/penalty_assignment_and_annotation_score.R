source("C:/Users/duany/OneDrive/Desktop/Srivastava Lab/sRNA-seq/Isana.2023.11/scripts/B00.head.R")

# 0 input ----
input.df0 <- name.short %>% select(PBid, pre_miRNA)

# 1 Mature & Star Length ----
length.df <- read.csv("PB_sRNA/File_Export/length.csv", header = T)

penalty.df1 <- input.df0 %>%
  left_join(length.df %>% select(PBid, Lm, Ls), by = "PBid") %>%
  mutate(Pen.Lm = case_when(
    Lm >= 20 & Lm <= 26 ~ 0,
    Lm == 19 ~ 5,
    Lm == 18 ~ 10,
    TRUE ~ NA_real_   
    )) %>%
  mutate(Pen.Ls = case_when(
    Ls >= 20 & Ls <= 26 ~ 0,
    Ls <= 19 ~ 5,
    TRUE ~ NA_real_   
  ))

# 2 MFE  ----
penalty.df2 <- penalty.df1 %>%
  left_join(read.csv("PB_sRNA/File_Export/hmi_MFE.csv", header = T) %>%  
              select(PBid, MFE), 
            by = "PBid") %>%
  mutate(Pen.MFE = ifelse(MFE <= -17.4, 0,
                          ifelse(MFE <= -14.7, round((14.81*MFE + 257.78),0), 
                                 40)))


# 3 Mature Base Pairng ratio ----
penalty.df3 <- penalty.df2 %>%
  left_join(filter.df10 %>% 
              select(PBid, Mature_BP_ratio_ziv),
            by = "PBid") %>%
  mutate(Pen.Mature_BP_ratio = ifelse(Mature_BP_ratio_ziv >= 0.69, 0,
                                      ifelse(Mature_BP_ratio_ziv >= 0.62, round((-416.67*Mature_BP_ratio_ziv + 288.33),0), 
                                             40))) %>%
  dplyr::rename(Mature_BP_ratio = Mature_BP_ratio_ziv)


# 4 Max bulge asymmetry ----
penalty.df4 <- penalty.df3 %>%
  left_join(filter.df10 %>% 
              select(PBid, Mature_max_bulge_asymmetry_ziv),
            by = "PBid") %>%
  mutate(Pen.max.Bulge.symmetry = ifelse(Mature_max_bulge_asymmetry_ziv <= 2, 0,
                                       ifelse(Mature_max_bulge_asymmetry_ziv == 3, 20, 
                                              ifelse(Mature_max_bulge_asymmetry_ziv == 4, 30, 
                                                     40)) )) %>%
  dplyr::rename(max_bulge_asymmetry = Mature_max_bulge_asymmetry_ziv)

# 5 3' overhang ----
penalty.df5 <- penalty.df4 %>%
  left_join(filter.df10 %>% 
              select(PBid, X3p_overhang_ziv),
            by = "PBid") %>%
  dplyr::rename(X3p_overhang = X3p_overhang_ziv) %>%
  mutate(
    Pen.3p.overhang = case_when(
      X3p_overhang == 2                     ~ 0,
      X3p_overhang == 3                     ~ 5,
      X3p_overhang == 1                     ~ 10,
      X3p_overhang %in% c(0, -1) |
        X3p_overhang >= 4                   ~ 15,
      X3p_overhang <= -2                    ~ 40,
      TRUE                                  ~ NA_real_
    )
  )

# 6 loop size ----
hmia_loop.df <- read.csv("PB_sRNA/File_Export/hmia_loop.df.csv", header = T)

penalty.df6 <- penalty.df5 %>%
  left_join(filter.df10 %>% 
              select(PBid, loop_size),
            by = "PBid") %>%
  dplyr::rename(loop_size_sequential = loop_size) %>%
  left_join(hmia_loop.df %>%
              select(PBid, loop_size_structural),
            by = "PBid") %>%
  mutate(
    Pen.loop_size_structural = case_when(
      loop_size_structural >= 6                     ~ 0,
      loop_size_structural == 5                     ~ 3,
      loop_size_structural == 4                     ~ 5,
      loop_size_structural == 3                     ~ 10,

      TRUE                                  ~ NA_real_
    )
  ) %>%
  mutate(
    Pen.loop_size_sequential = case_when(
      loop_size_sequential >= 6                     ~ 0,
      loop_size_sequential == 5                     ~ 5,
      loop_size_sequential <= 4                     ~ 10,
      TRUE                                  ~ NA_real_
    )
  )
  
# 7 Total ms reads ----
penalty.df7 <- penalty.df6 %>%
  left_join(count.df_iso %>% 
              select(PBid, total_mature, total_star, total_ms, ms_rank_percentile),
            by = "PBid") %>%
  mutate(
    Pen.total_ms = case_when(
      total_ms > 1995                     ~ 0,
      total_ms <= 1995 & total_ms > 270   ~ round((-0.0139*total_ms+ 28.757), 0),
      total_ms <= 270                    ~ 25,
      TRUE                                  ~ NA_real_
    ) ) 
  

# 8 Max mature per-base cov. ----
penalty.df8 <- penalty.df7 %>%
  left_join(count_max.df_iso %>% dplyr::select(PBid, max_mature, m.pBC), by = "PBid") %>%
  dplyr::mutate(
    m.pBC = as.numeric(m.pBC),   # <- crucial if m.pBC is character
    Pen.m_perBase_cov = dplyr::case_when(
      m.pBC > 6.1                     ~ 0,
      m.pBC <= 6.1 & m.pBC > 0.22   ~ round(-4.0816 * m.pBC + 25.898, 0),
      m.pBC <= 0.22                    ~ 25,
      TRUE                             ~ NA_real_
    )
  ) %>%
  dplyr::rename(max_mature_perBase_cov = m.pBC)


# 9 5prime heterogeneity ----
miRGE.df0 <- read_tsv("PB_sRNA/File_Export/isomiR_custom/global_5p_heterogeneity.tsv")
iso_revised <- read.csv("PB_sRNA/Key_files/revise_coord.df.csv", header = T)

miRGE.df1 <- miRGE.df0 %>%
  filter(PBid %notin% iso_revised$PBid)

miRGE.df2 <- miRGE.df0 %>%
  filter(PBid %in% iso_revised$PBid) %>%
  left_join(iso_revised %>%
              select(PBid, global_5p_pct),
            by = "PBid") %>%
  mutate(heterogeneity = (100 - global_5p_pct) / 100) %>%
  select(-global_5p_pct)

miRGE.df <- rbind(miRGE.df1, miRGE.df2)

penalty.df9 <- penalty.df8 %>%
  left_join(miRGE.df %>% select(PBid, heterogeneity), by = "PBid") %>%
  dplyr::rename(mature_5p_heterogeneity = heterogeneity) %>%
  dplyr::mutate(
    mature_5p_heterogeneity = as.numeric(round(mature_5p_heterogeneity * 100,4)),
    Pen.5_hg = dplyr::case_when(
      mature_5p_heterogeneity <= 10 ~ 0,
      mature_5p_heterogeneity > 10 & mature_5p_heterogeneity <= 50 ~ round(0.975 * mature_5p_heterogeneity - 8.75, 0),
      mature_5p_heterogeneity > 50 ~ 40,
      TRUE ~ 10
    )
  )

# 10 in-cluster ratio of m and s ----
icr.df <- read.csv("PB_sRNA/File_Export/icr.df.csv", header = T)

penalty.df10 <- penalty.df9 %>% 
  left_join(icr.df %>% 
              select(PBid, in_cluster_ratio_Mature, in_cluster_ratio_3p_plus_5p), 
            by = "PBid") %>%
  mutate(
    in_cluster_ratio_3p_plus_5p = as.numeric(in_cluster_ratio_3p_plus_5p),
         Pen.in_cluster_ratio = dplyr::case_when(
           in_cluster_ratio_3p_plus_5p >= 0.99 ~ 0,
           in_cluster_ratio_3p_plus_5p < 0.99 & in_cluster_ratio_3p_plus_5p >= 0.75 ~ round(-162.5 * in_cluster_ratio_3p_plus_5p + 161.875, 0),
           in_cluster_ratio_3p_plus_5p < 0.75 ~ 40,
           TRUE ~ NA_real_
         ))

# 11 Detection ----
detect.df <- full.df %>% 
  mutate(detection_algorithm = case_when(
    Description_mirdeep == "." ~ "sRNAbench",
    Description_sRNAbench == "." ~ "miRDeep2",
    Description_mirdeep != "." & Description_sRNAbench != "." ~ "Both",
    TRUE ~ NA_character_
  )) %>%
  select(PBid, detection_algorithm)

penalty.df11 <- penalty.df10 %>%
  left_join(detect.df, by = "PBid") %>%
  mutate(Pen.detection = ifelse(detection_algorithm == "Both", -20, 0))

## 12 summarise ----
penalty.df12 <- penalty.df11 %>%
  mutate(Penalty = rowSums(dplyr::select(., dplyr::starts_with("Pen.")), na.rm = TRUE)) %>%
  mutate(Score = pmax(pmin(100 - Penalty, 100), 0))

write.csv(penalty.df12, "PB_sRNA/Key_files/penalty_full_iso.csv", row.names = F)

## 13 plot big----

plot.df1 <- penalty.df12 %>%
  select(PBid, total_ms, Penalty, Score) %>%
  left_join(fam.origin, by = "PBid") %>%
  mutate(hi_conf = ifelse(PBid %in% hi_conf.list, TRUE, FALSE))

p1_label.df <- plot.df1 %>%
  dplyr::filter(species_label == TRUE) %>%
  dplyr::arrange(dplyr::desc(hi_conf), dplyr::desc(Score), dplyr::desc(total_ms)) %>%
  dplyr::slice_head(n = 18)

write.csv(plot.df1, "PB_sRNA/Key_files/penalty_short_iso.csv", row.names = F)

p1 <- ggplot(
  plot.df1,
  aes(x = log10(total_ms), y = Score)
) +
  geom_point(
    data = dplyr::filter(plot.df1, species_label == TRUE),
    shape = 23,
    fill = "white",
    color = "blue3",
    stroke = 0.85,
    size = 3.5,
    alpha = 0.75
  ) +
  geom_point(
    aes(color = hi_conf),
    shape = 16,
    size = 2.5,
    alpha = 0.85
  ) +
  # ggrepel::geom_text_repel(
  #   data = p1_label.df,
  #   aes(label = family),
  #   size = 4.2,
  #   box.padding = 0.8,
  #   point.padding = 0.5,
  #   min.segment.length = 0,
  #   max.overlaps = Inf,
  #   seed = 123,
  #   segment.alpha = 0.55,
  #   segment.color = "grey55",
  #   color = "black"
  # ) +
  scale_color_manual(
    values = c("TRUE" = "red2", "FALSE" = "grey55"),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(1, 7),
    breaks = 1:7,
    oob = scales::squish
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 10)
  ) +
  theme_bw() +
  labs(
    x = expression(log[10](total_ms)),
    y = "Score"
  ) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )

p2 <- p1 +
  annotate("point", x = 5.8, y = 10, color = "red2", size = 4.5) +
  annotate(
    "text",
    x = 5.9,
    y = 10,
    label = "High-confident",
    hjust = 0,
    size = 5.5
  ) +
  annotate("point", x = 5.8, y = 4, shape = 23, fill = "white", color = "blue3", stroke = 1, size = 4.5) +
  annotate(
    "text",
    x = 5.9,
    y = 4,
    label = "Bilateria + X/A",
    hjust = 0,
    size = 5.5
  )

p2

ggsave(p2, filename = "PB_sRNA/Key_files/key_pic1/Score_all_iso.jpeg", width = 275, height = 150, dpi = 450, unit= "mm")
ggsave(p2, filename = "PB_sRNA/Key_files/key_pic1/PDF_pic/Score_all_iso.pdf", width = 275, height = 150, dpi = 450, unit= "mm")


##  14 plot small ----
plot.df2 <- plot.df1 %>%
  mutate(species_group = dplyr::if_else(species_label, "Bilateria + X/A", "Other"))

cum_total_ms.df <- plot.df2 %>%
  dplyr::select(PBid, species_group, total_ms) %>%
  dplyr::filter(!is.na(total_ms)) %>%
  dplyr::group_by(species_group) %>%
  dplyr::arrange(dplyr::desc(total_ms), .by_group = TRUE) %>%
  dplyr::mutate(
    rank = dplyr::row_number(),
    rank_percentile = rank / dplyr::n(),
    cum_value = cumsum(total_ms),
    cum_fraction = cum_value / sum(total_ms)
  ) %>%
  dplyr::ungroup()

p_total_ms <- ggplot(
  cum_total_ms.df,
  aes(x = rank_percentile, y = cum_fraction, color = species_group)
) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c("Bilateria + X/A" = "blue3", "Other" = "grey45"),
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_y_continuous(
    limits = c(0.2, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  theme_test() +
  labs(
    x = " ",
    y = "Cumulative fraction",
    title = "Total reads"
  ) +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = c(0.65, 0.2),
    plot.title = element_text(size = 15, face = "bold")
  )

p_total_ms

cum_score.df <- plot.df2 %>%
  dplyr::select(PBid, species_group, Score) %>%
  dplyr::filter(!is.na(Score)) %>%
  dplyr::group_by(species_group) %>%
  dplyr::arrange(dplyr::desc(Score), .by_group = TRUE) %>%
  dplyr::mutate(
    rank = dplyr::row_number(),
    rank_percentile = rank / dplyr::n(),
    cum_value = cumsum(Score),
    cum_fraction = cum_value / sum(Score)
  ) %>%
  dplyr::ungroup()

# p_score <- ggplot(
#   cum_score.df,
#   aes(x = rank_percentile, y = cum_fraction, color = species_group)
# ) +
#   geom_line(linewidth = 1.2) +
#   scale_color_manual(
#     values = c("Bilateria + X/A" = "skyblue3", "Other" = "grey45"),
#     name = NULL
#   ) +
#   scale_x_continuous(
#     limits = c(0, 1),
#     breaks = seq(0, 1, 0.2)
#   ) +
#   scale_y_continuous(
#     limits = c(0, 1),
#     breaks = seq(0, 1, 0.25)
#   ) +
#   theme_test() +
#   labs(
#     x = "Rank percentile",
#     y = "Cumulative fraction",
#     title = "Annotation Score"
#   ) +
#   theme(
#     axis.title = element_text(size = 12, face = "bold"),
#     axis.text = element_text(size = 12),
#     legend.text = element_text(size = 12),
#     legend.position = c(0.65, 0.2),
#     plot.title = element_text(size = 15, face = "bold")
#   )
# 
# p_score
# 
# p3 <- p_total_ms / p_score
# p3
# 
# ggsave(p3, filename = "PB_sRNA/Key_files/key_pic1/Score_cum.jpeg", width = 80, height = 155, dpi = 300, unit= "mm")
# ggsave(p3, filename = "PB_sRNA/Key_files/key_pic1/PDF_pic/Score_cum.pdf", width = 80, height = 155, dpi = 300, unit= "mm")

violin_score <- ggplot(
  plot.df1,
  aes(
    x = ifelse(species_label, "Bilateria + X/A", "Other"),
    y = Score,
    fill = ifelse(species_label, "Bilateria + X/A", "Other")
  )
) +
  geom_violin(alpha = 0.35, trim = FALSE) +
  # geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.1, size = 0.5, alpha = 0.35) +
  scale_fill_manual(
    values = c("Bilateria + X/A" = "blue3", "Other" = "grey75"),
    guide = "none"
  ) +
  theme_bw() +
  scale_y_continuous(
    limits = c(0, 120),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    x = "",
    y = "Score",
    title = "Annotation Score"
  ) +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = c(0.65, 0.2),
    plot.title = element_text(size = 15, face = "bold")
  )


p4 <- patchwork::wrap_plots(violin_score, p_total_ms, ncol = 1)
p4

ggsave(p4, filename = "PB_sRNA/Key_files/key_pic1/Score_cum_iso.jpeg", width = 80, height = 155, dpi = 300, unit= "mm")
ggsave(p4, filename = "PB_sRNA/Key_files/key_pic1/PDF_pic/Score_cum_iso.pdf", width = 80, height = 155, dpi = 300, unit= "mm")


## 14 penalty term plots ----
# Edit the text on the right-hand side if you want custom axis labels for each term.
penalty_term_labels <- c(
  "Pen.Lm" = "Mature Strand Length",
  "Pen.Ls" = "Star Strand Length",
  "Pen.MFE" = "Minimum Free Energy",
  "Pen.Mature_BP_ratio" = "Mature BP ratio",
  "Pen.max.Bulge.symmetry" = "Max Bulge assymmetry",
  "Pen.3p.overhang" = "3' Pverhang",
  "Pen.loop_size_structural" = "Loop Size Structural",
  "Pen.loop_size_sequential" = "Loop Size Sequential",
  "Pen.total_ms" = "Total Mature + Star Reads",
  "Pen.m_perBase_cov" = "Max perBase Coverate",
  "Pen.5_hg" = "5' Heterogeneity",
  "Pen.in_cluster_ratio" = "In-cluster Ratio",
  "Pen.detection" = "Detected Algorithm"
)

penalty_term_order <- names(penalty_term_labels)

penalty_long.df <- penalty.df12 %>%
  dplyr::select(PBid, pre_miRNA, Penalty, Score, dplyr::starts_with("Pen.")) %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("Pen."),
    names_to = "Penalty_term",
    values_to = "Penalty_value"
  ) %>%
  dplyr::mutate(
    Penalty_value_plot = dplyr::if_else(
      Penalty_term == "Pen.detection",
      -Penalty_value,
      Penalty_value
    ),
    Penalty_value_heatmap = dplyr::if_else(
      Penalty_term == "Pen.detection",
      -Penalty_value_plot,
      Penalty_value_plot
    ),
    Penalty_term = factor(
      Penalty_term,
      levels = penalty_term_order
    ),
    Penalty_term_label = factor(
      penalty_term_labels[as.character(Penalty_term)],
      levels = rev(penalty_term_labels)
    ),
    Penalty_fill = dplyr::if_else(
      Penalty_term == "Pen.detection",
      "Detection",
      "Penalty"
    )
  )

heatmap.df <- penalty_long.df %>%
  dplyr::group_by(PBid, pre_miRNA, Penalty, Score) %>%
  dplyr::summarise(max_term_penalty = max(Penalty_value_plot, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(Penalty), dplyr::desc(max_term_penalty), pre_miRNA) %>%
  dplyr::mutate(pre_miRNA = factor(pre_miRNA, levels = unique(pre_miRNA)))

penalty_heatmap.df <- penalty_long.df %>%
  dplyr::left_join(
    heatmap.df %>% dplyr::select(PBid, pre_miRNA),
    by = c("PBid", "pre_miRNA")
  ) %>%
  dplyr::mutate(pre_miRNA = factor(pre_miRNA, levels = levels(heatmap.df$pre_miRNA)))

p_penalty_heatmap <- ggplot(
  penalty_heatmap.df,
  aes(x = pre_miRNA, y = Penalty_term_label, fill = Penalty_value_heatmap)
) +
  geom_tile() +
  scale_fill_gradientn(
    colours = c("springgreen1", "white", "khaki", "goldenrod3", "lightpink", "tomato", "firebrick4"),
    values = scales::rescale(c(-20, 0, 5, 10, 20, 30, 40)),
    limits = c(-20, 40),
    breaks = c(-20, 0, 10, 20, 30, 40),
    name = "Penalty"
  ) +
  theme_bw() +
  labs(
    x = "miRNA",
    y = ""
    # title = "Penalty term heatmap"
  ) +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 15, face = "bold"),
    panel.grid = element_blank(),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      barwidth = unit(70, "mm"),
      barheight = unit(4, "mm")
    )
  )

p_penalty_violin <- ggplot(
  penalty_long.df,
  aes(x = Penalty_value_plot, y = Penalty_term_label, color = Penalty_fill)
) +
  geom_violin(alpha = 0.5, trim = FALSE, color = NA) +
  geom_jitter(height = 0.12, width = 0, size = 0.85, alpha = 0.2) +
  theme_bw() +
  scale_color_manual(
    values = c("Penalty" = "orchid2", "Detection" = "springgreen3"),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 40),
    breaks = seq(0, 40, 10)
  ) +
  labs(
    x = "Penalty",
    y = ""
    # title = "Penalty term distributions"
  ) +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.ticks.x = element_line(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 15, face = "bold")
  )

p_penalty_pair <- patchwork::wrap_plots(
  p_penalty_heatmap,
  p_penalty_violin,
  nrow = 1,
  widths = c(4.5, 1.2)
)

p_penalty_pair

ggsave(
  p_penalty_pair,
  filename = "PB_sRNA/Key_files/key_pic1/penalty_term_pair.png",
  width = 325,
  height = 95,
  dpi = 300,
  unit = "mm"
)

ggsave(
  p_penalty_pair,
  filename = "PB_sRNA/Key_files/key_pic1/PDF_pic/penalty_term_pair.pdf",
  width = 325,
  height = 95,
  dpi = 300,
  unit = "mm"
)

# penalty distribution
p1_slide <- ggplot(
  plot.df1,
  aes(x = log10(total_ms), y = Score)
) +
  geom_point(
    data = plot.df1,
    shape = 16,
    color = "grey55",
    stroke = 0.85,
    size = 2.5,
    alpha = 0.75
  ) +
  scale_x_continuous(
    limits = c(1, 7),
    breaks = 1:7,
    oob = scales::squish
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 10)
  ) +
  theme_test() +
  labs(
    x = expression(log[10](total_ms)),
    y = "Score"
  ) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )


p2_slide <- ggplot(
  plot.df1,
  aes(x = log10(total_ms), y = Score)
) +
  geom_point(
    aes(color = hi_conf),
    shape = 16,
    size = 2.5,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = c("TRUE" = "red2", "FALSE" = "grey55"),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(1, 7),
    breaks = 1:7,
    oob = scales::squish
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 10)
  ) +
  theme_test() +
  labs(
    x = expression(log[10](total_ms)),
    y = "Score"
  ) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )

p2_slide

slide_exp <- p1_slide + p2_slide
slide_exp


p3_205_154 <- ggplot(
  plot.df1,
  aes(x = log10(total_ms), y = Score)
) +
  geom_point(
    color = "grey70",
    shape = 16,
    size = 1.5,
    alpha = 0.5
  ) +
  geom_point(
    data = dplyr::filter(
      plot.df1,
      PBid %in% hi_conf.list &
      family %in% c(
        "mir-205", "mir-154", "mir-15", "mir-142", "mir-335",
        "mir-929", "mir-1271", "mir-991", "mir-2056", "xbo-novel-1"
      )
    ),
    aes(color = family),
    shape = 16,
    size = 3,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = c(
      "mir-205" = "blue",
      "mir-154" = "magenta",
      "mir-15" = "springgreen2",
      "mir-142" = "darkorange3",
      "mir-335" = "goldenrod2",
      "mir-929" = "forestgreen",
      "mir-1271" = "cyan3",
      "mir-991" = "steelblue4",
      "mir-2056" = "purple3",
      "xbo-novel-1" = "red2"
    ),
    breaks = c(
      "mir-205", "mir-154", "mir-15", "mir-142", "mir-335",
      "mir-929", "mir-1271", "mir-991", "mir-2056", "xbo-novel-1"
    ),
    labels = c(
      "mir-205", "mir-154", "mir-15", "mir-142", "mir-335",
      "mir-929", "mir-1271", "mir-991", "mir-2056", "xbo-novel-1"
    ),
    name = "Family"
  ) +
  scale_x_continuous(
    limits = c(1, 7),
    breaks = 1:7,
    oob = scales::squish
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 10)
  ) +
  theme_test() +
  labs(
    x = expression(log[10](total_ms)),
    y = "Score"
  ) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14)
  )

p3_205_154
ggsave(p3_205_154, filename = "PB_sRNA/Key_files/key_pic1/PDF_pic/score_revise_miR.pdf", width = 140, height = 125, dpi = 300, unit= "mm")

