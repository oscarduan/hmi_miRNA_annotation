source("C:/Users/duany/OneDrive/Desktop/Srivastava Lab/sRNA-seq/Isana.2023.11/scripts/B00.head.R")
library("Mfuzz")
library("pheatmap")
library("paletteer")

xlsx_path <- "PB_sRNA/dev.count_ms.xlsx"

m.rpm <- read_xlsx(xlsx_path, sheet = 4)
s.rpm <- read_xlsx(xlsx_path, sheet = 5)
ms.rpm <-  read_xlsx(xlsx_path, sheet = 6)
# stage_cols <- names(m.rpm)[-1]
# 
# ms.rpm.joined <- m.rpm %>%
#   left_join(s.rpm, by = names(m.rpm)[1], suffix = c(".m", ".s"))
# 
# ms.rpm <- ms.rpm.joined %>%
#   transmute(
#     !!names(m.rpm)[1] := .data[[names(m.rpm)[1]]],
#     !!!setNames(
#       lapply(stage_cols, function(col) {
#         trunc((ms.rpm.joined[[paste0(col, ".m")]] + ms.rpm.joined[[paste0(col, ".s")]]) * 100) / 100
#       }),
#       stage_cols
#     )
#   )

# stages <- c("EC","GA", "DI","PDi","PDii", "PL", "PH", "HL", "IST", "AMP", "SMA")
stages <- c("GA", "DI","PDi","PDii", "PL", "PH", "HL", "IST", "AMP", "SMA")

miR.input0 <- ms.rpm %>% select(-EC) %>%
  filter(Geneid %in% hi_conf.list)

colnames(miR.input0)[1:11] <- c("GENE_ID", stages)

time_row2 <- miR.input0[1,]
time_row2[1,1] <- "TIME"
time_row2[1, 2:11] <- as.list(c(1, 2, 3, 4, 5, 6, 7, 8,9,10))

miR.input <-rbind(time_row2, miR.input0) 

write_delim(miR.input, "PB_sRNA/File_Export/fuzz/miR.dev2_renamed.txt", delim = "\t")


# Start mFuzz
your_data_eset1 <- table2eset(file = "PB_sRNA/File_Export/fuzz/miR.dev2_renamed.txt")

your_data_eset1.r <- filter.NA(your_data_eset1, thres= 1.5)

your_data_eset1.f <- fill.NA(your_data_eset1.r,mode="mean")

tmp1 <- filter.std(your_data_eset1.f,min.std=0)  

your_data_eset1.s <- standardise(your_data_eset1.f)

#this is the object that contains clustering information in it
#here, c is the number of clusters, which you can vary, and m is the fuzzifier value, which also can be varied. 

cl1 <- mfuzz(your_data_eset1.s, c = 4, m = 1.35)

mfuzz.plot2(your_data_eset1.s,cl=cl1,mfrow=c(1,5),colo="fancy",
            ax.col="black",bg = "white",col.axis="black",col.lab="black",
            col.main="black", col.sub="blue",col="blue", cex.main=3, cex.lab=2, 
            centre=T, centre.col = "grey15", centre.lwd = 6.5,
            Xwidth = 22, Xheight = 4
)

# Extract cluster assignments from the fclust object
cluster_assignments1 <- cl1$cluster


miR_cl1_index <- data.frame(Gene = rownames(your_data_eset1.s), Cluster = cluster_assignments1, row.names = NULL)

cl1.df <- miR_cl1_index %>%
  group_by(Cluster) %>%
  summarise(n = n()) %>%
  ungroup()

write.csv(miR_cl1_index, "PB_sRNA/File_Export/fuzz/miR_cl1_index.csv", row.names = F)

cluster_assignments1 <- cl1$cluster

miR_cl1_index <- data.frame(Gene = rownames(your_data_eset1.s), Cluster = cluster_assignments1, row.names = NULL)


centers_miR <- cl1$centers

colors <- paletteer::paletteer_c("grDevices::Tropic", n = 30)

p1 <- pheatmap(centers_miR, 
               main = "miRNA_clustering", 
               cluster_cols=FALSE, 
               color = colors,
               angle_col = 45)

heatmap_grob1 <- p1$gtable

ggsave(plot = heatmap_grob1, file = "PB_sRNA/Key_files/key_pic1/PDF_pic/mFuzz_heatmap.pdf", width = 135, height = 50, units = "mm", dpi = 300)




