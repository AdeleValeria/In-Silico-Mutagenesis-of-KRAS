#-----Import packages-----
library("data.table")
library("ComplexHeatmap")
library("reshape2")
library("ggplot2")
library("circlize")
library("ggrepel")
library("tidyverse")
library("RColorBrewer")
library("cowplot")
library("ggpubr")
library("vioplot")
library("gridExtra")

# DATA CLEANING
#------Import data-------
# Calculated by FoldX
# all possible mutations in 6BOF background (dimer)
dat <- fread("6bof_data.csv")
# all possible mutations in one monomer only (6BOF background)
dat.mono <- fread("6bof_mono_data.csv")
# all possible mutations in WT background (dimer)
dat.WT <-fread("WT_data.csv")
# all possible mutations in one monomer only (WT background)
dat.WT.mono <- fread("WT_mono_data.csv")

# Calculated by FreeSASA
# SASA of each residue in monomer A (6bof background)
sasa.seq.A <- fread("SASA_seq_monoA_6bof.txt", header = FALSE)
# SASA of each residue in monomer B (6bof background)
sasa.seq.B <- fread("SASA_seq_monoB_6bof.txt", header = FALSE)
# SASA of each residue in dimer form (6bof background)
sasa.seq <- fread("SASA_seq_6bof.txt", header = FALSE)
sasa.seq <- data.frame(sasa.seq[1:168,]$V3, sasa.seq[1:168,]$V4, sasa.seq[1:168,]$V6, 
                   sasa.seq[169:336,]$V6, sasa.seq.A$V6, sasa.seq.B$V6)

# SASA of each residue in monomer A ( WT background)
sasa.seq.A.WT <- fread("SASA_seq_monoA_WT.txt", header = FALSE)
# SASA of each residue in monomer B (WT background)
sasa.seq.B.WT <- fread("SASA_seq_monoB_WT.txt", header = FALSE)
# SASA of each residue in dimer form (WT background)
sasa.seq.WT <- fread("SASA_seq_WT.txt", header = FALSE)
sasa.seq.WT <- data.frame(sasa.seq.WT[1:168,]$V3, sasa.seq.WT[1:168,]$V4, sasa.seq.WT[1:168,]$V6, 
                       sasa.seq.WT[169:336,]$V6, sasa.seq.A.WT$V6, sasa.seq.B.WT$V6)

# add column names
# add SASA and SASA strand A
colnames(dat) <- c("Group","Pos","AA","Binding","dG WT",
                   "dG","ddG","SASA","SASA A", "Clashes1","Clashes2")
colnames(dat.WT) <- c("Group", "Pos","AA","Binding","dG WT",
                      "dG","ddG","SASA","SASA A", "Clashes1","Clashes2")
colnames(dat.mono) <- c("Pos", "AA", "dG WT", "dG", "ddG Mono","SASA Mono","Clashes1","Clashes2")
colnames(dat.WT.mono) <- c("Pos", "AA", "dG WT", "dG", "ddG Mono","SASA Mono","Clashes1","Clashes2")
colnames(sasa.seq) <- c("Pos", "AA", "A", "B", "Mono A", "Mono B")
colnames(sasa.seq.WT) <- c("Pos","AA", "A", "B", "Mono A", "Mono B")


# change starting Pos for WT data 
dat.WT$Pos <- dat.WT$Pos + 1
dat.WT.mono$Pos <- dat.WT.mono$Pos + 1
sasa.seq.WT$Pos <- sasa.seq.WT$Pos + 1

# subset mono datasets
dat.mono2 <- dat.mono[,c(1,2,5,6)]
dat.WT.mono2 <- dat.WT.mono[,c(1,2,5,6)]

# sort data
dat <- dat[order(dat$Pos),]
dat.WT <- dat.WT[order(dat.WT$Pos),]

# join the other datasets to the two main dataframe
dat <- left_join(dat, dat.mono2, by = c("Pos", "AA"))
dat.WT <- left_join(dat.WT, dat.WT.mono2, by = c("Pos", "AA"))


#-------AA Information--------------0
AA.type <- data.frame("AA" = c("LYS", "ARG", "HIS", "ASP", "GLU", "SER", "THR",
                               "CYS", "TYR", "ASN", "GLN", "GLY", "ALA", "VAL", 
                               "LEU", "ILE", "MET", "PHE", "TRP", "PRO"),
                      "Type" = c(rep("Basic", 3), rep("Acidic", 2), 
                                 rep("Polar", 6), rep("Non-Polar", 9)))
AA.oneletter <- data.frame("AA" = c("LYS", "ARG", "HIS", "ASP", "GLU", "SER", "THR",
                                    "CYS", "TYR", "ASN", "GLN", "GLY", "ALA", "VAL", 
                                    "LEU", "ILE", "MET", "PHE", "TRP", "PRO"),
                           "Letter" = c("K", "R", "H", "D", "E", "S", "T", "C",
                                        "Y", "N", "Q", "G", "A", "V", "L", "I",
                                        "M", "F", "W", "P"))
AA.type2 <- left_join(AA.type, AA.oneletter, by = "AA")[-1]
colnames(AA.type2) <- c("AA Property", "AA")

#--------Organize 6bof data frame--------------
# compute the mean of SASA of each residue
sasa.seq$`MeanRes SASA` <- rowMeans(sasa.seq[,c("A", "B")], na.rm=TRUE)

# compute the SASA difference between same strands in dimer vs. monomer form 
sasa.seq$`Diff A` <- sasa.seq$A - sasa.seq$`Mono A`
sasa.seq$`Diff B` <- sasa.seq$B - sasa.seq$`Mono B`

# identify buried residues (those with SASA below 0.15)
buried <- sasa.seq[sasa.seq$A <= 0.15 | sasa.seq$B <= 0.15 | 
                 sasa.seq$`Mono A` <= 0.15 | sasa.seq$`Mono B` <= 0.15,]
buried$Location <- "Buried"

# identify interacting residues 
# SASA difference between monomer vs bound form must be at least 15
interacting <- sasa.seq[sasa.seq$`Diff A` <= -15 | sasa.seq$`Diff B` <= -15,]
interacting$Location <- "Interacting"

# combine both buried and interacting data frames by rows
buried_interacting <- rbind(buried, interacting)
buried_interacting <- buried_interacting %>%
  select(Pos, Location)

# left join with the sasa data frame
sasa.seq <- left_join(sasa.seq, buried_interacting, by = "Pos")

# assign exposed for rows that are not categorized as buried or interacting
sasa.seq[is.na(sasa.seq)] <- "Exposed"


# Assign alpha helix and beta sheets region
sasa.seq$Secondary <- NA
sasa.seq[c(3:8, 27:31, 39:45, 48:55, 76:82, 110:115, 140:143),]$Secondary <- "Beta"
sasa.seq[c(15:26, 60:63, 66:73, 86:103, 126:136, 151:167),]$Secondary <- "Alpha"

# Assign KRAS main regions
sasa.seq$Region <- NA
sasa.seq[8:16,]$Region <- "P-Loop"
sasa.seq[27:33,]$Region <- "Switch I"
sasa.seq[57:58,]$Region <- "Switch II"
sasa.seq[114:117,]$Region <- "Base-Binding Loop"

# add amino acid information
sasa.seq <- left_join(sasa.seq, AA.type, by = "AA")

# add sasa data to main 6bof dataframe
sasa.seq <- sasa.seq %>%
  select(Pos, Type, `MeanRes SASA`, Location, Secondary, Region)
colnames(sasa.seq) <- c("Pos", "AA Property WT", "MeanRes SASA", 
                       "Location", "Secondary", "Region")
dat <- left_join(dat, sasa.seq, by = "Pos")
dat <- left_join(dat, AA.type2, by = "AA")
dat$`ddG Binding` <- dat$Binding - (-4.15)

# add stability categories
dat$`TE Category` <- cut(
  dat$ddG, c(-Inf, -1.84, -0.92, -0.46, 0.46, 0.92, 1.84, Inf),
  c("H. Stabilizing", "Stabilizing", "S. Stabilizing",
    "Neutral", "S. Destabilizing", "Destabilizing","H. Destabilizing"))

dat$`Dimerization Category` <- cut(
  dat$`ddG Binding`, 
  c(-Inf, -1.84, -0.92, -0.46, 0.46, 0.92, 1.84, Inf),
  c("H. Stabilizing", "Stabilizing", "S. Stabilizing","Neutral", 
    "S. Destabilizing", "Destabilizing", "H. Destabilizing"))

dat$`TE Mono Category` <- cut(
  dat$`ddG Mono`, c(-Inf, -1.84, -0.92, -0.46, 0.46, 0.92, 1.84, Inf),
  c("H. Stabilizing", "Stabilizing", "S. Stabilizing",
    "Neutral", "S. Destabilizing", "Destabilizing","H. Destabilizing"))

# ONCOGENIC MUTATIONS
# top 40: https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=KRAS
# with search keyword: "substitution"
onco.mut <- rbind(
  # pos 12
  dat[dat$Pos == 12 & dat$AA == "A",], #5
  dat[dat$Pos == 12 & dat$AA == "C",], #4
  dat[dat$Pos == 12 & dat$AA == "D",], #1
  dat[dat$Pos == 12 & dat$AA == "F",], #17
  dat[dat$Pos == 12 & dat$AA == "L",], #31
  dat[dat$Pos == 12 & dat$AA == "R",], #7
  dat[dat$Pos == 12 & dat$AA == "S",], #6
  dat[dat$Pos == 12 & dat$AA == "V",], #2
  # pos 13
  dat[dat$Pos == 13 & dat$AA == "A",], #19
  dat[dat$Pos == 13 & dat$AA == "C",], #8
  dat[dat$Pos == 13 & dat$AA == "D",], #3 
  dat[dat$Pos == 13 & dat$AA == "R",], #16
  dat[dat$Pos == 13 & dat$AA == "S",], #13
  dat[dat$Pos == 13 & dat$AA == "V",], #20
  # pos 14
  dat[dat$Pos == 14 & dat$AA == "I",], #22
  # pos 18
  dat[dat$Pos == 18 & dat$AA == "D",], #29
  # pos 19
  dat[dat$Pos == 19 & dat$AA == "F",], #26
  # pos 22
  dat[dat$Pos == 22 & dat$AA == "K",], #23
  #pos 58
  dat[dat$Pos == 58 & dat$AA == "I",], #28
  # pos 59
  dat[dat$Pos == 59 & dat$AA == "G",], #27
  dat[dat$Pos == 59 & dat$AA == "T",], #21
  # pos 60
  dat[dat$Pos == 60 & dat$AA == "D",], #30
  # pos 61
  dat[dat$Pos == 61 & dat$AA == "E",], #33
  dat[dat$Pos == 61 & dat$AA == "H",], #9
  dat[dat$Pos == 61 & dat$AA == "K",], #15
  dat[dat$Pos == 61 & dat$AA == "L",], #12
  dat[dat$Pos == 61 & dat$AA == "P",], #24
  dat[dat$Pos == 61 & dat$AA == "R",], #11
  # pos 63
  dat[dat$Pos == 63 & dat$AA == "K",], #32
  # pos 117
  dat[dat$Pos == 117 & dat$AA == "N",] #18
)

onco.mut$Status <- "Oncogenic"
onco.mut <- onco.mut %>%
  select(Group, Pos, AA, Status)
dat <- left_join(dat, onco.mut, by = c("Group", "Pos", "AA"))


# subset homozygous rows from the sorted data
homo <- dat[dat$Group == "homo",]
hetero <- dat[dat$Group == "hetero",]

# add TE Change column
homo2 <- homo %>%
  select(Pos, AA, SASA, `dG WT`, dG, ddG, `ddG Binding`, `AA Property WT`, `AA Property`, `TE Category`, `Dimerization Category`)
names(homo2)[names(homo2) == "TE Category"] <- "TE Category Homo"
names(homo2)[names(homo2) == "Dimerization Category"] <- "Dimerization Category Homo"
names(homo2)[names(homo2) == "dG WT"] <- "dG WT Homo"
names(homo2)[names(homo2) == "dG"] <- "dG Homo"
names(homo2)[names(homo2) == "ddG"] <- "ddG Homo"
names(homo2)[names(homo2) == "SASA"] <- "SASA Homo"
names(homo2)[names(homo2) == "ddG Binding"] <- "ddG Binding Homo"


hetero2 <- hetero %>%
  select(Pos, AA, SASA, `dG WT`, dG, ddG, `ddG Binding`, `ddG Mono`, `AA Property WT`, `AA Property`, `TE Category`, `Dimerization Category`)
names(hetero2)[names(hetero2) == "TE Category"] <- "TE Category Hetero"
names(hetero2)[names(hetero2) == "Dimerization Category"] <- "Dimerization Category Hetero"
names(hetero2)[names(hetero2) == "dG WT"] <- "dG WT Hetero"
names(hetero2)[names(hetero2) == "dG"] <- "dG Hetero"
names(hetero2)[names(hetero2) == "ddG"] <- "ddG Hetero"
names(hetero2)[names(hetero2) == "ddG Binding"] <- "ddG Binding Hetero"
names(hetero2)[names(hetero2) == "SASA"] <- "SASA Hetero"
homo_hetero_TE <- left_join(homo2, hetero2, by = c("Pos", "AA", "AA Property WT", "AA Property"))


homo_hetero_TE$`TE Status` <- "Unknown"
homo_hetero_TE$`TE Change` <- "Unknown"
homo_hetero_TE$`DA Status` <- "Unknown"
homo_hetero_TE$`DA Change` <- "Unknown"
homo_hetero_TE$`AA Status` <- "Unknown"
homo_hetero_TE$`AA Change` <- "Unknown"



for (i in 1:nrow(homo_hetero_TE)){
  if (homo_hetero_TE[i,"TE Category Homo"] == homo_hetero_TE[i,"TE Category Hetero"]){
    homo_hetero_TE[i,"TE Status"] <- "Same"
    homo_hetero_TE[i,"TE Change"] <- paste(unlist(homo_hetero_TE[i,"TE Category Hetero"]), " to ", unlist(homo_hetero_TE[i,"TE Category Homo"]))
  } else {
    homo_hetero_TE[i,"TE Status"] <- "Different"
    homo_hetero_TE[i,"TE Change"] <- paste(unlist(homo_hetero_TE[i,"TE Category Hetero"])," to ", unlist(homo_hetero_TE[i,"TE Category Homo"]))
  }
}

for (i in 1:nrow(homo_hetero_TE)){
  if (homo_hetero_TE[i,"Dimerization Category Homo"] == homo_hetero_TE[i,"Dimerization Category Hetero"]){
    homo_hetero_TE[i,"DA Status"] <- "Same"
    homo_hetero_TE[i,"DA Change"] <- paste(unlist(homo_hetero_TE[i,"Dimerization Category Hetero"]), " to ", unlist(homo_hetero_TE[i,"Dimerization Category Homo"]))
  } else {
    homo_hetero_TE[i,"DA Status"] <- "Different"
    homo_hetero_TE[i,"DA Change"] <- paste(unlist(homo_hetero_TE[i,"Dimerization Category Hetero"])," to ", unlist(homo_hetero_TE[i,"Dimerization Category Homo"]))
  }
}

for (i in 1:nrow(homo_hetero_TE)){
  if (homo_hetero_TE[i,"AA Property WT"] == homo_hetero_TE[i,"AA Property"]){
    homo_hetero_TE[i,"AA Status"] <- "Same"
    homo_hetero_TE[i,"AA Change"] <- paste(unlist(homo_hetero_TE[i,"AA Property WT"]), " to ", unlist(homo_hetero_TE[i,"AA Property"]))
  } else {
    homo_hetero_TE[i,"AA Status"] <- "Different"
    homo_hetero_TE[i,"AA Change"] <- paste(unlist(homo_hetero_TE[i,"AA Property WT"])," to ", unlist(homo_hetero_TE[i,"AA Property"]))
  }
}

homo_hetero_TE <- homo_hetero_TE %>%
  select(Pos, AA, `ddG Binding Homo`, `ddG Binding Hetero`, `dG WT Hetero`, `dG WT Homo`, `dG Hetero`, `dG Homo`, 
         `ddG Homo`, `ddG Hetero`, `ddG Mono`, `TE Status`, `TE Change`, `DA Status`, `DA Change`, 
         `AA Status`, `AA Change`, `SASA Homo`, `SASA Hetero`)
homo_hetero_TE2 <- homo_hetero_TE %>%
  select(Pos, AA, `TE Status`, `TE Change`, `DA Status`, `DA Change`, `AA Status`, `AA Change`)

dat <- left_join(dat, homo_hetero_TE2, by = c("Pos", "AA"))
dat$`TE Mono Dimer Status` <- as.factor(dat$`TE Mono Dimer Status`)
dat$`TE Mono Dimer` <- as.factor(dat$`TE Mono Dimer`)
dat$`TE Change` <- as.factor(dat$`TE Change`)
dat$`TE Status` <- as.factor(dat$`TE Status`)
dat$`DA Change` <- as.factor(dat$`DA Change`)
dat$`DA Status` <- as.factor(dat$`DA Status`)
dat$`AA Change` <- as.factor(dat$`AA Change`)
dat$`AA Status` <- as.factor(dat$`AA Status`)
dat$Status <- as.factor(dat$Status)
dat$Change <- as.factor(dat$Change)

for (i in 1:nrow(dat)){
  if (dat[i,"TE Mono Category"] == dat[i,"TE Category"]){
    dat[i,"TE Mono Dimer Status"] <- "Same"
    dat[i,"TE Mono Dimer"] <- paste(unlist(dat[i,"TE Mono Category"]), " to ", unlist(dat[i,"TE Category"]))
  } else {
    dat[i,"TE Mono Dimer Status"] <- "Different"
    dat[i,"TE Mono Dimer"] <- paste(unlist(dat[i,"TE Mono Category"])," to ", unlist(dat[i,"TE Category"]))
  }
}

for (i in 1:nrow(dat)){
  dat[i,"Change"] <- paste(unlist(dat[i,"TE Mono Category"])," to ", unlist(dat[i,"TE Change"]))
}

dat$`AA Property` <- as.factor(dat$`AA Property`)
dat$Region <- as.factor(dat$Region)
dat$Secondary <- as.factor(dat$Secondary)
dat$Location <- as.factor(dat$Location)
dat$`AA Property WT` <- as.factor(dat$`AA Property WT`)
                                         

#------------------PLOTTING (HEATMAP)-----------------------
# heatmap 6bof homozygous
homo.ddG <- homo[,c("AA","Pos","ddG")]

# transform data from long to wide to identify missing values
homo.wide <- reshape(homo.ddG, idvar = "AA", timevar = "Pos", direction = "wide")
homo.wide[is.na(homo.wide)] <- 0

# transform data back to long format
homo.long <- reshape(homo.wide, 
                     direction = "long",
                     varying = list(names(homo.wide)[2:169]),
                     v.names = "ddG",
                     idvar = "AA",
                     timevar = "Pos",
                     times = 2:169)

# transform to matrix (preferred input by ComplexHeatmap package)
homo.matrix <- acast(homo.long, AA~Pos, value.var = "ddG")

col_fun <- colorRamp2(c(-1.84, -0.92, -0.46, 0.46, 0.92, 1.84), 
                      RColorBrewer::brewer.pal(name = "RdBu", n = 6))

# stacked barplot annotation (below the heatmap)
annot <- HeatmapAnnotation(
  stability = anno_barplot(cbind(
    colSums(homo.matrix < -1.84), # H. stabilizing
    colSums(homo.matrix > -1.84 & homo.matrix < -0.92), # stabilizing
    colSums(homo.matrix > -0.92 & homo.matrix < -0.46), # S. stabilizing
    colSums(homo.matrix > -0.46 & homo.matrix < 0.46), #neutral
    colSums(homo.matrix > 0.46 & homo.matrix < 0.92), # S. destabilizing
    colSums(homo.matrix > 0.92 & homo.matrix < 1.84), #destabilizing
    colSums(homo.matrix > 1.84)), #H. destabilizing
    gp = gpar(fill = RColorBrewer::brewer.pal(name = "Paired", n = 7), 
              col = RColorBrewer::brewer.pal(name = "Paired", n = 7)),
    height = unit(4, "cm")),
  show_annotation_name = FALSE
)

# legend for barplot
lgd <- Legend(
  labels = c("H. Stabilizing", "Stabilizing", "S. Stabilizing", 
             "Neutral", "S. Stabilizing", 
             "Destabilizing", "H. Destabilizing"),
  title = "Category",
  legend_gp = gpar(fill = RColorBrewer::brewer.pal(name = "Paired", n = 7))
)

# simple bars annotation (above the heatmap)
bar.annotation <- data.frame(sasa.seq$Location, sasa.seq$Secondary, sasa.seq$Region)
colnames(bar.annotation) <- c("Solvent Exposure", "Secondary Structure", "Region")

# amino acid annotation 
ann.annotation <- data.frame(
  "Property" = c("Non-Polar", "Polar", "Acidic", "Acidic", "Non-Polar", 
              "Non-Polar", "Basic", "Non-Polar", "Basic", "Non-Polar",
              "Non-Polar", "Polar", "Non-Polar", "Polar", "Basic", "Polar",
              "Polar", "Non-Polar", "Non-Polar", "Polar")
)

# specify colors of above bars
colours <- list(
  "Solvent Exposure" = c("Exposed" = "gray", 
                 "Buried" = "#FB9A99", "Interacting" = "#E31A1C"),
  "Secondary Structure" = c("Alpha" = "#B2DF8A", "Beta" = "#33A02C"),
  "Region" = c("P-Loop" = "#FF7F00", "Switch I" = "#FDBF6F", 
                    "Switch II" = "#E31A1C","Base-Binding Loop" = "#FB9A99")
)

# specify colors of amino acid bar
colours2 <- list(
  "Property" = c("Non-Polar" = "#FB9A99", "Polar" = "#E31A1C", 
              "Basic" = "#FDBF6F", "Acidic" = "#FF7F00")
)

# build bars
annot2 <- HeatmapAnnotation(df = bar.annotation,
                            which = "col",
                            col = colours,
                            annotation_width = unit(c(1, 4), "cm"),
                            gap = unit(1, "mm"))
annot3 <- rowAnnotation(df = ann.annotation, col = colours2,
                        show_annotation_name = FALSE)


# build heatmap
hm <- Heatmap(
  homo.matrix, col = col_fun, name = "????G (kcal/mol)", 
  row_title = "Amino Acid Substitution", row_names_gp = gpar(fontsize = 7), 
  row_names_side = "left", show_row_dend = FALSE, 
  row_order = c("G", "A", "V", "L", "I", "M", "F", "W", "P", "S", 
                "T", "C", "Y", "N", "Q", "K", "R", "H", "D", "E"),
  column_title = "Mutational Effects on Stability of Homozygous Group",
  column_names_gp = gpar(fontsize = 6), column_names_rot = 45,
  column_order = order(as.numeric(colnames(homo.matrix))),
  width = ncol(homo.matrix)*unit(2.5, "mm"), 
  height = nrow(homo.matrix)*unit(4, "mm"),
  bottom_annotation = annot, top_annotation = annot2, left_annotation = annot3
)

draw(hm, ht_gap = unit(5, "mm"), heatmap_legend_list = list(lgd),
     annotation_legend_side = "bottom")

##################################################
# heatmap homozygous binding
homo.dimer <- homo[,c("AA","Pos","ddG Binding")]

# transform data from long to wide to identify missing values
homo.dimer.wide <- reshape(homo.dimer, idvar = "AA", timevar = "Pos", direction = "wide")
homo.dimer.wide[is.na(homo.dimer.wide)] <- 0

# transform data back to long format
homo.dimer.long <- reshape(homo.dimer.wide, 
                     direction = "long",
                     varying = list(names(homo.dimer.wide)[2:169]),
                     v.names = "ddG Binding",
                     idvar = "AA",
                     timevar = "Pos",
                     times = 2:169)

# transform to matrix (preferred input by ComplexHeatmap package)
homo.dimer.matrix <- acast(homo.dimer.long, AA~Pos, value.var = "ddG Binding")

# stacked barplot annotation (below the heatmap)
annot <- HeatmapAnnotation(
  stability = anno_barplot(cbind(
    colSums(homo.dimer.matrix < -1.84), # H. stabilizing
    colSums(homo.dimer.matrix > -1.84 & homo.dimer.matrix < -0.92), # stabilizing
    colSums(homo.dimer.matrix > -0.92 & homo.dimer.matrix < -0.46), # S. stabilizing
    colSums(homo.dimer.matrix > -0.46 & homo.dimer.matrix < 0.46), #neutral
    colSums(homo.dimer.matrix > 0.46 & homo.dimer.matrix < 0.92), # S. destabilizing
    colSums(homo.dimer.matrix > 0.92 & homo.dimer.matrix < 1.84), #destabilizing
    colSums(homo.dimer.matrix > 1.84)), #H. destabilizing
    gp = gpar(fill = RColorBrewer::brewer.pal(name = "Paired", n = 7), 
              col = RColorBrewer::brewer.pal(name = "Paired", n = 7)),
    height = unit(4, "cm")),
  show_annotation_name = FALSE
)


# build heatmap
hm <- Heatmap(
  homo.dimer.matrix, col = col_fun, name = "????G (kcal/mol)", 
  row_title = "Amino Acid Substitution", row_names_gp = gpar(fontsize = 7), 
  row_names_side = "left", show_row_dend = FALSE, 
  row_order = c("G", "A", "V", "L", "I", "M", "F", "W", "P", "S", 
                "T", "C", "Y", "N", "Q", "K", "R", "H", "D", "E"),
  column_title = "Mutational Effects on Dimerization Affinity of Homozygous Group",
  column_names_gp = gpar(fontsize = 6), column_names_rot = 45,
  column_order = order(as.numeric(colnames(homo.dimer.matrix))),
  width = ncol(homo.dimer.matrix)*unit(2.5, "mm"), 
  height = nrow(homo.dimer.matrix)*unit(4, "mm"),
  bottom_annotation = annot, top_annotation = annot2, left_annotation = annot3
)

draw(hm, ht_gap = unit(10, "mm"), heatmap_legend_list = list(lgd),
     annotation_legend_side = "bottom")


###############################################
# heatmap heterozygous 6bof
hetero.ddG <- hetero[,c("AA","Pos","ddG")]

# transform data from long to wide to identify missing values
hetero.wide <- reshape(hetero.ddG, idvar = "AA", timevar = "Pos", direction = "wide")
hetero.wide[is.na(hetero.wide)] <- 0

# transform data back to long format
hetero.long <- reshape(hetero.wide, 
                        direction = "long",
                        varying = list(names(hetero.wide)[2:169]),
                        v.names = "ddG",
                        idvar = "AA",
                        timevar = "Pos",
                        times = 2:169)

# transform to matrix (preferred input by ComplexHeatmap package)
hetero.matrix <- acast(hetero.long, AA~Pos, value.var = "ddG")


# stacked barplot annotation (below the heatmap)
annot <- HeatmapAnnotation(
  stability = anno_barplot(cbind(
    colSums(hetero.matrix < -1.84), # H. stabilizing
    colSums(hetero.matrix > -1.84 & hetero.matrix < -0.92), # stabilizing
    colSums(hetero.matrix > -0.92 & hetero.matrix < -0.46), # S. stabilizing
    colSums(hetero.matrix > -0.46 & hetero.matrix < 0.46), #neutral
    colSums(hetero.matrix > 0.46 & hetero.matrix < 0.92), # S. destabilizing
    colSums(hetero.matrix > 0.92 & hetero.matrix < 1.84), #destabilizing
    colSums(hetero.matrix > 1.84)), #H. destabilizing
    gp = gpar(fill = RColorBrewer::brewer.pal(name = "Paired", n = 7), 
              col = RColorBrewer::brewer.pal(name = "Paired", n = 7)),
    height = unit(4, "cm")),
  show_annotation_name = FALSE
)

# simple bars annotation (above the heatmap)
bar.annotation <- data.frame(sasa.seq$Location, sasa.seq$Secondary, sasa.seq$Region)
colnames(bar.annotation) <- c("Solvent Exposure", "Secondary Structure", "Region")

# build bars
annot2 <- HeatmapAnnotation(df = bar.annotation,
                            which = "col",
                            col = colours,
                            annotation_width = unit(c(1, 4), "cm"),
                            gap = unit(1, "mm"))
annot3 <- rowAnnotation(df = ann.annotation, col = colours2,
                        show_annotation_name = FALSE)


# build heatmap
hm <- Heatmap(
  hetero.matrix, col = col_fun, name = "????G (kcal/mol)", 
  row_title = "Amino Acid Substitution", row_names_gp = gpar(fontsize = 7), 
  row_names_side = "left", show_row_dend = FALSE, 
  row_order = c("G", "A", "V", "L", "I", "M", "F", "W", "P", "S", 
                "T", "C", "Y", "N", "Q", "K", "R", "H", "D", "E"),
  column_title = "Mutational Effects on Stability of Heterozygous Group",
  column_names_gp = gpar(fontsize = 6), column_names_rot = 45,
  column_order = order(as.numeric(colnames(hetero.matrix))),
  width = ncol(hetero.matrix)*unit(2.5, "mm"), 
  height = nrow(hetero.matrix)*unit(4, "mm"),
  bottom_annotation = annot, top_annotation = annot2, left_annotation = annot3
)

draw(hm, ht_gap = unit(5, "mm"), heatmap_legend_list = list(lgd),
     annotation_legend_side = "bottom")

##################################################
# heatmap heterozygous binding
hetero.dimer <- hetero[,c("AA","Pos","ddG Binding")]

# transform data from long to wide to identify missing values
hetero.dimer.wide <- reshape(hetero.dimer, idvar = "AA", timevar = "Pos", direction = "wide")
hetero.dimer.wide[is.na(hetero.dimer.wide)] <- 0

# transform data back to long format
hetero.dimer.long <- reshape(hetero.dimer.wide, 
                              direction = "long",
                              varying = list(names(hetero.dimer.wide)[2:169]),
                              v.names = "ddG Binding",
                              idvar = "AA",
                              timevar = "Pos",
                              times = 2:169)

# transform to matrix (preferred input by ComplexHeatmap package)
hetero.dimer.matrix <- acast(hetero.dimer.long, AA~Pos, value.var = "ddG Binding")

# stacked barplot annotation (below the heatmap)
annot <- HeatmapAnnotation(
  stability = anno_barplot(cbind(
    colSums(hetero.dimer.matrix < -1.84), # H. stabilizing
    colSums(hetero.dimer.matrix > -1.84 & hetero.dimer.matrix < -0.92), # stabilizing
    colSums(hetero.dimer.matrix > -0.92 & hetero.dimer.matrix < -0.46), # S. stabilizing
    colSums(hetero.dimer.matrix > -0.46 & hetero.dimer.matrix < 0.46), #neutral
    colSums(hetero.dimer.matrix > 0.46 & hetero.dimer.matrix < 0.92), # S. destabilizing
    colSums(hetero.dimer.matrix > 0.92 & hetero.dimer.matrix < 1.84), #destabilizing
    colSums(hetero.dimer.matrix > 1.84)), #H. destabilizing
    gp = gpar(fill = RColorBrewer::brewer.pal(name = "Paired", n = 7), 
              col = RColorBrewer::brewer.pal(name = "Paired", n = 7)),
    height = unit(4, "cm")),
  show_annotation_name = FALSE
)


# build heatmap
hm <- Heatmap(
  hetero.dimer.matrix, col = col_fun, name = "????G (kcal/mol)", 
  row_title = "Amino Acid Substitution", row_names_gp = gpar(fontsize = 7), 
  row_names_side = "left", show_row_dend = FALSE, 
  row_order = c("G", "A", "V", "L", "I", "M", "F", "W", "P", "S", 
                "T", "C", "Y", "N", "Q", "K", "R", "H", "D", "E"),
  column_title = "Mutational Effects on Dimerization Affinity of Heterozygous Group",
  column_names_gp = gpar(fontsize = 6), column_names_rot = 45,
  column_order = order(as.numeric(colnames(hetero.dimer.matrix))),
  width = ncol(hetero.dimer.matrix)*unit(2.5, "mm"), 
  height = nrow(hetero.dimer.matrix)*unit(4, "mm"),
  bottom_annotation = annot, top_annotation = annot2, left_annotation = annot3
)

draw(hm, ht_gap = unit(10, "mm"), heatmap_legend_list = list(lgd),
     annotation_legend_side = "bottom")

#--------------PLOTTING (PIE CHART)-------------------
# homo TE 6BOF
pie.homo <- data.frame(
  group = c("H. Stabilizing", "Stabilizing", "S. Stabilizing", "Neutral",
            "S. Destabilizing", "Destabilizing", "H. Destabilizing"),
  value = c(nrow(homo[homo$ddG < -1.84]),
            nrow(homo[homo$ddG > -1.84 & homo$ddG < -0.92]),
            nrow(homo[homo$ddG > -0.92 & homo$ddG < -0.46]),
            nrow(homo[homo$ddG > -0.46 & homo$ddG < 0.46]),
            nrow(homo[homo$ddG > 0.46 & homo$ddG < 0.92]),
            nrow(homo[homo$ddG > 0.92 & homo$ddG < 1.84]),
            nrow(homo[homo$ddG > 1.84]))
)
pie.homo <- pie.homo[order(-pie.homo$value),]

pie.homo <- pie.homo %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/1.6 + lead(csum, 1),
         pos = if_else(is.na(pos), value/1.6, pos))


ggplot(pie.homo, aes(x = "" , y = value, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Paired") +
  geom_label_repel(data = pie.homo,
                   aes(y = pos, label = paste0(value, " (", scales::percent(value / sum(value)), ")")),
                   size = 3.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Category")) +
  theme_void()

# homo binding 6BOF
pie.homo.binding <- data.frame(
  group = c("H. Stabilizing", "Stabilizing", "S. Stabilizing", "Neutral",
            "S. Destabilizing", "Destabilizing", "H. Destabilizing"),
  value = c(nrow(homo.dimer[homo.dimer$`ddG Binding` < -1.84]),
            nrow(homo.dimer[homo.dimer$`ddG Binding` > -1.84 & homo.dimer$`ddG Binding` < -0.92]),
            nrow(homo.dimer[homo.dimer$`ddG Binding` > -0.92 & homo.dimer$`ddG Binding` < -0.46]),
            nrow(homo.dimer[homo.dimer$`ddG Binding` > -0.46 & homo.dimer$`ddG Binding` < 0.46]),
            nrow(homo.dimer[homo.dimer$`ddG Binding` > 0.46 & homo.dimer$`ddG Binding` < 0.92]),
            nrow(homo.dimer[homo.dimer$`ddG Binding` > 0.92 & homo.dimer$`ddG Binding` < 1.84]),
            nrow(homo.dimer[homo.dimer$`ddG Binding` > 1.84]))
)
pie.homo.binding <- pie.homo.binding[order(-pie.homo.binding$value),]

pie.homo.binding <- pie.homo.binding %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/1.6 + lead(csum, 1),
         pos = if_else(is.na(pos), value/1.6, pos))


ggplot(pie.homo.binding, aes(x = "" , y = value, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Paired") +
  geom_label_repel(data = pie.homo.binding,
                   aes(y = pos, label = paste0(value, " (", scales::percent(value / sum(value)), ")")),
                   size = 3.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Category")) +
  theme_void()

# hetero TE 6BOF
pie.hetero <- data.frame(
  group = c("H. Stabilizing", "Stabilizing", "S. Stabilizing", "Neutral",
            "S. Destabilizing", "Destabilizing", "H. Destabilizing"),
  value = c(nrow(hetero[hetero$ddG < -1.84]),
            nrow(hetero[hetero$ddG > -1.84 & hetero$ddG < -0.92]),
            nrow(hetero[hetero$ddG > -0.92 & hetero$ddG < -0.46]),
            nrow(hetero[hetero$ddG > -0.46 & hetero$ddG < 0.46]),
            nrow(hetero[hetero$ddG > 0.46 & hetero$ddG < 0.92]),
            nrow(hetero[hetero$ddG > 0.92 & hetero$ddG < 1.84]),
            nrow(hetero[hetero$ddG > 1.84]))
)
pie.hetero <- pie.hetero[order(-pie.hetero$value),]

pie.hetero <- pie.hetero %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/1.6 + lead(csum, 1),
         pos = if_else(is.na(pos), value/1.6, pos))


ggplot(pie.hetero, aes(x = "" , y = value, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Paired") +
  geom_label_repel(data = pie.hetero,
                   aes(y = pos, label = paste0(value, " (", scales::percent(value / sum(value)), ")")),
                   size = 3.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Category")) +
  theme_void()

# hetero binding 6BOF
pie.hetero.binding <- data.frame(
  group = c("H. Stabilizing", "Stabilizing", "S. Stabilizing", "Neutral",
            "S. Destabilizing", "Destabilizing", "H. Destabilizing"),
  value = c(nrow(hetero.dimer[hetero.dimer$`ddG Binding` < -1.84]),
            nrow(hetero.dimer[hetero.dimer$`ddG Binding` > -1.84 & hetero.dimer$`ddG Binding` < -0.92]),
            nrow(hetero.dimer[hetero.dimer$`ddG Binding` > -0.92 & hetero.dimer$`ddG Binding` < -0.46]),
            nrow(hetero.dimer[hetero.dimer$`ddG Binding` > -0.46 & hetero.dimer$`ddG Binding` < 0.46]),
            nrow(hetero.dimer[hetero.dimer$`ddG Binding` > 0.46 & hetero.dimer$`ddG Binding` < 0.92]),
            nrow(hetero.dimer[hetero.dimer$`ddG Binding` > 0.92 & hetero.dimer$`ddG Binding` < 1.84]),
            nrow(hetero.dimer[hetero.dimer$`ddG Binding` > 1.84]))
)
pie.hetero.binding <- pie.hetero.binding[order(-pie.hetero.binding$value),]

pie.hetero.binding <- pie.hetero.binding %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/1.6 + lead(csum, 1),
         pos = if_else(is.na(pos), value/1.6, pos))


ggplot(pie.hetero.binding, aes(x = "" , y = value, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Paired") +
  geom_label_repel(data = pie.hetero.binding,
                   aes(y = pos, label = paste0(value, " (", scales::percent(value / sum(value)), ")")),
                   size = 3.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Category")) +
  theme_void()


dat.hetero <- dat[dat$Group == "hetero",]
dat.homo <- dat[dat$Group == "homo",]
#----------------------PLOTTING (CORRELATION)-----------------------
# hetero 6BOF


stability.binding1 <- ggplot(dat[dat$Group == "hetero"], aes(x = SASA, y = `ddG Binding`, fill = Allosteric)) +
  geom_point(size = 3.5, color = "black", 
             aes(fill = Allosteric), shape = 21) +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99", "#B2DF8A", "#E31A1C"),
                    name = "Solvent Exposure") +
  geom_smooth(method = lm, color = "black") +
  labs(x="\n SASA (Å)", y = "??????G Binding (kcal/mol) \n", tag = "A") +
  #ylim(min(dat.hetero$`ddG Binding`), max(dat.hetero$`ddG Binding`)) +
  # x = 30, y = 6
  annotate("text", x = 17600, y = 6, col = "black",
           label = paste("Spearman's rho\nAllosteric Buried = ", 
                         signif(cor(dat[dat$Group == "hetero" & dat$Allosteric == "Allosteric Buried",]$SASA, dat[dat$Group == "hetero" & dat$Allosteric == "Allosteric Buried",]$`ddG Binding`, method = "spearman"),2),
                         "\nAllosteric Exposed = ",  signif(cor(dat[dat$Group == "hetero" & dat$Allosteric == "Allosteric Exposed",]$SASA, dat[dat$Group == "hetero" & dat$Allosteric == "Allosteric Exposed",]$`ddG Binding`, method = "spearman"),2),
                         "\nInteracting = ",  signif(cor(dat[dat$Group == "hetero" & dat$Allosteric == "Interacting",]$SASA, dat[dat$Group == "hetero" & dat$Allosteric == "Interacting",]$`ddG Binding`, method = "spearman"),2),
                         "\nNon-Allosteric and Non-Interacting = ",  signif(cor(dat[dat$Group == "hetero" & dat$Allosteric == "Non-Allosteric & Non-Interacting",]$SASA, dat[dat$Group == "hetero" & dat$Allosteric == "Non-Allosteric & Non-Interacting",]$`ddG Binding`, method = "spearman"),2))) +
  theme_classic()

# homo 6bof
stability.binding2 <- ggplot(dat[dat$Group == "homo"], aes(x = SASA, y = `ddG Binding`, fill = Allosteric)) +
  geom_point(size = 3.5, color = "black", 
             aes(fill = Allosteric), shape = 21) +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99", "#B2DF8A", "#E31A1C"),
                    name = "Solvent Exposure") +
  geom_smooth(method = lm, color = "black") +
  labs(x="\n SASA (Å)", y = "??????G Binding (kcal/mol) \n", tag = "B") +
  #ylim(min(dat.homo$`ddG Binding`), max(dat.homo$`ddG Binding`)) +
  # x = 60, y = 10
  annotate("text", x = 17700, y = 8, col = "black",
           label = paste("Spearman's rho\nAllosteric Buried = ", 
                         signif(cor(dat[dat$Group == "homo" & dat$Allosteric == "Allosteric Buried",]$SASA, dat[dat$Group == "homo" & dat$Allosteric == "Allosteric Buried",]$`ddG Binding`, method = "spearman"),2),
                         "\nAllosteric Exposed = ",  signif(cor(dat[dat$Group == "homo" & dat$Allosteric == "Allosteric Exposed",]$SASA, dat[dat$Group == "homo" & dat$Allosteric == "Allosteric Exposed",]$`ddG Binding`, method = "spearman"),2),
                         "\nInteracting = ",  signif(cor(dat[dat$Group == "homo" & dat$Allosteric == "Interacting",]$SASA, dat[dat$Group == "homo" & dat$Allosteric == "Interacting",]$`ddG Binding`, method = "spearman"),2),
                         "\nNon-Allosteric and Non-Interacting = ",  signif(cor(dat[dat$Group == "homo" & dat$Allosteric == "Non-Allosteric & Non-Interacting",]$SASA, dat[dat$Group == "homo" & dat$Allosteric == "Non-Allosteric & Non-Interacting",]$`ddG Binding`, method = "spearman"),2))) +
  theme_classic()

grid.arrange(stability.binding1, stability.binding2)

# homo 6BOF 
cor2 <- ggplot(dat.homo, aes(x = ddG, y = `ddG Binding`, fill = Location)) +
  geom_point(size = 3.5, color = "black", 
             aes(fill = Location), shape = 21) +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99", "#B2DF8A", "#E31A1C"),
                    name = "Solvent Exposure") +
  geom_smooth(method = lm, color = "black") +
  labs(x="\n Total Energy", y = "Binding Energy \n", tags = "B") +
  ylim(min(dat.homo$`ddG Binding`), max(dat.homo$`ddG Binding`)) +
  #geom_text_repel(data = . %>% filter(`ddG Binding` > 0.92), aes(label = Pos)) +
  annotate("text", x = 90, y = 10, col = "black",
           label = paste("Spearman's rho\nOverall = ", 
                         signif(cor(dat.homo$ddG, dat.homo$`ddG Binding`, method = "spearman"),2),
                         "\nInteracting = ",  signif(cor(dat.homo[dat.homo$Location == "Interacting",]$ddG, dat.homo[dat.homo$Location == "Interacting",]$`ddG Binding`, method = "spearman"),2))) +
  theme_classic()

grid.arrange(cor1, cor2)

#------------DIMERIZATION AFFINITY-------------------
# number of mutations that do not affect binding (2713 mutations)
nrow(dat.hetero[dat.hetero$`DA Change` == "Neutral  to  Neutral",])

# mutations that affect binding
affect.binding <- dat.hetero[!(dat.hetero$`DA Change` == "Neutral  to  Neutral"),]
ggplot(affect.binding, aes(x = SASA, y = `ddG Binding`, color = Location)) +
  geom_point() +
  geom_smooth(method = lm)


# allosteric mutations that need both strands
allosteric <- affect.binding[affect.binding$Location != "Interacting",]

# allosteric sites in protein interior
buried.allosteric <- allosteric[allosteric$Location =="Buried",]
buried.allosteric.rec <- buried.allosteric[buried.allosteric$`Dimerization Category` == "Neutral",]
buried.allosteric.non.rec <- buried.allosteric[buried.allosteric$`Dimerization Category` != "Neutral",]
buried.allosteric.overdominant <- dat.homo[dat.homo$Location == "Buried" & dat.homo$`DA Change` != "Neutral  to  Neutral" &
                                        dat.homo$`Dimerization Category` == "Neutral",]



exposed.allosteric <- allosteric[allosteric$Location == "Exposed"]
exposed.allosteric.rec <- exposed.allosteric[exposed.allosteric$`Dimerization Category` == "Neutral",]
exposed.allosteric.non.rec <- exposed.allosteric[exposed.allosteric$`Dimerization Category` != "Neutral",]
exposed.allosteric.overdominant <- dat.homo[dat.homo$Location == "Exposed" & dat.homo$`DA Change` != "Neutral  to  Neutral" &
                                             dat.homo$`Dimerization Category` == "Neutral",]

interacting <- affect.binding[affect.binding$Location == "Interacting",]
interacting.rec <- interacting[interacting$`Dimerization Category` == "Neutral",]
interacting.non.rec <- interacting[interacting$`Dimerization Category` != "Neutral",]

dat.hetero$Allosteric <- "Non-Allosteric & Non-Interacting"
dat.hetero[dat.hetero$Location == "Buried" & dat.hetero$`Dimerization Category` != "Neutral",]$Allosteric <- "Allosteric Buried"
dat.hetero[dat.hetero$Location == "Exposed" & dat.hetero$`Dimerization Category` != "Neutral",]$Allosteric <- "Allosteric Exposed"
dat.hetero[dat.hetero$Location == "Interacting" & dat.hetero$`Dimerization Category` != "Neutral",]$Allosteric <- "Interacting"

dat.homo$Allosteric <- "Non-Allosteric & Non-Interacting"
dat.homo[dat.homo$Location == "Buried" & dat.homo$`Dimerization Category` != "Neutral",]$Allosteric <- "Allosteric Buried"
dat.homo[dat.homo$Location == "Exposed" & dat.homo$`Dimerization Category` != "Neutral",]$Allosteric <- "Allosteric Exposed"
dat.homo[dat.homo$Location == "Interacting" & dat.homo$`Dimerization Category` != "Neutral",]$Allosteric <- "Interacting"

dat$Allosteric <- "Non-Allosteric & Non-Interacting"
dat[dat$Location == "Buried" & dat$`DA Change` != "Neutral  to  Neutral",]$Allosteric <- "Allosteric Buried"
dat[dat$Location == "Exposed" & dat$`DA Change` != "Neutral  to  Neutral",]$Allosteric <- "Allosteric Exposed"
dat[dat$Location == "Interacting" & dat$`DA Change` != "Neutral  to  Neutral",]$Allosteric <- "Interacting"

itr.comp.opp <- rbind(
  interacting.non.rec[interacting.non.rec$`DA Change` == "Destabilizing  to  S. Stabilizing",],
  interacting.non.rec[interacting.non.rec$`DA Change` == "S. Destabilizing  to  S. Stabilizing",],
  interacting.non.rec[interacting.non.rec$`DA Change` == "S. Stabilizing  to  Destabilizing",],
  interacting.non.rec[interacting.non.rec$`DA Change` == "S. Stabilizing  to  S. Destabilizing",]
)


pai <- ggplot(dat.hetero[!(dat.hetero$Allosteric == "Non-Allosteric & Non-Interacting"),], 
       aes(x = Allosteric, y = SASA)) +
  geom_violin(fill = "#A6CEE3") +
  geom_boxplot(fill = "gray", width = 0.1) +
  stat_compare_means(comparisons = list(c("Allosteric Buried", "Allosteric Exposed"), c("Allosteric Buried", "Interacting"), c("Allosteric Exposed", "Interacting")), method = "wilcox.test") +
  labs(x = "", y = "??????G Binding\n", tag = "A") +
  theme_classic()

pai2 <- ggplot(dat.homo[!(dat.homo$Allosteric == "Non-Allosteric & Non-Interacting"),], 
       aes(x = Allosteric, y = SASA)) +
  geom_violin(fill = "#A6CEE3") +
  geom_boxplot(fill = "gray", width = 0.1) +
  stat_compare_means(comparisons = list(c("Allosteric Buried", "Allosteric Exposed"),
                                        c("Allosteric Buried", "Interacting"),
                                        c("Allosteric Exposed", "Interacting")), method = "wilcox.test") +
  labs(x = "", y = "??????G Binding\n", tag = "B") +
  theme_classic()

grid.arrange(pai, pai2)

ggplot(dat, aes(x = `AA Change`, y = `ddG Binding`, fill = Group)) +
  geom_boxplot(position = position_dodge(1), colour = "black", outlier.shape = NA) +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99"),
                    labels = c("hetero" = "Heterozygous", "homo" = "Homozygous")) + 
  labs(x = "", y = "??????G Binding (kcal/mol)\n") +
  facet_wrap(vars(Allosteric)) +
  coord_cartesian(ylim = c(-5,5)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_classic() 


x1 <- ggplot(dat[dat$Group == "hetero",], aes(x = ddG, y = `ddG Binding`, fill = `AA Property WT`)) +
  geom_point(color = "gray", size = 2) +
  geom_smooth(method = lm, color = "black") +
  facet_wrap(vars(Allosteric )) +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99", "#B2DF8A", "#E31A1C")) +
  labs(x = "\nSASA (Å)", y = "??????G Binding (kcal/mol)\n", tag = "A") +
  theme_light()

x2 <- ggplot(dat[dat$Group == "homo",], aes(x = ddG, y = `ddG Binding`, fill = `AA Property WT`)) +
  geom_point(color = "gray", size = 2) +
  geom_smooth(method = lm, color = "black") +
  facet_wrap(vars(Allosteric )) +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99", "#B2DF8A", "#E31A1C")) +
  labs(x = "\nSASA (Å)", y = "??????G Binding (kcal/mol)\n", tag = "B") +
  theme_light()

grid.arrange(x1, x2)

# check summary of allosteric residues
summary(dat.homo[dat.homo$Allosteric == "Allosteric" & dat.homo$`Dimerization Category` != "Neutral",])

# check summary of interacting residues
summary(dat.homo[dat.homo$Allosteric == "Interacting" & dat.homo$`Dimerization Category` != "Neutral",])

# Number of mutations located in allosteric buried residues that affect binding
unique(dat.homo[dat.homo$Location == "Buried" & dat.homo$`Dimerization Category` != "Neutral",]$Pos)

table(dat.homo[dat.homo$Location == "Buried" & dat.homo$`Dimerization Category` != "Neutral",]$`AA`)

nrow(dat.homo[dat.homo$Location == "Buried" & dat.homo$`Dimerization Category` != "Neutral" & dat.homo$`AA Property WT` == "Non-Polar",])


ggplot(buried.allosteric[buried.allosteric$`AA Property WT` == "Non-Polar",],
       aes(x = `AA Change`, y = `ddG Binding`)) +
  geom_violin(fill = "#A6CEE3") +
  geom_boxplot(fill = "gray", width = 0.1) +
  theme_classic()


        

ggplot(dat.homo[dat.homo$Location == "Buried" & dat.homo$`Dimerization Category` != "Neutral" & dat.homo$`AA Property WT` == "Non-Polar",],
       aes(x = `AA Change`, fill = `DA Change`)) +
  geom_bar(position = dodge)



ba <- ggplot(buried.allosteric, aes(x = fct_rev(fct_infreq(`DA Change`)))) +
  geom_bar(fill = "#A6CEE3", colour = "black") +
  labs(x = "", y = "\n Number of Mutations \n") +
  geom_text(stat='count', aes(label=..count..), hjust = -0.3) +
  labs(tag = "A", x = "") +
  coord_flip() +
  theme_classic() 

ea <- ggplot(exposed.allosteric, aes(x = fct_rev(fct_infreq(`DA Change`)))) +
  geom_bar(fill = "#A6CEE3", colour = "black") +
  labs(x = "", y = "\n Number of Mutations \n") +
  geom_text(stat='count', aes(label=..count..), hjust = -0.3) +
  labs(tag = "B", x = "") +
  coord_flip() +
  theme_classic() 

itr <- ggplot(interacting, aes(x = fct_rev(fct_infreq(`DA Change`)))) +
  geom_bar(fill = "#A6CEE3", colour = "black") +
  labs(x = "", y = "\n Number of Mutations \n") +
  geom_text(stat='count', aes(label=..count..), hjust = -0.3) +
  labs(tag = "C") +
  coord_flip() +
  theme_classic() 

grid.arrange(ba, ea, itr)

ggplot(dat.hetero[!(dat.hetero$`DA Change` == "Neutral  to  Neutral"),], aes(x = fct_rev(fct_infreq(`DA Change`)))) +
  geom_bar(fill = "#A6CEE3", colour = "black") +
  labs(x = "", y = "\n Number of Mutations \n") +
  geom_text(stat='count', aes(label=..count..), hjust = -0.3) +
  coord_flip() +
  theme_classic() 



#-------------------STABILITY EFFECTS------------------
same.effects <- rbind(
  dat.hetero[dat.hetero$`TE Change` == "H. Destabilizing  to  H. Destabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "Destabilizing  to  Destabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "S. Destabilizing  to  S. Destabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "S. Stabilizing  to  S. Stabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "Stabilizing  to  Stabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "H. Stabilizing  to  H. Stabilizing",]
)

stronger.effects <- rbind(
  dat.hetero[dat.hetero$`TE Change` == "Destabilizing  to  H. Destabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "S. Destabilizing  to  Destabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "S. Stabilizing  to  Stabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "Stabilizing  to  H. Stabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "S. Destabilizing  to  H. Destabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "S. Stabilizing  to  H. Stabilizing",]
)

weaker.effects <- rbind(
  dat.hetero[dat.hetero$`TE Change` == "S. Destabilizing  to  Neutral",],
  dat.hetero[dat.hetero$`TE Change` == "S. Stabilizing  to  Neutral",],
  dat.hetero[dat.hetero$`TE Change` == "Destabilizing  to  S. Destabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "Stabilizing  to  S. Stabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "Destabilizing  to  Neutral",],
  dat.hetero[dat.hetero$`TE Change` == "H. Destabilizing  to  Destabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "H. Stabilizing  to  Neutral",],
  dat.hetero[dat.hetero$`TE Change` == "Stabilizing  to  Neutral",]
)

neutral.stronger <- rbind(
  dat.hetero[dat.hetero$`TE Change` == "Neutral  to  S. Destabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "Neutral  to  S. Stabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "Neutral  to  Destabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "Neutral  to  Stabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "Neutral  to  H. Destabilizing",],
  dat.hetero[dat.hetero$`TE Change` == "Neutral  to  H. Stabilizing",]
)

neutral <- dat.hetero[dat.hetero$`TE Change` == "Neutral  to  Neutral",]

mono.dimer.diff <- dat.hetero[dat.hetero$`TE Mono Dimer Status` == "Different",]
dodge <- position_dodge(width = 0.6)

# monomer vs dimer
ggplot(mono.dimer.diff, aes(x = fct_rev(fct_infreq(`TE Mono Dimer`)))) +
  geom_bar(fill = "#A6CEE3", colour = "black") +
  labs(x = "", y = "\n Number of Mutations \n") +
  geom_text(stat='count', aes(label=..count..), hjust = -0.3) +
  coord_flip() +
  theme_classic() 

ggplot(dat, aes(x = `TE Mono Dimer Status`, y = ddG)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.09, fill = "gray",
               outlier.shape = NA, position = dodge) +
  #coord_cartesian(ylim = quantile(dat.violin.WT$ddG, c(0, 0.985))) +
  facet_wrap(vars(Group)) +
  labs(x = "", y = "????G (kcal/mol)\n") +
  stat_compare_means() +
  theme_classic()


# homo hetero
mono.hetero <- dat.violin[dat.violin$Category == "Mutant" | dat.violin$Category == "Mutant-WT",]
mono.hetero.annot <- dat.hetero[,c("Pos", "AA", "TE Mono Dimer Status"),]
mono.hetero <- left_join(mono.hetero, mono.hetero.annot, by = c("Pos", "AA"))

mono.homo <- dat.violin[dat.violin$Category == "Mutant" | dat.violin$Category == "Mutant-Mutant",]
mono.homo.annot <- dat.homo[,c("Pos", "AA", "TE Mono Dimer Status"),]
mono.homo <- left_join(mono.homo, mono.homo.annot, by = c("Pos", "AA"))

pmn <- ggplot(mono.hetero, aes(x = Category, y = ddG)) +
          geom_violin(position = dodge, fill = "#A6CEE3") +
          geom_boxplot(outlier.shape = NA, width = 0.1, position = dodge,
                       show.legend = FALSE, fill = "gray") +
          labs(x = "", y = "????G (kcal/mol)\n", tag = "A") +
          facet_wrap(vars(`TE Mono Dimer Status`)) +
          #coord_cartesian(ylim = quantile(dat.violin.WT$ddG, c(0.01, 0.95))) +
          stat_compare_means(comparisons = list(c("Mutant", "Mutant-WT")), method = "wilcox.test", paired = TRUE) +
          theme_classic()

pmn2 <- ggplot(mono.homo, aes(x = Category, y = ddG)) +
          geom_violin(position = dodge, fill = "#A6CEE3") +
          geom_boxplot(outlier.shape = NA, width = 0.1, position = dodge,
                       show.legend = FALSE, fill = "gray") +
          labs(x = "", y = "????G (kcal/mol)\n", tag = "B") +
          facet_wrap(vars(`TE Mono Dimer Status`)) +
          #coord_cartesian(ylim = quantile(dat.violin.WT$ddG, c(0.005, 0.984))) +
          stat_compare_means(comparisons = list(c("Mutant", "Mutant-Mutant")), method = "wilcox.test") +
          theme_classic()

pmn3 <- ggplot(mono.hetero, aes(x = `TE Mono Dimer Status`, y = ddG)) +
  geom_violin(position = dodge, fill = "#A6CEE3") +
  geom_boxplot(outlier.shape = NA, width = 0.1, position = dodge,
               show.legend = FALSE, fill = "gray") +
  labs(x = "", y = "????G (kcal/mol)\n", tag = "B") +
  facet_wrap(vars(Category)) +
  #coord_cartesian(ylim = quantile(dat.violin.WT$ddG, c(0.005, 0.984))) +
  stat_compare_means(comparisons = list(c("Different", "Same")), method = "wilcox.test") +
  theme_classic()

grid.arrange(pmn, pmn2)

#------------------BOXPLOT--------------------
# monomer vs dimer
dat.violin <- dat %>%
  select(Group, Pos, AA, ddG, `ddG Mono`)
dat.violin <- reshape(dat.violin, idvar = c("Pos", "AA", "ddG Mono"), timevar = "Group", direction = "wide")
dat.violin <- gather(dat.violin, Category, ddG, `ddG Mono`:ddG.hetero, factor_key = TRUE)
levels(dat.violin$Category)[1] <- "Mutant"
levels(dat.violin$Category)[2] <- "Mutant-Mutant"
levels(dat.violin$Category)[3] <- "Mutant-WT"
dat.violin$Category <- factor(dat.violin$Category, levels = c("Mutant", "Mutant-WT", "Mutant-Mutant"))

dat.violin.WT <- dat.WT %>%
  select(Group, Pos, AA, ddG, `ddG Mono`)
dat.violin.WT <- reshape(dat.violin.WT, idvar = c("Pos", "AA", "ddG Mono"), timevar = "Group", direction = "wide")
dat.violin.WT <- gather(dat.violin.WT, Category, ddG, `ddG Mono`:ddG.hetero, factor_key = TRUE)
levels(dat.violin.WT$Category)[1] <- "Mutant"
levels(dat.violin.WT$Category)[2] <- "Mutant-Mutant"
levels(dat.violin.WT$Category)[3] <- "Mutant-WT"
dat.violin.WT$Category <- factor(dat.violin.WT$Category, levels = c("Mutant", "Mutant-WT", "Mutant-Mutant"))

dat2 <- dat.hetero %>%
  select (Pos, AA, Location, Region, `AA Change`, Secondary, Status)
dat.violin <- left_join(dat.violin, dat2, by = c("Pos", "AA"))

dat2.WT <- dat.hetero.WT %>%
  select (Pos, AA, Location, Region, `AA Change`, Secondary, Status)
dat.violin.WT <- left_join(dat.violin.WT, dat2.WT, by = c("Pos", "AA"))
dodge <- position_dodge(width = 0.6)

# monomer vs dimer
ggplot(dat.hetero, aes(x = fct_rev(fct_infreq(`Change`)))) +
  geom_bar(fill = "#A6CEE3", colour = "black") +
  labs(x = "", y = "\n Number of Mutations \n") +
  geom_text(stat='count', aes(label=..count..), hjust = -0.3) +
  coord_flip() +
  theme_classic() 

ggplot(dat.hetero, aes(x = fct_rev(fct_infreq()))) +
  geom_bar(fill = "#A6CEE3", colour = "black") +
  labs(x = "", y = "\n Number of Mutations \n") +
  geom_text(stat='count', aes(label=..count..), hjust = -0.3) +
  coord_flip() +
  theme_classic() 


ggplot(dat, aes(x = `TE Mono Dimer Status`, fill = Group)) +
  geom_bar(position = "dodge") +
  geom_text(stat='count', aes(label=..count..), vjust = -1)

# monomer vs dimer (WT)
ggplot(dat.violin[dat.violin$Category == "Mutant",], aes(x = ddG)) +
  geom_histogram(binwidth = 1) 

hist(dat.violin[dat.violin$Category == "Mutant",]$ddG, xlim = c(-5,20))
hist(dat.violin[dat.violin$Category == "Mutant-WT",]$ddG, xlim = c(-5,50))
hist(dat.violin[dat.violin$Category == "Mutant-Mutant",]$ddG, xlim = c(-10,50))




ggplot(dat.violin, aes(x = Category, y = ddG, fill = Category)) + 
  geom_violin(trim = TRUE) +
  labs(x = "", y = "????G (kcal/mol)\n") +
  geom_boxplot(width = 0.09, fill = "gray",
               outlier.shape = NA) +
  stat_compare_means(comparisons = list(c("Mutant", "Mutant-WT"), 
                                        c("Mutant-WT", "Mutant-Mutant"),
                                        c("Mutant", "Mutant-Mutant")),
                     method = "wilcox.test", paired = TRUE) +
  #coord_cartesian(ylim = quantile(dat.violin.WT$ddG, c(0, 0.985))) +
  theme_classic()

ggplot(dat, aes(x = `TE Status`, y = ddG, fill = Group)) +
  geom_boxplot()


# Correlate SASA with ddG (AA Property WT)
ggplot(dat, aes(x = SASA, y = ddG, fill = `AA Property WT`)) +
  geom_point(size = 3.5, color = "black", 
             aes(fill = `AA Property WT`), shape = 21) +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99", "#B2DF8A", "#E31A1C"),
                    name = "Property") +
  geom_smooth(method = lm, color = "black") +
  labs(x="\n SASA (Angstrom)", y = "????G Total Energy (kcal/mol) \n") +
  annotate("text", x = 17750, y = 100, col = "black",
           label = paste("Spearman's rho\nOverall = ", 
                         signif(cor(dat$SASA, dat$ddG, method = "spearman"),2),
                         "\nNon-Polar = ",  signif(cor(dat[dat$`AA Property WT` == "Non-Polar",]$SASA, dat[dat$`AA Property WT` == "Non-Polar",]$ddG, method = "spearman"),2))) +
  theme_classic()

# Correlate SASA with ddG (Solvent Exposure)
ggplot(dat, aes(x = SASA, y = ddG, fill = Location)) +
  geom_point(size = 3.5, color = "black", 
             aes(fill = Location), shape = 21) +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99", "#B2DF8A", "#E31A1C"),
                    name = "Property") +
  geom_smooth(method = lm, color = "black") +
  labs(x="\n SASA (Angstrom)", y = "????G Total Energy (kcal/mol) \n") +
  annotate("text", x = 17750, y = 100, col = "black",
           label = paste("Spearman's rho\nOverall = ", 
                         signif(cor(dat$SASA, dat$ddG, method = "spearman"),2),
                         "\nNon-Polar = ",  signif(cor(dat[dat$Location == "Buried",]$SASA, dat[dat$Location == "Buried",]$ddG, method = "spearman"),2))) +
  theme_classic()


# Heterozygous vs Homozygous
ggplot(dat, aes(x = `AA Change`, y = ddG, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  coord_cartesian(ylim = c(min(dat$ddG), 27)) +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99"),
                    labels = c("hetero" = "Heterozygous", "homo" = "Homozygous")) + 
  labs(x = "", y = "????G (kcal/mol)\n") +
  theme_classic()

ggplot(dat.violin, aes(x = Category, y = ddG)) +
  geom_violin(position = dodge, fill = "#A6CEE3") +
  geom_boxplot(outlier.shape = NA, width = 0.1, position = dodge,
               show.legend = FALSE, fill = "gray") +
  labs(x = "", y = "????G (kcal/mol)\n") +
  facet_wrap(vars(Location)) +
  stat_compare_means(comparisons = list(c("Mutant", "Mutant-WT"),
                                        c("Mutant-WT", "Mutant-Mutant"),
                                        c("Mutant", "Mutant-Mutant")),
                     method = "wilcox.test") +
  theme_classic()


ggplot(dat.violin, aes(x = Category, y = ddG)) +
  geom_violin(position = dodge, fill = "#A6CEE3") +
  geom_boxplot(outlier.shape = NA, width = 0.1, position = dodge,
               show.legend = FALSE, fill = "gray") +
  labs(x = "", y = "????G Folding (kcal/mol)\n") +
  stat_compare_means(comparisons = list(c("Mutant", "Mutant-WT"),
                                        c("Mutant-WT", "Mutant-Mutant"),
                                        c("Mutant", "Mutant-Mutant")),
                     method = "wilcox.test") +
  coord_cartesian(ylim = quantile(dat.violin.WT$ddG, c(0, 0.985))) +
  theme_classic()


# compare oncogenic and normal
onco <- dat[dat$Status == "Oncogenic",]
ggplot(dat, aes(x = Status, y = SASA)) +
  geom_violin(fill = "#A6CEE3") +
  geom_boxplot(fill = "gray", width = 0.2, outlier.shape = NA) +
  facet_wrap(vars(Group), labeller = labeller(Group = c("hetero" = "Heterozygous", "homo" = "Homozygous"))) +
  stat_compare_means(method = "wilcox.test") +
  labs(x = "", y = "SASA (Å)\n") +
  theme_classic()

ggplot(dat.violin.WT, aes(x = Status, y = ddG)) +
  geom_violin(position = dodge) +
  geom_boxplot(width = 0.1, outlier.shape = NA, position = dodge) +
  coord_cartesian(ylim = c(min(dat.violin.WT$ddG), 18)) +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99", "#B2DF8A", "#E31A1C")) +
  facet_wrap(vars(Category)) +
  theme_classic()


dat.anova <- dat[,c(1,6,4)]
wt <- data.frame(Group = "WT",
                 dG = -32.63,
                 Binding = -4.15)
dat.anova <- rbind(dat.anova, wt)
my_comparisons <- list( c("hetero", "homo"), c("homo", "WT"), c("hetero", "WT") )
ggplot(dat.anova, aes(x = Group, y = Binding)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons) +
  theme_classic()




#----------------------Compare Homozygous and Heterozygous--------------------
homo_hetero_TE$`Mean dG WT` <- (homo_hetero_TE$`dG WT Hetero` + homo_hetero_TE$`dG WT Homo`)/2
homo_hetero_TE$`Half Pheno` <- (homo_hetero_TE$`dG Homo` + homo_hetero_TE$`Mean dG WT`)/2
homo_hetero_TE$`Half Homo` <- homo_hetero_TE$`dG Homo`/2
homo_hetero_TE$`Half ddG Homo` <- homo_hetero_TE$`ddG Homo`/2
homo_hetero_TE$`Half ddG Binding Homo` <- homo_hetero_TE$`ddG Binding Homo`/2
homo_hetero_TE <- homo_hetero_TE %>%
  select(Pos, AA, `SASA Homo`, `SASA Hetero`, `Mean dG WT`, `Half Pheno`, `Half Homo`, `dG Hetero`, `dG Homo`, `ddG Mono`, `ddG Homo`, `Half ddG Homo`, `ddG Hetero`, `ddG Binding Hetero`, `ddG Binding Homo`, `Half ddG Binding Homo`)
homo_hetero_TE$`Diff Hetero HPheno`<- homo_hetero_TE$`dG Hetero` - homo_hetero_TE$`Half Pheno`
homo_hetero_TE$`Diff Hetero Homo` <- homo_hetero_TE$`dG Hetero` - homo_hetero_TE$`dG Homo`
dat.hetero2 <- dat.hetero[,c("Pos", "AA", "TE Change", "DA Change", "Location", "Region", "AA Change", "Secondary", "Status")]
homo_hetero_TE <- left_join(homo_hetero_TE, dat.hetero2, by = c("Pos", "AA"))


strong <- homo_hetero_TE[homo_hetero_TE$`TE Change` == "H. Destabilizing  to  H. Destabilizing" | 
                         homo_hetero_TE$`TE Change` == "H. Stabilizing  to  H. Stabilizing" | 
                         homo_hetero_TE$`TE Change` == "Destabilizing  to  H. Destabilizing" |
                         homo_hetero_TE$`TE Change` == "Stabilizing  to  H. Stabilizing" |
                         homo_hetero_TE$`TE Change` == "S. Destabilizing  to  Destabilizing" | 
                         homo_hetero_TE$`TE Change` == "S. Stabilizing  to  Stabilizing",]

moderate <- homo_hetero_TE[!(homo_hetero_TE$`TE Change` == "H. Destabilizing  to  H. Destabilizing" | 
                           homo_hetero_TE$`TE Change` == "H. Stabilizing  to  H. Stabilizing" | 
                           homo_hetero_TE$`TE Change` == "Destabilizing  to  H. Destabilizing" |
                           homo_hetero_TE$`TE Change` == "Stabilizing  to  H. Stabilizing" |
                           homo_hetero_TE$`TE Change` == "S. Destabilizing  to  Destabilizing" | 
                           homo_hetero_TE$`TE Change` == "S. Stabilizing  to  Stabilizing" |
                           homo_hetero_TE$`TE Change` == "Neutral  to  Neutral"),]

complete.opp <- homo_hetero_TE[homo_hetero_TE$`TE Change` == "Destabilizing  to  Stabilizing" | 
                                   homo_hetero_TE$`TE Change` == "S. Destabilizing  to  S. Stabilizing" | 
                                   homo_hetero_TE$`TE Change` == "S. Destabilizing  to  Stabilizing" |
                                   homo_hetero_TE$`TE Change` == "S. Stabilizing  to  H. Destabilizing" |
                                   homo_hetero_TE$`TE Change` == "Stabilizing  to  Destabilizing" | 
                                   homo_hetero_TE$`TE Change` == "Stabilizing  to  S. Destabilizing",]



# correlate half pheno with dG Hetero
h1 <- ggplot(strong, aes(x = `ddG Binding Hetero`, y = `Half ddG Binding Homo`, fill = `TE Change`)) + 
  geom_point(size = 3.5, color = "black", 
             shape = 21) +
  geom_smooth(method = lm, color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.direction = "horizontal", legend.position = "bottom") +
  guides(fill=guide_legend(title="")) +
  labs(x = "Heterozygous ????G Folding", y = "Homozygous 1/2 ????G Folding\n", tag = "A") 

h2 <- ggplot(moderate, aes(x = `ddG Binding Hetero`, y = `Half ddG Binding Homo`, fill = `TE Change`)) + 
  geom_point(size = 3.5, color = "black", 
             shape = 21) +
  geom_smooth(method = lm, color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.direction = "horizontal", legend.position = "bottom") +
  guides(fill=guide_legend(title="")) +
  labs(x = "\nHeterozygous ????G Folding", y = "Homozygous  1/2 ????G Folding\n", tag = "B") 

grid.arrange(h1, h2)


dimerization <- homo_hetero_TE[!(homo_hetero_TE$`DA Change` == "Neutral  to  Neutral"),]
dimerization.allosteric <- dimerization[!(dimerization$Location == "Interacting"),]
dimerization.interacting <- dimerization[dimerization$Location == "Interacting",]

# correlate hetero binding and half homo binding
h1 <- ggplot(dimerization[dimerization$`DA Change` == "H. Stabilizing  to  H. Stabilizing" | dimerization$`DA Change` == "Destabilizing  to  H. Destabilizing" | dimerization$`DA Change` == "S. Destabilizing  to  H. Destabilizing",], aes(x = `ddG Binding Hetero`, y = `Half ddG Binding Homo`, fill = `DA Change`)) + 
  geom_point(size = 3.5, color = "black", 
             shape = 21) +
  geom_smooth(method = lm, color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.direction = "horizontal", legend.position = "bottom") +
  guides(fill = guide_legend(title="")) +
  facet_grid(vars(Location)) +
  labs(x = "Heterozygous ????G Binding", y = "Homozygous 1/2 ????G Binding\n", tag = "A") 

h2 <- ggplot(dimerization, aes(x = `ddG Binding Hetero`, y = `Half ddG Binding Homo`)) + 
  geom_point(size = 3.5, color = "black", 
             shape = 21) +
  geom_smooth(method = lm, color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.direction = "horizontal", legend.position = "bottom") +
  guides(fill = guide_legend(title="")) +
  #facet_grid(vars(Location)) +
  labs(x = "\nHeterozygous ????G Binding (kcal/mol)", y = "Homozygous  1/2 ????G Binding (kcal/mol)\n") 

grid.arrange(h1, h2)



ggplot(strong, aes(x = fct_rev(fct_infreq(`AA Change`)), fill = Location)) +
  geom_bar(colour = "black", position = dodge) +
  labs(x = "", y = "\n Number of Mutations \n") +
  geom_text(stat='count', aes(label=..count..), hjust = -0.7, size = 3) +
  coord_flip() + 
  theme_classic() 


# number of amino acid types per group
ggplot(dat.hetero[!duplicated(dat.hetero$Pos), ], aes(x = fct_rev(fct_infreq(`AA Property WT`)))) +
  geom_bar(colour = "black", position = dodge, fill ="#A6CEE3") +
  labs(x = "", y = "\n Number of Residues \n") +
  facet_wrap(vars(Region)) +
  geom_text(stat='count', aes(label=..count..), vjust = -0.5, size = 4) +
  theme_classic() 

dodge <- position_dodge(width = 0.6)
ggplot(dat, aes(x = Location, y = ddG)) +
  geom_violin(aes(fill = Group), position = dodge) +
  geom_boxplot(aes(group = interaction(Group, Location)), width = 0.07, outlier.shape = NA, position = dodge, fill = "gray") +
  labs(x = "", y = "\n ????G Folding (kcal/mol) \n") +
  coord_cartesian(ylim = quantile(dat.violin.WT$ddG, c(0, 0.982))) +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99"),
                    labels = c("hetero" = "Heterozygous", "homo" = "Homozygous")) + 
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_classic() 

# ddG Location
new_labels <- c("hetero" = "Heterozygous", "homo" = "Homozygous")
ggplot(dat, aes(x = Location, y = `ddG Binding`)) +
  geom_violin(position = dodge, fill = "#A6CEE3", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, position = dodge, fill = "gray") +
  labs(x = "", y = "\n ????G Folding (kcal/mol) \n") +
  facet_wrap(vars(Group), labeller = labeller(Group = c("hetero" = "Heterozygous", "homo" = "Homozygous"))) +
  stat_compare_means(comparisons = list(c("Buried", "Exposed"), c("Exposed", "Interacting"), c("Buried", "Interacting"))) +
  theme_classic()

ggplot(dat, aes(x = `AA Change`, y = `ddG Binding`, fill = Group)) +
  geom_boxplot(position = position_dodge(1), colour = "black", outlier.shape = NA) +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99"),
                    labels = c("hetero" = "Heterozygous", "homo" = "Homozygous")) + 
  labs(x = "", y = "????G Binding (kcal/mol)\n") +
  facet_grid(vars(Location)) +
  coord_cartesian(ylim = c(-3, 4)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_classic() 


sasa <- ggplot(dat.hetero, aes(x = SASA, y = `ddG Binding`, fill = `AA Property WT`)) +
  geom_point(size = 3.5, color = "black", 
             aes(fill = `AA Property WT`), shape = 21) +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99", "#B2DF8A", "#E31A1C"),
                    name = "Original Amino Acid Group") +
  geom_smooth(method = lm, color = "black") +
  labs(x = "\nSASA (Å)", y = "????G Folding (kcal/mol)\n", tag = "A") +
  facet_wrap(vars(Location)) +
  theme_classic()

sasa2 <- ggplot(dat.homo, aes(x = SASA, y = `ddG Binding`, fill = `AA Property WT`)) +
  geom_point(size = 3.5, color = "black", 
             aes(fill = `AA Property WT`), shape = 21) +
  scale_fill_manual(values = c("#A6CEE3", "#FB9A99", "#B2DF8A", "#E31A1C"),
                    name = "Original Amin Acid Group") +
  geom_smooth(method = lm, color = "black") +
  labs(x = "\nSASA (Å)", y = "????G Folding (kcal/mol)\n", tag = "B") +
  facet_wrap(vars(Location)) +
  theme_classic()


grid.arrange(sasa, sasa2)


# number of residues
ggplot(dat.hetero, aes(x = fct_rev(fct_infreq(`TE Change`)))) +
  geom_bar(fill = "#A6CEE3", colour = "black") +
  labs(x = "", y = "\n Number of Mutations \n") +
  geom_text(stat='count', aes(label=..count..), hjust = -0.3) +
  coord_flip() +
  theme_classic() 

ggplot(dat.homo, aes(x = `Location`, y = `ddG Binding`)) +
  geom_boxplot()


#-------------PLOTTING (BARPLOT)--------------
ggplot(dat, aes(x = Group, y = ddG)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#A6CEE3", "#B2DF8A"),
                    labels = c("hetero" = "Heterozygous", "homo" = "Homozygous")) +
  stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p.format") +
  theme_classic()


ggplot(dat.hetero, aes(x = fct_infreq(`TE Status`))) +
  geom_bar(color = "black", fill = "#A6CEE3") +
  labs(x = "Total Energy Change (Heterozygous to Homozygous) \n", y = "\n Number of Residues \n") +
  theme_classic() 
  

ggplot(dat.hetero, aes(x = fct_infreq(`TE Change`))) +
  geom_bar(color = "black", fill = "#A6CEE3") +
  labs(x = "Stability Change (Heterozygous to Homozygous) \n", y = "\n Number of Residues \n") +
  coord_flip() +
  theme_classic() 
