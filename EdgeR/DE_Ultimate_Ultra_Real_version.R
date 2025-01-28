# Function to check BiocManager version and update if needed
check_and_update_biocmanager <- function() {
  # Check if BiocManager is installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("BiocManager is not installed. Would you like to install it? (yes/no)")
    response <- tolower(readline())
    if (response == "yes") {
      install.packages("BiocManager")
      message("BiocManager has been installed.")
    } else {
      stop("BiocManager installation is required to proceed. Exiting.")
    }
  }
  
  # Check the installed version of BiocManager
  bioc_version <- packageVersion("BiocManager")
  required_version <- "1.30.0"
  
  # Compare versions and update if necessary
  if (bioc_version < required_version) {
    message("Your BiocManager version is ", bioc_version, ". A version of ", required_version, " or higher is required.")
    message("Would you like to update BiocManager? (yes/no)")
    response <- tolower(readline())
    if (response == "yes") {
      install.packages("BiocManager")
      message("BiocManager has been successfully updated to version ", packageVersion("BiocManager"), ".")
    } else {
      stop("BiocManager update is required to proceed. Exiting.")
    }
  } else {
    message("Your BiocManager version (", bioc_version, ") is up-to-date.")
  }
}	
# Run the function
check_and_update_biocmanager()

using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  n<-length(need)
  if(n>0){
    libsmsg<-if(n>2) paste(paste(need[1:(n-1)],collapse=", "),",",sep="") else need[1]
    print(libsmsg)
    if(n>1){
      libsmsg<-paste(libsmsg," and ", need[n],sep="")
    }
    libsmsg<-paste("The following packages could not be found: ",libsmsg,"\n\r\n\rInstall missing packages?",collapse="")
    if(winDialog(type = c("yesno"), libsmsg)=="YES"){       
      BiocManager::install(need)
      lapply(need,require,character.only=TRUE)
    }
  }
}

using("edgeR","limma","viridis","mixOmics","RColorBrewer","biomaRt","dplyr",
      "GOfuncR","tidyr","topGO","org.Hs.eg.db")

library(edgeR)
library(limma)
library(data.table)
library(viridis)
library(mixOmics)
library(RColorBrewer)
library(biomaRt)
library(dplyr)
library(GOfuncR) # brak wersji dla BioCmanagera 3.20
library(tidyr)
library(topGO)
library(org.Hs.eg.db)

# Pobranie argumentów z terminala
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
num_experimental <- as.integer(args[2])  # Liczba grup eksperymentalnych
num_control <- as.integer(args[3])       # Liczba grup kontrolnych

# Wczytanie danych (zakładamy, że dane są w formacie CSV)
rawdata <- read.table(input_file, header = TRUE, sep = ";", stringsAsFactors = FALSE, row.names = 1)

# Pomijanie pierwszych dwóch kolumn
data_counts <- rawdata[, -c(0, 1)]

# Sprawdzenie struktury danych
print(head(data_counts))
print(dim(data_counts))  # Liczba genów i próbek

# Upewnij się, że liczba kolumn zgadza się z sumą grup eksperymentalnych i kontrolnych
if (ncol(data_counts) != (num_experimental + num_control)) {
  stop("Suma grup eksperymentalnych i kontrolnych musi odpowiadać liczbie kolumn danych.")
}

# Przypisanie grup
group <- c(rep("UD", num_experimental), rep("Ctrl", num_control))

# Utworzenie obiektu DGEList
y <- DGEList(
  counts = data_counts,
  group = group
)

# Normalizacja
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)

# Wyświetlenie informacji o obiekcie DGEList
print(y$samples)

keep <- filterByExpr(y, min.count = 10, min.prop = 0.8) # można zmienić poziomy
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

# Ponowne przeliczenie wielkości bibliotek i czynników normalizujących
y$samples$lib.size <- colSums(y$counts)
y = calcNormFactors(y)
y$samples

cpm_y <- cpm(y, log=T, prior.count = 1)
head(cpm_y)
png("MDS_plot.png", width = 800, height = 600)
plotMDS(y)
dev.off()
group=c(rep("UD", num_experimental), rep("Ctrl", num_control)) #grupy jak wyżej
data.frame(Sample=colnames(y), group)

design = model.matrix(~0+group)
rownames(design) = colnames(y)
design

contrast <- makeContrasts(groupUD-groupCtrl, levels = design) #groupUD-groupCtrl - zmienić według swoich nazw
contrast
y = estimateDisp(y, design, robust=TRUE)
y$common.dispersion

# BCV plot
jpeg("BCV_plot.jpeg", width = 800, height = 600)
plotBCV(y, cex=1)
dev.off()

fit = glmQLFit(y, design)
lrt = glmQLFTest(fit, contrast=contrast)

colnames(design)

summary(decideTests(lrt))
png("MD_plot.png", width = 800, height = 600)
plotMD(lrt)
abline(h=c(-1.5, 1.5), col="blue")
dev.off()

# Korekta p-value na wielokrotne testowanie
#lrt$table$padj <- p.adjust(lrt$table$PValue, method = "BH")

all.res <- topTags(lrt, n=Inf)
all.res = as.data.frame(all.res)

# Filtracja genów z korektą p-value
up.res <- all.res[all.res$logFC > 0 & all.res$padj <= 0.05,]
up.res = up.res[order(up.res$logFC, decreasing = T),]

down.res <- all.res[all.res$logFC < 0 & all.res$padj <= 0.05,]
down.res = down.res[order(down.res$logFC, decreasing = F),]
up.res[1:10,]
down.res[1:10,]

table = topTags(lrt, n=Inf)
table = as.data.frame(table)
selY <- cpm_y[rownames(table)[table$FDR<=0.05 & abs(table$logFC)>=2],]

cimColor <- turbo(255, begin = 0.1, end = 0.9, direction = 1)
png("Heatmap.png", width = 800, height = 600)
cim(t(selY), color = cimColor, symkey=T, transpose = T)
dev.off()

# Ustaw większy timeout
httr::set_config(httr::config(timeout = 360))

# Pobranie danych z Ensembl
mart <- useMart("ensembl", host = "https://www.ensembl.org")
dataset <- useDataset("hsapiens_gene_ensembl", mart = mart)


url <- "https://www.ensembl.org"
if (!httr::http_error(url)) {
  print("Połączenie z Ensembl działa poprawnie!")
} else {
  print("Nie można połączyć się z Ensembl!")
}

###########################################################################

library(dplyr)

# Pobranie wyników analizy
table = topTags(lrt, n=Inf)
table = as.data.frame(table)
table$Geneid = row.names(table)

# Sprawdzenie kolumn
print(colnames(table))
print(head(table))

# Pobranie dodatkowych informacji o genach z Ensembl
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl", host = "https://www.ensembl.org"))
genes = table$Geneid
G_list = getBM(filters = "ensembl_gene_id", 
               attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"), 
               values = genes, mart = mart)

# Sprawdzenie wyników z Ensembl
print(head(G_list))

# Łączenie tabel
table = full_join(table, G_list, by = c("Geneid" = "ensembl_gene_id"))

# Sprawdzenie wyników po połączeniu
print(head(table))

# Wybór istniejących kolumn
if (all(c("Geneid", "external_gene_name", "logFC", "logCPM", "F", "PValue", "FDR") %in% colnames(table))) {
  table = table[,c("Geneid", "external_gene_name", "logFC", "logCPM", "F", "PValue", "FDR")]
  colnames(table)[2] = "gene_name"
} else {
  stop("Nie można wybrać kolumn, ponieważ nie istnieją w tabeli.")
}

# Obliczanie CPM dla grup eksperymentalnych i kontrolnych
cpm_y2 = as.data.frame(cpm_y)
cpm_y2$Geneid = row.names(cpm_y2)
cpm_y2$mean_Ctrl = rowMeans(cpm_y2[,c(6:8)], na.rm = TRUE)
cpm_y2$mean_Exp = rowMeans(cpm_y2[,c(1:5)], na.rm = TRUE)

# Łączenie tabeli z CPM
table1 <- merge(table, cpm_y2, by = "Geneid", all.x = TRUE)  # Utrzymuje wszystkie wiersze z 'table'

# Add a new column using case_when
table1 <- table1 %>%
	mutate(DE = case_when(
		FDR <= 0.05 & logFC > 0 ~ "UP",
		FDR <= 0.05 & logFC < 0 ~ "DN",
		TRUE ~ "NC"  # Default case
	))
# Sprawdzenie wyników po połączeniu
print(head(table1))

# Konwersja kolumn typu factor na character
table1 %>% mutate_if(is.factor, as.character) -> table1

# Uzupełnianie nazw genów, jeśli są puste
table1$gene_name = ifelse(is.na(table1$gene_name) | table1$gene_name == "", table1$Geneid, table1$gene_name)

# Zapis wyników do pliku CSV z separatorem tabulatora
write.csv(format(table1, scientific = TRUE), 
            file = "DE_hs_min1_ult.csv", 
            row.names = FALSE, 
            quote = FALSE)
