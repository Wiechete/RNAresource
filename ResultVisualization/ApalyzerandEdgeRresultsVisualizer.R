# Ustawienie lustrzanego serwera CRAN
options(repos = c(CRAN = "https://cran.rstudio.com"))

# Załaduj wymagane pakiety
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")

# Załaduj pakiety
library(ggplot2)
library(dplyr)
library(readr)
library(VennDiagram)

# Wczytaj dane
final_result <- read_csv("final_result.csv", na = c("NA", ""))

# Usunięcie wartości NA
final_result <- na.omit(final_result)

# Tworzenie katalogu na wizualizacje
output_dir <- "visualization_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Definiowanie kolorów dla kategorii
category_colors <- c(
  "UP/UP" = "red",
  "UP/DN" = "blue",
  "DN/UP" = "green",
  "DN/DN" = "purple",
  "NC" = "grey"
)

# Tworzenie nowej kolumny dla kategorii
final_result <- final_result %>%
  mutate(Category = case_when(
    DE == "UP" & APAreg == "UP" ~ "UP/UP",
    DE == "UP" & APAreg == "DN" ~ "UP/DN",
    DE == "DN" & APAreg == "UP" ~ "DN/UP",
    DE == "DN" & APAreg == "DN" ~ "DN/DN",
    TRUE ~ "NC"
  ))

# 🔹 **Wizualizacja 1: Oddzielne histogramy dla DN, NC, UP**
for (de_status in c("DN", "NC", "UP")) {
  sub_data <- final_result %>% filter(DE == de_status)
  
  hist_plot <- ggplot(sub_data, aes(x = APAreg, fill = Category)) +
    geom_bar(position = "dodge") +
    scale_fill_manual(values = category_colors) +
    labs(
      title = paste("Liczba genów dla DE =", de_status),
      x = "Deregulacja długości 3'UTR (APAreg)",
      y = "Liczba genów",
      fill = "Kategoria (DE/APAreg)"
    ) +
    theme_classic() +
    labs(caption = "Legenda: Pierwsze UP/DN oznacza ekspresję różnicową genu, drugie UP/DN oznacza zmianę długości 3'UTR.")
  
  filename <- paste0("DE_", de_status, "_histogram.png")
  ggsave(file.path(output_dir, filename), hist_plot, width = 8, height = 6)
}

# 🔹 Sprawdzenie liczby wyników APAlyzer per gen
gene_counts <- final_result %>%
  group_by(gene_symbol) %>%
  summarise(n = n())

# 🔹 Odfiltrowanie genów z więcej niż 3 wpisami
filtered_genes <- gene_counts %>% filter(n <= 3) %>% pull(gene_symbol)
filtered_result <- final_result %>% filter(gene_symbol %in% filtered_genes)

# 🔹 Wizualizacja 2: Nowy wykres wulkaniczny z filtrowanymi genami
for (de_status in c("DN", "NC", "UP")) {
  if (de_status == "NC") {
    sub_data <- filtered_result %>% filter(DE == de_status & abs(logFC) > 0.5) # Zmieniony próg
  } else {
    sub_data <- filtered_result %>% filter(DE == de_status & FDR < 0.05 & abs(logFC) > 1)
  }
  
  volcano_plot <- ggplot(sub_data, aes(x = logFC, y = -log10(FDR), color = Category, label = gene_symbol)) +
    geom_point(size = 3, alpha = 0.7) +  # Usunięcie legendy dla rozmiaru
    geom_text(vjust = -1, size = 3) +
    scale_color_manual(values = category_colors) +
    labs(
      title = paste("Wykres wulkaniczny dla próbek o niezmienionym poziomie ekspresji DE =", de_status),
      x = "logFC",
      y = "-log10(FDR)"
    ) +
    theme_classic()
  
  filename <- paste0("DE_", de_status, "_volcano_filtered.png")
  ggsave(file.path(output_dir, filename), volcano_plot, width = 10, height = 6)
}


# 🔹 **Wizualizacja 3: Diagram Venna (UP/DN vs UP/DN)**
venn_filename1 <- file.path(output_dir, "venn_DE_UP_DN_vs_APA_UP_DN.png")

DE_genes_UP_DN <- final_result$gene_symbol[final_result$DE %in% c("UP", "DN")]
APAreg_genes_UP_DN <- final_result$gene_symbol[final_result$APAreg %in% c("UP", "DN")]

venn.plot1 <- draw.pairwise.venn(
  area1 = length(DE_genes_UP_DN),
  area2 = length(APAreg_genes_UP_DN),
  cross.area = length(intersect(DE_genes_UP_DN, APAreg_genes_UP_DN)),
  category = c("Ekspresja różnicowa (UP/DN)", "Zmiana długości 3'UTR (UP/DN)"),
  fill = c("lightblue", "pink"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-30, 30),
  cat.dist = 0.05
)

png(venn_filename1, width = 800, height = 600)
grid.draw(venn.plot1)
grid.text("Legenda: Geny UP/DN w ekspresji vs UP/DN w APAreg", x = 0.5, y = 0.1, gp = gpar(fontsize = 12, col = "black"))
dev.off()

# 🔹 **Wizualizacja 4: Diagram Venna (NC vs UP/DN)**
venn_filename2 <- file.path(output_dir, "venn_DE_NC_vs_APA_UP_DN.png")

DE_genes_NC <- final_result$gene_symbol[final_result$DE == "NC"]
APAreg_genes_UP_DN <- final_result$gene_symbol[final_result$APAreg %in% c("UP", "DN")]

venn.plot2 <- draw.pairwise.venn(
  area1 = length(DE_genes_NC),
  area2 = length(APAreg_genes_UP_DN),
  cross.area = length(intersect(DE_genes_NC, APAreg_genes_UP_DN)),
  category = c("Brak deregulacji ekspresji (NC)", "Zmiana długości 3'UTR (UP/DN)"),
  fill = c("lightgray", "pink"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-30, 30),
  cat.dist = 0.05
)

png(venn_filename2, width = 800, height = 600)
grid.draw(venn.plot2)
grid.text("Legenda: Geny bez zmiany ekspresji (NC) vs UP/DN w APAreg", x = 0.5, y = 0.1, gp = gpar(fontsize = 12, col = "black"))
dev.off()

cat("Wizualizacje zapisano w katalogu 'visualization_results'.\n")
