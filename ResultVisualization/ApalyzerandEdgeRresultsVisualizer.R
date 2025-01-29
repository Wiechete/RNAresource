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
  "UP/UP" = "red",     # UP w ekspresji + UP w długości 3'UTR
  "UP/DN" = "blue",    # UP w ekspresji + DN w długości 3'UTR
  "DN/UP" = "green",   # DN w ekspresji + UP w długości 3'UTR
  "DN/DN" = "purple",  # DN w ekspresji + DN w długości 3'UTR
  "NC" = "gray"        # NC oznacza brak deregulacji w jednej z kategorii
)

# Tworzenie nowej kolumny dla kategorii
final_result <- final_result %>%
  mutate(Category = case_when(
    DE == "UP" & APAreg == "UP" ~ "UP/UP",
    DE == "UP" & APAreg == "DN" ~ "UP/DN",
    DE == "DN" & APAreg == "UP" ~ "DN/UP",
    DE == "DN" & APAreg == "DN" ~ "DN/DN",
    TRUE ~ "NC"  # Jeśli NC pojawia się w DE lub APAreg, cała kategoria jest NC
  ))

# 🔹 **Wizualizacja 1: Wykres słupkowy**
bar_plot <- ggplot(final_result, aes(x = DE, fill = Category)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = category_colors) +
  labs(
    title = "Liczba genów z deregulacją ekspresji (DE) i długości 3'UTR (APAreg)",
    x = "Deregulacja ekspresji (DE)",
    y = "Liczba genów",
    fill = "Kategoria (DE/APAreg)"
  ) +
  theme_classic()  # Białe tło

ggsave(file.path(output_dir, "DE_APA_comparison.png"), bar_plot, width = 8, height = 6)

# 🔹 **Wizualizacja 2: Wykres punktowy logFC vs FDR**
scatter_plot <- ggplot(final_result, aes(x = logFC, y = FDR, color = Category, label = gene_symbol)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text(vjust = -1, size = 3) +
  scale_color_manual(values = category_colors) +
  labs(
    title = "Wykres logFC vs FDR dla genów z APAreg i DE",
    x = "logFC",
    y = "FDR",
    color = "Kategoria (DE/APAreg)"
  ) +
  theme_classic()  # Białe tło

ggsave(file.path(output_dir, "logFC_vs_FDR_selected_genes.png"), scatter_plot, width = 10, height = 6)

# 🔹 **Wizualizacja 3: Diagram Venna dla DE i APAreg**
venn_filename <- file.path(output_dir, "venn_DE_APAreg.png")

# Geny unikalne dla DE i APAreg
DE_genes <- final_result$gene_symbol[final_result$DE != "NC"]
APAreg_genes <- final_result$gene_symbol[final_result$APAreg != "NC"]

# Tworzenie diagramu Venna
venn.plot <- draw.pairwise.venn(
  area1 = length(DE_genes),
  area2 = length(APAreg_genes),
  cross.area = length(intersect(DE_genes, APAreg_genes)),
  category = c("Różnicowa ekspresja (DE)", "Zmiana długości 3'UTR (APAreg)"),
  fill = c("lightblue", "pink"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-30, 30),
  cat.dist = 0.05
)

# Zapis diagramu Venna na białym tle
png(venn_filename, width = 800, height = 600)
grid.draw(venn.plot)
dev.off()

cat("Wizualizacje zapisano w katalogu 'visualization_results'.\n")
