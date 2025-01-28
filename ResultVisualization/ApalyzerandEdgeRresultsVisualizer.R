# Ustawienie lustrzanego serwera CRAN
options(repos = c(CRAN = "https://cran.rstudio.com"))

# Załaduj wymagane pakiety
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(VennDiagram)


# Sprawdź argumenty przekazane do skryptu
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Użycie: Rscript <skrypt.R> <plik_apalazer.csv> <plik_edger.csv> <liczba_grup_badawczych> <liczba_grup_kontrolnych>")
}

# Wczytaj argumenty
apalazer_file <- args[1]
edger_file <- args[2]
num_exp <- as.integer(args[3])
num_ctrl <- as.integer(args[4])

# Wczytaj dane
apalazer_data <- read_tsv(apalazer_file)
edger_data <- read_tsv(edger_file)

# Debugowanie - wyświetl nazwy kolumn
cat("Kolumny w pliku Apalazer:\n")
print(colnames(apalazer_data))

cat("Kolumny w pliku EdgeR:\n")
print(colnames(edger_data))

# Sprawdź poprawność danych
if (!"Geneid" %in% colnames(apalazer_data) | !"Geneid" %in% colnames(edger_data)) {
  stop("Pliki muszą zawierać kolumnę 'Geneid'.")
}

if (!"logFC" %in% colnames(apalazer_data) | !"logFC" %in% colnames(edger_data)) {
  stop("Pliki muszą zawierać kolumnę 'logFC'.")
}

# Połącz dane na podstawie 'Geneid'
merged_data <- apalazer_data %>%
  inner_join(edger_data, by = "Geneid", suffix = c("_apalazer", "_edger"))

# Debugowanie - wyświetl kolumny po połączeniu
cat("Kolumny w merged_data:\n")
print(colnames(merged_data))

# Oblicz średnie wartości dla grup kontrolnych i badawczych
merged_data <- merged_data %>%
  mutate(
    mean_Ctrl = rowMeans(select(., starts_with("Control_")), na.rm = TRUE),
    mean_Exp = rowMeans(
      select(., where(is.numeric)), # Uwzględnia wyłącznie kolumny liczbowe
      na.rm = TRUE
    )
  )


# Wizualizacja różnic logFC
plot <- ggplot(merged_data, aes(x = logFC_apalazer, y = logFC_edger)) +
  geom_point(aes(color = FDR_apalazer < 0.05 & FDR_edger < 0.05), size = 2) +
  labs(
    title = "Porównanie logFC pomiędzy wynikami Apalazera i EdgeR",
    x = "logFC Apalazer",
    y = "logFC EdgeR",
    color = "Istotne FDR (<0.05)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )

# Stwórz katalog 'visualization_results', jeśli nie istnieje
output_dir <- "visualization_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Zapisz wykres do pliku w folderze 'visualization_results'
output_file <- file.path(output_dir, "logFC_comparison_plot.png")
ggsave(output_file, plot, width = 10, height = 6)

histogram_logFC <- ggplot(merged_data, aes(x = logFC_apalazer)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = "Rozkład logFC (Apalazer)",
    x = "logFC Apalazer",
    y = "Częstość"
  ) +
  theme_classic()

# Zapisz histogram
ggsave(file.path(output_dir, "logFC_histogram_apalazer.png"), histogram_logFC, width = 8, height = 6)

boxplot_FDR <- ggplot(merged_data, aes(x = factor(1), y = FDR_apalazer)) +
  geom_boxplot(fill = "orange", alpha = 0.6) +
  labs(
    title = "Boxplot FDR (Apalazer)",
    x = "",
    y = "FDR"
  ) +
  theme_classic()

# Zapisz wykres pudełkowy
ggsave(file.path(output_dir, "FDR_boxplot_apalazer.png"), boxplot_FDR, width = 8, height = 6)

scatter_FDR_logFC <- ggplot(merged_data, aes(x = logFC_apalazer, y = FDR_apalazer)) +
  geom_point(alpha = 0.5, color = "red") +
  labs(
    title = "Zależność FDR i logFC (Apalazer)",
    x = "logFC Apalazer",
    y = "FDR Apalazer"
  ) +
  theme_classic()

# Zapisz wykres punktowy
ggsave(file.path(output_dir, "FDR_vs_logFC_apalazer.png"), scatter_FDR_logFC, width = 8, height = 6)

combined_histogram <- ggplot(merged_data) +
  geom_histogram(aes(x = logFC_apalazer, fill = "Apalazer"), binwidth = 0.5, alpha = 0.6, position = "identity") +
  geom_histogram(aes(x = logFC_edger, fill = "EdgeR"), binwidth = 0.5, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("Apalazer" = "blue", "EdgeR" = "green")) +
  labs(
    title = "Porównanie rozkładów logFC (Apalazer vs EdgeR)",
    x = "logFC",
    y = "Częstość",
    fill = "Metoda"
  ) +
  theme_classic()

# Zapisz wykres
ggsave(file.path(output_dir, "logFC_histogram_comparison.png"), combined_histogram, width = 8, height = 6)

heatmap_data <- merged_data %>%
  mutate(diff_logFC = abs(logFC_apalazer - logFC_edger)) %>%
  select(Geneid, diff_logFC) %>%
  arrange(desc(diff_logFC)) %>%
  head(50)  # Pokaż tylko 50 genów z największymi różnicami

heatmap <- ggplot(heatmap_data, aes(x = Geneid, y = diff_logFC, fill = diff_logFC)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(
    title = "Heatmapa różnic logFC między Apalazer i EdgeR",
    x = "Gen",
    y = "Różnica logFC"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Zapisz heatmapę
ggsave(file.path(output_dir, "logFC_diff_heatmap.png"), heatmap, width = 12, height = 6)

library(VennDiagram)

significant_genes <- list(
  Apalazer = merged_data %>% filter(FDR_apalazer < 0.05) %>% pull(Geneid),
  EdgeR = merged_data %>% filter(FDR_edger < 0.05) %>% pull(Geneid)
)

venn <- venn.diagram(
  significant_genes,
  filename = file.path(output_dir, "venn_diagram.png"),
  fill = c("blue", "green"),
  alpha = 0.5,
  cat.col = c("blue", "green"),
  cat.cex = 1.5,
  main = "Porównanie genów istotnych (FDR < 0.05)"
)


cat(paste0("Wykresy zapisano w katalogu '", output_dir, "' jako 'logFC_comparison_plot.png'.\n"))
