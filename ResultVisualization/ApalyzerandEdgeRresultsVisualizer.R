# Ustawienie lustrzanego serwera CRAN
options(repos = c(CRAN = "https://cran.rstudio.com"))

# Załaduj wymagane pakiety
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")

# Załaduj pakiety
# Załaduj pakiety
library(ggplot2)
library(dplyr)
library(readr)

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
  "NC" = "gray"
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

# Filtrowanie do istotnych wartości (bez NC)
selected_genes <- final_result %>%
  filter(Category != "NC")

# Wykres słupkowy
bar_plot <- ggplot(selected_genes, aes(x = DE, fill = Category)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = category_colors) +
  labs(
    title = "Liczba genów z deregulacją ekspresji (DE) i APA",
    x = "Deregulacja ekspresji (DE)",
    y = "Liczba genów",
    fill = "Kategoria"
  ) +
  theme_classic()  # Białe tło

ggsave(file.path(output_dir, "DE_APA_comparison.png"), bar_plot, width = 8, height = 6)

# Wykres punktowy logFC vs FDR
scatter_plot <- ggplot(selected_genes, aes(x = logFC, y = FDR, color = Category, label = gene_symbol)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text(vjust = -1, size = 3) +
  scale_color_manual(values = category_colors) +
  labs(
    title = "Wykres logFC vs FDR dla genów z APAreg i DE",
    x = "logFC",
    y = "FDR",
    color = "Kategoria"
  ) +
  theme_classic()  # Białe tło

ggsave(file.path(output_dir, "logFC_vs_FDR_selected_genes.png"), scatter_plot, width = 10, height = 6)

cat("Wizualizacje zapisano w katalogu 'visualization_results'.\n")
