1. Uruchomienie skryptu
Aby uruchomić skrypt, użyj poniższego polecenia w terminalu (zakładając, że plik ApalyzerandEdgeRresultsVisualizer.R jest w tym samym katalogu co pliki wejściowe):

bash
Kopiuj
Edytuj
Rscript ApalyzerandEdgeRresultsVisualizer.R <plik_apalazer.csv> <plik_edger.csv> <liczba_grup_badawczych> <liczba_grup_kontrolnych>
Argumenty wejściowe:

- plik_apalazer.csv>: Ścieżka do pliku wynikowego Apalazer (format CSV).
- plik_edger.csv>: Ścieżka do pliku wynikowego EdgeR (format CSV).
- liczba_grup_badawczych>: Liczba grup badawczych.
- liczba_grup_kontrolnych>: Liczba grup kontrolnych.

2. Generowanie wykresów
Skrypt generuje następujące wykresy w katalogu visualization_results:

Porównanie logFC (log fold change) pomiędzy wynikami Apalazera i EdgeR:

Typ wykresu: Scatter plot
Opis: Wykres punktowy przedstawiający zależność logFC pomiędzy wynikami Apalazera i EdgeR. Punkty są kolorowane na podstawie istotności statystycznej (FDR < 0.05).
Plik wynikowy: logFC_comparison_plot.png
Rozkład logFC (Apalazer):

Typ wykresu: Histogram
Opis: Histogram pokazujący rozkład logFC z wyników Apalazera.
Plik wynikowy: logFC_histogram_apalazer.png
Boxplot FDR (Apalazer):

Typ wykresu: Wykres pudełkowy (boxplot)
Opis: Boxplot przedstawiający rozkład wartości FDR dla wyników Apalazera.
Plik wynikowy: FDR_boxplot_apalazer.png
Zależność FDR i logFC (Apalazer):

Typ wykresu: Scatter plot
Opis: Wykres punktowy przedstawiający zależność wartości FDR i logFC dla wyników Apalazera.
Plik wynikowy: FDR_vs_logFC_apalazer.png
Porównanie rozkładów logFC (Apalazer vs EdgeR):

Typ wykresu: Histogram
Opis: Porównanie rozkładów logFC pomiędzy wynikami Apalazera i EdgeR.
Plik wynikowy: logFC_histogram_comparison.png
Heatmapa różnic logFC:

Typ wykresu: Heatmapa
Opis: Heatmapa pokazująca 50 genów z największymi różnicami w logFC pomiędzy Apalazerem a EdgeR.
Plik wynikowy: logFC_diff_heatmap.png
Wykres Venn'a dla genów istotnych (FDR < 0.05):

Typ wykresu: Diagram Venn'a
Opis: Diagram Venn'a pokazujący geny, które są istotne w obu analizach (Apalazer i EdgeR) z FDR < 0.05.
Plik wynikowy: venn_diagram.png

3. Gdzie znaleźć wyniki
Wszystkie wykresy są zapisywane w folderze visualization_results, który jest tworzony automatycznie, jeśli jeszcze nie istnieje. Pliki wykresów będą miały odpowiednie nazwy, np. logFC_comparison_plot.png, logFC_histogram_apalazer.png itd.