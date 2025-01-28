Aby uruchomić skrypt, użyj poniższego polecenia w terminalu (zakładając, że skrypt jest w tym samym katalogu co pliki wejściowe):

Rscript Apalyzer_BF.R [liczba próbek eksperymentalnych] [liczba próbek kontrolnych] [folder z plikami BAM i GTF]

Jeżeli folder z plikami BAM i GTF nie zotanie podany na wejściu, skrypt automatycznie zacznie pracować na plikach w tym samym folderze, w którym się znajduje.
W folderze ze skryptem MUSI istnieć plik sample_names_metadata.txt, w którym po przecinku bez spacji są wymieniane nazwy próbek.
