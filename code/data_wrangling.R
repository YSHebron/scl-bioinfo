options(repos = "https://cloud.r-project.org/")

# Installs pacman ("package manager") if needed
if (!require("pacman"))
  install.packages("pacman")
library(pacman)

# Install xml2 in advance to prep for tidyverse installation in Linux.
if (pacman::p_detectOS() == "Linux" && !pacman::p_exists(xml2, local = TRUE)) {
  install.packages("xml2", dependencies = TRUE, INSTALL_opts = c("--no-lock"))
  pacman::p_load(xml2)
}

# Use pacman to load add-on packages as desired.
pacman::p_load(plyr, GGally, ggthemes, ggvis, plotly, psych,
               htmlwidgets, rio, markdown, shiny, tidyverse,
               ape, seqinr, kmer, validate, gsubfn,
               Rtsne, tsne, umap, factoextra, scales,
               RColorBrewer, ggfortify, devtools,
               ggdendro, dendextend, cluster, colorspace,
               microbenchmark, data.table, highcharter,
               wordcloud, tm)

test <- readr::read_table("data_yeast.txt",
                          col_names = c("p1", "p2", "feature", "score")) %>%
  dplyr::filter(feature %in% list("PPIREL", "STRING", "PUBMED", "PPIREL_CD")) %>%
  dplyr::select(p1, p2, score) %>%
  dplyr::distinct(p1, p2, .keep_all = TRUE)

df <- readr::read_table("data_yeast.txt",
                          col_names = c("p1", "p2", "feature", "score")) %>%
  dplyr::filter(feature %in% list("PPIREL")) %>%
  dplyr::select(p1, p2, score) %>%
  dplyr::distinct(p1, p2, .keep_all = TRUE)

rand_df <- df[sample(nrow(df), size = 10000),]

rand_df %>% readr::write_csv("data_yeast_rand.csv")

test %>% readr::write_csv("data_yeast_scl_noheader.csv", col_names = FALSE)
