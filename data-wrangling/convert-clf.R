library(tidyverse)
library(purrr)
library(here)

extract_data <- function(elem, idx) {
  elem$spot |>
    as_tibble(rownames="slide_id") |>
    bind_cols(
      elem$xy |> as_tibble() |> rename(slide_rep = id, x = V2, y = V3)
    ) |>
    bind_cols(
      elem$pr |> as_tibble()
    ) |>
    mutate(
      tnbc_id = idx
    ) |>
    relocate(tnbc_id, .before = 1)
}


classif_all <- readRDS(here("data/ST/classification/classifAll.RDS"))

df_classifications <- classif_all |>
  purrr::imap(extract_data) |>
  purrr::reduce(bind_rows)

df_classifications |>
  write_csv(here("data/all_classifications.csv"))
