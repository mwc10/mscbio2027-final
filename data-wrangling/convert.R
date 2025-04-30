library(reticulate)
library(tidyverse)
library(stringr)
library(readxl)
library(here)

py_require("numpy")
py_require("anndata")
py_require("pandas")
py_require("scipy")

np <- import("numpy")
ad <- import("anndata")
pd <- import("pandas")
csr_matrix <- import("scipy")$sparse$csr_matrix

###### Usage #################################################
# This script converts the matrices and dataframes
# into python anndata files for later use
#
### Rscript convert.R <path-to-downloaded-data>
#
##############################################################


### Conversion Functions #####
get_patient_id <- function(fname) {
  fname |>
    str_extract("\\d+") |>
    as.numeric()
}

output_counts <- here("data/anndata/counts/")
output_annos <- here("data/anndata/annotations/")

convert_counts <- function(file) {
  # read in counts
  id <- basename(file) |> get_patient_id()
  cnts <- readRDS(file)
  
  # they didn't label their counts correctly,
  # as they left out the slide in the rownames for the counts...
  # but it seems that they should match the associated "spots" df
  if (!all(str_ends(rownames(cnts$spots), rownames(cnts$cnts)))) {
    stop(paste0("Rownames in spots did not match rownames in counts matrix for ", id))
  }
  
  # create anndata object based on counts
  adata <- ad$AnnData(csr_matrix(np$array(cnts$cnts)))
  adata$obs_names = rownames(cnts$spots)
  adata$var_names = colnames(cnts$cnts)
  
  # add in count spot info 
  #adata$obs = pd$concat(c(adata$obs, r_to_py(cnts$spots)))
  adata$obs = r_to_py(cnts$spots)
  
  # output file
  out_fname = paste0(output_counts,'cnt-',id,'.h5ad')
  adata$write_h5ad(out_fname)
  message(out_fname)
}

convert_anno <- function(file) {
  # get sample id
  id <- basename(file) |> get_patient_id()
  
  # read in RDS
  anno = readRDS(file)
  
  # create anndata object 
  adata <- ad$AnnData(np$array(anno$annots))
  adata$obs_names = rownames(anno$annots)
  adata$var_names = colnames(anno$annots)
  
  # add spot metadata
  df = r_to_py(anno$spots)
  df = df$drop(columns=c("slide"))
  adata$obs = df
  
  # write out file
  out_fname = paste0(output_annos,'annotation-',id,'.h5ad')
  adata$write_h5ad(out_fname)
  
  message(out_fname)
}

###################################################
############# MAIN ################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  data_path <- here("data/ST/")
} else {
  data_path <- args[1]
}

message(paste0('Reading files from: ', data_path))

message(paste0('outputting converted files to: ', here('data/anndata')))

dir.create(output_counts, recursive = TRUE, showWarnings = FALSE)
dir.create(output_annos, recursive = TRUE, showWarnings = FALSE)

## Output the sample clinical info as a CSV
readRDS(paste0(data_path, "/Clinical/Clinical.RDS")) |>
  as_tibble() |>
  select(-c(annotations, RFS, iBCFS, iDFS, DRFS, OS)) |>
  write_csv(here('data/clinical-info.csv'))

## clinical annotation ids to map to spots 
paste0(data_path, "/Clinical/ids.RDS") |>
  readRDS() |>
  as_tibble(rownames = 'name') |>
  write_csv(here('data/clinical-ids.csv'))


### Deal with Counts and Spot Annotation
COUNT_PATH <- "/Robjects/counts/"
ANNO_PATH <- "/Robjects/annotsBySpot/"

counts_files <- list.files(paste0(data_path, COUNT_PATH), full.names = TRUE)
anno_files <- list.files(paste0(data_path, ANNO_PATH), full.names = TRUE)

counts_files |>
  purrr::walk(convert_counts)

anno_files |>
  purrr::walk(convert_anno)
