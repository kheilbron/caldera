

# Set the path to the CALDERA repo
caldera_path <- "~/repos/caldera/"

# Set paths
ex_in_path   <- file.path( caldera_path, "example", "inputs" )
ex_out_path  <- file.path( caldera_path, "example", "output" )

# Load libraries
library(data.table)
library(tidyverse)
source( file.path( caldera_path, "z_caldera.R" ) )

# Specify inputs
cs_file   <- file.path( ex_in_path, "lange2025_parkinsons_disease_cred_sets.tsv" )
pops_file <- file.path( ex_in_path, "lange2025_parkinsons_disease_pops.preds" )

# Run CALDERA
out <- caldera( pops_file    = pops_file, 
                cs_file      = cs_file, 
                assembly     = 37,
                caldera_path = caldera_path )

# Write output to file
outfile <- file.path( ex_out_path, "lange2025_parkinsons_disease_caldera_output.tsv" )
fwrite( x=out, file=outfile, sep="\t" )

