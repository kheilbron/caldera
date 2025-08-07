
# Input arguments
maindir      <- "~/projects/causal_genes/"
caldera_path <- "~/repos/caldera/"

# Load libraries
library(data.table)
library(tidyverse)
source( file.path( caldera_path, "z_caldera.R" ) )

# mem_used: Memory used by an R session
mem_used <- function(){
  out <- gc()
  currentMb <- sum(out[,2])
  maxMb <- sum(out[,7])
  r_mem <- paste0('current ', currentMb, 'Mb, max ', maxMb, 'Mb')
  r_mem
}

# message_header: Nicely-formatted text to break up your log files into readable chunks
message_header <- function(...){
  message( "\n\n#--------------------------------------------------" )
  message( "#   ", ...)
  message( "#--------------------------------------------------" )
}

# Find all CS and V2G files
indir    <- file.path( maindir, "run_caldera_on_all_traits/cs_and_v2g_files" )
cs_paths   <- list.files( path=indir, pattern="_cs.tsv", full.names=TRUE )
pops_paths <- list.files( path=indir, pattern="_pops.tsv", full.names=TRUE )

# Loop through traits
# Run CALDERA
traits <- sub( pattern="^.*cs_and_v2g_files/(.*)_cs.tsv$", replacement="\\1", x=cs_paths )
out <- list()
mem_used()
t0 <- proc.time()
for( i in seq_along(traits) ){
  
  # Message
  message_header( "Starting trait ", i, "/", length(traits), ": ", traits[i] )
  
  # Specify inputs
  cs_file   <- cs_paths[i]
  pops_file <- pops_paths[i]
  
  # Run CALDERA
  out[[i]] <- caldera( pops_file    = pops_file, 
                       cs_file      = cs_file, 
                       assembly     = 37,
                       caldera_path = caldera_path )
  out[[i]]$trait <- traits[i]
}
t1 <- proc.time()

# How long did it take?
t1 - t0
mem_used()

# Format
out <- do.call( rbind, out )

# Add a trait description
td1 <- fread("~/projects/causal_genes/weeks2023_table_s1.csv")
names(td1) <- tolower( names(td1) )
td2 <- td1 %>% select( trait, description ) %>%
  distinct( trait, description, .keep_all=TRUE )
out2 <- left_join( x=out, y=td2, by="trait" )

# Write output to file
outfile <- file.path( maindir, "run_caldera_on_all_traits/caldera_results_93_ukb_traits.tsv" )
fwrite( x=out2, file=outfile, sep="\t" )

















