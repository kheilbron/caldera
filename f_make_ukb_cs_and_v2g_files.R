
# Input arguments
maindir <- "~/projects/causal_genes/"

# Load libraries
library(data.table)
library(tidyverse)

# z_to_p: Convert a z-score into a P value
z_to_p <- function( z, log.p=FALSE ){
  if(log.p){
    log(2) + pnorm( -abs(z), log.p=TRUE )
  }else{
    2*pnorm( -abs(z) )
  }
}  

# Read in UKB V2G
vg_infile <- file.path( maindir, "v2g_file.tsv" )
vg <- fread( input=vg_infile )
tapply( vg$pops_score, INDEX=vg$trait, FUN=function(x) mean( !is.na(x) ) ) #%>% table()

# Read in CS column names
# Read in CSs, add column names
cs_colnames <- fread( file.path( maindir, "release1.1", "UKBB_94traits_release1.cols" ), 
                      header=FALSE )
cs <- fread( file.path( maindir, "release1.1", "UKBB_94traits_release1.bed.gz" ) )
names(cs) <- cs_colnames$V1

# Remove SNPs that are not in CSs, or are in CSs with the 5th+ worst ABF
# Subset to SUSIE CSs only
# Add P value and TCP columns
n_css <- 5
cs2 <- cs[  cs$cs_id > 0 & 
              cs$cs_id <= n_css , ]
cs3 <- cs2[ cs2$method == "SUSIE" , ]
cs3$p <- z_to_p( z=sqrt(cs3$chisq_marginal) )
cs3$tcp <- paste( cs3$trait, cs3$region, cs3$cs_id, sep="_" )
head( cs3, 1 )
tapply( cs3$tcp, INDEX=cs3$trait, FUN=function(x) length( unique(x) ) ) %>% sort( decreasing=TRUE )

# Ensure that the same traits are present in both CSs and V2G
shared_traits <- intersect( vg$trait, cs3$trait )
vg  <- vg[  vg$trait  %in% shared_traits , ]
cs3 <- cs3[ cs3$trait %in% shared_traits , ]

# Make a directory to store CS and V2G files
outdir <- file.path( maindir, "run_caldera_on_all_traits/cs_and_v2g_files" )
dir.create( path=outdir, showWarnings=FALSE, recursive=TRUE )

# Loop through traits
# Write trait-specific CS and V2G files
for( i in unique(cs3$trait) ){
  
  # Subset to focal CS and required columns
  sub_cs <- cs3 %>%
    filter( trait == i ) %>%
    select( tcp, chromosome, end, pip, rsid ) %>%
    rename( locus = tcp,
            chr   = chromosome,
            bp    = end ) %>%
    mutate( chr = sub( pattern="^chr", replacement="", x=chr ) ) %>%
    mutate( chr = as.integer(chr) )
  
  # Subset to focal POPS
  sub_vg <- vg %>%
    filter( trait == i ) %>%
    select( ensgid, pops_score ) %>%
    rename( ENSGID     = ensgid,
            PoPS_Score = pops_score )
  
  # Write
  cs_outfile <- paste0( outdir, "/", i, "_cs.tsv" )
  vg_outfile <- paste0( outdir, "/", i, "_pops.tsv" )
  fwrite( x=sub_cs, file=cs_outfile, sep="\t" )
  fwrite( x=sub_vg, file=vg_outfile, sep="\t" )
}



























