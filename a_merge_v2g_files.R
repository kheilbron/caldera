
# Input arguments
maindir <- "~/projects/causal_genes/"

# Load libraries
suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(dplyr) )

# Read in gencode gene boundaries
gencode  <- fread("~/projects/causal_genes/gene_locations.tsv")

# Read in main V2G
vg_file1 <- file.path( maindir, "UKB_AllMethods_GenePrioritization.txt.gz" )
vg1 <- fread(vg_file1)
vg1 <- vg1[ order( vg1$trait, vg1$region, vg1$cs_id, vg1$ensgid ) , ]

# Read in MAGMA and SMR V2G
# One row per credible set variant. We extract gene-level results for the top variant.
vg_file2 <- file.path( maindir, "PoPS_UKBB_noncoding_validation_1348CSs_v2.txt.gz" )
vg2 <- fread(vg_file2)
vg2 <- vg2[ order( vg2$trait, vg2$region, vg2$cs_id, vg2$ensgid, -vg2$pip ) , ]
idx <- duplicated( vg2[ , c( "trait", "region", "cs_id", "ensgid" ) ] )
vg2 <- vg2[ !idx , ]
gcols_vg2 <- c( "trait", "region", "cs_id", "ensgid", "smr_p", "smr_rank", 
                "magma_score", "magma_rank" )

# Read in DEPICT and NetWAS V2G
vg_file3 <- file.path( maindir, "netwas_depict.txt.gz" )
vg3 <- fread(vg_file3)
vg3 <- vg3[ !grepl( pattern="^PASS_", x=vg3$trait ) , ]
names(vg3)[ names(vg3) == "ENSGID" ] <- "ensgid"
vg3$trait <- sub( pattern="^UKB_", replacement="", x=vg3$trait )

# Merge
j1 <- left_join( x=vg1, y=vg2[,..gcols_vg2] )
j2 <- left_join( x=j1, y=vg3 )

# Remove genes >300kb from the lead variant and non-canonical ENSGIDs
# j3 <- j2[ j2$distance_genebody < 300e3 & j2$ensgid %in% gencode$ENSGID , ]
j3 <- j2[ j2$ensgid %in% gencode$ENSGID , ]
NROW(vg1); NROW(vg2); NROW(vg3)
NROW(j1); NROW(j3)

# Write
vg_outfile <- file.path( maindir, "v2g_file.tsv" )
fwrite( x=j3, file=vg_outfile, sep="\t" )










