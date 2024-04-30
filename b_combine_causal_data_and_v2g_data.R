

#-------------------------------------------------------------------------------
#   simple_functions
#-------------------------------------------------------------------------------

# logit10: Convert a probability to a log (base 10) odds ratio
logit10 <- function(p) log10( p / (1-p) )

# p_to_z: Convert a P value into a z-score
p_to_z <- function( p, direction=NULL, limit=.Machine$double.xmin, log.p=FALSE ){
  
  # Set P value lower limit to avoid Inf/-Inf
  if( !is.null( limit ) )  p[ which( p < limit ) ] <- limit
  
  # Get z
  if(log.p){
    z <- -qnorm( p - log(2), log.p=TRUE )
  }else{
    z <- -qnorm(p/2)
  }
  
  # Correct sign, return
  if ( !is.null( direction) )  z <-  z * sign(direction)
  z
}

# prop_to_odds: Converts a proportion to an odds
prop_to_odds <- function(proportion){
  proportion / (1 - proportion)
}

# z_to_p: Convert a z-score into a P value
z_to_p <- function( z, log.p=FALSE ){
  if(log.p){
    log(2) + pnorm( -abs(z), log.p=TRUE )
  }else{
    2*pnorm( -abs(z) )
  }
}  



#-------------------------------------------------------------------------------
#   Merge causal genes with V2G
#-------------------------------------------------------------------------------

# Input arguments
maindir <- "~/projects/causal_genes/"

# Load libraries
library(data.table)
library(dplyr)

# Read in causal/non-causal trait-gene pairs
cnc_file <- file.path( maindir, "causal_noncausal_trait_gene_pairs", 
                       "causal_noncausal_trait_gene_pairs_300kb.tsv" )
cnc <- fread(cnc_file)

# Read in V2G
vg_file1 <- file.path( maindir, "UKB_AllMethods_GenePrioritization.txt.gz" )
vg1 <- fread(vg_file1)
vg1 <- vg1[ order( vg1$trait, vg1$region, vg1$cs_id, vg1$ensgid ) , ]

# Read in extended V2G (MAGMA and SMR)
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

# Remove genes >300kb from the lead variant
j3 <- j2[ j2$distance_genebody < 300e3 , ]
NROW(vg1); NROW(vg2); NROW(vg3)
NROW(j1); NROW(j3) 

# Subset both to shared trait-ENSGID-CS triplets (TECs)
# Note: yields about half as many missed causal/non-causal TGPs than using gene
# Note: that match() will only select the first instance for the small...
# ...number of cases where there are multiple ENSGIDs per gene
tec_cnc <- paste( cnc$trait,  cnc$ensgid,  cnc$cs_id,  sep="_" )
tec_vg <- paste( j3$trait, j3$ensgid, j3$cs_id, sep="_" )
tec_shared <- intersect( tec_cnc, tec_vg )
cnc2 <- cnc[  match( tec_shared, tec_cnc ) , ]
j4 <- j3[ match( tec_shared, tec_vg ) , ]

# Explore: What percentage of triplets were shared?
length(tec_shared) /  length(tec_cnc) #93% are found in the causal gene dataset
length(tec_shared) /  length(tec_vg)  #8% are found in the V2G dataset

# Merge the two datasets
gcols_vg <- c( "trait", "gene", "ensgid", "region", "cs_id", "chromosome",
               "start", "end", "tss", "pops_score", "pops_rank", "abc_score",
               "abc_rank", "pchic_jung_score", "pchic_jung_rank", 
               "pchic_javierre_score", "pchic_javierre_rank", 
               "corr_andersson_score", "corr_andersson_rank", "corr_liu_score",
               "corr_liu_rank", "corr_ulirsch_score", "corr_ulirsch_rank",
               "twas_z", "twas_p", "twas_rank", "clpp_prob", "clpp_rank", 
               "coding_prob", "coding_rank", "distance_tss", 
               "distance_tss_rank", "distance_genebody", 
               "distance_genebody_rank", "smr_p", "smr_rank", 
               "magma_score", "magma_rank", 
               "full_loco_netwas_score", "full_loco_netwas_bon_score",
               "no_loco_nc_netwas_score", "no_loco_nc_netwas_bon_score", 
               "no_loco_nc_depict_p", "full_loco_depict_p" )
gcols_cnc <- c( "n_cs_in_region", "n_cs_snps", "cs_width", "ngenes_nearby" )
m <- cbind( causal=cnc2$causal, j4[ , ..gcols_vg ], cnc2[ , ..gcols_cnc ] )

# Add a column for the number of genes in the locus
m$prior_n_genes_locus <- logit10( 1/m$ngenes_nearby )
tcp_m <- paste( m$trait, m$region, m$cs_id, sep="_" )


#-------------------------------------------------------------------------------
#   Add columns for global features
#-------------------------------------------------------------------------------

# POPS
m$pops_glo <- m$pops_score

# Distance to gene body 
m$dist_gene_raw_l2g <- log( m$distance_genebody + 1 )
m$dist_gene_glo     <- log10( m$distance_genebody + 1e3 )

# Distance to TSS
m$dist_tss_raw_l2g  <- log( m$distance_tss + 1 )
m$dist_tss_glo      <- log10( m$distance_tss + 1e3 )

# TWAS
m$twas_p <- ifelse( m$twas_p == 0, .Machine$double.xmin, m$twas_p )
m$twas_logp_glo <- ifelse( is.na(m$twas_p), -log10(1), -log10(m$twas_p) )
m$twas_glo <- ifelse( is.na(m$twas_z), 0,            abs(m$twas_z) )

# E-P correlation
m$corr_liu_raw_l2g <- ifelse( is.na(m$corr_liu_score), 0, m$corr_liu_score )
m$corr_liu_glo  <- log10( ifelse( m$corr_liu_raw_l2g < 1, 1, m$corr_liu_raw_l2g ) )
m$corr_and_glo  <- ifelse( is.na(m$corr_andersson_score), 0, m$corr_andersson_score )
m$corr_uli_glo  <- ifelse( is.na(m$corr_ulirsch_score), 0, m$corr_ulirsch_score )

# PCHi-C
m$pchic_jung_glo <- ifelse( is.na(m$pchic_jung_score), 0, m$pchic_jung_score )
m$pchic_jav_glo  <- ifelse( is.na(m$pchic_javierre_score), 0, m$pchic_javierre_score )

# CLPP
m$clpp_raw_l2g <- ifelse( is.na(m$clpp_prob), 0, m$clpp_prob )
m$clpp_glo  <- ifelse( is.na(m$clpp_prob), log10(0.0001), log10(m$clpp_prob) )

# ABC
m$abc_raw_l2g <- ifelse( is.na(m$abc_score), 0, m$abc_score )
m$abc_glo  <- ifelse( is.na(m$abc_score) | m$abc_score<0.01, 
                      log10(0.01), log10(m$abc_score) )

# MAGMA
m$magma_glo <- ifelse( is.na(m$magma_score), 
                       median( abs(m$magma_score), na.rm=TRUE), 
                       abs(m$magma_score) )
m$magma_glo[ m$magma_glo > 10 ] <- 10

# SMR
m$smr_glo <- ifelse( is.na(m$smr_p),
                     median( -log10(m$smr_p), na.rm=TRUE ), 
                     -log10(m$smr_p) )
m$smr_glo[ m$smr_glo > 10 ] <- 10

# Coding
m$coding_glo <- ifelse( is.na(m$coding_prob) | logit10(m$coding_prob) < log10(10^-3), 
                        log10(10^-3), 
                        logit10(m$coding_prob) )
m$coding     <- ifelse( is.na(m$coding_prob), 0, 1 )

# DEPICT
m$depict_z_glo <- p_to_z(m$no_loco_nc_depict_p)
m$depict_z_glo[ is.na(m$depict_z_glo) ] <- 0

# NetWAS
m$netwas_score_glo     <- m$no_loco_nc_netwas_score
m$netwas_bon_score_glo <- m$no_loco_nc_netwas_bon_score
m$netwas_score_glo[     is.na(m$netwas_score_glo) ]     <- 0
m$netwas_bon_score_glo[ is.na(m$netwas_bon_score_glo) ] <- 0


#-------------------------------------------------------------------------------
#   Add columns for relative features
#-------------------------------------------------------------------------------

m$pops_rel                <- as.numeric(NA)
m$dist_gene_rel           <- as.numeric(NA)
m$dist_tss_rel            <- as.numeric(NA)
m$twas_rel                <- as.numeric(NA)
m$corr_liu_rel            <- as.numeric(NA)
m$corr_and_rel            <- as.numeric(NA)
m$corr_uli_rel            <- as.numeric(NA)
m$pchic_jung_rel          <- as.numeric(NA)
m$pchic_jav_rel           <- as.numeric(NA)
m$clpp_rel                <- as.numeric(NA)
m$abc_rel                 <- as.numeric(NA)
m$magma_rel               <- as.numeric(NA)
m$smr_rel                 <- as.numeric(NA)
m$coding_rel              <- as.numeric(NA)
m$depict_z_rel            <- as.numeric(NA)
m$netwas_score_rel        <- as.numeric(NA)
m$netwas_bon_score_rel    <- as.numeric(NA)
tcp_m <- paste( m$trait, m$region, m$cs_id, sep="_" )
for( i in unique(tcp_m) ){
  
  # Subset
  sub <- m[ tcp_m == i , ]
  
  # Compute L2G-like values
  pops_vals     <- sub$pops_glo             - max(sub$pops_glo)
  gene_vals     <- log10( sub$distance_genebody - min(sub$distance_genebody) + 1e3 )
  tss_vals      <- log10( sub$distance_tss      - min(sub$distance_tss) + 1e3 )
  twas_vals     <- sub$twas_glo             - max(sub$twas_glo)
  liu_vals      <- sub$corr_liu_glo         - max(sub$corr_liu_glo)
  and_vals      <- sub$corr_and_glo         - max(sub$corr_and_glo)
  uli_vals      <- sub$corr_uli_glo         - max(sub$corr_uli_glo)
  jung_vals     <- sub$pchic_jung_glo       - max(sub$pchic_jung_glo)
  jav_vals      <- sub$pchic_jav_glo        - max(sub$pchic_jav_glo)
  clpp_vals     <- sub$clpp_glo             - max(sub$clpp_glo)
  abc_vals      <- sub$abc_glo              - max(sub$abc_glo)
  magma_vals    <- sub$magma_glo            - max(sub$magma_glo)
  smr_vals      <- sub$smr_glo              - max(sub$smr_glo)
  coding_vals   <- sub$coding_glo           - max(sub$coding_glo)
  depict_vals   <- sub$depict_z_glo         - max(sub$depict_z_glo)
  netwas_vals   <- sub$netwas_score_glo     - max(sub$netwas_score_glo)
  netwas_b_vals <- sub$netwas_bon_score_glo - max(sub$netwas_bon_score_glo)
  
  # Replace: POPS
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "pops_rel", 
       value = pops_vals )
  
  # Replace: gene body
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "dist_gene_rel", 
       value = gene_vals )
  
  # Replace: TSS
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "dist_tss_rel", 
       value = tss_vals )
  
  # Replace: TWAS
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "twas_rel", 
       value = twas_vals )
  
  # Replace: E-P Liu
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "corr_liu_rel", 
       value = liu_vals )
  
  # Replace: E-P Andersson
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "corr_and_rel", 
       value = and_vals )
  
  # Replace: E-P Ulirsch
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "corr_uli_rel", 
       value = uli_vals )
  
  # Replace: Jung PCHi-C
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "pchic_jung_rel", 
       value = jung_vals )
  
  # Replace: Javierre PCHi-C
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "pchic_jav_rel", 
       value = jav_vals )
  
  # Replace: CLPP
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "clpp_rel", 
       value = clpp_vals )
  
  # Replace: ABC
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "abc_rel", 
       value = abc_vals )
  
  # Replace: MAGMA z
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "magma_rel", 
       value = magma_vals )
  
  # Replace: SMR
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "smr_rel", 
       value = smr_vals )
  
  # Replace: coding
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "coding_rel", 
       value = coding_vals )
  
  # Replace: DEPICT
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "depict_z_rel", 
       value = depict_vals )
  
  # Replace: NetWAS score
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "netwas_score_rel", 
       value = netwas_vals )
  
  # Replace: NetWAS bon score
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "netwas_bon_score_rel", 
       value = netwas_b_vals )
}

# Add BIL
m$depict_z_bil         <- m$depict_z_rel         == 0
m$netwas_score_bil     <- m$netwas_score_rel     == 0
m$netwas_bon_score_bil <- m$netwas_bon_score_rel == 0

# Look at the first TCP
m[ tcp_m == unique(tcp_m)[1] , ]


#-------------------------------------------------------------------------------
#   Modify rank variables to reflect what was used in the paper
#-------------------------------------------------------------------------------

# If TWAS P is not Bonferroni-significant, set rank to NA
bonf <- 0.05/2e4
set( x = m, 
     i = which( m$twas_p > bonf ), 
     j = "twas_rank", 
     value = NA )

# If CLPP <= 0.1, set rank to NA
set( x = m, 
     i = which( m$clpp_prob <= 0.1 ), 
     j = "clpp_rank", 
     value = NA )

# Create new columns for whether a given TGP has best-in-locus support from a given method
rank_cols <- grep( pattern="rank$", x=names(m), value=TRUE )
bil_colnames <- sub( pattern="_rank$", replacement="_bil", x=rank_cols )
bil_colnames <- sub( pattern="javierre", replacement="jav", x=bil_colnames )
bil_colnames <- sub( pattern="andersson", replacement="and", x=bil_colnames )
bil_colnames <- sub( pattern="ulirsch", replacement="uli", x=bil_colnames )
bil_colnames <- sub( pattern="distance", replacement="dist", x=bil_colnames )
bil_colnames <- sub( pattern="genebody", replacement="gene", x=bil_colnames )
empty_bil <- as.data.table( matrix( nrow=NROW(m), ncol=length(bil_colnames) ) )
names(empty_bil) <- bil_colnames
m2 <- cbind( m, empty_bil )

# For each V2G method,
# Assign rank 1 TGP as having V2G support and all other TGPs as not
for( j in seq_along(bil_colnames) ){
  
  # Subset the data
  vg_col  <- bil_colnames[j]
  rank_col <- rank_cols[j]
  
  # Assign the lowest-ranked TGP as having V2G support and all other TGPs as not
  has_bil <- m2[[rank_col]] == 1 & !is.na(m2[[rank_col]])
  set( x=m2, j=vg_col, value=has_bil )
}

# Add columns for any V2G across E-P correlations
# And again for PC-HiC
ep_cols    <- grep( pattern="^corr_.*_bil$",  x=names(m2), value=TRUE )
pchic_cols <- grep( pattern="^pchic_.*_bil$", x=names(m2), value=TRUE )
m2$corr_any_bil  <- apply( X=m2[ , ..ep_cols ],    MARGIN=1, FUN=any )
m2$pchic_any_bil <- apply( X=m2[ , ..pchic_cols ], MARGIN=1, FUN=any )

# Add columns for POPS+local and POPS+nearest_gene, like in the original paper
m2$pops_and_local <- m2$pops_bil & 
  ( m2$dist_gene_bil | m2$magma_bil | m2$twas_bil | m2$clpp_bil | m2$abc_bil | 
      m2$corr_any_bil | m2$pchic_any_bil | m2$smr_bil )
m2$pops_and_nearest <- m2$pops_bil & m2$dist_gene_bil


#-------------------------------------------------------------------------------
#   Add protein attenuation level as a covariate
#-------------------------------------------------------------------------------

# Read in the data
pa <- fread( file.path( maindir, "protein_attenuation.csv" ) )

# Merge
pa_cols <- c( "gene", "prot_att", "prot_att_class" )
m3 <- left_join( x=m2, y=pa[,..pa_cols] )
m3$prot_att_miss <- ifelse( is.na(m3$prot_att_class), 1, 0 )
set( x     = m3, 
     i     = which( is.na(m3$prot_att) ), 
     j     = "prot_att", 
     # value = 0.093843 )
     value = median( m3$prot_att, na.rm=TRUE ) )
m3$prot_att_poly2 <- m3$prot_att * m3$prot_att
set( x     = m3, 
     i     = which( is.na(m3$prot_att_class) ), 
     j     = "prot_att_class", 
     value = "d_miss" )


#-------------------------------------------------------------------------------
#   Add King et al. 2019 gene-level covariates
#-------------------------------------------------------------------------------

# Read in the data
king <- fread( file.path( maindir, "king_2019_gene_covariates.tsv" ) )

# Merge
king_cols <- c( "ensgid", "rvis" )
m4 <- left_join( x=m3, y=king[,..king_cols] )
m4$rvis_miss <- ifelse( is.na(m4$rvis), 1, 0 )
set( x     = m4, 
     i     = which( is.na(m4$rvis) ), 
     j     = "rvis", 
     # value = -0.603598 )
     value = median( m4$rvis, na.rm=TRUE ) )
m4$rvis4 <- ifelse( abs(m4$rvis) > 4, 4*sign(m4$rvis), m4$rvis )
m4$rvis4_poly2 <- m4$rvis4^2


#-------------------------------------------------------------------------------
#   Write to file
#-------------------------------------------------------------------------------

# Write
merged_outfile <- file.path( maindir, "causal_noncausal_trait_gene_pairs", 
                             "causal_tgp_and_gene_mapping_data_300kb.tsv" )
fwrite( x=m4, file=merged_outfile, sep="\t" )


#-------------------------------------------------------------------------------
#   Done
#-------------------------------------------------------------------------------




