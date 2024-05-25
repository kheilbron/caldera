

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
suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(dplyr) )

# Read in gencode gene boundaries
gencode  <- fread("~/projects/causal_genes/gene_locations.tsv")

# Read in causal/non-causal trait-gene pairs
cnc_file <- file.path( maindir, "causal_noncausal_trait_gene_pairs", 
                       "causal_noncausal_trait_gene_pairs_300kb.tsv" )
cnc <- fread(cnc_file)

# Read in V2G
vg_file <- file.path( maindir, "v2g_file.tsv" )
vg <- fread(vg_file)

# Subset both to shared trait-ENSGID-CS triplets (TECs)
# Note: yields about half as many missed causal/non-causal TGPs than using gene
tec_cnc <- paste( cnc$trait, cnc$ensgid, cnc$cs_id, sep="_" )
tec_vg  <- paste( vg$trait,  vg$ensgid,  vg$cs_id,  sep="_" )
tec_shared <- intersect( tec_cnc, tec_vg )
cnc2 <- cnc[ match( tec_shared, tec_cnc ) , ]
vg2  <- vg[  match( tec_shared, tec_vg ) , ]

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
gcols_cnc <- c( "n_cs_snps", "cs_width", "ngenes_nearby" )
m <- cbind( causal=cnc2$causal, vg2[ , ..gcols_vg ], cnc2[ , ..gcols_cnc ] )

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
m$dist_gene_glo     <- -log10( m$distance_genebody + 1e3 )

# Distance to TSS
m$dist_tss_raw_l2g  <- log( m$distance_tss + 1 )
m$dist_tss_glo      <- -log10( m$distance_tss + 1e3 )

# TWAS
m$twas_p <- ifelse( m$twas_p == 0, .Machine$double.xmin, m$twas_p )
m$twas_logp_glo <- ifelse( is.na(m$twas_p), -log10(1), -log10(m$twas_p) )
m$twas_glo <- ifelse( is.na(m$twas_z), 0, abs(m$twas_z) )

# E-P correlation
m$corr_liu_raw_l2g <- ifelse( is.na(m$corr_liu_score), 
                              min( m$corr_liu_score, na.rm=TRUE ), 
                              m$corr_liu_score )
m$corr_liu_glo  <- log10(m$corr_liu_raw_l2g)
m$corr_and_glo  <- ifelse( is.na(m$corr_andersson_score), 0, m$corr_andersson_score )
m$corr_uli_glo  <- ifelse( is.na(m$corr_ulirsch_score), 0, m$corr_ulirsch_score )

# PCHi-C
m$pchic_jung_glo <- ifelse( is.na(m$pchic_jung_score), 0, m$pchic_jung_score )
m$pchic_jav_glo  <- ifelse( is.na(m$pchic_javierre_score), 0, m$pchic_javierre_score )

# CLPP
m$clpp_raw_l2g <- ifelse( is.na(m$clpp_prob), 0, m$clpp_prob )
m$clpp_glo     <- ifelse( is.na(m$clpp_prob), 
                          log10( min( m$clpp_prob, na.rm=TRUE ) ), 
                          log10(m$clpp_prob) )

# ABC
m$abc_raw_l2g <- ifelse( is.na(m$abc_score), 0, m$abc_score )
m$abc_glo  <- ifelse( is.na(m$abc_score) | m$abc_score == 0, 
                      log10( min( m$abc_score[ m$abc_score > 0 ], na.rm=TRUE ) ), 
                      log10(m$abc_score) )

# MAGMA
m$magma_glo <- ifelse( is.na(m$magma_score), 
                       median( m$magma_score, na.rm=TRUE ), 
                       m$magma_score )

# SMR
m$smr_glo <- ifelse( is.na(m$smr_p),
                     -log10( max( m$smr_p, na.rm=TRUE ) ), 
                     -log10(m$smr_p) )

# Coding
m$coding_glo <- ifelse( is.na(m$coding_prob), 
                        logit10( min( m$coding_prob, na.rm=TRUE ) ), 
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
  gene_vals     <- -log10( sub$distance_genebody - min(sub$distance_genebody) + 1e3 )
  tss_vals      <- -log10( sub$distance_tss      - min(sub$distance_tss) + 1e3 )
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
#   Add Mostafavi et al. 2023 gene-level covariates
#-------------------------------------------------------------------------------

# Read in the data, merge
mo_file <- "~/projects/causal_genes/mostafavi2023_gene_annots/pc_genes.txt"
mo <- fread(mo_file)
names(mo)[ names(mo) == "GeneSymbol" ] <- "ensgid"
mo$gene <- NULL
m3 <- left_join( x=m2, y=mo, by="ensgid" )

# Missingness
m3$pritchard_miss <- ifelse( is.na(m3$length), 1, 0 )

# Gene length
m3$length[ m3$pritchard_miss == 1 ] <- min( m3$length, na.rm=TRUE )
m3$gene_bp_log10 <- log10( m3$length*1e3 )

# CDS length
m3$CDS_length[ m3$pritchard_miss == 1 ] <- min( m3$CDS_length, na.rm=TRUE )
m3$cds_bp_log10 <- log10( m3$CDS_length*1e3 )

# pLI
m3$pLI[ is.na(m3$pLI) ] <- 0.5
m3$pLI[ m3$pLI == 1 ] <- max( m3$pLI[ m3$pLI < 1 ])
m3$pLI_gt_0.9 <- ifelse( m3$pLI > 0.9, 1, 0 )
m3$pLI_lt_0.1 <- ifelse( m3$pLI < 0.1, 1, 0 )
m3$pLI_log10OR <- logit10(m3$pLI)
m3$pLI_log10OR <- ifelse( m3$pLI_log10OR < -20, -20, m3$pLI_log10OR )

# LOEUF
m3$LOEUF[ is.na(m3$LOEUF) ] <- 0

# ABC_count
m3$ABC_count[ m3$pritchard_miss == 1 ] <- min( m3$ABC_count, na.rm=TRUE )

# ABC_length_per_type
m3$ABC_length_per_type[ m3$pritchard_miss == 1 ] <- min( m3$ABC_length_per_type, na.rm=TRUE )
m3$abc_bp_log10 <- ifelse( m3$ABC_length_per_type == 0, 
                           log10( 0.2 * 1e3 ),
                           log10( m3$ABC_length_per_type * 1e3 ) )

# Roadmap_count
m3$Roadmap_count[ m3$pritchard_miss == 1 ] <- min( m3$Roadmap_count, na.rm=TRUE )

# Roadmap_length_per_type
m3$Roadmap_length_per_type[ m3$pritchard_miss == 1 ] <- min( m3$Roadmap_length_per_type, na.rm=TRUE )
m3$roadmap_bp_log10 <- ifelse( m3$Roadmap_length_per_type == 0, 
                               log10( min( m3$Roadmap_length_per_type[ 
                                 m3$Roadmap_length_per_type > 0 ] ) * 1e3 ),
                               log10( m3$Roadmap_length_per_type * 1e3 ) )

# promoter_count
m3$promoter_count[ m3$pritchard_miss == 1 ] <- min( m3$promoter_count, na.rm=TRUE )
m3$promoter_count_log10 <- log10( m3$promoter_count + 1 )

# connect_decile
m3$connect_decile[ m3$pritchard_miss == 1 ] <- min( m3$connect_decile, na.rm=TRUE )

# connect_quantile
m3$connect_quantile[ m3$pritchard_miss == 1 ] <- min( m3$connect_quantile, na.rm=TRUE )

# connectedness
m3$connectedness[ m3$pritchard_miss == 1 ] <- min( m3$connectedness, na.rm=TRUE )

# PPI_degree_decile
m3$PPI_degree_decile[ m3$pritchard_miss == 1 ] <- min( m3$PPI_degree_decile, na.rm=TRUE )

# PPI_degree_quantile
m3$PPI_degree_quantile[ m3$pritchard_miss == 1 ] <- min( m3$PPI_degree_quantile, na.rm=TRUE )

# PPI_degree_cat
m3$PPI_degree_cat[ m3$pritchard_miss == 1 ] <- min( m3$PPI_degree_cat, na.rm=TRUE )

# TF
m3$TF[ m3$pritchard_miss == 1 ] <- min( m3$TF, na.rm=TRUE )

# hs
m3$hs[ is.na(m3$hs) ] <- max( m3$hs, na.rm=TRUE )
m3$hs_log10 <- ifelse( m3$hs < 1e-5, log10(1e-5), log10(m3$hs) )


#-------------------------------------------------------------------------------
#   Write to file
#-------------------------------------------------------------------------------

# Write
merged_outfile <- file.path( maindir, "causal_noncausal_trait_gene_pairs", 
                             "causal_tgp_and_gene_mapping_data_300kb.tsv" )
fwrite( x=m3, file=merged_outfile, sep="\t" )


#-------------------------------------------------------------------------------
#   Done
#-------------------------------------------------------------------------------





















