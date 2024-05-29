

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


#-------------------------------------------------------------------------------
#   Add columns for global features
#-------------------------------------------------------------------------------

# n_genes
m$prior_n_genes_locus <- logit10( 1/m$ngenes_nearby )
tcp_m <- paste( m$trait, m$region, m$cs_id, sep="_" )

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
#   Add columns for relative and best-in-locus features
#-------------------------------------------------------------------------------

# Set up global feature columns and TCPs
glo_cols <- grep( "_glo$", x=names(m), value=T)
tcp_m    <- paste( m$trait, m$region, m$cs_id, sep="_" )

# Loop through global columns
for( j in glo_cols ){
  
  # Initialize relative and best-in-locus columns
  rel_col <- sub( pattern="_glo$", replacement="_rel", x=j )
  bil_col <- sub( pattern="_glo$", replacement="_bil", x=j )
  m[[rel_col]] <- as.numeric(NA)
  m[[bil_col]] <- as.logical(NA)
  
  # Loop through TCPs
  for( i in unique(tcp_m) ){
    
    # Subset
    sub <- m[ tcp_m == i , ]
    
    # Compute relative and best-in-locus values
    r_vals <- sub[[j]] - max( sub[[j]] )
    b_vals <- r_vals == 0 & sum( r_vals == 0 ) == 1
    
    # Insert relative values into the main data table
    set( x     = m, 
         i     = which( tcp_m == i ), 
         j     = rel_col, 
         value = r_vals )
    
    # Insert best-in-locus values into the main data table
    set( x     = m, 
         i     = which( tcp_m == i ), 
         j     = bil_col, 
         value = b_vals ) 
  }
}


#-------------------------------------------------------------------------------
#   Add some final V2G columns
#-------------------------------------------------------------------------------

# Add BIL
rel_cols <- grep( pattern="_rel$", x=names(m), value=TRUE )
for( i in rel_cols ){
  bil_col <- sub( pattern="_rel$", replacement="_bil", x=i )
  m[[bil_col]] <- m[[i]] == 0
}

# Add columns for any V2G across E-P correlations
# And again for PC-HiC
ep_cols    <- grep( pattern="^corr_.*_bil$",  x=names(m), value=TRUE )
pchic_cols <- grep( pattern="^pchic_.*_bil$", x=names(m), value=TRUE )
m$corr_any_bil  <- apply( X=m[ , ..ep_cols ],    MARGIN=1, FUN=any )
m$pchic_any_bil <- apply( X=m[ , ..pchic_cols ], MARGIN=1, FUN=any )

# Add columns for POPS+local and POPS+nearest_gene, like in the original paper
m$pops_and_local <- m$pops_bil & 
  ( m$dist_gene_bil | m$magma_bil | m$twas_bil | m$clpp_bil | m$abc_bil | 
      m$corr_any_bil | m$pchic_any_bil | m$smr_bil )
m$pops_and_nearest <- m$pops_bil & m$dist_gene_bil

# Look at the first TCP
m[ tcp_m == unique(tcp_m)[1] , ]


#-------------------------------------------------------------------------------
#   Add Mostafavi et al. 2023 gene-level covariates
#-------------------------------------------------------------------------------

# Read in the data, merge
mo_file <- "~/projects/causal_genes/mostafavi2023_gene_annots/pc_genes.txt"
mo <- fread(mo_file)
names(mo)[ names(mo) == "GeneSymbol" ] <- "ensgid"
mo$gene <- NULL
m2 <- left_join( x=m, y=mo, by="ensgid" )

# Missingness
m2$pritchard_miss <- ifelse( is.na(m2$length), 1, 0 )

# pLI
m2$pLI[ is.na(m2$pLI) ] <- 0.5
m2$pLI_lt_0.9 <- ifelse( m2$pLI < 0.9, 1, 0 )
m2$pLI_lt_0.1 <- ifelse( m2$pLI < 0.1, 1, 0 )
m2$pLI_log10 <- -log10(m2$pLI)
pLI99 <- quantile( m2$pLI_log10, probs=0.99 )
m2$pLI_log10 <- ifelse( m2$pLI_log10 > pLI99, pLI99, m2$pLI_log10 )

# LOEUF
m2$LOEUF[ is.na(m2$LOEUF) ] <- 0

# hs
m2$hs[ is.na(m2$hs) ] <- max( m2$hs, na.rm=TRUE )
m2$hs_log10 <- -log10(m2$hs)

# Gene length
m2$length[ m2$pritchard_miss == 1 ] <- min( m2$length, na.rm=TRUE )
m2$gene_bp_log10 <- log10( m2$length*1e3 )

# CDS length
m2$CDS_length[ m2$pritchard_miss == 1 ] <- min( m2$CDS_length, na.rm=TRUE )
m2$cds_bp_log10 <- log10( m2$CDS_length*1e3 )

# ABC_count
m2$ABC_count[ m2$pritchard_miss == 1 ] <- min( m2$ABC_count, na.rm=TRUE )

# ABC_length_per_type
m2$ABC_length_per_type[ m2$pritchard_miss == 1 ] <- min( m2$ABC_length_per_type, na.rm=TRUE )
m2$abc_bp_log10 <- ifelse( m2$ABC_length_per_type == 0, 
                           log10( 0.2 * 1e3 ),
                           log10( m2$ABC_length_per_type * 1e3 ) )

# Roadmap_count
m2$Roadmap_count[ m2$pritchard_miss == 1 ] <- min( m2$Roadmap_count, na.rm=TRUE )

# Roadmap_length_per_type
m2$Roadmap_length_per_type[ m2$pritchard_miss == 1 ] <- min( m2$Roadmap_length_per_type, na.rm=TRUE )
m2$roadmap_bp_log10 <- ifelse( m2$Roadmap_length_per_type == 0, 
                               log10( min( m2$Roadmap_length_per_type[ 
                                 m2$Roadmap_length_per_type > 0 ] ) * 1e3 ),
                               log10( m2$Roadmap_length_per_type * 1e3 ) )

# promoter_count
m2$promoter_count[ m2$pritchard_miss == 1 ] <- min( m2$promoter_count, na.rm=TRUE )
m2$promoter_count_log10 <- log10( m2$promoter_count + 1 )

# connect_decile
m2$connect_decile[ m2$pritchard_miss == 1 ] <- min( m2$connect_decile, na.rm=TRUE )

# connect_quantile
m2$connect_quantile[ m2$pritchard_miss == 1 ] <- min( m2$connect_quantile, na.rm=TRUE )

# connectedness
m2$connectedness[ m2$pritchard_miss == 1 ] <- min( m2$connectedness, na.rm=TRUE )

# PPI_degree_decile
m2$PPI_degree_decile[ m2$pritchard_miss == 1 ] <- min( m2$PPI_degree_decile, na.rm=TRUE )

# PPI_degree_quantile
m2$PPI_degree_quantile[ m2$pritchard_miss == 1 ] <- min( m2$PPI_degree_quantile, na.rm=TRUE )

# PPI_degree_cat
m2$PPI_degree_cat[ m2$pritchard_miss == 1 ] <- min( m2$PPI_degree_cat, na.rm=TRUE )

# TF
m2$TF[ m2$pritchard_miss == 1 ] <- min( m2$TF, na.rm=TRUE )


#-------------------------------------------------------------------------------
#   Write to file
#-------------------------------------------------------------------------------

# Write
merged_outfile <- file.path( maindir, "causal_noncausal_trait_gene_pairs", 
                             "causal_tgp_and_gene_mapping_data_300kb.tsv" )
fwrite( x=m2, file=merged_outfile, sep="\t" )


#-------------------------------------------------------------------------------
#   Done
#-------------------------------------------------------------------------------





















