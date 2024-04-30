

#-------------------------------------------------------------------------------
#   simple_functions
#-------------------------------------------------------------------------------

# irnt: quantile-normalize a variable
irnt <- function(x){
  qnorm( ( rank( x, na.last="keep" ) -0.5 ) / sum( !is.na(x) ) )
}

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
                       "causal_noncausal_trait_gene_pairs_300kb_pubmed.tsv" )
cnc <- fread(cnc_file)

# Read in V2G
vg_file1 <- file.path( maindir, "UKB_AllMethods_GenePrioritization.txt.gz" )
vg1 <- fread(vg_file1)
vg1 <- vg1[ order( vg1$trait, vg1$region, vg1$cs_id, vg1$ensgid ) , ]

# Read in extended V2G
vg_file2 <- file.path( maindir, "PoPS_UKBB_noncoding_validation_1348CSs_v2.txt.gz" )
vg2 <- fread(vg_file2)
vg2 <- vg2[ order( vg2$trait, vg2$region, vg2$cs_id, vg2$ensgid, -vg2$pip ) , ]
idx <- duplicated( vg2[ , c( "trait", "region", "cs_id", "ensgid" ) ] )
vg2 <- vg2[ !idx , ]
gcols_vg2 <- c( "trait", "region", "cs_id", "ensgid", "smr_p", "smr_rank", 
                 "magma_score", "magma_rank" )

# Read in promoter SNP V2G
vg_file3 <- file.path( maindir, "promoter_v2g.tsv" )
vg3 <- fread(vg_file3)

# Read in rare variant burden test V2G
vg_file4 <- file.path( maindir, "rare_variant_burden.tsv" )
vg4 <- fread(vg_file4)
names(vg4) <- tolower( names(vg4) )
idx <- duplicated( vg4[ , c( "trait", "gene" ) ] )
table(idx)
vg4 <- vg4[ !idx , ]
gcols_vg4 <- c( "trait", "gene", "burden_ms", "burden_plof", 
                 "skato_ms", "skato_plof" )

# Read in DEPICT and NetWAS V2G
vg_file5 <- file.path( maindir, "netwas_depict.txt.gz" )
vg5 <- fread(vg_file5)
vg5 <- vg5[ !grepl( pattern="^PASS_", x=vg5$trait ) , ]
names(vg5)[ names(vg5) == "ENSGID" ] <- "ensgid"
vg5$trait <- sub( pattern="^UKB_", replacement="", x=vg5$trait )
vg5

# Merge
j1 <- left_join( x=vg1, y=vg2[,..gcols_vg2] )
j2 <- left_join( x=j1, y=vg3 )
j3 <- left_join( x=j2, y=vg4[,..gcols_vg4] )
j4 <- left_join( x=j3, y=vg5 )

# Remove genes >300kb from the lead variant
j5 <- j4[ j4$distance_genebody < 300e3 , ]
NROW(vg1); NROW(vg2); NROW(vg3); NROW(vg4); NROW(vg5) 
NROW(j1); NROW(j5) 

# Subset both to shared trait-ENSGID-CS triplets
# Note: yields about half as many missed causal/non-causal TGPs than using gene
# Note: that match() will only select the first instance for the small...
# ...number of cases where there are multiple ENSGIDs per gene
tec_cnc <- paste( cnc$trait,  cnc$ensgid,  cnc$cs_id,  sep="_" )
tec_vg <- paste( j5$trait, j5$ensgid, j5$cs_id, sep="_" )
tec_shared <- intersect( tec_cnc, tec_vg )
cnc2 <- cnc[  match( tec_shared, tec_cnc ) , ]
j6 <- j5[ match( tec_shared, tec_vg ) , ]

# Explore: What percentage of triplets were shared?
length(tec_shared) /  length(tec_cnc)
length(tec_shared) /  length(tec_vg)

# Explore: Do the datasets have the same starts, and ends?
# About half of the time they differ (different references?)
table( cnc2$gene_start  == j6$start )
table( cnc2$gene_end    == j6$end )

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
               "magma_score", "magma_rank", "promoter_pip", 
               "burden_ms", "burden_plof", "skato_ms", "skato_plof",
               "full_loco_netwas_score", "full_loco_netwas_bon_score",
               "no_loco_nc_netwas_score", "no_loco_nc_netwas_bon_score", 
               "no_loco_nc_depict_p", "full_loco_depict_p" )
gcols_cnc <- c( "pubmed", "pubmed_trait", "pubmed_gene",
                "n_cs_in_region", "n_cs_snps", "cs_width", "ngenes_nearby" )
m <- cbind( causal=cnc2$causal, j6[ , ..gcols_vg ], cnc2[ , ..gcols_cnc ] )

# Explore: Loop through TCPs and check whether locus-based ranks are still good
rank_cols <- grep( pattern="rank$", x=names(m), value=TRUE )
tcp_m <- paste( m$trait, m$region, m$cs_id, sep="_" )
mat <- matrix( nrow = length( unique(tcp_m) ),
               ncol = length(rank_cols) )
for( i in seq_along( unique(tcp_m) ) ){
  sub <- m[ tcp_m == tcp_m[i] , ..rank_cols ]
  for( j in seq_along(rank_cols) ){
    sub2 <- na.omit( sub[[j]] )
    mat[i,j] <- all( sort(sub2) == seq_along(sub2) )
  }
}
table(mat)

# Add a column for the number of genes in the locus
logit10 <- function(p) log10( p / (1-p) )
m$prior_n_genes_locus <- logit10( 1/m$ngenes_nearby )


#-------------------------------------------------------------------------------
#   PIP in gene, PIP-weighted distance to TSS (Weibull) converted to P(causal)
#-------------------------------------------------------------------------------

# # Read in eQTLGen TSS distances, extract density
# dte <- fread( file.path("~/projects/causal_genes/fauman_and_hyde_2022_supp_files/",
#                         "fauman_and_hyde_2022_supp_file3.txt") )
# fit.w  <- fitdist( data=dte2$dist, distr="weibull" )

# Read in CS column names
# Read in CSs, add column names
cs_colnames <- fread( file.path( maindir, "release1.1", "UKBB_94traits_release1.cols" ),
                      header=FALSE )
cs <- fread( file.path( maindir, "release1.1", "UKBB_94traits_release1.bed.gz" ) )
names(cs) <- cs_colnames$V1

# Remove SNPs that are not in CSs
# Subset to SUSIE CSs only
# Add P value and TCP columns
cs2 <- cs[  cs$cs_id > 0 , ]
cs3 <- cs2[ cs2$method == "SUSIE" , ]
cs3$p <- z_to_p( z=sqrt(cs3$chisq_marginal) )
cs3$tcp <- paste( cs3$trait, cs3$region, cs3$cs_id, sep="_" )
head( cs3, 1 )

# For each TGP
# m$dist_tss_wei <- as.numeric(NA)
m$dist_gene_pip <- as.numeric(NA)
for( i in unique(tcp_m) ){

  # Subset to the correct CS
  cs4 <- cs3[ cs3$tcp == i , c( "rsid", "end", "pip" ) ]

  # Subset to the genes in the correct locus
  loc <- m[ tcp_m == i , c( "gene", "start", "end", "tss" ) ]

  # # Compute the distance to the TSS of all genes in the locus
  # d <- outer( X=loc$tss, Y=cs4$end, FUN=function(i,j) abs(i-j)  )
  # dimnames(d)[[1]] <- loc$gene
  # dimnames(d)[[2]] <- cs4$rsid
  # 
  # # Replace distance with the Weibull density fitted to eQTLGen data
  # # d2 <- dweibull( x=log10(d), shape=7.375, scale=4.7 )
  # d2 <- dweibull( x=d, shape=fit.w$estimate["shape"],
  #                 scale=fit.w$estimate["scale"] )
  # d2[ d2 > 2.5e-5 ] <- 2.5e-5
  # 
  # # Multiply Weibull densities by PIP and sum (matrix multiplication)
  # d3 <- d2 %*% cs4$pip
  # 
  # # Set
  # set( x     = m,
  #      i     = which( tcp_m == i ),
  #      j     = "dist_tss_wei",
  #      value = d3[,1] )

  # Is each SNP inside the boundaries of each gene?
  ge <- outer( X=loc$start, Y=cs4$end, FUN=function(i,j) ifelse( j >= i, 1, 0 )  )
  le <- outer( X=loc$end,   Y=cs4$end, FUN=function(i,j) ifelse( j <= i, 1, 0 )  )
  dimnames(ge)[[1]] <- dimnames(le)[[1]] <- loc$gene
  dimnames(ge)[[2]] <- dimnames(le)[[2]] <- cs4$rsid
  b <- ge*le

  # Multiply Weibull densities by PIP and sum (matrix multiplication)
  b2 <- b %*% cs4$pip

  # Set
  set( x     = m,
       i     = which( tcp_m == i ),
       j     = "dist_gene_pip",
       value = b2[,1] )
}
m$whole_cs_in_gene <- m$dist_gene_pip >= 0.95


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
# m$dist_tss_wei      <- dweibull( x=log10(m$distance_tss), shape=7.375, scale=4.7 )

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

# Promoter
m$promoter_glo <- logit10(m$promoter_pip)
m$promoter_glo[ is.na(m$promoter_glo) ] <- -3
m$promoter_glo[ m$promoter_glo > 0 ] <- 0

# PubMed
m$pubmed_glo <- ifelse( m$pubmed == 0, -1, log10(m$pubmed) )
m$pubmed0    <- m$pubmed == 0

# DEPICT
m$depict_z_glo <- p_to_z(m$no_loco_nc_depict_p)
m$depict_z_glo[ is.na(m$depict_z_glo) ] <- 0

# NetWAS
m$netwas_score_glo     <- m$no_loco_nc_netwas_score
m$netwas_bon_score_glo <- m$no_loco_nc_netwas_bon_score
m$netwas_score_glo[     is.na(m$netwas_score_glo) ]     <- 0
m$netwas_bon_score_glo[ is.na(m$netwas_bon_score_glo) ] <- 0


#-------------------------------------------------------------------------------
#   Add columns for IRNT features
#-------------------------------------------------------------------------------

m$pops_int           <- as.numeric(NA)
m$dist_gene_int      <- as.numeric(NA)
m$dist_tss_int       <- as.numeric(NA)
m$twas_int           <- as.numeric(NA)
m$magma_int          <- as.numeric(NA)
m$smr_int            <- as.numeric(NA)
m$corr_liu_int       <- as.numeric(NA)
m$coding_int         <- as.numeric(NA)
for( i in unique(m$trait) ){

  # Subset
  sub <- m[ m$trait == i , ]

  # Compute L2G-like values
  pops_vals   <- irnt(sub$pops_glo)
  gene_vals   <- irnt(sub$dist_gene_glo)
  tss_vals    <- irnt(sub$dist_tss_glo)
  twas_vals   <- irnt(sub$twas_glo)
  magma_vals  <- irnt(sub$magma_glo)
  smr_vals    <- irnt(sub$smr_glo)
  liu_vals    <- irnt(sub$corr_liu_glo)
  coding_vals <- irnt(sub$coding_glo)

  # Replace: POPS
  set( x     = m,
       i     = which( m$trait == i ),
       j     = "pops_int",
       value = pops_vals )

  # Replace: gene body
  set( x     = m,
       i     = which( m$trait == i ),
       j     = "dist_gene_int",
       value = gene_vals )

  # Replace: TSS
  set( x     = m,
       i     = which( m$trait == i ),
       j     = "dist_tss_int",
       value = tss_vals )

  # Replace: TWAS
  set( x     = m,
       i     = which( m$trait == i ),
       j     = "twas_int",
       value = twas_vals )

  # Replace: MAGMA
  set( x     = m,
       i     = which( m$trait == i ),
       j     = "magma_int",
       value = magma_vals )

  # Replace: SMR
  set( x     = m,
       i     = which( m$trait == i ),
       j     = "smr_int",
       value = smr_vals )

  # Replace: E-P Liu
  set( x     = m,
       i     = which( m$trait == i ),
       j     = "corr_liu_int",
       value = liu_vals )

  # Replace: coding
  set( x     = m,
       i     = which( m$trait == i ),
       j     = "coding_int",
       value = coding_vals )
}

# Look at the first 5 traits
plot(  density( m$pops_int[ m$trait == unique(m$trait)[1] ] ) )
lines( density( m$pops_int[ m$trait == unique(m$trait)[2] ] ) )
lines( density( m$pops_int[ m$trait == unique(m$trait)[3] ] ) )
lines( density( m$pops_int[ m$trait == unique(m$trait)[4] ] ) )
lines( density( m$pops_int[ m$trait == unique(m$trait)[5] ] ) )
lines( density( m$pops_glo[ m$trait == unique(m$trait)[1] ] ), col="steelblue" )
lines( density( m$pops_glo[ m$trait == unique(m$trait)[2] ] ), col="steelblue" )
lines( density( m$pops_glo[ m$trait == unique(m$trait)[3] ] ), col="steelblue" )
lines( density( m$pops_glo[ m$trait == unique(m$trait)[4] ] ), col="steelblue" )
lines( density( m$pops_glo[ m$trait == unique(m$trait)[5] ] ), col="steelblue" )


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
m$promoter_rel            <- as.numeric(NA)
m$pubmed_rel              <- as.numeric(NA)
m$depict_z_rel            <- as.numeric(NA)
m$netwas_score_rel        <- as.numeric(NA)
m$netwas_bon_score_rel    <- as.numeric(NA)
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
  promoter_vals <- sub$promoter_glo         - max(sub$promoter_glo)
  pubmed_vals   <- sub$pubmed_glo           - max(sub$pubmed_glo)
  depict_vals   <- sub$depict_z_glo         - max(sub$depict_z_glo)
  netwas_vals   <- sub$netwas_score_glo     - max(sub$netwas_score_glo)
  netwas_b_vals <- sub$netwas_bon_score_glo - max(sub$netwas_bon_score_glo)
  
  # Replace: POPS
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "pops_rel", 
       value = pops_vals )
  
  # Replace: gene body linear
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "dist_gene_rel", 
       value = gene_vals )
  
  # Replace: TSS linear
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "dist_tss_rel", 
       value = tss_vals )
  
  # Replace: TWAS linear
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "twas_rel", 
       value = twas_vals )
  
  # Replace: E-P Liu linear
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
  
  # Replace: CLPP linear
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "clpp_rel", 
       value = clpp_vals )
  
  # Replace: ABC linear
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
  
  # Replace: promoter
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "promoter_rel", 
       value = promoter_vals )
  
  # Replace: PubMed
  set( x     = m, 
       i     = which( tcp_m == i ), 
       j     = "pubmed_rel", 
       value = pubmed_vals )
  
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
m$promoter_bil         <- m$promoter_rel         == 0
m$pubmed_bil           <- m$pubmed_rel           == 0
m$depict_z_bil         <- m$depict_z_rel         == 0
m$netwas_score_bil     <- m$netwas_score_rel     == 0
m$netwas_bon_score_bil <- m$netwas_bon_score_rel == 0

# Look at the first CS
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
bil_colnames <- sub( pattern="_rank$", replacement="_bil", x=rank_cols )
bil_colnames <- sub( pattern="javierre", replacement="jav", x=bil_colnames )
bil_colnames <- sub( pattern="andersson", replacement="and", x=bil_colnames )
bil_colnames <- sub( pattern="ulirsch", replacement="uli", x=bil_colnames )
bil_colnames <- sub( pattern="distance", replacement="dist", x=bil_colnames )
bil_colnames <- sub( pattern="genebody", replacement="gene", x=bil_colnames )
empty_bil <- as.data.table( matrix( nrow=NROW(m), ncol=length(bil_colnames) ) )
names(empty_bil) <- bil_colnames
m2 <- cbind( m, empty_bil )

# For each vg method,
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
# king_cols <- names(king)
king_cols <- c( "ensgid", "rvis" )
m4 <- left_join( x=m3, y=king[,..king_cols] )
m4$rvis_miss <- ifelse( is.na(m4$rvis), 1, 0 )
set( x     = m4, 
     i     = which( is.na(m4$rvis) ), 
     j     = "rvis", 
     # value = -0.603598 )
     value = median( m4$rvis, na.rm=TRUE ) )
m4$rvis7 <- ifelse( abs(m4$rvis) > 7, 7*sign(m4$rvis), m4$rvis )
m4$rvis4 <- ifelse( abs(m4$rvis) > 4, 4*sign(m4$rvis), m4$rvis )
m4$rvis4_poly2 <- m4$rvis4^2


#-------------------------------------------------------------------------------
#   Add pLI
#-------------------------------------------------------------------------------

# Read in the data
psyops_file <- "~/repos/PsyOPS/constraint_scores.tsv"
po <- fread(psyops_file)
names(po) <- c( "gene", "pLI", "oe", "loeuf", "cds_length", 
                "chrom", "start", "end" )

# Merge
po_cols <- c( "gene", "pLI" )
po2 <- po[ order(-po$pLI) , ..po_cols ]
po3 <- po2[ !duplicated(po2$gene) , ]
m5 <- left_join( x=m4, y=po3, by="gene" )
m5$pLI0.99 <- m5$pLI > 0.99 & !is.na(m5$pLI)


#-------------------------------------------------------------------------------
#   Add burden results and covariate
#-------------------------------------------------------------------------------

m5$burden <- ifelse( m5$burden_plof < 0.05/2e5 & !is.na(m5$burden_plof), TRUE, FALSE )
burden_prop  <- tapply( m5$burden, m5$trait, mean )
burden_prop2 <- data.frame( trait=names(burden_prop), burden_prop=burden_prop )
m6           <- left_join( x=m5, y=burden_prop2, by="trait" )


#-------------------------------------------------------------------------------
#   Add PubMed evidence and covariates
#-------------------------------------------------------------------------------

# Add a binary feature for having any PubMed hits
# Add a continuous feature for log10(PubMed hits) with NAs converted to 0's
# (to preserve the meaning of the intercept)
m6$pubmed_any <- m6$pubmed > 0
m6$pubmed_log <- ifelse( m6$pubmed_any, log10(m6$pubmed), 0 )

# Add a trait-level covariate for the proportion of genes with PubMed evidence
pm_trait_prop  <- tapply( m6$pubmed_any, m6$trait, mean )
pm_trait_prop2 <- data.frame( trait             = names(pm_trait_prop), 
                              pubmed_trait_prop = pm_trait_prop )
m7             <- left_join( x=m6, y=pm_trait_prop2, by="trait" )

# Add a trait-level covariate for the mean log10(PubMed hits)
pm_log_mean  <- tapply( m7$pubmed_log[m7$pubmed_any], 
                         m7$trait[m7$pubmed_any], mean )
pm_log_mean2 <- data.frame( trait           = names(pm_log_mean), 
                            pubmed_log_mean = pm_log_mean )
m8           <- left_join( x=m7, y=pm_log_mean2, by="trait" )
m8$pubmed_log_mean[ is.na(m8$pubmed_log_mean) ] <- 0


#-------------------------------------------------------------------------------
#   Write to file
#-------------------------------------------------------------------------------

# Write
merged_outfile <- file.path( maindir, "causal_noncausal_trait_gene_pairs", 
                             "causal_tgp_and_gene_mapping_data_300kb.tsv" )
fwrite( x=m8, file=merged_outfile, sep="\t" )


#-------------------------------------------------------------------------------
#   Done
#-------------------------------------------------------------------------------




