

#-------------------------------------------------------------------------------
#   Find coding TGPs
#-------------------------------------------------------------------------------

# Load libraries
suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(dplyr) )

# Read in V2G
maindir <- file.path( Sys.getenv("HOME"), "projects/causal_genes/" )
v2g_file <- file.path( maindir, "UKB_AllMethods_GenePrioritization.txt.gz" )
v2g <- fread(v2g_file)
v2g$tcp <- paste( v2g$trait, v2g$region, v2g$cs_id, sep="_" )
NROW(v2g)
head( v2g, 2 )

# Extract coding TCPs and TGPs (coding PIP > 50%)
idx <- !is.na(v2g$coding_prob) & v2g$coding_prob > 0.5
tcp_c <- unique( v2g$tcp[idx] )
tgp_c <- v2g[ idx , ]
length(tcp_c)
head( tgp_c, 2 )


#-------------------------------------------------------------------------------
#   Format credible sets
#-------------------------------------------------------------------------------

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

# Add "number of CSs in region" to each CS to be used as a covariate later
trp <- paste( cs3$trait, cs3$region, sep="_" )
cs3$n_cs_in_region <- as.numeric(NA)
for( i in unique(trp) ){
  val <- length( unique( cs3$tcp[ trp == i ] ) )
  set( x     = cs3, 
       i     = which( trp == i ), 
       j     = "n_cs_in_region", 
       value = val )
}


#-------------------------------------------------------------------------------
#   Find non-coding TCPs near coding TGPs
#-------------------------------------------------------------------------------

# Extract non-coding TCPs
tcp_n <- setdiff( cs3$tcp, tcp_c )
length(tcp_n)

# Subset coding and non-coding TCPs to those in shared TRPs
# (to speed up computation)
trp_c <- sub( pattern="^(.*_.*)_([[:digit:]]+)$", replacement="\\1", x=tcp_c )
trp_n <- sub( pattern="^(.*_.*)_([[:digit:]]+)$", replacement="\\1", x=tcp_n )
trp_b <- intersect( trp_c, trp_n )
tcp_c2 <- tcp_c[ trp_c %in% trp_b ]
tcp_n2 <- tcp_n[ trp_n %in% trp_b ]
length(tcp_c); length(tcp_c2)
length(tcp_n); length(tcp_n2)

# Remove CSs >300kb from causal genes
# Loop through non-coding TCPs and flag those that are near a coding TGP
window <- 3e5
tcp_n_near0 <- list()
for( i in seq_along(tcp_n2) ){
  
  # Compute start and end boundaries
  sub      <- cs3[ cs3$tcp == tcp_n2[i] , ]
  trait    <- sub$trait[1]
  region   <- sub$region[1]
  cs_id    <- sub$cs_id[1]
  chr      <- as.integer( sub( "chr", "", sub$chromosome[1] ) )
  cs_start <- min(sub$start)
  cs_end   <- max(sub$end)
  n_cs     <- sub$n_cs_in_region[1]
  
  # Do any coding TGPs have?: same trait, same chromosome, 
  #   causal start < CS end + 500kb, causal end > CS start - 500kb
  idx <- tgp_c$trait == trait &
    tgp_c$chromosome == chr &
    tgp_c$start <= cs_end   + window &
    tgp_c$end   >= cs_start - window
  
  # If so, store information about the CS
  if( any(idx) ){
    dt <- data.table( tcp=tcp_n2[i], trait=trait, region=region, cs_id=cs_id, 
                      chr=chr, cs_start=cs_start, cs_end=cs_end, 
                      n_cs_snps=NROW(sub), cs_width=cs_end-cs_start,
                      n_cs_in_region=n_cs )
    tcp_n_near0[[i]] <- dt
  }
}
tcp_n_near <- do.call( rbind, tcp_n_near0 )
tcp_n_near <- tcp_n_near[ order( tcp_n_near$trait, 
                                 tcp_n_near$chr, 
                                 tcp_n_near$cs_start ) , ]
cs_n_near  <- cs3[ cs3$tcp %in% tcp_n_near$tcp , ]
NROW(tcp_n_near); length(tcp_n2); NROW(cs_n_near)
head(tcp_n_near)


#-------------------------------------------------------------------------------
#   Extract all genes near non-coding TCPs (near causal TGPs)
#-------------------------------------------------------------------------------

# Read in gencode gene boundaries
gencode  <- fread("~/projects/causal_genes/gene_locations.tsv")

# For each non-coding CS near causal TGPs, extract all nearby genes
all_genes_near_nc_cs0 <- list()
for( i in seq_along(tcp_n_near$tcp) ){
  
  # Extract all genes between:
  # The leftmost CS SNP minus "window" and the rightmost CS SNP plus "window"
  chr    <- tcp_n_near$chr[i]
  lmost  <- tcp_n_near$cs_start[i] - window
  rmost  <- tcp_n_near$cs_end[i]   + window
  gc <- gencode[ gencode$CHR == chr & 
                           gencode$END   >= lmost &
                           gencode$START <= rmost , ]
  
  dt <- cbind( tcp_n_near[i], gene=gc$NAME, ensgid=gc$ENSGID, 
               gene_start=gc$START, gene_end=gc$END, ngenes_nearby=NROW(gc) )
  all_genes_near_nc_cs0[[i]] <- dt
}
all_genes_near_nc_cs <- do.call( rbind, all_genes_near_nc_cs0 )
colorder <- c( "trait", "gene", "ensgid", "chr", "gene_start", "gene_end",
               "ngenes_nearby", "tcp", "region", "cs_id", "cs_start", "cs_end",
               "n_cs_snps", "cs_width", "n_cs_in_region" )
all_genes_near_nc_cs <- setcolorder( all_genes_near_nc_cs, colorder )
NROW(all_genes_near_nc_cs)


#-------------------------------------------------------------------------------
#   Assign TGPs as causal or non-causal using coding TGPs
#-------------------------------------------------------------------------------

# Assign all TGPs near non-coding CSs (near coding TGPs) as non-causal,...
cnc <- all_genes_near_nc_cs
cnc$causal <- FALSE

# Convert to causal for TEPs in the coding dataset
tep_a <- paste0( cnc$trait,   cnc$ensgid,   sep="_" )
tep_c <- paste0( tgp_c$trait, tgp_c$ensgid, sep="_" )
cnc$causal[ tep_a %in% tep_c ] <- TRUE
table(cnc$causal)
table( duplicated( cnc[ cnc$causal  , c( "trait", "ensgid" ) ] ) )
table( duplicated( cnc[ !cnc$causal , c( "trait", "ensgid" ) ] ) )

# Removed loci with < 2 genes in them
cnc2 <- cnc[ cnc$ngenes_nearby >= 2 , ]
NROW(cnc); NROW(cnc2)

# Write causal/non-causal trait-gene pairs to file
causal_gene_file <- file.path( maindir, "causal_noncausal_trait_gene_pairs",
                               "causal_noncausal_trait_gene_pairs_300kb.tsv" )
fwrite( x=cnc2, file=causal_gene_file, sep="\t" )


#-------------------------------------------------------------------------------
#   Write out coding genes, distal nearest genes, and distal NG+TP genes
#-------------------------------------------------------------------------------

# Non-coding TCPs that are far from fine-mapped (PIP > 50%) coding TCPs
tcp_n_far <- setdiff( tcp_n, tcp_n_near$tcp )
length(tcp_n); NROW(tcp_n_near); length(tcp_n_far)

# All (coding and non-coding) fine-mapped TCPs
tcp_a_fine <- unique( cs3$tcp[ cs3$pip > 0.5 ] )

# Non-coding fine-mapped TCPs that are far from fine-mapped coding TCPs
tcp_n_far_fine <- intersect( tcp_n_far, tcp_a_fine )
# tcp_n_far_fine <- tcp_n_far
length(tcp_n_far); length(tcp_a_fine); length(tcp_n_far_fine); length(tcp_c2)

# Find the nearest gene to each distal non-coding fine-mapped TCP
# If multiple genes are equal, take the one with the higher POPS score
ng_n_far_fine0 <- v2g[ v2g$tcp %in% tcp_n_far_fine & 
                       v2g$distance_genebody_rank == 1 , ]
ng_n_far_fine0 <- ng_n_far_fine0[ order( ng_n_far_fine0$tcp, -ng_n_far_fine0$pops_score ) , ]
ng_n_far_fine  <- ng_n_far_fine0[ !duplicated(ng_n_far_fine0$tcp) , ]
NROW(ng_n_far_fine0); NROW(ng_n_far_fine)

# Is there any overlap between these distal non-coding causal genes...
# ...and the coding causal genes (for other traits)?
length( unique(ng_n_far_fine$gene) )
length( unique( cnc2$gene[ cnc2$causal ] ) )
length( intersect( ng_n_far_fine$gene, cnc2$gene[cnc2$causal] ) )

# Establish the set of unique high-PIP coding genes
coding_genes0 <- cnc2[ cnc2$causal , c( "gene", "ensgid" ) ]
coding_genes <- coding_genes0[ !duplicated(coding_genes0) , ]
coding_genes <- coding_genes[ order(coding_genes$gene) , ]

# Establish the set of non-causal genes in the same loci as causal genes
noncausal_near_coding0 <- cnc2[ !cnc2$causal , c( "gene", "ensgid" ) ]
noncausal_near_coding1 <- noncausal_near_coding0[ !duplicated(noncausal_near_coding0) , ]
noncausal_near_coding  <- noncausal_near_coding1[ !( noncausal_near_coding1$ensgid %in% coding_genes$ensgid ) , ]

# Remove coding genes from distal genes
distal_nearest_genes0 <- ng_n_far_fine[ !( ng_n_far_fine$ensgid %in% coding_genes$ensgid ) , c( "gene", "ensgid" ) ]
distal_nearest_genes  <- distal_nearest_genes0[ !duplicated(distal_nearest_genes0) , ]
distal_nearest_genes <- distal_nearest_genes[ order(distal_nearest_genes$gene) , ]
NROW(coding_genes); NROW(distal_nearest_genes)

# Find noncausal genes near distal nearest genes
noncausal_near_ng0 <- list()
for( i in seq_along(distal_nearest_genes$ensgid) ){
  chr   <- gencode$CHR[   gencode$ENSGID == distal_nearest_genes$ensgid[i] ]
  start <- gencode$START[ gencode$ENSGID == distal_nearest_genes$ensgid[i] ]
  end   <- gencode$END[   gencode$ENSGID == distal_nearest_genes$ensgid[i] ]
  local <- gencode[ gencode$CHR == chr & 
                      gencode$END > start - 3e5 & 
                      gencode$START < end + 3e5 &
                      gencode$ENSGID != distal_nearest_genes$ensgid[i], ]
  noncausal_near_ng0[[i]] <- local[ , c( "NAME", "ENSGID" ) ]
}
noncausal_near_ng1 <- do.call( rbind, noncausal_near_ng0 )
noncausal_near_ng <- noncausal_near_ng1[ !duplicated(noncausal_near_ng1) , ]
names(noncausal_near_ng) <- c( "gene", "ensgid" )
noncausal_near_ng <- noncausal_near_ng[ order(noncausal_near_ng$gene) , ]

# Summary
NROW(coding_genes); NROW(noncausal_near_coding)
NROW(distal_nearest_genes); NROW(noncausal_near_ng)

# Write to file
coding_genes_file          <- file.path( maindir, 
                                         "causal_noncausal_trait_gene_pairs",
                                         "coding_genes.tsv" )
noncausal_near_coding_file <- file.path( maindir, 
                                         "causal_noncausal_trait_gene_pairs",
                                         "noncausal_genes_near_coding.tsv" )
distal_genes_file          <- file.path( maindir, 
                                         "causal_noncausal_trait_gene_pairs",
                                         "distal_genes.tsv" )
noncausal_near_distal_file <- file.path( maindir, 
                                         "causal_noncausal_trait_gene_pairs",
                                         "noncausal_genes_near_distal.tsv" )
fwrite( x=coding_genes,           sep="\t", file=coding_genes_file )
fwrite( x=noncausal_near_coding,  sep="\t", file=noncausal_near_coding_file )
fwrite( x=distal_ng_and_tp_genes, sep="\t", file=distal_genes_file )
fwrite( x=noncausal_near_ng,    sep="\t", file=noncausal_near_distal_file )


#-------------------------------------------------------------------------------
#   Done
#-------------------------------------------------------------------------------













