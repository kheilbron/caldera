

#-------------------------------------------------------------------------------
#   Format credible sets
#-------------------------------------------------------------------------------

# Load libraries
suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(dplyr) )

# Read in gencode gene boundaries
gencode <- fread("~/projects/causal_genes/gene_locations.tsv")

# Independent traits
it0 <- fread("~/projects/causal_genes/weeks2023_table_s1.csv")
it <- it0$Trait[ it0$Cohort == "UKB" & it0$`Independent traits (for meta-analysis)` ]
length(it)

# Read in CS column names
# Read in CSs, add column names
maindir <- file.path( Sys.getenv("HOME"), "projects/causal_genes/" )
cs_colnames <- fread( file.path( maindir, "release1.1", "UKBB_94traits_release1.cols" ), 
                      header=FALSE )
cs <- fread( file.path( maindir, "release1.1", "UKBB_94traits_release1.bed.gz" ) )
names(cs) <- cs_colnames$V1

# z_to_p: Convert a z-score into a P value
z_to_p <- function( z, log.p=FALSE ){
  if(log.p){
    log(2) + pnorm( -abs(z), log.p=TRUE )
  }else{
    2*pnorm( -abs(z) )
  }
}  

# Remove SNPs that are not in CSs, or are in CSs with the 6th+ worst ABF
# Subset to SUSIE CSs only
# Add P value and TCP columns
cs2 <- cs[  cs$cs_id > 0 & cs$cs_id < 6 & cs$trait %in% it , ]
cs3 <- cs2[ cs2$method == "SUSIE" , ]
cs3$p <- z_to_p( z=sqrt(cs3$chisq_marginal) )
cs3$tcp <- paste( cs3$trait, cs3$region, cs3$cs_id, sep="_" )
head( cs3, 1 )

# Add pip_max and pip_ratio
pip_max <- tapply( X=cs3$pip, INDEX=cs3$tcp, FUN=max )
cs3$pip_max <- pip_max[ match( cs3$tcp, names(pip_max) ) ]
cs3$pip_ratio <- cs3$pip / cs3$pip_max
summary(cs3$pip_ratio)

# Subset to CSs where PIP > 25% and pip_ratio = 1
cs4 <- cs3[ cs3$pip > 0.25 & cs3$pip_ratio == 1 , ]
NROW(cs3); NROW(cs4)


#-------------------------------------------------------------------------------
#   Find coding TGPs
#-------------------------------------------------------------------------------

# Read in V2G
vg_file <- file.path( maindir, "v2g_file.tsv" )
vg0 <- fread(vg_file)
vg <- vg0[ vg0$trait %in% it , ]
vg$tcp <- paste( vg$trait, vg$region, vg$cs_id, sep="_" )
NROW(vg0); NROW(vg)
head( vg, 2 )

# Extract coding TGPs and TCPs (coding PIP > 25% + maxPIP)
# idx1 <- !is.na(vg$coding_prob)
# idx2 <- vg$coding_prob > 0.5
# idx3 <- vg$coding_prob > 0.25 & vg$tcp %in% cs4$tcp
# idx <- idx1 & ( idx2 | idx3 )
# tgp_c0 <- vg[ idx , ]
# NROW(tgp_c0)
# tcp_c_all <- unique(tgp_c0$tcp)

# Extract coding TGPs and TCPs (coding PIP > 50%)
idx <- !is.na(vg$coding_prob) & vg$coding_prob > 0.5
tgp_c0 <- vg[ idx , ]
tcp_c_all <- unique(tgp_c0$tcp)

# Subset to coding TGPs that 1) use canonical ENSGIDs, 2) have cs_id <= 5
# (5th strongest CS in the region or better), and removing duplicated TGPs
# (preferentially retaining low cs_id)
tgp_c0 <- tgp_c0[ order( tgp_c0$trait, tgp_c0$region, tgp_c0$cs_id ) , ]
tgp_c  <- tgp_c0[ tgp_c0$ensgid %in% gencode$ENSGID & 
                  tgp_c0$cs_id <= 5 & 
                  !duplicated( tgp_c0[ , c( "trait", "ensgid" ) ] ) , ]
tcp_c <- unique(tgp_c$tcp)
NROW(tgp_c0); NROW(tgp_c); length(tcp_c)
head( tgp_c, 2 )


#-------------------------------------------------------------------------------
#   Find non-coding TCPs near coding TGPs
#-------------------------------------------------------------------------------

# Extract non-coding TCPs
tcp_n <- setdiff( cs3$tcp, tcp_c_all )
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

# Remove non-coding TCPs >300kb from causal genes
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
  #   causal start < CS end + window, causal end > CS start - window
  idx <- tgp_c$trait == trait &
    tgp_c$chromosome == chr &
    tgp_c$start <= cs_end   + window &
    tgp_c$end   >= cs_start - window &
    cs_end - cs_start < 4e5
  
  # If a non-coding TCP is near a single coding TGP, store TCP information
  if( sum(idx) == 1 ){
    dt <- data.table( tcp=tcp_n2[i], trait=trait, region=region, cs_id=cs_id,
                      chr=chr, cs_start=cs_start, cs_end=cs_end,
                      n_cs_snps=NROW(sub), cs_width=cs_end-cs_start,
                      n_cs_in_region=n_cs, tcp_c=tgp_c$tcp[idx],
                      ensgid_c=tgp_c$ensgid[idx] )
    tcp_n_near0[[i]] <- dt
  }
}
tcp_n_near0 <- do.call( rbind, tcp_n_near0 )

# If a coding TCP is associated with multiple non-coding TCPs, only keep the best cs_id (ABF)
# tcp_n_near0 <- tcp_n_near0[ order( tcp_n_near0$trait, 
#                                    tcp_n_near0$region, 
#                                    tcp_n_near0$cs_id ) , ]
# tcp_n_near  <- tcp_n_near0[ !duplicated(tcp_n_near0$tcp_c) , ]
# NROW(tcp_n_near0); NROW(tcp_n_near)
tcp_n_near <- tcp_n_near0

# Extract full CS data for non-coding TCPs near coding TGPs
cs_n_near  <- cs3[ cs3$tcp %in% tcp_n_near$tcp , ]
NROW(cs_n_near)


#-------------------------------------------------------------------------------
#   Extract all genes near non-coding TCPs (near causal TGPs)
#-------------------------------------------------------------------------------

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
               "n_cs_snps", "cs_width" )
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

# Remove loci with < 2 genes in them
cnc2 <- cnc[ cnc$ngenes_nearby >= 2 , ]

# Remove traits with <5 causal genes
n_causal_per_trait <- sort( table( cnc2$trait[ cnc2$causal ] ), decreasing=TRUE )
ge_5_causal_per_trait <- names(n_causal_per_trait)[ n_causal_per_trait >= 5 ]
cnc3 <- cnc2[ cnc2$trait %in% ge_5_causal_per_trait , ]
table(cnc$causal); table(cnc2$causal); table(cnc3$causal)

# Write causal/non-causal trait-gene pairs to file
causal_gene_file <- file.path( maindir, "causal_noncausal_trait_gene_pairs",
                               "causal_noncausal_trait_gene_pairs_300kb.tsv" )
fwrite( x=cnc3, file=causal_gene_file, sep="\t" )


#-------------------------------------------------------------------------------
#   Done
#-------------------------------------------------------------------------------




























