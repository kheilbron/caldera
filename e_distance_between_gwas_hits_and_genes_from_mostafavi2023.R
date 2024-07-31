

#-------------------------------------------------------------------------------
#   Read in data
#-------------------------------------------------------------------------------

# Read in gencode gene boundaries
gencode <- fread("~/projects/causal_genes/gene_locations.tsv")

# Read in the data, merge
gw_file <- "~/projects/causal_genes/mostafavi2023_gene_annots/gwas_hits.csv"
gw <- fread(gw_file)
pattern <- "^([[:digit:]]+):([[:digit:]]+):.*$"
gw$chr <- as.integer( sub( pattern=pattern, replacement="\\1", x=gw$Variant ) )
gw$bp  <- as.integer( sub( pattern=pattern, replacement="\\2", x=gw$Variant ) )

# Read in the data, merge
eq_file <- "~/projects/causal_genes/mostafavi2023_gene_annots/eqtl_hits.csv"
eq0 <- fread(eq_file)
eq <- eq0[ eq0$eGene_Symbol %in% gencode$ENSGID , ]
eq$chr <- as.integer( sub( pattern=pattern, replacement="\\1", x=eq$Variant ) )
eq$bp  <- as.integer( sub( pattern=pattern, replacement="\\2", x=eq$Variant ) )
NROW(eq0); NROW(eq)


#-------------------------------------------------------------------------------
#   Extract nearest gene/TSS for GWAS
#-------------------------------------------------------------------------------

# For each GWAS hit
# Find the distance to the nearest gene (and TSS)
gw_dist0 <- list()
for( i in seq_along(gw$chr) ){
  
  # Subset to the right chromosome
  sub <- gencode[ gencode$CHR == gw$chr[i] , ]
  
  # Find distances between focal variant and all genes
  diff_start <- abs( sub$START - gw$bp[i] )
  diff_end   <- abs( sub$END   - gw$bp[i] )
  diff_tss   <- abs( sub$TSS   - gw$bp[i] )
  
  # Find minimum distance
  min_tss    <- min(diff_tss)
  min_body   <- min( c( diff_end, diff_start ) )
  if( any( gw$bp[i] > sub$START & gw$bp[i] < sub$END ) ) min_body <- 0
  
  # Return
  gw_dist0[[i]] <- data.table( body=min_body, tss=min_tss )
}
gw_dist <- do.call( rbind, gw_dist0 )
head(gw_dist)


#-------------------------------------------------------------------------------
#   Extract nearest/actual gene/TSS for eQTLs
#-------------------------------------------------------------------------------

# For each eQTL hit
# Find the distance to the nearest gene (and TSS), as well as the actual gene
eq_dist0 <- list()
for( i in seq_along(eq$chr) ){
  
  # Subset to the right chromosome
  sub  <- gencode[ gencode$CHR == eq$chr[i] , ]
  sub2 <- gencode[ gencode$ENSGID == eq$eGene_Symbol[i] , ]
  
  # Find distances between focal variant and all genes
  diff_start <- abs( sub$START - eq$bp[i] )
  diff_end   <- abs( sub$END   - eq$bp[i] )
  diff_tss   <- abs( sub$TSS   - eq$bp[i] )
  
  # Find distance to nearest gene
  min_tss    <- min(diff_tss)
  min_body   <- min( c( diff_end, diff_start ) )
  if( any( eq$bp[i] > sub$START & eq$bp[i] < sub$END ) ) min_body <- 0
  
  # Find distance to causal gene
  real_tss    <- abs( sub2$TSS   - eq$bp[i] )
  real_body   <- min( c( abs( sub2$START - eq$bp[i] ),
                         abs( sub2$END   - eq$bp[i] ) ) )
  if( eq$bp[i] > sub2$START & eq$bp[i] < sub2$END ) real_body <- 0
  
  # Return
  eq_dist0[[i]] <- data.table( nearest_body=min_body, nearest_tss=min_tss,
                               actual_body=real_body, actual_tss=real_tss )
}
eq_dist <- do.call( rbind, eq_dist0 )
head(eq_dist)


#-------------------------------------------------------------------------------
#   Plot number of genes in various distance bins
#-------------------------------------------------------------------------------

# Get the percentage in each bin
dists <- seq( 0, 5e5, 1e5 )
perc0 <- list()
for( i in seq_len( length(dists) - 1 ) ){
  p_eq <- mean( eq_dist$actual_body >= dists[i] & eq_dist$actual_body < dists[i+1] )
  p_gw <- mean( gw_dist$body        >= dists[i] & gw_dist$body        < dists[i+1] )
  perc0[[i]] <- c( p_eq, p_gw )
}

# Format
df <- as.data.frame( do.call( rbind, perc0 ) )
names(df) <- c( "eQTL", "GWAS" )
row.names(df) <- c( "0-100kb", "100-200kb", "200-300kb", 
                    "300-400kb", "400-500kb" )

# Plot
maindir     <- "~/projects/causal_genes/"
fig_dir     <- file.path( maindir, "figures" )
fig_s4_file <- file.path( fig_dir, "figure_s4.jpg" )
jpeg( filename=fig_s4_file, width=600*4, height=300*4, res=75*4 )
par( mfrow=c(1,2) )
par( mar=c(6,5,1,1) )
bp1 <- barplot( as.matrix(df), beside=TRUE, las=1, ylab="Proportion of genes",
         col=brewer.pal( n=6, name="Greens" )[5:1] )
legend( x="topright", legend=row.names(df), cex=0.7,
        fill=brewer.pal( n=6, name="Greens" )[5:1] )
barplot( t(df), beside=TRUE, las=2, ylab="Proportion of genes",
         col=brewer.pal( n=6, name="Greens" )[5:4] )
legend( x="topright", legend=names(df), 
        fill=brewer.pal( n=6, name="Greens" )[5:4] )
par( mfrow=c(1,1) )
dev.off()


#-------------------------------------------------------------------------------
#   Done
#-------------------------------------------------------------------------------
























