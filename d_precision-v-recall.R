
#-------------------------------------------------------------------------------
#   Functions
#-------------------------------------------------------------------------------

# agresti.95ci.prop: get the Agresti-Coull 95% confidence interval for a proportion
agresti.95ci.prop <- function(X, n, lo.or.hi, Z=qnorm(.025, lower.tail=FALSE) ){
  #X = number of successes
  #n = total number of trials
  #Z = z-score for CI threshold
  p <- (X + 2) / (n + 4)
  num <- p * (1-p)
  denom <- n+4
  diff <- Z * sqrt(num/denom)
  if( lo.or.hi == "lo" ){
    ans <- p - diff
  }else if( lo.or.hi == "hi" ){
    ans <- p + diff
  }else{
    stop("lo.or.hi must be either lo or hi")
  }
  return(ans)
}


#-------------------------------------------------------------------------------
#   Read in the data
#-------------------------------------------------------------------------------

# Arguments
maindir <- "~/projects/causal_genes/"

# Load libraries and sources
library(data.table)
library(psych)

# Read in Weeks et al. 2023 naive precision-recall results
pub_file <- file.path( maindir, "weeks_2023_precision_recall_by_method.csv" )
pub <- fread(pub_file)

# Read in causal/non-causal TGPs annotated with gene mapping evidence
cnc_map_file <- file.path( maindir, "causal_noncausal_trait_gene_pairs", 
                           "causal_tgp_and_gene_mapping_data_300kb.tsv" )
dt <- fread(file=cnc_map_file)


#-------------------------------------------------------------------------------
#   Calculate precision and recall
#-------------------------------------------------------------------------------

# Calculate naive precision and recall
# For each V2G method
# Calculate TP, FP, TN, and FN
# Calculate precision and recall
v2g_colnames <- sort( grep( pattern="_bil$", x=names(dt), value=TRUE ) )
naive0 <- list()
for( i in v2g_colnames ){
  
  # Get set up
  method <- sub( "_bil", "", i )
  focal_cols <- c( "causal", i )
  sub <- dt[ , ..focal_cols ]
  
  # Calculate values
  tp  <- sum(  sub[[i]] &  sub[["causal"]] )
  fp  <- sum(  sub[[i]] & !sub[["causal"]] )
  tn  <- sum( !sub[[i]] & !sub[["causal"]] )
  fn  <- sum( !sub[[i]] &  sub[["causal"]] )
  p   <- tp / ( tp + fp )
  r   <- tp / ( tp + fn )
  
  # Return
  naive0[[i]] <- data.table( method=method, precision=p, recall=r, 
                             tp=tp, fp=fp, tn=tn, fn=fn )
}
naive <- do.call( rbind, naive0 )

# Subset to methods found in the published results
naive$method[ naive$method == "corr_and" ] <- "corr_andersson"
naive$method[ naive$method == "corr_uli" ] <- "corr_ulirsch"
naive$method[ naive$method == "dist_gene" ] <- "distance_genebody"
naive$method[ naive$method == "dist_tss" ] <- "distance_tss"
naive$method[ naive$method == "pchic_jav" ] <- "pchic_javierre"
bmethods <- intersect( pub$method, naive$method )
n2 <- naive[ match( bmethods, naive$method ) , ]
n2


#-------------------------------------------------------------------------------
#   Plot published v. internal results
#-------------------------------------------------------------------------------

# Add error bars to published data
pub$p_lo <- agresti.95ci.prop( X=pub$tp, n=pub$tp+pub$fp, lo.or.hi="lo" )
pub$p_hi <- agresti.95ci.prop( X=pub$tp, n=pub$tp+pub$fp, lo.or.hi="hi" )
pub$r_lo <- agresti.95ci.prop( X=pub$tp, n=pub$tp+pub$fn, lo.or.hi="lo" )
pub$r_hi <- agresti.95ci.prop( X=pub$tp, n=pub$tp+pub$fn, lo.or.hi="hi" )

# Add error bars to internal data
n2$p_lo <- agresti.95ci.prop( X=n2$tp, n=n2$tp+n2$fp, lo.or.hi="lo" )
n2$p_hi <- agresti.95ci.prop( X=n2$tp, n=n2$tp+n2$fp, lo.or.hi="hi" )
n2$r_lo <- agresti.95ci.prop( X=n2$tp, n=n2$tp+n2$fn, lo.or.hi="lo" )
n2$r_hi <- agresti.95ci.prop( X=n2$tp, n=n2$tp+n2$fn, lo.or.hi="hi" )

# Add error bars to naive data
naive$p_lo <- agresti.95ci.prop( X=naive$tp, n=naive$tp+naive$fp, lo.or.hi="lo" )
naive$p_hi <- agresti.95ci.prop( X=naive$tp, n=naive$tp+naive$fp, lo.or.hi="hi" )
naive$r_lo <- agresti.95ci.prop( X=naive$tp, n=naive$tp+naive$fn, lo.or.hi="lo" )
naive$r_hi <- agresti.95ci.prop( X=naive$tp, n=naive$tp+naive$fn, lo.or.hi="hi" )

# Plot: published metrics v. internal metrics
colnames <- names(pub)[2:7]
par( mfrow=c(2,3) )
for( i in colnames ){
  xymax <- max( c( pub[[i]], n2[[i]] ) )
  plot( x=pub[[i]], y=n2[[i]], las=1, cex=2, pch=19, col="#70AD47", 
        xlim=c(0,xymax), ylim=c(0,xymax),
        xlab="Published", ylab="Internal", main=i)
  abline( a=0, b=1, lty=2, lwd=2, col="grey" )
  for( j in seq_along(pub$method) ){
    if( i == colnames[1] ){
      lines( x=c( pub$p_lo[j], pub$p_hi[j] ),
             y=c( n2$precision[j], n2$precision[j] ),
             lwd=2, col="#70AD47" )
      lines( x=c( pub$precision[j], pub$precision[j] ),
             y=c( n2$p_lo[j], n2$p_hi[j] ),
             lwd=2, col="#70AD47" )
    }
    if( i == colnames[2] ){
      lines( x=c( pub$r_lo[j], pub$r_hi[j] ),
             y=c( n2$recall[j], n2$recall[j] ),
             lwd=2, col="#70AD47" )
      lines( x=c( pub$recall[j], pub$recall[j] ),
             y=c( n2$r_lo[j], n2$r_hi[j] ),
             lwd=2, col="#70AD47" )
    }
  }
}
par( mfrow=c(1,1) )

# Plot: precision v. recall for both datasets
pr_file <- file.path( maindir, "plots", "precision_recall_plots.jpg" )
jpeg( filename=pr_file, width=10, height=5, units="in", res=300 )
par( mfrow=c(1,2) )
par( mar=c(4,4,2,1) )
xymax <- max( c( pub$recall, pub$precision, n2$recall, n2$precision ) )
plot( x=pub$recall, y=pub$precision, las=1, cex=2, pch=19, col=pub$col,
      xlim=c( 0, xymax ), ylim=c( 0, xymax ),
      xlab="Recall", ylab="Precision", main="Published")
legend( "bottomright", legend=pub$method, fill=pub$col, cex=0.7 )
plot( x=n2$recall, y=n2$precision, las=1, cex=2, pch=19, col=pub$col,
      xlim=c( 0, xymax ), ylim=c( 0, xymax ),
      xlab="Recall", ylab="Precision", main="CALDERA")
legend( "bottomright", legend=pub$method, fill=pub$col, cex=0.7 )
par( mfrow=c(1,1) )
dev.off()

# Plot: precision v. recall: naive
library(colorspace)
plot( x=naive$recall, y=naive$precision, las=1, cex=2, pch=19,
      xlim=c( 0, 1 ), ylim=c( 0, 1 ),
      col=rainbow_hcl( n=NROW(naive) ),
      xlab="Recall", ylab="Precision", main="Internal")
legend( "topright", legend=naive$method, fill=rainbow_hcl( n=NROW(naive) ), cex=0.8 )


#-------------------------------------------------------------------------------
#   Plot precision v. recall for each POPS quantile method
#-------------------------------------------------------------------------------

# Set up the data
n3 <- naive[ grepl( "^pops", naive$method ) , ]
n3$col <- viridis( NROW(n3) )
n3$f1 <- ( 2 * n3$precision * n3$recall ) / ( n3$precision + n3$recall )

# F1 contour line
f1 <- n3$f1[1]
rs <- seq( 0.25, 0.5, 0.01 )
ps <- (f1*rs) / (2*rs - f1)

# Plot
par( mar=c(4.5,4.5,1,1) ) 
plot( x=n3$recall, y=n3$precision, las=1, cex=2, pch=19, col=n3$col,
      xlim=c( 0, max(n3$recall) ), ylim=c( 0, max(n3$precision) ),
      xlab="Recall", ylab="Precision" )
lines( x=rs, y=ps, lty=2 )
legend( "bottomright", legend=n3$method, fill=n3$col )


#-------------------------------------------------------------------------------
#   Done
#-------------------------------------------------------------------------------








