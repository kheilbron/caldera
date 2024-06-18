

#-------------------------------------------------------------------------------
#   simple_functions
#-------------------------------------------------------------------------------

# logit10: Convert a probability into a log10 odds ratio
logit10 <- function(p) log10( p / (1-p) )

# prop_to_odds: Converts a proportion to an odds
prop_to_odds <- function(proportion){
  proportion / (1 - proportion)
}


#-------------------------------------------------------------------------------
#   largest_rg_trait: Find the largest genetic correlation for each phenotype
#-------------------------------------------------------------------------------

largest_rg_trait <- function(trait){
  
  # Input arguments
  maindir <- "~/projects/causal_genes/"
  
  # Load libraries
  suppressPackageStartupMessages( library(data.table) )
  suppressPackageStartupMessages( library(dplyr) )
  
  # Read in rg
  rg <- fread( file.path( maindir, "LDSC/UKBB.gencor.txt" ) )
  rg$p1 <- sub( pattern="UKBB.", replacement="", x=rg$p1 )
  rg$p2 <- sub( pattern="UKBB.", replacement="", x=rg$p2 )
  
  # Read in CALDERA LOTO model
  bgl_las       <- readRDS( file=file.path( maindir, "bgl_las.rds" ) )
  tr_traits     <- bgl_las$trait
  non_tr_traits <- setdiff( sort( unique( c( rg$p1, rg$p2 ) ) ), tr_traits )
  if( trait %in% tr_traits )  return(trait)
  
  # Subset to rows with only one independent and one non-independent trait
  rg$p1_in_tr <- rg$p1 %in% tr_traits
  rg$p2_in_tr <- rg$p2 %in% tr_traits
  rg$n_tr     <- rg$p1_in_tr  + rg$p2_in_tr
  idx         <- which( rg$n_tr == 1 )
  rg2         <- rg[ idx, ]
  
  # Sort by nominally-significant P and abs(rg)
  rg2 <- rg2[ order( rg2$p > 0.05, -abs(rg2$rg) ) , ]
  
  # Select top rg
  idx <- min( which( rg2$p1 == trait | rg2$p2 == trait ) )
  top_indep_trait <- ifelse( rg2$p1_in_tr[idx], rg2$p1[idx], rg2$p2[idx] )
  return(top_indep_trait)
}


#-------------------------------------------------------------------------------
#   caldera_internal: Run CALDERA on all 94 UKB traits from original PoPS pub.
#-------------------------------------------------------------------------------

caldera_internal <- function(){
  
  #-----------------------------------------------------------------------------
  #   Read in data
  #-----------------------------------------------------------------------------
  
  # Input arguments
  maindir <- "~/projects/causal_genes/"
  
  # Load libraries
  suppressPackageStartupMessages( library(data.table) )
  suppressPackageStartupMessages( library(dplyr) )
  suppressPackageStartupMessages( library(glmnet) )
  
  # Read in gencode gene boundaries
  gencode  <- fread("~/projects/causal_genes/gene_locations.tsv")
  
  # Read in P(causal) prediction and recalibration models
  bgl_las   <- readRDS( file=file.path( maindir, "bgl_las.rds" ) )
  
  # Read in V2G
  vg_file <- file.path( maindir, "v2g_file.tsv" )
  vg0 <- fread(vg_file)
  vg0$tcp <- paste0( vg0$trait, "_", vg0$region, "_", vg0$cs_id )
  vg0 <- vg0[ order( vg0$trait, vg0$region, vg0$cs_id ) , ]
  
  # Read in the data, merge
  mo_file <- file.path( maindir, "mostafavi2023_gene_annots/pc_genes.txt" )
  mo <- fread(mo_file)
  names(mo)[ names(mo) == "GeneSymbol" ] <- "ensgid"
  mo$gene <- NULL
  vg <- left_join( x=vg0, y=mo, by="ensgid" )
  
  # Read in credible sets
  cs_colnames <- fread( file.path( maindir, "release1.1", "UKBB_94traits_release1.cols" ), 
                        header=FALSE )
  cs <- fread( file.path( maindir, "release1.1", "UKBB_94traits_release1.bed.gz" ) )
  names(cs) <- cs_colnames$V1
  cs2 <- cs[  cs$cs_id > 0 , ]
  cs3 <- cs2[ cs2$method == "SUSIE" , ]
  cs3$tcp <- paste( cs3$trait, cs3$region, cs3$cs_id, sep="_" )
  
  
  #-----------------------------------------------------------------------------
  #   Subset to genes within locus boundaries
  #-----------------------------------------------------------------------------
  
  lbounds <- cbind( tapply( X=cs3$start, INDEX=cs3$tcp, FUN=min ),
                    tapply( X=cs3$end,   INDEX=cs3$tcp, FUN=max ) )
  lbounds <- as.data.frame(lbounds)
  names(lbounds) <- c( "locus_start", "locus_end" )
  
  # Subset V2G to genes that are within the locus boundaries
  idx <- vg$start < lbounds[ vg$tcp , "locus_end" ] & 
         vg$end   > lbounds[ vg$tcp , "locus_start" ]
  vg <- vg[ idx & !is.na(idx) , ]
  
  
  #-----------------------------------------------------------------------------
  #   Construct global features
  #-----------------------------------------------------------------------------
  
  # n_genes
  n_genes_per_tcp0 <- list()
  count <- 0
  for( i in seq_along( unique(vg$tcp) ) ){
    
    # Message
    count <- count + 1
    if( count %% 1000 == 0 ) 
      message( date(), "   TCP ", count, "/", length( unique(vg$tcp) ) )
    
    # Do
    tcp     <- unique(vg$tcp)[i]
    lstart  <- lbounds[ tcp, "locus_start" ]
    lend    <- lbounds[ tcp, "locus_end" ]
    lead    <- vg$lead_variant[ vg$tcp == unique(vg$tcp)[i] ][1]
    pattern <- "chr([[:digit:]]+):([[:digit:]]+):.*"
    chr     <- as.integer( sub( pattern=pattern, replacement="\\1", x=lead ) )
    n_genes_per_tcp0[[i]] <- sum( gencode$CHR   == chr & 
                                  gencode$START  < lend   + 3e5 & 
                                  gencode$END    > lstart - 3e5 )
  }
  names(n_genes_per_tcp0) <- unique(vg$tcp)
  n_genes_per_tcp <- unlist(n_genes_per_tcp0)
  vg$ngenes_nearby <- n_genes_per_tcp[ match( vg$tcp, names(n_genes_per_tcp) ) ]
  vg$prior_n_genes_locus <- logit10( 1 / vg$ngenes_nearby )
  vg$prior_n_genes_locus[ vg$prior_n_genes_locus == Inf ] <- logit10(0.75)
  vg <- vg[ !is.na(vg$prior_n_genes_locus) , ]
  
  # POPS
  vg$pops_glo <- vg$pops_score
  
  # Distance to gene body 
  vg$dist_gene_glo <- -log10( vg$distance_genebody + 1e3 )
  
  # MAGMA
  vg$magma_glo <- ifelse( is.na(vg$magma_score), 
                          median( vg$magma_score, na.rm=TRUE ), 
                          vg$magma_score )
  
  # Coding
  vg$coding_glo <- ifelse( is.na(vg$coding_prob), 
                           logit10( min( vg$coding_prob, na.rm=TRUE ) ), 
                           logit10(vg$coding_prob) )
  vg$coding_glo[ vg$coding_glo == Inf ] <- logit10(0.99)
  
  
  #-----------------------------------------------------------------------------
  #   Construct relative and best-in-locus features
  #-----------------------------------------------------------------------------
  
  # Set up global feature columns and TCPs
  glo_cols <- grep( "_glo$", x=names(vg), value=TRUE )
  
  # Loop through global columns
  for( j in glo_cols ){
    
    # Message
    message( date(), "   Making relative/BIL feature ", which( glo_cols == j ), 
             "/", length(glo_cols), ": ", j )
    
    # Initialize relative and best-in-locus columns
    rel_col <- sub( pattern="_glo$", replacement="_rel", x=j )
    bil_col <- sub( pattern="_glo$", replacement="_bil", x=j )
    vg[[rel_col]] <- as.numeric(NA)
    vg[[bil_col]] <- as.logical(NA)
    
    # Loop through TCPs
    for( i in unique(vg$tcp) ){
      
      # Subset
      sub <- vg[ vg$tcp == i , ]
      
      # Compute relative and best-in-locus values
      r_vals <- sub[[j]] - max( sub[[j]] )
      b_vals <- r_vals == 0 & sum( r_vals == 0 ) == 1
      
      # Insert relative values into the main data table
      set( x     = vg, 
           i     = which( vg$tcp == i ), 
           j     = rel_col, 
           value = r_vals )
      
      # Insert best-in-locus values into the main data table
      set( x     = vg, 
           i     = which( vg$tcp == i ), 
           j     = bil_col, 
           value = b_vals ) 
    }
  }
  
  
  #-----------------------------------------------------------------------------
  #   Construct GLC features
  #-----------------------------------------------------------------------------
  
  # Missingness
  vg$pritchard_miss <- FALSE
  
  # pLI
  vg$pLI[ is.na(vg$pLI) ] <- 0.5
  vg$pLI_lt_0.9 <- ifelse( vg$pLI < 0.9, 1, 0 )
  vg$pLI_lt_0.1 <- ifelse( vg$pLI < 0.1, 1, 0 )
  vg$pLI_log10 <- -log10(vg$pLI)
  pLI99 <- 34.3
  vg$pLI_log10 <- ifelse( vg$pLI_log10 > pLI99, pLI99, vg$pLI_log10 )
  
  # LOEUF
  vg$LOEUF[ is.na(vg$LOEUF) ] <- 0
  
  # hs
  vg$hs[ is.na(vg$hs) ] <- max( vg$hs, na.rm=TRUE )
  vg$hs_log10 <- -log10(vg$hs)
  
  # Gene length
  vg$length[ is.na(vg$length) ] <- min( vg$length, na.rm=TRUE )
  vg$gene_bp_log10 <- log10( vg$length*1e3 )
  
  # CDS length
  vg$CDS_length[ is.na(vg$CDS_length) ] <- min( vg$CDS_length, na.rm=TRUE )
  vg$cds_bp_log10 <- log10( vg$CDS_length*1e3 )
  
  # ABC_length_per_type
  vg$ABC_length_per_type[ is.na(vg$ABC_length_per_type) ] <- 
    min( vg$ABC_length_per_type, na.rm=TRUE )
  vg$abc_bp_log10 <- ifelse( vg$ABC_length_per_type == 0, 
                             log10( 0.2 * 1e3 ),
                             log10( vg$ABC_length_per_type * 1e3 ) )
  
  # Roadmap_length_per_type
  vg$Roadmap_length_per_type[ is.na(vg$Roadmap_length_per_type) ] <- 
    min( vg$Roadmap_length_per_type, na.rm=TRUE )
  vg$roadmap_bp_log10 <- ifelse( vg$Roadmap_length_per_type == 0, 
                                 log10( min( vg$Roadmap_length_per_type[ 
                                   vg$Roadmap_length_per_type > 0 ] ) * 1e3 ),
                                 log10( vg$Roadmap_length_per_type * 1e3 ) )
  
  
  #-----------------------------------------------------------------------------
  #   For each trait, get predictions and recalibrate
  #-----------------------------------------------------------------------------
  
  ukb_preds0 <- list()
  traits <- unique(vg$trait)
  for( trait in traits ){
    
    # Message
    message( date(), "   Starting trait ", which( traits == trait), "/",
             length(traits), ", ", trait )
    
    # Select the right LOTO and recalibration models
    mod_trait <- largest_rg_trait(trait)
    main_mod  <- bgl_las$model[[mod_trait]]
    rc_mod    <- bgl_las$recal_model[[mod_trait]]
    
    # Subset to necessary columns
    bias_cols <- c( "pLI_lt_0.9", "pLI_lt_0.1", "pLI_log10", "hs_log10", 
                    "gene_bp_log10", "cds_bp_log10", "abc_bp_log10", 
                    "roadmap_bp_log10", "LOEUF", "pritchard_miss" )
    gcols <- c( "trait", "gene", "ensgid", "tcp",
                grep( pattern="_glo$", x=names(vg), value=TRUE ),
                grep( pattern="_rel$", x=names(vg), value=TRUE ),
                grep( pattern="_bil$", x=names(vg), value=TRUE ),
                "prior_n_genes_locus",
                bias_cols )
    row_idx <- vg$trait == trait
    vg2 <- vg[ row_idx , ..gcols ]
    
    # Fix potentially-biased covariates to their mean value
    for( i in bias_cols ){
      vg2[[i]] <- mean( vg2[[i]] )
    }
    
    # Get predictions
    feat_cols <- main_mod$glmnet.fit$beta@Dimnames[[1]]
    vg2$causal_p <- as.vector( 
      predict( object=main_mod, s="lambda.min", type="response",
               newx=as.matrix( vg2[ , ..feat_cols ] ) ) )
    
    # Convert to global and relative scores
    vg2$global   <- log( prop_to_odds(vg2$causal_p) )
    vg2$relative <- as.numeric(NA)
    
    # Loop through TCPs
    for( i in unique(vg2$tcp) ){
      
      # Subset
      sub <- vg2[ vg2$tcp == i , ]
      
      # Compute relative and best-in-locus values
      r_vals <- sub[["global"]] - max( sub[["global"]] )
      
      # Insert relative values into the main data table
      set( x     = vg2, 
           i     = which( vg2$tcp == i ), 
           j     = "relative", 
           value = r_vals )
    }
    
    # Recalibrate
    vg2$causal_r <- as.vector( 
      predict( object=rc_mod, s="lambda.min", type="response",
               newx=as.matrix( vg2[ , c( "global", "relative" ) ] ) ) )
    
    # Return
    out_cols <- c( "trait", "gene", "causal_p", "causal_r", "ensgid", "tcp" )
    ukb_preds0[[trait]] <- vg2[ , ..out_cols ]
  }
  ukb_preds <- do.call( rbind, ukb_preds0 )
  
  
  #-----------------------------------------------------------------------------
  #   Return original and recalibrated predictions
  #-----------------------------------------------------------------------------
  outfile <- file.path( maindir, "ukb_preds.tsv" )
  fwrite( x=ukb_preds, file=outfile, sep="\t" )
  return(ukb_preds)
}
sum( ukb_preds$causal_r > 0.589 )
sum( ukb_preds$causal_r > 0.75 )
a <- ukb_preds[ ukb_preds$causal_r > 0.589 , ]
a <- ukb_preds[ ukb_preds$causal_r > 0.75 , ]
table( duplicated( a[ , c( "trait", "ensgid" ) ] ) )
























