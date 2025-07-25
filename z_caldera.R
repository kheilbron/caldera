

#-------------------------------------------------------------------------------
#   Simple functions
#-------------------------------------------------------------------------------

# logit(10): Convert a probability into a log or log10 odds ratio
logit   <- function(p) log( p / (1-p) )
logit10 <- function(p) log10( p / (1-p) )

# message2: Just like message, but with the date and a gap pasted in
message2 <- function(...) message(date(), "     ", ...)


#-------------------------------------------------------------------------------
#   CALDERA
#-------------------------------------------------------------------------------

caldera <- function( pops_file, cs_file, caldera_path ){
  
  #-----------------------------------------------------------------------------
  #   Inputs
  #-----------------------------------------------------------------------------
  
  # pops_file:    PoPS output file. Must contain columns: ENSGID and PoPS_Score.
  # cs_file:      Credible set (CS) file containing all CS variants across all
  #               loci. Used to define loci, the genes therein, and coding PIPs 
  #               for these genes. Must contain columns: locus, chr (chromosome),
  #               bp (GRCh37 position), and pip (posterior inclusion probability).
  #               Optionally, you can include a c_gene column containing the
  #               ENSGID of the gene whose coding sequence is affected by the 
  #               given variant (and NA if the variant does not affect any coding
  #               sequences). Alternatively, an "rsid" column containing rsIDs can be
  #               provided and variants will automatically be annotated using
  #               a list of ~76k UK Biobank coding variants. If "rsid" is not
  #               provided, variants will be automatically annotated based on 
  #               chr and (GRCh37) bp.
  # caldera_path: Path to local CALDERA GitHub repository.
  
  
  #-----------------------------------------------------------------------------
  #   Read in data
  #-----------------------------------------------------------------------------
  
  # Libraries
  library(data.table)
  library(dplyr)
  
  # Read in CS file
  cs <- fread(cs_file)
  
  # Read in POPS file
  pop <- fread(pops_file)
  
  # Read in coding variants
  data_dir    <- file.path( caldera_path, "data" )
  coding_file <- file.path( data_dir, "high_h2_coding_SNPs.tsv" )
  cod <- fread(coding_file)
  
  # Read in GENCODE genes
  gene_file <- file.path( data_dir, "gene_locations.tsv" )
  gen <- fread(gene_file)
  
  # Read in CALDERA model
  mod_dir  <- file.path( caldera_path, "trained_models" )
  cal_file <- file.path( mod_dir, "caldera_model.rds" )
  cal_mod  <- readRDS(cal_file)
  
  # Read in covariate means
  cov_file <- file.path( mod_dir, "cov_means.rds" )
  cov_means  <- readRDS(cov_file)
  
  
  #-----------------------------------------------------------------------------
  #   Check inputs
  #-----------------------------------------------------------------------------
  
  # pops_file
  pop_cols <- setdiff( c( "ENSGID", "PoPS_Score" ), names(pop) )
  if( length(pop_cols) > 0 ){
    stop( "PoPS file must contain columns: ENSGID and PoPS_Score" )
  }
  
  # cs_file: required columns
  cs_cols <- setdiff( c( "locus", "chr", "bp", "pip" ), names(cs) )
  if( length(cs_cols) > 0 ){
    stop( "CS file must contain columns: locus, chr, bp, pip" )
  }
  
  # cs_file: report what kind of coding variant annotation will be used
  if( "c_gene" %in% names(cs) ){
    message2( "A c_gene column has been provided and will be used to",
              " determine genes affected by coding variants")
  }else if( "rsid" %in% names(cs) ){
    message2( "An rsID column has been provided and will be used to determine",
              " genes affected by a list of ~76k UK Biobank coding variants")
  }else{
    message2( "Chromosome and (assumed GRCh37) position will be used to determine",
              " genes affected by a list of ~76k UK Biobank coding variants")
  }
  
  
  #-----------------------------------------------------------------------------
  #   Define loci, extract genes in loci
  #-----------------------------------------------------------------------------
  
  # Define locus boundaries
  loci0 <- list()
  window <- 3e5
  for( i in unique(cs$locus) ){
    sub    <- cs[ cs$locus == i , ]
    chr    <- sub$chr[1]
    bp     <- mean( sub$bp[ sub$pip == max(sub$pip) ] )
    bp_min <- min(sub$bp) - window
    if( bp_min < 1 )  bp_min <- 1
    bp_max <- max(sub$bp) + window
    id     <- paste0( chr, "_", bp_min, "_", bp_max )
    loci0[[i]] <- data.table( locus=i, locus_pos=id, chr=chr, bp=bp,
                             bp_min=bp_min, bp_max=bp_max )
  }
  loci <- do.call( rbind, loci0 )
  
  # Find all genes in each GWAS locus and their distance to the GWAS hit
  genes0 <- list()
  for( i in seq_along(loci$locus) ){
    sub <- gen[ gen$CHR  == loci$chr[i] & 
                gen$END   > loci$bp_min[i] & 
                gen$START < loci$bp_max[i] , ]
    sub$dist_s <- abs( sub$START - loci$bp[i] )
    sub$dist_e <- abs( sub$END   - loci$bp[i] )
    sub$dist   <- ifelse( sub$dist_s < sub$dist_e, sub$dist_s, sub$dist_e)
    sub$dist[ loci$bp[i] > sub$START & loci$bp[i] < sub$END ] <- 0
    if( NROW(sub) > 0 ){
      genes0[[i]] <- data.table( locus=loci$locus[i], locus_pos=loci$locus_pos[ loci$locus == i ], 
                                 gene=sub$NAME, ensgid=sub$ENSGID, 
                                 n_genes=NROW(sub), dist=sub$dist )
    }
  }
  genes <- do.call( rbind, genes0 )
  
  
  #-----------------------------------------------------------------------------
  #   Decorate genes with transformed distance, POPS, and coding PIP
  #-----------------------------------------------------------------------------
  
  # Transform distance
  genes$dist_gene_glo <- -log10( genes$dist + 1e3 )
  
  # Decorate these genes with their POPS scores
  genes$pops_glo <- pop$PoPS_Score[ match( genes$ensgid, pop$ENSGID ) ]
  genes$pops_glo[ is.na(genes$pops_glo) ] <- 0
  
  # If not already done, annotate coding CS variants
  if( !( "c_gene" %in% names(cs) ) ){
    
    if( "rsid" %in% names(cs) ){                            # Using rsID
      cs$c_gene <- cod$ensgid[ match( cs$rsid, cod$SNP ) ]
    }else{                                                 # Using chr and pos
      bcols <- c( "CHR", "BP", "ensgid" )
      ccols <- c( "chr", "bp", "c_gene" )
      cod2 <- cod[ , ..bcols ]
      names(cod2) <- ccols
      cs <- left_join( x=cs, y=cod2, by=c( "chr", "bp" ) )
    }
  }
  
  # Aggregate PIPs across loci and genes
  # Decorate with coding PIP
  cs2 <- cs[ !is.na(cs$c_gene) & cs$c_gene != "" , ]
  if( NROW(cs2) > 0 ){
    cs3 <- aggregate( x=cs2$pip, by=list( locus=cs2$locus, ensgid=cs2$c_gene ), FUN=sum )
    names(cs3)[3] <- "coding_glo"
    genes <- left_join( x=genes, y=cs3, by=c( "locus", "ensgid" ) )
    genes$coding_glo[ is.na(genes$coding_glo) ] <- 0
  }else{
    genes$coding_glo <- 0
  }
  
  # Add prior
  genes$prior_n_genes_locus <- logit10( 1/genes$n_genes )
  genes$prior_n_genes_locus[ genes$n_genes < 2 ] <- logit10(0.75)
  
  
  #-----------------------------------------------------------------------------
  #   Add covariates
  #-----------------------------------------------------------------------------
  
  for( i in names(cov_means) ){
    genes[[i]] <- cov_means[i]
  }
  
  
  #-----------------------------------------------------------------------------
  #   Get raw and normalized predictions
  #-----------------------------------------------------------------------------
  
  # Loop through loci, assign values needed for recalibration
  genes$raw <- predict( object=cal_mod, newdata=genes, type="response" )
  genes$caldera <- as.numeric(NA)
  for( i in unique(genes$locus) ){
    
    # Subset, scale to 1
    sub  <- genes[ genes$locus == i , ]
    scaled <- sub[["raw"]] / sum( sub[["raw"]] )
    
    # Assign values: scaled
    set( x     = genes, 
         i     = which( genes$locus == i ), 
         j     = "caldera", 
         value = scaled )
  }
  
  
  #-----------------------------------------------------------------------------
  #   Format and return
  #-----------------------------------------------------------------------------
  
  # Subset to interesting columns
  bcols <- c( "locus", "locus_pos", "gene", "caldera", "n_genes", 
              "dist", "pops_glo", "coding_glo", "ensgid" )
  gcols <- c( "locus", "locus_pos", "gene", "caldera", "n_genes", 
              "dist", "pops", "coding", "ensgid" )
  genes2 <- genes[ order( genes$locus, -genes$caldera ) , ..bcols ]
  names(genes2) <- gcols
  return(genes2)
  
  
  #-----------------------------------------------------------------------------
  #   Done
  #-----------------------------------------------------------------------------
  
}


#-------------------------------------------------------------------------------
#   Done
#-------------------------------------------------------------------------------























