#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Functions
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   simple_functions
#-------------------------------------------------------------------------------

# aucpr.conf.int.expit: Calculates a logit 95% CI for an AUPRC
aucpr.conf.int.expit <- function( estimate, num.pos, conf.level=0.95 ){
  
  ## convert to logit scale
  est.logit = log(estimate/(1-estimate))
  
  ## standard error (from Kevin Eng)
  se.logit = sqrt( estimate * (1-estimate) / num.pos ) * ( 1/estimate + 1/(1-estimate) )
  
  ## confidence interval in logit
  ci.logit = est.logit + qnorm( c( (1-conf.level)/2, (1+conf.level)/2) ) * se.logit
  
  ## back to original scale
  ci = exp(ci.logit)/(1+exp(ci.logit))
  attr(ci,"conf.level") = conf.level
  attr(ci,"method") = "expit"
  return(ci)
}

# better_pr: Return best precision for a given recall value (and vice versa)
better_pr <- function( p, r, t ){
  df <- as.data.frame(prc$curve)
  names(df) <- c( "recall", "precision", "threshold" )
  idx_r <- max( which( df$recall    > r ) )
  idx_p <- min( which( df$precision > p ) )
  idx_t <- min( which( df$threshold > t ) )
  df2 <- df[ c( idx_r, idx_p, idx_t ) , ]
  row.names(df2) <- c( "same_recall", "same_precision", paste0( "threshold_", t ) )
  df2
}

# ci95_lo and ci95_hi: Return the lower/upper 95% confidence interval
ci95_lo <- function( beta, se, zcrit=qnorm(.025, lower.tail=FALSE) )  beta - zcrit * se
ci95_hi <- function( beta, se, zcrit=qnorm(.025, lower.tail=FALSE) )  beta + zcrit * se

# dev_exp: Calculates the deviance explained by a model
dev_exp <- function( model, in_sample=TRUE ){
  if(in_sample){
    de <- ( model$null.deviance - model$deviance ) / model$null.deviance
  }else{
    de <- ( model$null.deviance - AIC(model) ) / model$null.deviance
  }
  return(de)
}

# logistic10: Convert a log10 odds ratio into a probability 
logistic10 <- function(x) ( 1 / ( 1 + 10^-x ) )

# logit10: Convert a probability into a log10 odds ratio
logit10 <- function(p) log10( p / (1-p) )

# logistic: Convert a logOR into a probability
logistic <- function(x) ( 1 / ( 1 + exp(-x) ) )

# make_formula: Constructs a formula for use in regression
make_formula <- function(lhs, rhs, formula.or.character="formula"){
  
  # Inputs:
  #   lhs: the left-hand side of the regression, should be a vector with one element
  #   rhs: the right-hand side of the regression, should be a vector with as many elements as you like
  rhs_string  <- paste( rhs, collapse=" + ")
  full_string <- paste(lhs, "~", rhs_string )
  if( formula.or.character == "formula" )    out <- as.formula(full_string)
  if( formula.or.character == "character" )  out <- full_string
  if( !( formula.or.character %in% c("formula", "character" ) ) )  stop("Must specify formula or character")
  out
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

# quad_miss_imp: Find the optimal value to impute NAs to for a quadratic relationship
quad_miss_imp <- function( data, model, miss_var ){
  
  # Quadratic equation
  quad <- function(a, b, c){
    a <- as.complex(a)
    answer <- c((-b + sqrt(b^2 - 4 * a * c)) / (2 * a),
                (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
    if(all(Im(answer) == 0)) answer <- Re(answer)
    if(answer[1] == answer[2]) return(answer[1])
    answer
  }
  
  # Solve it
  y <- log( prop_to_odds( mean( data$causal[ data[[miss_var]] == 1 ] ) ) )
  a <- summary(model)$coef[ 3, "Estimate" ]
  b <- summary(model)$coef[ 2, "Estimate" ]
  c <- summary(model)$coef[ 1, "Estimate" ] - y
  quad( a, b, c )
}

# scientific_notation: Takes a numeric vector and returns it in scientific notation with one decimal
scientific_notation <- function(vector) formatC(x=vector, digits=1, format="e")

# v2g_devexp: Deviance explained by BIL/global/relative features for a given method
v2g_devexp <- function( data, prefix ){
  
  forms <- paste0( "causal ~ ", prefix, "_", c( "bil", "glo", "rel" ) )
  mod1  <- glm( as.formula(forms[1]), data=data, family="binomial" )
  mod2  <- glm( as.formula(forms[2]), data=data, family="binomial" )
  mod3  <- glm( as.formula(forms[3]), data=data, family="binomial" )
  de <- c( dev_exp(mod1),
           dev_exp(mod2),
           dev_exp(mod3) )
  names(de) <- paste0( prefix, "_", c( "bil", "glo", "rel" ) )
  de
}

# Function: correlation between BIL/global/relative features for a given method
v2g_cor <- function( data, prefix ){
  cols <- paste0( prefix, c( "_bil", "_glo", "_rel" ) )
  cor( data[,..cols])
}

#-------------------------------------------------------------------------------
#   recalibrate_preds: Get model predictions, plus several re-calibrations
#-------------------------------------------------------------------------------

recalibrate_preds <- function( df, model, method, 
                               bias_cols=NULL, glc_cols=NULL, rm_glc=TRUE ){
  
  #-----------------------------------------------------------------------------
  #   Average covariate values that shouldn't be used for prediction
  #-----------------------------------------------------------------------------
  
  # Fix potentially-biased covariates to their mean value in the test dataset
  for( i in bias_cols ){
    df[[i]] <- mean( df[[i]] )
  }
  
  # Fix gene-level covariates to their mean value in the test dataset
  if(rm_glc){
    for( i in glc_cols ){
      df[[i]] <- mean( df[[i]] )
    }
  }
  
  
  #-----------------------------------------------------------------------------
  #   Extract model predictions
  #-----------------------------------------------------------------------------
  
  # GLM
  if( method == "glm" ){
    df$pred <- predict( object=model, newdata=df, type="response" )
  }
  
  # GLMM
  if( method == "glmm" ){
    df$pred <- predict( object=model, newdata=df, type="response", re.form=NA )
  }
  
  # LASSO or ElNet: lambda.1se
  if( method %in% c( "lasso1", "elnet1" ) ){
    feat_cols <- model$glmnet$beta@Dimnames[[1]]
    df$pred <- predict( model, newx=as.matrix( df[ , ..feat_cols ] ), 
                        s="lambda.1se", type="response" )
  }
  
  # LASSO or ElNet: lambda.min
  if( method %in% c( "lasso2", "elnet2" ) ){
    feat_cols <- model$glmnet$beta@Dimnames[[1]]
    df$pred <- predict( model, newx=as.matrix( df[ , ..feat_cols ] ), 
                        s="lambda.min", type="response" )
  }
  
  # XGBoost
  if( method == "xgb" ){
    feat_cols <- sub( pattern="TRUE$", replacement="", x=model$features[-1] )
    te_data <- data.frame( causal=as.factor( as.integer( df[["causal"]] ) ), 
                           model.matrix( ~.+1, data=df[ , ..feat_cols ] )  )
    te_task <- makeClassifTask( data=te_data, target="causal" )
    df$pred <- predict( model, te_task )$data$prob.1
  }
  
  
  #-----------------------------------------------------------------------------
  #   Loop through TCPs and generate columns needed for recalibration
  #-----------------------------------------------------------------------------
  
  # Initialize output
  out <- data.table( trait=df[["trait"]], causal=df[["causal"]], 
                     original=as.vector(df$pred) )
  out$scaled <- out$best <- out$relative <- out$global <- as.numeric(NA)
  
  # Loop through test set loci
  df$tcp <- paste( df$trait, df$region, df$cs_id, sep="_" )
  for( i in unique(df$tcp) ){
    
    # Subset to locus
    sub <- df[ df$tcp == i , ]
    idx <- order( -sub[["pred"]] )
    
    # Re-calibrate: scaled
    scaled <- sub[["pred"]] / sum( sub[["pred"]] )
    
    # Make variables for global, relative, and BIL scores
    glo <- log( prop_to_odds( sub[["pred"]] ) )
    rel <- glo - max(glo)
    bil <- ifelse( glo == max(glo), 1, 0 )
    
    # Assign values: scaled
    set( x     = out, 
         i     = which( df$tcp == i ), 
         j     = "scaled", 
         value = scaled )
    
    # Assign values: global score
    set( x     = out, 
         i     = which( df$tcp == i ), 
         j     = "global", 
         value = glo )
    
    # Assign values: relative score
    set( x     = out, 
         i     = which( df$tcp == i ), 
         j     = "relative", 
         value = rel )
    
    # Assign values: BIL
    set( x     = out, 
         i     = which( df$tcp == i ), 
         j     = "best", 
         value = bil )
  }
  
  
  #-----------------------------------------------------------------------------
  #   Return
  #-----------------------------------------------------------------------------
  return(out)
}


#-------------------------------------------------------------------------------
#   recalibration_model: Train a LASSO to re-calibrate model outputs
#-------------------------------------------------------------------------------

recalibration_model <- function(data){
  
  # Build the model to predict causal genes based on step 1 predictions
  data$foldid <- as.integer( as.factor(data$trait) )
  pred_cols <- c( "global", "relative" )
  mod_las <- cv.glmnet( x=as.matrix( data[ , ..pred_cols ] ), 
                        y=data[["causal"]], family="binomial", foldid=data$foldid )
  
  # Fit a GAM
  mod_gam <- gam( causal ~ s(original), data=data )
  
  # Return
  out <- list()
  out$lasso <- mod_las
  out$gam   <- mod_gam
  return(out)
}


#-------------------------------------------------------------------------------
#   loto:            Leave-one-trait-out models, PR curves, and predictions
#-------------------------------------------------------------------------------

loto <- function( data, method, feat_cols, 
                  bias_cols=NULL, glc_cols=NULL, rm_glc=FALSE,
                  backwards_selection=FALSE, maxit=10 ){
  
  # Format input
  data <- data[ order( data$trait, data$region, data$cs_id, data$start ) , ]
  
  # Initialize output
  out <- list( trait       = unique(data$trait), 
               model       = list(), 
               recal_model = list(), 
               preds       = data.table( trait  = data$trait,
                                         causal = data$causal,
                                         pred   = NA,
                                         recal  = NA,
                                         local  = NA,
                                         scaled = NA ) )
  
  # Loop through traits (to leave out)
  for( i in out$trait ){
    
    #---------------------------------------------------------------------------
    #   Get set up
    #---------------------------------------------------------------------------
    
    # Message
    message( date(), "   Starting trait ", which( out$trait == i ), "/", 
             length(out$trait), ": ", i )
    
    # Split into training and testing sets
    trn <- data[ data$trait != i , ]
    tst <- data[ data$trait == i , ]
    
    # Ensure that each trait is its own fold for CV (to prevent leakage)
    trn$foldid <- as.integer( as.factor(trn$trait) )
    
    # Fix potentially-biased covariates to their mean value in the test dataset
    for( j in bias_cols ){
      tst[[i]] <- mean( tst[[j]] )
    }
    xgb_cols <- setdiff( feat_cols, bias_cols )
    
    # Fix gene-level covariates to their mean value in the test dataset
    if(rm_glc){
      for( j in glc_cols ){
        tst[[j]] <- mean( tst[[j]] )
      }
      xgb_cols <- setdiff( xgb_cols, glc_cols )
    }
    
    
    #---------------------------------------------------------------------------
    #   Train models
    #---------------------------------------------------------------------------
    
    # GLM
    if( method == "glm" ){
      l_form <- make_formula( lhs=feat_cols[1], rhs=feat_cols[-1] )
      l_glm0 <- glm( formula=l_form, data=trn, family="binomial" )
      if(backwards_selection){
        l_mod <- stats::step( object=l_glm0, k=qchisq( 1-0.05/5, df=1 ), trace=0 )
      }else{
        l_mod <- l_glm0
      }
    }
    
    # GLMM
    if( method == "glmm" ){
      library(lme4)
      l_form <- make_formula( lhs=feat_cols[1], rhs=feat_cols[-1] )
      l_form <- update( l_form, . ~ . + ( 1 | ensgid_c ) )
      l_mod  <- glmer( formula=l_form, data=trn, family="binomial" )
    }
    
    # LASSO
    if( method %in% c( "lasso1", "lasso2" ) ){
      l_mod <- cv.glmnet( x=as.matrix( trn[ , ..feat_cols ] )[,-1],
                          y=trn[["causal"]], family="binomial", 
                          foldid=trn$foldid, alpha=1 )
    }
    
    # Elastic net
    if( method %in% c( "elnet1", "elnet2" ) ){
      l_mod <- cv.glmnet( x=as.matrix( trn[ , ..feat_cols ] )[,-1], 
                          y=trn[["causal"]], family="binomial", 
                          foldid=trn$foldid, alpha=0.5 )
    }
    
    # XGB
    if( method == "xgb" ){
      l_mod <- suppressMessages( suppressWarnings(
        train_xgb( data=trn, feat_cols=xgb_cols[-1], maxit=maxit ) ) )
    }
    
    
    #---------------------------------------------------------------------------
    #   Get test set predictions and various calibrations
    #---------------------------------------------------------------------------
    
    # Get model predictions (training and test), plus info for calibration
    preds_tr <- recalibrate_preds( df=trn, model=l_mod, method=method, 
                                   bias_cols=bias_cols )
    preds_te <- recalibrate_preds( df=tst, model=l_mod, method=method, 
                                   bias_cols=bias_cols )
    
    # Train LASSO and GAM smodel to recalibrate model outputs
    cal_mod <- recalibration_model( data=preds_tr )
    
    # Apply GAM model: training
    preds_tr$gam <- predict( object=cal_mod$gam, newdata=preds_tr )
    preds_tr$gam[ preds_tr$gam < 0.001 ] <- 0.001
    preds_tr$gam[ preds_tr$gam > 0.999 ] <- 0.999
    
    # Apply GAM model: testing
    preds_te$gam <- predict( object=cal_mod$gam, newdata=preds_te )
    preds_te$gam[ preds_te$gam < 0.001 ] <- 0.001
    preds_te$gam[ preds_te$gam > 0.999 ] <- 0.999
    
    # Apply LASSO model
    pred_cols <- c( "global", "relative" )
    preds_tr$modeled <- predict( object=cal_mod$lasso, s="lambda.min", type="response",
                                 newx=as.matrix( preds_tr[ , ..pred_cols ] ) )
    preds_te$modeled <- predict( object=cal_mod$lasso, s="lambda.min", type="response",
                                 newx=as.matrix( preds_te[ , ..pred_cols ] ) )
    
    
    #---------------------------------------------------------------------------
    #   Save outputs
    #---------------------------------------------------------------------------
    
    # Models
    out$model[[i]]       <- l_mod
    out$recal_model[[i]] <- cal_mod
    
    # Predicted probabilities
    out$preds$pred[   out$preds$trait == i ] <- preds_te$original
    out$preds$scaled[ out$preds$trait == i ] <- preds_te$scaled
    out$preds$recal[  out$preds$trait == i ] <- preds_te$gam
    out$preds$local[  out$preds$trait == i ] <- preds_te$modeled
    out$preds$bil[    out$preds$trait == i ] <- preds_te$best
  }
  
  # Return
  return(out)
}


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Read in data
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------


# Arguments
maindir <- "~/projects/causal_genes/"

# Load libraries and sources
suppressPackageStartupMessages( library(corrplot) )
suppressPackageStartupMessages( library(cowplot) )
suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(e1071) )
suppressPackageStartupMessages( library(glmnet) )
suppressPackageStartupMessages( library(ggrepel) )
suppressPackageStartupMessages( library(gridExtra) )
suppressPackageStartupMessages( library(lme4) )
suppressPackageStartupMessages( library(mgcv) )
suppressPackageStartupMessages( library(mlr) )
suppressPackageStartupMessages( library(PRROC) )
suppressPackageStartupMessages( library(parallel) )
suppressPackageStartupMessages( library(parallelMap) )
suppressPackageStartupMessages( library(patchwork) )
suppressPackageStartupMessages( library(probably) )
suppressPackageStartupMessages( library(RColorBrewer) )
suppressPackageStartupMessages( library(rpart.plot) )
suppressPackageStartupMessages( library(tidymodels) )
suppressPackageStartupMessages( library(tidyr) )
suppressPackageStartupMessages( library(tidyverse) )
suppressPackageStartupMessages( library(viridis) )
suppressPackageStartupMessages( library(xgboost) )

# Read in causal/non-causal TGPs annotated with gene mapping evidence
cnc_map_file <- file.path( maindir, "causal_noncausal_trait_gene_pairs",
                           "causal_tgp_and_gene_mapping_data_300kb.tsv" )
dt <- fread(file=cnc_map_file)
dt <- dt[ order( dt$trait, dt$region, dt$cs_id, dt$start ) , ]

# Split training and testing traits 
# Sort traits by number of causal genes and pick every fifth one
all_traits  <- sort( table( dt$trait[ dt$causal ] ), decreasing=TRUE )
te_traits <- names(all_traits)[ seq( 4, length(all_traits), by=4 ) ]
tr_traits <- setdiff( names(all_traits), te_traits )

#read data and make training datatable
dt <- fread(file=cnc_map_file)
dt <- dt[ order( dt$trait, dt$region, dt$cs_id, dt$start ) , ]
tr <- dt

# Create an empty list to store results
saved_matches <- list()

# Get unique entries in the ensg_c column
unique_ensg_c <- unique(tr$ensgid_c)
seedno <- 0
# Loop through each unique ensg_c
for (ensg in unique_ensg_c) {
  # Increment the seed
  seedno <- seedno + 1
  
  # Subset the dataframe for the current ensg_c
  subset_ensg <- tr[tr$ensgid_c == ensg, ]
  
  # Set the seed for reproducibility
  set.seed(seedno)
  
  # Pick a unique row using the seed
  picked_row <- subset_ensg[sample(nrow(subset_ensg), 1), ]
  
  # Save trait, region, and cs from the picked row
  saved_matches[[ensg]] <- picked_row[, c("trait", "region", "cs_id")]
}

# Combine saved matches into a dataframe
saved_matches_df <- do.call(rbind, saved_matches)


# Convert saved_matches_df to a data frame
saved_matches_df <- as.data.frame(saved_matches_df)

# Perform a join to subset the original dataframe
subset_tr <- tr %>%
  inner_join(saved_matches_df, by = c("trait", "region", "cs_id"))


# Run adjusted LOTO LASSO training
bl_las <- loto( data=subset_tr, method="lasso2", feat_cols=b_cols )
bgl_las_all  <- loto( data=tr, method="lasso2", feat_cols=bg_cols )

#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Figures
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Figure S2: AUPRC barplots all vs unique genes
#-------------------------------------------------------------------------------

# Extract AUPRC and 95% CI for all models
bgl_las_pr <- plot_loto_pr(loto_obj = bgl_las, color = viridis(11), legend = TRUE)
bgl_las_all_pr <- plot_loto_pr(loto_obj = bgl_las_all, color = viridis(11), legend = TRUE)



fig1_file <- file.path( "figure_S2.jpg" )
jpeg( filename=fig1_file, width=300*4, height=300*4, res=75*4 )

par(mar = c(7, 5, 1, 1))

f_pr <- rbind(bgl_las_all_pr, bgl_las_pr)
dimnames(f_pr)[[1]] <- c("All genes", "Unique genes")
print(f_pr) # Ensure the data is as expected

# Plot barplot
bp <- barplot(height = f_pr[, "auprc"], las = 2, 
              ylab = "AUPRC (±95% CI)", ylim = c(0, 0.7),
              col = brewer.pal(n = 5, name = "Greens")[1:2])
for (i in seq_along(bp)) {
  lines(x = c(bp[i], bp[i]), lwd = 1.5,
        y = c(f_pr[i, "lo"], f_pr[i, "hi"]))
}

dev.off()

