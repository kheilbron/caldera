
#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Functions
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   simple_functions
#-------------------------------------------------------------------------------

# agresti_95ci_prop: The Agresti-Coull 95% confidence interval for a proportion
agresti_95ci_prop <- function( X, n, Z=qnorm(.025, lower.tail=FALSE ) ){
  #X = number of successes
  #n = total number of trials
  #Z = z-score for CI threshold
  p <- (X + 2) / (n + 4)
  num <- p * (1-p)
  denom <- n+4
  diff <- Z * sqrt(num/denom)
  lo <- p - diff
  hi <- p + diff
  out <- c( X/n, lo, hi)
  names(out) <- c( "prop", "lo", "hi" )
  return(out)
}

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
  
  # Return
  out <- c( estimate, ci, est.logit, se.logit )
  names(out) <- c( "auprc", "lo", "hi", "auprc_logit", "se_logit" )
  return(out)
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

# logistic(10): Convert a log(10)OR into a probability
logistic   <- function(x) ( 1 / ( 1 + exp(-x) ) )
logistic10 <- function(x) ( 1 / ( 1 + 10^-x ) )

# logit(10): Convert a probability into a log or log10 odds ratio
logit   <- function(p) log( p / (1-p) )
logit10 <- function(p) log10( p / (1-p) )

# make_formula: Constructs a formula for use in regression
make_formula <- function(lhs, rhs, formula.or.character="formula"){
  
  # Inputs:
  #   lhs: the left-hand side of the regression, a vector with one element
  #   rhs: the right-hand side of the regression, a vector with 1+ elements
  rhs_string  <- paste( rhs, collapse=" + ")
  full_string <- paste(lhs, "~", rhs_string )
  if( formula.or.character == "formula" )    out <- as.formula(full_string)
  if( formula.or.character == "character" )  out <- full_string
  if( !( formula.or.character %in% c("formula", "character" ) ) ){
    stop("Must specify formula or character")
  }
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

# scientific_notation: 
# Takes a numeric vector and returns it in scientific notation with one decimal
scientific_notation <- function(vector) formatC(x=vector, digits=1, format="e")

# two_sample_z_test: Return a P value from a two-sample z-test
two_sample_z_test <- function( mean1, mean2, se1, se2, alpha=0.05 ){
  z <- ( mean1 - mean2 ) / sqrt( se1^2 + se2^2 )
  p <- 2 * (1 - pnorm( abs(z) ) )
  return(p)
}

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

# v2g_cor: correlation between BIL/global/relative features for a given method
v2g_cor <- function( data, prefix ){
  cols <- paste0( prefix, c( "_bil", "_glo", "_rel" ) )
  cor( data[,..cols])
}


#-------------------------------------------------------------------------------
#   forest_plot:     Makes a forest plot
#-------------------------------------------------------------------------------

forest_plot <- function( df, value.col="effect", lo.col="lo", hi.col="hi", 
                         label.col=NULL, mtext.col="pvalue", colour.col=NULL, 
                         xlab="", xmin=NULL, xmax=NULL, exp.x.labels=NULL, 
                         vert.line.pos=c(-1,0,1), margins=c(4.5,11,1,5), 
                         main="" ){
  
  # If a data.table is provided, convert to a data.frame
  if( any( class(df) == "data.table" ) ){
    df <- as.data.frame(df)
  }
  
  # Prepare x- and y-axes
  ys <- 1:nrow(df)
  rev.ys <- rev(ys)
  if( is.null(lo.col) & is.null(hi.col) ){
    lo.col <- "lo"
    hi.col <- "hi"
    df[,lo.col] <- df[,hi.col] <- df[,value.col]
  }
  if( is.null(label.col) ){
    label.col <- "label"
    df[,label.col] <- row.names(df)
  }
  if( is.null(xmin) ) xmin <- min( c(0, df[,lo.col]) )
  if( is.null(xmax) ) xmax <- max( c(0, df[,hi.col]) )
  if( !is.null(exp.x.labels) ){
    x.labels <- log(exp.x.labels)
    xaxt <- "n"
  } else xaxt <- NULL
  if( is.null(colour.col) ){
    df$colour.col <- "black"
    colour.col <- "colour.col"
  } 
  
  # Plot
  par(mar=margins)
  plot(df[,value.col], rev.ys, pch=19, xlim=c(xmin, xmax), ylim=c( min(ys)-0.5, max(ys)+0.5 ), 
       xaxt=xaxt, yaxt="n", ylab="", xlab=xlab, col=df[,colour.col], main=main)
  if( !is.null(mtext.col)) axis( side=4, at=rev.ys, labels=df[,mtext.col], las=1 )
  axis( side=2, at=rev.ys, labels=df[,label.col], las=1 )
  if( !is.null(exp.x.labels) ) axis( side=1, at=x.labels, labels=exp.x.labels, las=1 )
  
  # Add lines
  for(i in vert.line.pos){
    abline(v=i, lty=2)
  }
  for(i in ys){
    lines( c(df[rev.ys[i],lo.col], df[rev.ys[i],hi.col]), c(i,i), lwd=2, col=df[rev.ys[i],colour.col] )
  }
}


#-------------------------------------------------------------------------------
#   glm_to_forest:   Make a forest plot from a GLM object
#-------------------------------------------------------------------------------

glm_to_forest <- function( mod, suffix="TRUE", mtext=TRUE, xmax=NULL, 
                           standardize=FALSE ){
  
  # Clean up row and column names
  mod2 <- as.data.frame( summary(mod)$coef[ -1, , drop=FALSE ] )
  names(mod2) <- c( "effect", "se", "z", "p" )
  pattern <- paste0( suffix, "$" )
  mod3 <- mod2[ grepl( pattern=pattern, row.names(mod2) ) , ]
  full_names  <- row.names(mod3)
  short_names <- sub( pattern, "", row.names(mod3) )
  row.names(mod3) <- short_names
  
  # If specified, standardize effects to reflect a 1 SD change
  if(standardize){
    stddevs <- apply( X=mod$model, MARGIN=2, FUN=sd )[-1]
    mod3$effect <- mod3$effect * stddevs
    mod3$se     <- mod3$se     * stddevs
  }
  
  # Add columns for plotting
  mod3$or <- exp(mod3$effect)
  mod3$or_lo <- exp( ci95_lo( beta=mod3$effect, se=mod3$se ) )
  mod3$or_hi <- exp( ci95_hi( beta=mod3$effect, se=mod3$se ) )
  mod3$lo <- ci95_lo( beta=mod3$effect, se=mod3$se )
  mod3$hi <- ci95_hi( beta=mod3$effect, se=mod3$se )
  mod3$col <- ifelse( mod3$p < 0.05, "#70AD47", "grey" )
  mod3$pvalue <- paste0( "P=", scientific_notation(mod3$p) )
  if(mtext) mtext.col="pvalue" else mtext.col=NULL
  
  # Plot
  forest_plot( df=mod3, colour.col="col", xlab="Odds ratio (±95% CI)", xmax=xmax, 
               margins=c(5,9,1,6), vert.line.pos=1, mtext.col=mtext.col,
               value.col="or", lo.col="or_lo", hi.col="or_hi" )
  par( mar=c(5,5,4,1) )
}


#-------------------------------------------------------------------------------
#   glm_to_forest_p: As above, but for P values
#-------------------------------------------------------------------------------

glm_to_forest_p <- function( mod, suffix="", mtext=TRUE, xmax=NULL ){
  
  # Clean up row and column names
  mod2 <- as.data.frame( summary(mod)$coef[-1,] )
  names(mod2) <- c( "effect", "se", "z", "p" )
  pattern <- paste0( suffix, "$" )
  mod3 <- mod2[ grepl( pattern=pattern, row.names(mod2) ) , ]
  full_names  <- row.names(mod3)
  short_names <- sub( pattern, "", row.names(mod3) )
  row.names(mod3) <- short_names
  
  # Add columns for plotting
  mod3$logp <- -log10(mod3$p)
  mod3$col <- ifelse( mod3$p < 0.05/NROW(mod3), "#70AD47", "grey" )
  mod3$pvalue <- paste0( "P=", scientific_notation(mod3$p) )
  if(mtext) mtext.col="pvalue" else mtext.col=NULL
  
  # Plot
  forest_plot( df=mod3, colour.col="col", xlab="-log10 P value", xmax=xmax, 
               margins=c(5,10,1,1), vert.line.pos=0, mtext.col=mtext.col,
               value.col="logp", lo.col=NULL, hi.col=NULL )
  par( mar=c(5,5,4,1) )
}


#-------------------------------------------------------------------------------
#   lasso_to_forest: Make a forest plot from a LASSO object
#-------------------------------------------------------------------------------

lasso_to_forest <- function( las, lambda_type="lambda.1se", standardize=TRUE, 
                             train=NULL, rm_feats=NULL, xmax=NULL ){
  
  # Extract coefficients
  mod1 <- as.data.frame( coef( object=las, s=lambda_type )[ -1, ] )
  names(mod1) <- "effect"
  
  # If specified, standardize coefficients to reflect a 1 SD change
  if(standardize){
    feat_cols <- las$glmnet.fit$beta@Dimnames[[1]]
    stddevs <- apply( X=train[ , ..feat_cols ], MARGIN=2, FUN=sd )
    mod1$effect <- mod1$effect * stddevs
  }
  
  # Remove dropped features
  mod2 <- mod1[ mod1$effect !=0 , , drop=FALSE ]
  
  # Drop specified (bias) features
  mod3 <- mod2[ setdiff( row.names(mod2), rm_feats ) , , drop=FALSE ]
  
  # Add columns for plotting
  mod3$or <- exp(mod3$effect)
  mod3$col <- brewer.pal( n=6, name="Greens" )[5]
  mod3 <- mod3[ order(-mod3$effect) , ]
  
  # Clean up labels
  clean_labels <- TRUE
  if(clean_labels){
    row.names(mod3) <- sub( x=row.names(mod3), pattern="_rel", 
                            replacement=", relative" )
    row.names(mod3) <- sub( x=row.names(mod3), pattern="_glo", 
                            replacement=", global" )
    row.names(mod3) <- sub( x=row.names(mod3), pattern="_bil", 
                            replacement=", best-in-locus" )
    row.names(mod3) <- sub( x=row.names(mod3), pattern="pops", 
                            replacement="PoPS" )
    row.names(mod3) <- sub( x=row.names(mod3), pattern="magma", 
                            replacement="MAGMA" )
    row.names(mod3) <- sub( x=row.names(mod3), pattern="dist_gene", 
                            replacement="Distance*" )
    row.names(mod3) <- sub( x=row.names(mod3), pattern="coding", 
                            replacement="Coding PIP" )
    row.names(mod3) <- sub( x=row.names(mod3), pattern="prior_n_genes_locus", 
                            replacement="# local genes*" )
  }
  
  # Plot
  xmin <- min( c( 1, min(mod3$or) ) )
  forest_plot( df=mod3, colour.col="col", xlab="Odds ratio", xmin=xmin, 
               xmax=xmax, margins=c(5,10,1,1), vert.line.pos=1, mtext.col=NULL,
               value.col="or", lo.col="or", hi.col="or" )
  par( mar=c(5,5,4,1) )
}


#-------------------------------------------------------------------------------
#   n_inc_las:       Plot % of LASSO LOTO models containing each feature
#-------------------------------------------------------------------------------

n_inc_las <- function(loto_obj){
  
  # Extract whether a feature is present or not
  out <- list()
  for( i in loto_obj$trait ){
    sub <- loto_obj$model[[i]]
    sub2 <- sub$glmnet.fit$beta[ , sub$index[1] ]
    out[[i]] <- sub2 != 0
  }
  out2 <- as.data.frame( do.call( rbind, out ) )
  
  # Aggregate across traits
  out3 <- sort( colSums(out2) / NROW(out2), decreasing=TRUE )
  out4 <- data.frame( row.names=names(out3), inc=out3 )
  out4$col <- ifelse( out4$inc >= 0.25, "lightgreen", "grey" )
  out4$col[ out4$inc >= 0.5 ] <- "#70AD47"
  forest_plot( df=out4, xlab="Proportion included in LASSO model", colour.col="col",
               margins=c(5,11,4,1), vert.line.pos=c(0,1), mtext.col=NULL, 
               value.col="inc", lo.col="inc", hi.col="inc", xmax=1, main="LASSO" )
  
  # Aggregate across groups of traits
  groups <- sub( pattern="_bil$|_glo$|_rel$", replacement="", x=names(out3) )
  out5 <- sort( tapply( X=out3, INDEX=groups, FUN=mean ), decreasing=TRUE )
  out6 <- data.frame( row.names=names(out5), inc=out5 )
  out6$col <- ifelse( out6$inc >= 0.25, "lightgreen", "grey" )
  out6$col[ out6$inc >= 0.5 ] <- "#70AD47"
  forest_plot( df=out6, xlab="Proportion included in LASSO model", colour.col="col",
               margins=c(5,11,4,1), vert.line.pos=c(0,1), mtext.col=NULL, 
               value.col="inc", lo.col="inc", hi.col="inc", xmax=1, main="LASSO" )
}


#-------------------------------------------------------------------------------
#   train_xgb:       Train an XGBoost model
#-------------------------------------------------------------------------------

train_xgb <- function( data, feat_cols, bias_cols=NULL, glc_cols=NULL, 
                       rm_glc=FALSE, maxit=20 ){
  
  # Remove potentially-biased covariates from the feature list
  feat_cols <- setdiff( feat_cols, bias_cols )
  
  # If specified, remove gene-level covariates
  if(rm_glc){
    feat_cols <- setdiff( feat_cols, glc_cols )
  }
  
  # XGBoost with basic hyperparameter tuning
  tr_data <- data.frame( causal=as.factor( as.integer( data[["causal"]] ) ), 
                         model.matrix( ~.+1, data=data[ , ..feat_cols ] ) )
  tr_task <- makeClassifTask( data=tr_data, target="causal", 
                              blocking=as.factor(data$trait) )
  lrn <- makeLearner( cl="classif.xgboost", predict.type="prob",
                      objective="binary:logistic", verbose=0 )
  rdesc <- makeResampleDesc( method="CV", blocking.cv=TRUE, iters=length( unique(data$trait) ) )
  parallelStartSocket( cpus=detectCores() )
  pars <- makeParamSet( makeDiscreteParam( "booster",          values = "gbtree" ), 
                        makeIntegerParam(  "nrounds",          lower = 100L,  upper = 250L ), 
                        makeNumericParam(  "eta",              lower = 0.01, upper = 0.15 ), 
                        makeIntegerParam(  "max_depth",        lower = 1L,   upper = 6L ), 
                        makeNumericParam(  "min_child_weight", lower = 0,    upper = 5 ), 
                        makeNumericParam(  "gamma",            lower = 0,    upper = 8 ), 
                        makeNumericParam(  "subsample",        lower = 0.5,  upper = 1 ), 
                        makeNumericParam(  "colsample_bytree", lower = 0.1,  upper = 0.6 ) )
  ctrl <- makeTuneControlRandom( maxit=maxit )
  mytune <- tuneParams( learner    = lrn, 
                        task       = tr_task, 
                        measures   = logloss,
                        resampling = rdesc, 
                        par.set    = pars, 
                        control    = ctrl, 
                        show.info  = TRUE )
  lrn_tune <- setHyperPars( lrn, par.vals=mytune$x )
  xgb_mod <- train( learner=lrn_tune, task=tr_task )
  return(xgb_mod)
}


#-------------------------------------------------------------------------------
#   xgb_to_forest:   Make a forest plot from an XGBoost model
#-------------------------------------------------------------------------------

xgb_to_forest <- function( xgb_mod, suffix="_bilTRUE", xmax=NULL ){
  
  # Format
  pattern <- paste0( suffix, "$" )
  fi  <- getFeatureImportance(xgb_mod)$res
  fi2 <- fi[ grepl( pattern=pattern, fi$variable ) , ]
  fi3 <- data.frame( row.names = sub( suffix, "", fi2$variable ),
                     fi        = fi2$importance,
                     col       = ifelse( fi2$importance == 0, "grey", "#70AD47" ) )
  fi3 <- fi3[ order(-fi3$fi) , ]
  
  # Plot
  forest_plot( df=fi3, colour.col="col", xlab="Feature importance", 
               margins=c(5,9,1,1), vert.line.pos=0, mtext.col=NULL, 
               value.col="fi", lo.col="fi", hi.col="fi", xmax=xmax )
  par( mar=c(5,5,4,1) )
}


#-------------------------------------------------------------------------------
#   plot_fi:         Plot grouped feature importance
#-------------------------------------------------------------------------------

plot_fi <- function(xgb_mod){
  fi <- getFeatureImportance(xgb_mod)$res[-1,]
  fi$group <- sub( pattern="_bilTRUE$|_glo$|_rel$", replacement="", x=fi$variable )
  fi2 <- sort( tapply( X=fi$importance, INDEX=fi$group, FUN=sum ), decreasing=TRUE )
  fi3 <- data.frame( row.names=names(fi2), fi=fi2, cumsum=cumsum(fi2) )
  fi3$col <- ifelse( fi3$cumsum < 0.99, "lightgreen", "grey" )
  fi3$col[ fi3$cumsum < 0.9 ] <- "#70AD47"
  forest_plot( df=fi3, xlab="Feature importance", colour.col="col", main="XGBoost",
               margins=c(5,9,4,1), vert.line.pos=0, mtext.col=NULL, 
               value.col="fi", lo.col="fi", hi.col="fi", xmax=NULL )
  par( mar=c(5,5,4,1) )
}

#-------------------------------------------------------------------------------
#   auprc_test_set:  Barplot of AUPRCs from a GLM, LASSO, and XGBoost model
#-------------------------------------------------------------------------------

auprc_test_set <- function( test_df, rand_mod=NULL, glm_mod=NULL, glmm_mod=NULL,
                            las_mod=NULL, xgb_mod=NULL, ymax=NULL, 
                            bias_cols=NULL, glc_cols=NULL, rm_glc=FALSE ){
  
  # Fix potentially-biased covariates to their mean value in the test dataset
  for( i in bias_cols ){
    test_df[[i]] <- mean( test_df[[i]] )
  }
  
  # Fix gene-level covariates to their mean value in the test dataset
  if(rm_glc){
    for( i in glc_cols ){
      test_df[[i]] <- mean( test_df[[i]] )
    }
  }
  
  # Random predictor
  if( !is.null(rand_mod) ){
    p_rand0 <- predict( object=rand_mod, newdata=test_df, type="response" )
    p_rand  <- pr.curve( scores.class0 = p_rand0[  test_df$causal ], 
                         scores.class1 = p_rand0[ !test_df$causal ], curve = TRUE )
  }else{
    p_rand  <- NULL
  }
  
  # GLM
  if( !is.null(glm_mod) ){
    p_glm0 <- predict( object=glm_mod, newdata=test_df, type="response" )
    p_glm  <- pr.curve( scores.class0 = p_glm0[  test_df$causal ], 
                        scores.class1 = p_glm0[ !test_df$causal ], curve = TRUE )
  }else{
    p_glm  <- NULL
  }
  
  # GLMM
  if( !is.null(glmm_mod) ){
    p_glmm0 <- predict( object=glmm_mod, newdata=test_df, type="response", re.form=NA )
    p_glmm  <- pr.curve( scores.class0 = p_glmm0[  test_df$causal ], 
                         scores.class1 = p_glmm0[ !test_df$causal ], curve = TRUE )
  }else{
    p_glmm  <- NULL
  }
  
  # LASSO
  if( !is.null(las_mod) ){
    feat_cols <- las_mod$glmnet$beta@Dimnames[[1]]
    p_las0 <- predict( las_mod, newx=as.matrix( test_df[ , ..feat_cols ] ), 
                       s="lambda.min", type="response" )
    p_las  <- pr.curve( scores.class0 = p_las0[  test_df$causal ], 
                        scores.class1 = p_las0[ !test_df$causal ], curve = TRUE )
  }else{
    p_las <- NULL
  }
  
  # XGBoost
  if( !is.null(xgb_mod) ){
    feat_cols <- sub( pattern="TRUE$", replacement="", x=xgb_mod$features[-1] )
    te_data <- data.frame( causal=as.factor( as.integer( test_df[["causal"]] ) ), 
                           model.matrix( ~.+1, data=test_df[ , ..feat_cols ] )  )
    te_task <- makeClassifTask( data=te_data, target="causal" )
    p_xgb0 <- predict( xgb_mod, te_task )
    p_xgb  <- pr.curve( scores.class0 = p_xgb0$data$prob.1[  test_df$causal ], 
                        scores.class1 = p_xgb0$data$prob.1[ !test_df$causal ], 
                        curve = TRUE )
  }else{
    p_xgb <- NULL
  }
  
  # Compare AUPRCs
  auprcs <- c( Random=p_rand$auc.integral,  GLM=p_glm$auc.integral, 
               GLMM=p_glmm$auc.integral, LASSO=p_las$auc.integral, 
               XGBoost=p_xgb$auc.integral )
  if( is.null(ymax) ) ymax <- max(auprcs)
  par( mar=c(5,9,1,1) )
  barplot( auprcs, las=2, ylab="AUPRC", ylim=c(0,ymax),
           col=brewer.pal( n=length(auprcs), name="RdBu" ) )
  par( mar=c(5,5,4,1) )
  return(auprcs)
}


#-------------------------------------------------------------------------------
#   plot_pr:         Plot precision-recall curve
#-------------------------------------------------------------------------------

plot_pr <- function( model, test_df, p=0.78, r=0.34, v=0.75, 
                     bias_cols=NULL, glc_cols=NULL, rm_glc=TRUE ){
  
  # Fix potentially-biased covariates to their mean value in the test dataset
  for( i in bias_cols ){
    test_df[[i]] <- mean( test_df[[i]] )
  }
  
  # Fix gene-level covariates to their mean value in the test dataset
  if(rm_glc){
    for( i in glc_cols ){
      test_df[[i]] <- mean( test_df[[i]] )
    }
  }
  
  # GLM
  if( any( class(model) == "glm" ) ){
    p_mod <- predict( object=model, newdata=test_df, type="response" )
    pr    <- pr.curve( scores.class0 = p_mod[  test_df$causal ], 
                       scores.class1 = p_mod[ !test_df$causal ], curve = TRUE )
    
  # LASSO
  }else if( class(model) == "cv.glmnet" ){
    feat_cols <- model$glmnet$beta@Dimnames[[1]]
    p_mod <- predict( model, newx=as.matrix( test_df[ , ..feat_cols ] ), 
                       s="lambda.min", type="response" )
    pr    <- pr.curve( scores.class0 = p_mod[  test_df$causal ], 
                       scores.class1 = p_mod[ !test_df$causal ], curve = TRUE )
    
  # XGBoost
  }else if( class(model) == "WrappedModel" ){
    feat_cols <- sub( pattern="TRUE$", replacement="", x=model$features[-1] )
    te_data <- data.frame( causal=as.factor( as.integer( test_df[["causal"]] ) ), 
                           model.matrix( ~.+1, data=test_df[ , ..feat_cols ] )  )
    te_task <- makeClassifTask( data=te_data, target="causal" )
    p_mod <- predict( model, te_task )
    pr    <- pr.curve( scores.class0 = p_mod$data$prob.1[  test_df$causal ], 
                       scores.class1 = p_mod$data$prob.1[ !test_df$causal ], 
                       curve = TRUE )
  }
  
  # Plot
  plot( pr, scale.color=viridis(11), las=1 )
  
  # Find the point on the curve that meets the specified criteria: precision
  if( !is.null(p) ){
    idx1 <- head( which( pr$curve[,2] > p ), 1 )
    point1 <- pr$curve[idx1,]
  }else{
    point1 <- NULL
  }
  
  # Find the point on the curve that meets the specified criteria: recall
  if( !is.null(r) ){
    idx2 <- tail( which( pr$curve[,1] > r ), 1 )
    point2 <- pr$curve[idx2,]
  }else{
    point2 <- NULL
  }
  
  # Find the point on the curve that meets the specified criteria: value
  if( !is.null(v) ){
    idx3 <- head( which( pr$curve[,3] > v ), 1 )
    point3 <- pr$curve[idx3,]
  }else{
    point3 <- NULL
  }
  
  # Combine results
  dt <- as.data.table( rbind( point1, point2, point3 ) )
  names(dt) <- c( "recall", "precision", "value" )
  setcolorder( x=dt, neworder=c( "precision", "recall", "value" ) )
  return(dt)
}


#-------------------------------------------------------------------------------
#   plot_mod_pr:        Plot PR curve and return AUPRC ± 95% CI
#-------------------------------------------------------------------------------

plot_mod_pr <- function( pred_col, data, make_plot=TRUE ){
  
  # Extract the curve
  prc <- list()
  out0 <- list()
  for( i in pred_col ){
    prc[[i]] <- pr.curve( scores.class0 = data[[i]][ data$causal == 1 ], 
                          scores.class1 = data[[i]][ data$causal == 0 ], 
                          curve = TRUE )
    out0[[i]] <- aucpr.conf.int.expit( estimate = prc[[i]]$auc.integral,
                                       num.pos  = sum(data$causal) )
  }
  out <- as.data.frame( do.call( rbind, out0 ) )
  
  # Prepare to plot depending on whether it is one trait or many
  if( length(pred_col) == 1 ){
    legend       <- TRUE
    scale.color  <- viridis(11)
    color        <- TRUE
    names(color) <- pred_col
    main         <- pred_col
    auc.main     <- TRUE
    legend_names <- NULL
    par( mar=c(5,5,4,1) )
  }else{
    legend       <- FALSE
    scale.color  <- NULL
    color        <- c( brewer.pal( n=6, name="Greens" )[5], 
                       brewer.pal( n=6, name="Reds" )[4], 
                       brewer.pal( n=6, name="Blues" )[4], 
                       brewer.pal( n=6, name="Oranges"  )[3] )
    names(color) <- c( "CALDERA", "FLAMES", "L2G", "cS2G" )
    main         <- ""
    auc.main     <- FALSE
    par( mar=c(5,5,1,1) )
    legend_names <- rev( paste0( row.names(out), " (", 
                                 round( out$auprc, digits=2 ), ")" ) )
  }
  
  # Plot
  if(make_plot){
    for( i in pred_col ){
      add <- pred_col[1] != i
      plot( prc[[i]], scale.color=scale.color, color=color[i], las=1, 
            legend=legend, add=add, main=main, auc.main=auc.main )
    }
    if( length(pred_col) > 1 ){
      legend( x=1, y=1, xjust=1, legend=legend_names,
              col=color, lty=1, lwd=3, cex=0.5 )
      # legend( x=0.05, y=0, yjust=0, legend=legend_names,
      #         col=color, lty=1, lwd=3, cex=0.7 )
    }
  }
  
  # Return AUPRC and 95% CI
  return(out)
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
  out$global <- logit(out$original)
  out$scaled <- out$best <- out$relative <- as.numeric(NA)
  
  # Loop through test set loci
  df$tcp <- paste( df$trait, df$region, df$cs_id, sep="_" )
  for( i in unique(df$tcp) ){
    
    # Subset to locus
    sub <- df[ df$tcp == i , ]
    
    # Re-calibrate: scaled
    scaled <- sub[["pred"]] / sum( sub[["pred"]] )
    
    # Make variables for global, relative, and BIL scores
    glo  <- logit( sub[["pred"]] )
    rel  <- glo - max(glo)
    bil  <- ifelse( glo == max(glo), 1, 0 )
    
    # Assign values: scaled
    set( x     = out, 
         i     = which( df$tcp == i ), 
         j     = "scaled", 
         value = scaled )
    
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
  
  # Fit GAMs
  mod_gam  <- gam( causal ~ s(original),               data=data )
  
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
                                         gam    = NA,
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
    #   Get test set predictions and various recalibrations
    #---------------------------------------------------------------------------
    
    # Get model predictions (training and test), plus info for recalibration
    preds_tr <- recalibrate_preds( df=trn, model=l_mod, method=method, 
                                   bias_cols=bias_cols )
    preds_te <- recalibrate_preds( df=tst, model=l_mod, method=method, 
                                   bias_cols=bias_cols )
    
    # Train LASSO and GAM models to recalibrate model outputs
    cal_mod <- recalibration_model( data=preds_tr )
    
    # Apply GAM model
    preds_te$gam <- predict( object=cal_mod$gam, newdata=preds_te )
    preds_te$gam[ preds_te$gam < 0.001 ] <- 0.001
    preds_te$gam[ preds_te$gam > 0.999 ] <- 0.999
    
    # Apply LASSO model
    pred_cols <- c( "global", "relative" )
    preds_te$recal <- predict( object=cal_mod$lasso, s="lambda.min", type="response",
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
    out$preds$recal[  out$preds$trait == i ] <- preds_te$recal
    out$preds$gam[    out$preds$trait == i ] <- preds_te$gam
    out$preds$bil[    out$preds$trait == i ] <- preds_te$best
  }
  
  # Return
  return(out)
}


#-------------------------------------------------------------------------------
#   plot_loto_pr:    Plot precision-recall curve for LOTO
#-------------------------------------------------------------------------------

plot_loto_pr <- function( loto_obj, type="original", show_plot=TRUE, 
                          color="black", legend=FALSE ){
  
  # Extract the right curve
  if( type == "original" ){
    prc <- pr.curve( scores.class0 = loto_obj$preds$pred[  loto_obj$preds$causal ], 
                     scores.class1 = loto_obj$preds$pred[ !loto_obj$preds$causal ], 
                     curve = TRUE )
  }else if( type == "recalibrated" ){
    prc <- pr.curve( scores.class0 = loto_obj$preds$recal[  loto_obj$preds$causal ], 
                     scores.class1 = loto_obj$preds$recal[ !loto_obj$preds$causal ], 
                     curve = TRUE )
  }else if( type == "gam" ){
    prc <- pr.curve( scores.class0 = loto_obj$preds$gam[  loto_obj$preds$causal ], 
                     scores.class1 = loto_obj$preds$gam[ !loto_obj$preds$causal ], 
                     curve = TRUE )
  }else if( type == "scaled" ){
    prc <- pr.curve( scores.class0 = loto_obj$preds$scaled[  loto_obj$preds$causal ], 
                     scores.class1 = loto_obj$preds$scaled[ !loto_obj$preds$causal ], 
                     curve = TRUE )
  }
  dimnames(prc$curve)[[2]] <- c( "recall", "precision", "threshold" )
  
  # Plot
  if(show_plot){
    par( mar=c(5,5,4,1) )
    plot( prc, scale.color=color, las=1, legend=legend )
  }
  
  # Return AUPRC and 95% CI
  out <- aucpr.conf.int.expit( estimate = prc$auc.integral, 
                               num.pos  = sum(loto_obj$preds$causal) )
  return(out)
}


#-------------------------------------------------------------------------------
#   plot_logOR_relationship_smooth: Smoothed P(causal) v. a variable
#-------------------------------------------------------------------------------

plot_logOR_relationship_smooth <- function( data, varname, add_v0=TRUE ){
  
  # Format data
  prop_miss  <- mean( is.na( data[[varname]] ) )
  pcausal_na <- mean( data$causal[ is.na( data[[varname]] ) ] )
  data <- as.data.frame(data)
  data <- data.frame( x=data[[varname]], y=data$causal )
  data <- data[ complete.cases(data) , ]
  
  # Run model, get predictions
  pmode <- max( table(data$x) ) / NROW(data)
  if( pmode < 0.5 ){
    model <- loess( y ~ x, data=data )
  }else{
    span  <- (1 + pmode) / 2 
    model <- loess( y ~ x, data=data, span=span )
  }
  xrange <- range( data$x )
  xseq <- seq( from=xrange[1], to=xrange[2], length=100 )
  pred <- predict( model, newdata=data.frame(x=xseq), se=TRUE )
  
  # Format fitted values
  y    <- pred$fit
  ci   <- pred$se.fit * qt( 0.95/2 + 0.5, pred$df )
  ymin <- y - ci
  ymax <- y + ci
  loess.DF <- data.frame(x = xseq, y, ymin, ymax, se = pred$se.fit)
  
  # Restrict y values to 0.01 to 0.99
  loess.DF$y[    loess.DF$y    > 0.99 ] <- 0.99
  loess.DF$ymin[ loess.DF$ymin > 0.99 ] <- 0.99
  loess.DF$ymax[ loess.DF$ymax > 0.99 ] <- 0.99
  loess.DF$y[    loess.DF$y    < 0.01 ] <- 0.01
  loess.DF$ymin[ loess.DF$ymin < 0.01 ] <- 0.01
  loess.DF$ymax[ loess.DF$ymax < 0.01 ] <- 0.01
  
  # Plot
  xmin  <- min( c( 0, data$x ) )
  xmax  <- max(data$x)
  label <- paste0( "Proportion missing: ", round( x=prop_miss, digits=2 ) )
  p1 <- ggplot( data, aes(x=x) ) +
    geom_histogram() +
    ggtitle(varname)
  if(add_v0){
    p1 <- p1 + geom_vline( xintercept=0, linetype=2 )
  }
  p2 <- ggplot( loess.DF, aes( x=x, y=y) ) +
    geom_smooth( aes_auto(loess.DF), data=loess.DF, stat="identity", col="#70AD47" ) +
    scale_y_continuous( trans="logit",
                        breaks=c( 0.01, 0.05, seq(0.1,0.9,0.1), 0.95, 0.99 ) ) +
    geom_hline( yintercept=0.01 ) +
    geom_hline( yintercept=0.99 ) +
    geom_hline( yintercept=pcausal_na, linetype=2 ) +
    geom_text( x=xmin + (xmax-xmin)*0.05, y=logit(0.98), label=label, hjust=0 ) +
    theme_bw()
  if(add_v0){
    p2 <- p2 + geom_vline( xintercept=0, linetype=2 )
  }
  p1 + p2 + plot_layout( nrow=2, heights=c(1,3) )
}


#-------------------------------------------------------------------------------
#   plot_prob_relationship_smooth: Smoothed P(causal) v. a variable
#-------------------------------------------------------------------------------

plot_prob_relationship_smooth <- function( data, varname ){
  
  # Format data
  data <- as.data.frame(data)
  data <- data.frame( x=data[[varname]], y=data$causal )
  
  # Run model, get predictions
  pmode <- max( table(data$x) ) / NROW(data)
  if( pmode < 0.5 ){
    model <- loess( y ~ x, data=data, span=0.2 )
  }else{
    span  <- (1 + pmode) / 2 
    model <- loess( y ~ x, data=data, span=span )
  }
  xrange <- range( data$x )
  xseq <- seq( from=xrange[1], to=xrange[2], length=100 )
  pred <- predict( model, newdata=data.frame(x=xseq), se=TRUE )
  
  # Format fitted values
  y    <- pred$fit
  ci   <- pred$se.fit * qt( 0.95/2 + 0.5, pred$df )
  ymin <- y - ci
  ymax <- y + ci
  loess.DF <- data.frame(x = xseq, y, ymin, ymax, se = pred$se.fit)
  
  # Restrict y values to 0 to 1
  loess.DF$y[    loess.DF$y    > 1 ] <- 1
  loess.DF$ymin[ loess.DF$ymin > 1 ] <- 1
  loess.DF$ymax[ loess.DF$ymax > 1 ] <- 1
  loess.DF$y[    loess.DF$y    < 0 ] <- 0
  loess.DF$ymin[ loess.DF$ymin < 0 ] <- 0
  loess.DF$ymax[ loess.DF$ymax < 0 ] <- 0
  
  # Plot
  p1 <- ggplot( data, aes(x=x) ) +
    geom_histogram() +
    labs( x="", y="Count" ) +
    ggtitle(varname)
  p2 <- ggplot( loess.DF, aes( x=x, y=y) ) +
    geom_abline( slope=1, intercept=0, lty=2, col="grey", lwd=1 ) +
    geom_smooth( aes_auto(loess.DF), data=loess.DF, stat="identity", col="#70AD47" ) +
    scale_y_continuous( breaks=seq( 0, 1, 0.1 ) ) +
    labs( x="Model prediction", y="Truth" ) +
    theme_bw()
  p1 + p2 + plot_layout( nrow=2, heights=c(1,3) )
}


#-------------------------------------------------------------------------------
#   plot_logOR_relationship_model:  Modeled P(causal) v. a variable
#-------------------------------------------------------------------------------

plot_logOR_relationship_model <- function( data, model, focal_var, 
                                           qmin=0.1, qmax=0.9, bins=30,
                                           g_trans_fun=NULL ){
  
  # Format data
  feat_cols <- names(model$model)[-1]
  xycols <- c( "causal", feat_cols )
  data <- data[ , ..xycols ]
  data <- as.data.frame( lapply( X=data[ , ..xycols ], FUN=as.numeric ) )
  
  # Create new, evenly-spaced values for the focal global feature
  n      <- 100
  xrange <- quantile( data[[focal_var]], probs=c( qmin, qmax ) )
  xseq   <- seq( from=xrange[1], to=xrange[2], length=n )
  
  # Construct a dataset to make predictions in
  # Set all global features to their mean from the training dataset
  feat_means <- apply( data[,-1], 2, mean )
  mat <- matrix( rep( feat_means, n ), 
                 nrow=n, ncol=length(feat_means), byrow=TRUE )
  nx  <- as.data.frame(mat)
  names(nx) <- names(feat_means)
  
  # Set the focal global feature to the evenly-spaced values
  nx[[focal_var]] <- xseq
  
  # Get predictions
  pred <- predict( model, newdata=nx[ , feat_cols ], type="response" )
  nx$pred <- as.vector(pred)
  df <- data.frame( x=nx[[focal_var]], y=nx$pred )
  
  # If necessary, transform the scale of the global feature
  if( !is.null(g_trans_fun) ){
    df$x <- g_trans_fun(df$x)
  }
  
  # Clean up labels
  clean_label <- focal_var
  clean_label <- sub( x=clean_label, pattern="_glo",      replacement="" )
  clean_label <- sub( x=clean_label, pattern="pops",      replacement="PoPS" )
  clean_label <- sub( x=clean_label, pattern="magma",     replacement="MAGMA z-score" )
  clean_label <- sub( x=clean_label, pattern="dist_gene", replacement="Distance (bp)" )
  clean_label <- sub( x=clean_label, pattern="coding",    replacement="Coding PIP" )
  clean_label <- sub( x=clean_label, pattern="prior_n_genes_locus",
                      replacement="Number of local genes" )
  
  # Prepare the histogram
  hdf <- data.frame( x=data[[focal_var]] )
  hdf <- hdf[ hdf$x >= xrange[1] & hdf$x <= xrange[2] , , drop=FALSE ]
  if( !is.null(g_trans_fun) ){
    hdf$x <- g_trans_fun(hdf$x)
  }
  
  # Plot
  p1 <- ggplot( hdf, aes(x=x) ) +
    geom_histogram(bins=bins) +
    labs( y="Count", x="" ) +
    theme_bw() + 
    ggtitle(clean_label)
  if( focal_var == "coding_glo" ){
    p1 <- p1 + coord_cartesian( ylim=c(0,20) )
  }
  ymax <- 0.6
  p2 <- ggplot( df, aes( x=x, y=y ) ) +
    geom_line( aes( y=y ), linewidth=1 ) + 
    scale_y_continuous( trans="logit",
                        breaks=c( seq(0.01,0.05,0.01), seq(0.1,ymax,0.05) ) ) +
    geom_hline( yintercept=0.01, col="lightgrey" ) +
    geom_hline( yintercept=ymax, col="lightgrey" ) +
    labs( y="Predicted causal probability", x=clean_label ) +
    scale_color_manual( values=brewer.pal( n=6, name="Greens" )[5] ) +
    theme_bw()
  p1 + p2 + plot_layout( nrow=2, heights=c(1,3) )
}


#-------------------------------------------------------------------------------
#   benchmarking_auprc: Get AUPRC (and associated values) for CALDERA and L2G in
#                       either the Open Targets or the ExWAS benchmarking set
#-------------------------------------------------------------------------------

benchmarking_auprc <- function(bench){
  
  #-------------------------------------------------------------------------------
  #   Read in data
  #-------------------------------------------------------------------------------
  
  # Libraries
  library(data.table)
  library(glmnet)
  maindir <- "~/projects/causal_genes/"
  
  # Read in CALDERA and recalibration models
  if( bench == "l2g" ){
    bg_las    <- readRDS( file=file.path( maindir, "models", "caldera_model.rds" ) )
    cov_means <- readRDS( file=file.path( maindir, "covariate_means", "cov_means.rds" ) )
    
  }else if( bench == "exwas" ){
    bg_las    <- readRDS( file=file.path( maindir, "models", "caldera_model_ex.rds" ) )
    cov_means <- readRDS( file=file.path( maindir, "covariate_means", "cov_means_ex.rds" ) )
    
  }else if( bench == "3mt" ){
    bg_las    <- readRDS( file=file.path( maindir, "models", "caldera_model_tmt.rds" ) )
    cov_means <- readRDS( file=file.path( maindir, "covariate_means", "cov_means_tmt.rds" ) )
    
  }else{
    stop( "bench must be either 'l2g', 'exwas', or '3mt'" )
  }
  
  # Read in L2G benchmarking dataset
  if( bench == "l2g" ){
    lg <- fread( file.path( maindir, "benchmarking_datasets", "L2G_benchmark.tsv" ) )
    bench <- "Open Targets"
  }else if( bench == "exwas" ){
    lg <- fread( file.path( maindir, "benchmarking_datasets", "ExWAS_benchmark.tsv" ) )
    bench <- "ExWAS"
  }else if( bench == "3mt" ){
    lg <- fread( file.path( maindir, "benchmarking_datasets", "3MT_benchmark.tsv" ) )
    bench <- "3MT"
  }else{
    stop( "bench must be either 'l2g', 'exwas', or '3mt'" )
  }
  
  # Update column names
  names(lg)[ names(lg) == "ensg" ]                <- "ensgid"
  names(lg)[ names(lg) == "filename" ]            <- "CS_ID"
  names(lg)[ names(lg) == "TP" ]                  <- "causal"
  names(lg)[ names(lg) == "weighted_distance" ]   <- "distance"
  names(lg)[ names(lg) == "Overall L2G score" ]   <- "L2G"
  names(lg)[ names(lg) == "weighted_cS2G_score" ] <- "cS2G"
  names(lg)[ names(lg) == "FLAMES_scaled" ]       <- "FLAMES"
  names(lg)[ names(lg) == "FLAMES_scaled_with_path" ] <- "FLAMES"
  
  # Read in mean bias column values, merge
  for( i in names(cov_means) ){
    lg[[i]] <- cov_means[i]
  }
  
  
  #-------------------------------------------------------------------------------
  #   Make global features
  #-------------------------------------------------------------------------------
  
  # Subset to non-missing genes
  lg2 <- lg[ !is.na(lg$L2G) & !is.na(lg$FLAMES) , ]
  
  # Process features: global
  lg2$pops_glo      <- lg2$PoPS_Score
  lg2$dist_gene_glo <- -log10( lg2$distance + 1e3 )
  lg2$coding_glo    <- lg2$coding
  
  # Make a column for number of local genes
  n_genes <- tapply( X     = lg2$distance, 
                     INDEX = lg2$CS_ID, 
                     FUN   = length )
  lg2$ngenes_nearby <- n_genes[ match( lg2$CS_ID, names(n_genes) ) ]
  lg2$prior_n_genes_locus <- logit10( 1/lg2$ngenes_nearby )
  lg2$prior_n_genes_locus[ lg2$ngenes_nearby <= 1 ] <- logit10(0.75)
  
  
  #-------------------------------------------------------------------------------
  #   Get CALDERA predictions, recalibrate
  #-------------------------------------------------------------------------------
  
  # Get CALDERA predictions
  if( any( class(bg_las) == "glm" ) ){
    lg2$CALDERA_m <- predict( object=bg_las, newdata=lg2, type="response" )
  }else if( any( class(bg_las) == "cv.glmnet" ) ){
    feat_cols <- bg_las$glmnet$beta@Dimnames[[1]]
    lg2$CALDERA_m <- predict( bg_las, newx=as.matrix( lg2[ , ..feat_cols ] ), 
                             s="lambda.min", type="response" )
  }
  
  # Make sure probabilities sum to 1 in each locus
  lg2$CALDERA <- as.numeric(NA)
  for( i in unique(lg2$CS_ID) ){
    
    # Subset, scale to 1
    sub  <- lg2[ lg2$CS_ID == i , ]
    if( NROW(sub) == 0 )  next
    scaled <- sub[["CALDERA_m"]] / sum( sub[["CALDERA_m"]] )
    
    # Assign values: scaled
    set( x     = lg2, 
         i     = which( lg2$CS_ID == i ), 
         j     = "CALDERA", 
         value = scaled )
  }
  
  
  #-------------------------------------------------------------------------------
  #   Extract AUPRCs, make plots, return
  #-------------------------------------------------------------------------------
  
  # Remove trait-gene pairs that are missing L2G information
  lg3 <- lg2[ !is.na(lg2$L2G) & !is.na(lg2$FLAMES) & !is.na(lg2$CALDERA) , ]
  
  # Make PR curves, extract AUPRC
  cal_pr  <- plot_mod_pr( pred_col="CALDERA",   data=lg3, make_plot=FALSE )
  mul_pr  <- plot_mod_pr( pred_col="CALDERA_m", data=lg3, make_plot=FALSE )
  l2g_pr  <- plot_mod_pr( pred_col="L2G",       data=lg3, make_plot=FALSE )
  fla_pr  <- plot_mod_pr( pred_col="FLAMES",    data=lg3, make_plot=FALSE )
  s2g_pr  <- plot_mod_pr( pred_col="cS2G",      data=lg3, make_plot=FALSE )
  
  # Make a single PR curve with all methods
  methods <- c( "cS2G", "L2G", "FLAMES", "CALDERA" )
  plot_mod_pr( pred_col=methods, data=lg3 )
  
  # Bind AUPRCs together
  b_pr <- rbind( cal_pr, mul_pr, fla_pr, l2g_pr, s2g_pr )
  dimnames(b_pr)[[1]] <- c( "CALDERA", "CALDERA_m", "FLAMES", "L2G", "cS2G" )
  
  # Plots: CALDERA single
  p1 <- cal_plot_logistic( .data=lg3, truth=causal, estimate=CALDERA,
                           include_rug=TRUE, conf_level=0.95 ) +
    ggtitle( paste0( bench, " dataset, CALDERA single" ) ) +
    xlab("Predicted probability") + 
    ylab("Ground truth probability")
  
  # Plots: CALDERA multi
  p2 <- cal_plot_logistic( .data=lg3, truth=causal, estimate=CALDERA_m,
                           include_rug=TRUE, conf_level=0.95 ) +
    ggtitle( paste0( bench, " dataset, CALDERA multi" ) ) +
    xlab("Predicted probability") + 
    ylab("Ground truth probability")
  
  # Plots: FLAMES
  p3 <- cal_plot_logistic( .data=lg3, truth=causal, estimate=FLAMES,
                           include_rug=TRUE, conf_level=0.95 ) +
    ggtitle( paste0( bench, " dataset, FLAMES" ) ) +
    xlab("Predicted probability") + 
    ylab("Ground truth probability")
  
  # Plots: L2G
  p4 <- cal_plot_logistic( .data=lg3, truth=causal, estimate=L2G,
                           include_rug=TRUE, conf_level=0.95 ) +
    ggtitle( paste0( bench, " dataset, L2G" ) ) +
    xlab("Predicted probability") + 
    ylab("Ground truth probability")
  
  # Return
  out <- list( pr=b_pr, p1=p1, p2=p2, p3=p3, p4=p4, data=lg3 )
  return(out)
}


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Read in data
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Read in data
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
suppressPackageStartupMessages( library(metafor) )
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

# Split training and testing sets based on trait
# tr <- dt[ dt$trait %in% tr_traits , ]
tr <- dt
te <- dt[ dt$trait %in% te_traits , ]
NROW(tr) / NROW(dt)

# If a gene is causal for multiple traits, remove it from the larger trait
n_causal_per_trait <- table( tr$trait[ tr$causal ] ) #|> sort()
unique_ensg_c <- unique(tr$ensgid_c)
uniq_tcps0 <- list()
for( i in seq_along(unique_ensg_c) ){
  
  # Subset to the focal causal gene, find unique TCPs
  gcols <- c( "trait", "region", "cs_id" )
  sub  <- tr[ tr$ensgid_c == unique_ensg_c[i] & tr$causal , ..gcols ]
  
  # Save the TCP for the trait with the fewest causal genes
  smallest_trait <- n_causal_per_trait[sub$trait] |> head(1) |> names()
  uniq_tcps0[[i]] <- sub[ sub$trait == smallest_trait , ]
}
uniq_tcps <- as.data.frame( do.call( rbind, uniq_tcps0 ) )
tr <- tr %>% inner_join(uniq_tcps)

# Create a subset that removes traits that overlap with the ExWAS dataset
exwas_traits <- c( "Ca", "eBMD", "Hb", "HbA1c", "Height", "LDLC", "TBil" )
tr_ex <- tr[ !( tr$trait %in% exwas_traits ) , ]

# Create a subset that removes traits that overlap with the 3MT dataset
tmt_traits <- "IGF1"
tr_tmt <- tr[ !( tr$trait %in% tmt_traits ) , ]


#-------------------------------------------------------------------------------
#   Random predictor
#-------------------------------------------------------------------------------

# Run a regression using a random predictor
set.seed(1)
tr$rand <- rnorm( n=NROW(tr) )
te$rand <- rnorm( n=NROW(te) )


#-------------------------------------------------------------------------------
#   Define feature vectors
#-------------------------------------------------------------------------------

# Global features
glo_cols <- c( "causal", "pops_glo", 
               "dist_gene_glo", 
               #"dist_tss_glo",
               "coding_glo",
               "twas_glo",
               "corr_liu_glo", "corr_and_glo", "corr_uli_glo", 
               "pchic_jung_glo", "pchic_jav_glo", 
               "clpp_glo", "smr_glo", "abc_glo",
               "depict_z_glo", "netwas_score_glo", "netwas_bon_score_glo" )

# Features that may capture bias in how causal genes were selected
bias_cols <- c( "pLI_lt_0.1",
                "cds_bp_log10",
                "roadmap_bp_log10" )

# Full feature set
f_cols <- unique( c( glo_cols, "prior_n_genes_locus" ) )

# Basic feature set
pattern <- "^causal$|^pops_|^dist_gene_|^magma_|^coding_|^prior"
b_cols  <- grep( pattern=pattern, x=f_cols, value=TRUE )

# Full + GLC feature set
fg_cols <- c( f_cols, bias_cols )

# Basic + GLC feature set
bg_cols <- c( b_cols, bias_cols )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Full model
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   LOTO
#-------------------------------------------------------------------------------

# Raw feature correlations
corrplot( cor( tr[ , ..f_cols ][,-1] ), order="hclust" )

# Run LOTO
fl_glm <- loto( data=tr, method="glm",    feat_cols=f_cols )
# fl_las <- loto( data=tr, method="lasso2", feat_cols=f_cols )
# fl_xgb <- loto( data=tr, method="xgb",    feat_cols=f_cols, maxit=100 )
# saveRDS( object=fl_xgb, file=file.path( maindir, "loto/loto_full_xgb.rds" ) )
fl_xgb  <- readRDS( file=file.path( maindir, "loto/loto_full_xgb.rds" ) )


#-------------------------------------------------------------------------------
#   Train on the full dataset (for external benchmarking)
#-------------------------------------------------------------------------------

# Make GLM
f_form <- make_formula( lhs="causal", rhs=f_cols[-1] )
f_glm <- glm( formula=f_form, data=tr, family="binomial" )

# Make LASSO model
f_las <- cv.glmnet( x=as.matrix( tr[ , ..f_cols ] )[,-1], 
                    foldid=as.integer( as.factor(tr$trait) ),
                    y=tr[["causal"]], family="binomial" )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Full model + GLCs
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   LOTO
#-------------------------------------------------------------------------------

# Raw feature correlations
corrplot( cor( tr[ , ..fg_cols ][,-1] ), order="hclust" )

# Run LOTO
fgl_glm <- loto( data=tr, method="glm",    feat_cols=fg_cols )
# fgl_las <- loto( data=tr, method="lasso2", feat_cols=fg_cols )


#-------------------------------------------------------------------------------
#   Train on the full dataset (for external benchmarking)
#-------------------------------------------------------------------------------

# Make GLM
fg_form <- make_formula( lhs="causal", rhs=fg_cols[-1] )
fg_glm <- glm( formula=fg_form, data=tr, family="binomial" )

# Make LASSO model
fg_las <- cv.glmnet( x=as.matrix( tr[ , ..fg_cols ] )[,-1], 
                     foldid=as.integer( as.factor(tr$trait) ),
                     y=tr[["causal"]], family="binomial" )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Basic model
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   LOTO
#-------------------------------------------------------------------------------

# Raw feature correlations
corrplot( cor( tr[ , ..b_cols ][,-1] ), order="hclust" )

# Run LOTO
bl_glm <- loto( data=tr, method="glm",    feat_cols=b_cols )
# bl_las <- loto( data=tr, method="lasso2", feat_cols=b_cols )


#-------------------------------------------------------------------------------
#   Train on the full dataset (for external benchmarking)
#-------------------------------------------------------------------------------

# Make GLM
b_form <- make_formula( lhs="causal", rhs=b_cols[-1] )
b_glm <- glm( formula=b_form, data=tr, family="binomial" )

# Removing ExWAS traits
b_glm_ex <- glm( formula=b_form, data=tr_ex, family="binomial" )

# Removing 3MT traits
b_glm_tmt <- glm( formula=b_form, data=tr_tmt, family="binomial" )

# # Make LASSO model
# b_las <- cv.glmnet( x=as.matrix( tr[ , ..b_cols ] )[,-1], 
#                      foldid=as.integer( as.factor(tr$trait) ),
#                      y=tr[["causal"]], family="binomial" )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Basic model + GLCs
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   LOTO
#-------------------------------------------------------------------------------

# Raw feature correlations
corrplot( cor( tr[ , ..bg_cols ][,-1] ), order="hclust" )

# Run LOTO
bgl_glm <- loto( data=tr, method="glm",    feat_cols=bg_cols, bias_cols=bias_cols )
# bgl_las <- loto( data=tr, method="lasso2", feat_cols=bg_cols, bias_cols=bias_cols )


#-------------------------------------------------------------------------------
#   Train on the full dataset (for external benchmarking)
#-------------------------------------------------------------------------------

# Make GLM
bg_form <- make_formula( lhs="causal", rhs=bg_cols[-1] )
bg_glm <- glm( formula=bg_form, data=tr, family="binomial" )

# Removing ExWAS traits
bg_glm_ex <- glm( formula=bg_form, data=tr_ex, family="binomial" )

# Removing 3MT traits
bg_glm_tmt <- glm( formula=bg_form, data=tr_tmt, family="binomial" )

# Make LASSO model
# bg_las <- cv.glmnet( x=as.matrix( tr[ , ..bg_cols ] )[,-1], 
#                      foldid=as.integer( as.factor(tr$trait) ),
#                      y=tr[["causal"]], family="binomial" )


#-------------------------------------------------------------------------------
#   Save CALDERA models to file
#-------------------------------------------------------------------------------

# Save covariate mean values
cov_means     <- colMeans( x=tr[ , ..bias_cols ] )
cov_means_ex  <- colMeans( x=tr_ex[ , ..bias_cols ] )
cov_means_tmt <- colMeans( x=tr_tmt[ , ..bias_cols ] )
# saveRDS( object=cov_means,     file=file.path( maindir, "covariate_means", "cov_means.rds" ) )
# saveRDS( object=cov_means_ex,  file=file.path( maindir, "covariate_means", "cov_means_ex.rds" ) )
# saveRDS( object=cov_means_tmt, file=file.path( maindir, "covariate_means", "cov_means_tmt.rds" ) )

# Save models: basic + GLCs
# saveRDS( object=bg_glm,     file=file.path( maindir, "models", "caldera_model.rds" ) )
# saveRDS( object=bg_glm_ex,  file=file.path( maindir, "models", "caldera_model_ex.rds" ) )
# saveRDS( object=bg_glm_tmt, file=file.path( maindir, "models", "caldera_model_tmt.rds" ) )

# Save models: basic
# saveRDS( object=b_glm,     file=file.path( maindir, "models", "caldera_model_no_covs.rds" ) )
# saveRDS( object=b_glm_ex,  file=file.path( maindir, "models", "caldera_model_no_covs_ex.rds" ) )
# saveRDS( object=b_glm_tmt, file=file.path( maindir, "models", "caldera_model_no_covs_tmt.rds" ) )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Figures
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Figure 1: AUPRC barplots, scatterplot with/without GLCs
#-------------------------------------------------------------------------------

# Get AUPRCs
fl_xgb_pr        <- plot_loto_pr( loto_obj=fl_xgb,   color=viridis(11), legend=TRUE )
fl_glm_pr        <- plot_loto_pr( loto_obj=fl_glm,   color=viridis(11), legend=TRUE )
fgl_glm_pr       <- plot_loto_pr( loto_obj=fgl_glm,  color=viridis(11), legend=TRUE )
bgl_glm_pr       <- plot_loto_pr( loto_obj=bgl_glm,  color=viridis(11), legend=TRUE )
bgl_glm_sca_pr   <- plot_loto_pr( loto_obj=bgl_glm,  color=viridis(11), legend=TRUE,
                                  type="scaled")

# Get the random AUPRC
rand_prc <- pr.curve( scores.class0 = tr$rand[  tr$causal ], 
                 scores.class1 = tr$rand[ !tr$causal ], 
                 curve = TRUE )
rand_pr  <- aucpr.conf.int.expit( estimate = rand_prc$auc.integral, 
                                  num.pos  = sum(tr$causal) )

# Figure 1: set up
fig_dir   <- file.path( maindir, "figures" )
fig1_file <- file.path( fig_dir, "Fig1.tif" )
dir.create( path=fig_dir, showWarnings=FALSE, recursive=TRUE )
tiff( filename=fig1_file, width=300*4, height=300*4, res=75*4 )
par( mfrow=c(1,1) )
par( mar=c(9,5,1,1) )

# Plot
all_pr <- rbind( rand_pr, fl_xgb_pr, fl_glm_pr, fgl_glm_pr,
                 bgl_glm_pr, bgl_glm_sca_pr)
dimnames(all_pr)[[1]] <- c( "Random", "XGBoost (full)", "GLM (full)",
                            "GLM (full + GLCs)", "GLM (basic + GLCs)",
                            "CALDERA")
bp1 <- barplot( height=all_pr[,"auprc"], las=2, 
                ylab="AUPRC (±95% CI)", ylim=c( 0, 0.75 ),
                col=brewer.pal( n=9, name="Greens" )[1:7] )
for( i in seq_along(bp1) ){
  lines( x = c( bp1[i],            bp1[i] ), lwd=1.5,
         y = c( all_pr[ i, "lo" ], all_pr[ i, "hi" ] ) )
}

# Figure 1: wrap up
par( mfrow=c(1,1) )
dev.off()

# Right panel
# plot( logit(bgl_glm$preds$scaled), logit(bl_glm$preds$scaled), las=1 )
plot( bgl_glm$preds$scaled, bl_glm$preds$scaled, las=1 )
abline( a=0, b=1, lty=2, lwd=2, col="grey" )
hist( logit(bgl_glm$preds$scaled) - logit(bl_glm$preds$scaled), 
      breaks=seq(-0.35,0.35,0.05) )
abline( v=0, lwd=2 )
abline( v=mean( logit(bgl_glm$preds$scaled) - logit(bl_glm$preds$scaled), na.rm=TRUE ),
        lty=2, lwd=2 )
t.test( x=logit(bgl_glm$preds$scaled) - logit(bl_glm$preds$scaled) )


#-------------------------------------------------------------------------------
#   Figure 2: Benchmarking AUPRC
#-------------------------------------------------------------------------------

# Plot
fig_dir   <- file.path( maindir, "figures" )
fig2_file <- file.path( fig_dir, "Fig2.tif" )
tiff( filename=fig2_file, width=750*4, height=250*4, res=75*4 )
par( mfrow=c(1,3) )
par( mar=c(7,5,1,1) )
lg_pr <- benchmarking_auprc( bench="l2g" )
ex_pr <- benchmarking_auprc( bench="exwas" )
mt_pr <- benchmarking_auprc( bench="3mt" )
par( mfrow=c(1,1) )
dev.off()

# FEMA
library(metafor)
bmark <- rbind( lg_pr$pr, ex_pr$pr, mt_pr$pr )
bmark$bench <- c( rep( "OT", 5 ),
                  rep( "Ex", 5 ),
                  rep( "MT", 5 ) )
bmark$method <- rep( row.names(lg_pr$pr), 3 )
ma0 <- list()
for( i in unique(bmark$method) ){
  sub <- bmark[ bmark$method == i , ]
  fit_fe <- rma.uni( yi=sub$auprc_logit, sei=sub$se_logit, method="REML" )
  ma0[[i]] <- data.table( method = i,
                          effect = fit_fe$b[1],
                          se     = fit_fe$se,
                          ci_lb  = fit_fe$ci.lb,
                          ci_ub  = fit_fe$ci.ub,
                          z      = fit_fe$zval,
                          p      = fit_fe$pval )
}
ma <- do.call( rbind, ma0 )
ma$auprc <- logistic(ma$effect)
ma$lo    <- logistic(ma$ci_lb)
ma$hi    <- logistic(ma$ci_ub)

# Two-sample t-tests: CALDERA vs. FLAMES
two_sample_z_test( mean1 = ma$effect[1], 
                   mean2 = ma$effect[3], 
                   se1   = ma$se[1], 
                   se2   = ma$se[3] )

# Two-sample t-tests: CALDERA vs. L2G
two_sample_z_test( mean1 = ma$effect[1], 
                   mean2 = ma$effect[4], 
                   se1   = ma$se[1], 
                   se2   = ma$se[4] )

# Two-sample t-tests: CALDERA vs. cS2G
two_sample_z_test( mean1 = ma$effect[1], 
                   mean2 = ma$effect[5], 
                   se1   = ma$se[1], 
                   se2   = ma$se[5] )


#-------------------------------------------------------------------------------
#   Figure 3: Calibration plots
#-------------------------------------------------------------------------------

# Figure 3
fig3_file <- file.path( fig_dir, "Fig3.tif" )
tiff( filename=fig3_file, width=300*4, height=300*4, res=75*4 )
cal_plot_logistic( .data=bgl_glm$preds, truth=causal, estimate=scaled, 
                   include_rug=FALSE, conf_level=0.95 ) + 
  ggtitle("CALDERA") +
  xlab("CALDERA predicted probability") + 
  ylab("Ground truth probability")
dev.off()

# LOESS
plot_prob_relationship_smooth( data=bgl_glm$preds, varname="pred" )
plot_prob_relationship_smooth( data=bgl_glm$preds, varname="scaled" )
plot_prob_relationship_smooth( data=bgl_glm$preds, varname="gam" )
plot_prob_relationship_smooth( data=bgl_glm$preds, varname="recal" )

# How do predictions compare for GLM and GLMM?
plot( bgl_glm$preds$scaled, bgl_glmm$preds$scaled, las=1, col="steelblue", lwd=1.5,
      xlab="GLM", ylab="GLMM" )
abline( a=0, b=1, lty=2, lwd=2 )
summary(  abs( bgl_glm$preds$scaled - bgl_glmm$preds$scaled ) )
quantile( abs( bgl_glm$preds$scaled - bgl_glmm$preds$scaled ), 
          probs=seq( 0.9, 1, 0.01 ) )

# Plot the largest 4 traits
ncausal <- sort( tapply(X=bgl_las_uniq$preds$causal, 
                        INDEX=bgl_las_uniq$preds$trait, 
                        FUN=sum ), decreasing=TRUE )
p1 <- cal_plot_logistic( .data=bgl_las_uniq$preds[ 
  bgl_las_uniq$preds$trait == names(ncausal)[1] ], 
                         truth=causal, estimate=recal, conf_level=0.95 ) +
  ggtitle( names(ncausal)[1] ) +
  xlab("CALDERA predicted probability") + 
  ylab("Ground truth probability")
p2 <- cal_plot_logistic( .data=bgl_las_uniq$preds[ 
  bgl_las_uniq$preds$trait == names(ncausal)[2] ], 
                         truth=causal, estimate=recal, conf_level=0.95 ) +
  ggtitle( names(ncausal)[2] ) +
  xlab("CALDERA predicted probability") + 
  ylab("Ground truth probability")
p3 <- cal_plot_logistic( .data=bgl_las_uniq$preds[ 
  bgl_las_uniq$preds$trait == names(ncausal)[3] ], 
                         truth=causal, estimate=recal, conf_level=0.95 ) +
  ggtitle( names(ncausal)[3] ) +
  xlab("CALDERA predicted probability") + 
  ylab("Ground truth probability")
p4 <- cal_plot_logistic( .data=bgl_las_uniq$preds[ 
  bgl_las_uniq$preds$trait == names(ncausal)[4] ], 
                         truth=causal, estimate=recal, conf_level=0.95 ) +
  ggtitle( names(ncausal)[4] ) +
  xlab("CALDERA predicted probability") + 
  ylab("Ground truth probability")
grid.arrange( p1, p2, p3, p4, ncol=2 )


#-------------------------------------------------------------------------------
#   Figure S1: Prediction comparisons across the evolution of CALDERA
#-------------------------------------------------------------------------------

# Get set up
fig_s1_file <- file.path( fig_dir, "FigS1.tif" )
tiff( filename=fig_s1_file, width=200*4, height=600*4, res=75*4 )
par( mar=c(5,5,1,1) )
par( mfrow=c(3,1) )

# How do predictions compare for XGBoost and LASSO using the full feature set?
plot( fl_xgb$preds$pred, fl_glm$preds$pred, las=1, lwd=1.5,
      col=ifelse( fl_glm$pred$causal,
                  brewer.pal( n=6, name="Greens" )[5], 
                  brewer.pal( n=6, name="Greens" )[5] ),
      xlab="XGBoost (full)", ylab="GLM (full)" )
abline( a=0, b=1, lty=2, lwd=2 )

# How do predictions compare for GLM full v. basic?
plot( fgl_glm$preds$pred, bgl_glm$preds$pred, las=1, lwd=1.5,
      col=ifelse( fl_glm$pred$causal,
                  brewer.pal( n=6, name="Greens" )[5], 
                  brewer.pal( n=6, name="Greens" )[5] ),
      xlab="GLM (full + GLCs)", ylab="GLM (basic + GLCs)" )
abline( a=0, b=1, lty=2, lwd=2 )

# How do predictions compare for GLM basic v. CALDERA?
plot( bgl_glm$preds$pred, bgl_glm$preds$recal, las=1, lwd=1.5,
      col=ifelse( fl_glm$pred$causal,
                  brewer.pal( n=6, name="Greens" )[5], 
                  brewer.pal( n=6, name="Greens" )[5] ),
      xlab="GLM (basic + GLCs)", ylab="CALDERA" )
abline( a=0, b=1, lty=2, lwd=2 )
dev.off()

# Large LASSO (full) predictions that are small XGBoost (full) predictions
idx <- fl_xgb$preds$pred < 0.7 & fl_glm$preds$pred > 0.9
tr[idx,]


#-------------------------------------------------------------------------------
#   Figure S2: P(causal) v. input feature values
#-------------------------------------------------------------------------------

# POPS
bins <- 30
fig_s2a_file <- file.path( fig_dir, "FigS2a.tif" )
tiff( filename=fig_s2a_file, width=300*4, height=400*4, res=75*4 )
plot_logOR_relationship_model( data=tr, model=bg_glm, focal_var="pops_glo", 
                               qmin=0.05, qmax=0.95, bins=bins, 
                               g_trans_fun=NULL )
dev.off()

# Distance
fig_s2b_file <- file.path( fig_dir, "FigS2b.tif" )
tiff( filename=fig_s2b_file, width=300*4, height=400*4, res=75*4 )
plot_logOR_relationship_model( data=tr, model=bg_glm, focal_var="dist_gene_glo", 
                               qmin=0.05, qmax=0.95, bins=bins,
                               g_trans_fun=function(x) 10^-x - 1e3 )
dev.off()

# Coding PIP
fig_s2c_file <- file.path( fig_dir, "FigS2c.tif" )
tiff( filename=fig_s2c_file, width=300*4, height=400*4, res=75*4 )
plot_logOR_relationship_model( data=tr, model=bg_glm, focal_var="coding_glo",
                               qmin=0, qmax=1, bins=bins, 
                               g_trans_fun=NULL )
dev.off()

# Number of genes
fig_s2d_file <- file.path( fig_dir, "FigS2d.tif" )
tiff( filename=fig_s2d_file, width=300*4, height=400*4, res=75*4 )
plot_logOR_relationship_model( data=tr, model=bg_glm, 
                               focal_var="prior_n_genes_locus",
                               qmin=0.05, qmax=0.95, bins=bins,
                               g_trans_fun=function(x) 1 / logistic10(x))
dev.off()

# Number of genes
df1 <- data.frame( x=bgl_glm$preds$scaled[  bgl_glm$preds$causal ] )
df2 <- data.frame( x=bgl_glm$preds$scaled[ !bgl_glm$preds$causal ] )
p1 <- ggplot( df1, aes(x=x) ) +
  geom_histogram(bins=30) +
  xlim( 0, 1 ) +
  labs( y="Count", x="CALDERA score" ) +
  theme_bw() +
  ggtitle("Causal genes")
p2 <- ggplot( df2, aes(x=x) ) +
  geom_histogram(bins=30) +
  xlim( 0, 1 ) +
  labs( y="Count", x="CALDERA score" ) +
  theme_bw() +
  ggtitle("Non-causal genes")
fig_s2e_file <- file.path( fig_dir, "FigS2e.tif" )
tiff( filename=fig_s2e_file, width=300*4, height=600*4, res=75*4 )
p1 + p2 + plot_layout( nrow=2, heights=c(1,1) )
dev.off()


#-------------------------------------------------------------------------------
#   Figure S3: AUPRC barplots
#-------------------------------------------------------------------------------

# Set up
fig_s3_file <- file.path( fig_dir, "FigS3.tif" )
tiff( filename=fig_s3_file, width=600*4, height=200*4, res=75*4 )
par( mfrow=c(1,3) )
par( mar=c(5,5,2,1) )

# Open Targets benchmarking set
bar_cols <- c( brewer.pal( n=6, name="Greens" )[5], 
               brewer.pal( n=6, name="Reds" )[4], 
               brewer.pal( n=6, name="Blues" )[4], 
               brewer.pal( n=6, name="Oranges"  )[3] )
data <- lg_pr$pr[ row.names(lg_pr$pr) != "CALDERA_m" , ]
bp <- barplot( height=data$auprc, las=2, names.arg=row.names(data),
               main="Open Targets", ylab="AUPRC (±95% CI)", ylim=c( 0, 1 ),
               col=bar_cols )
for( i in seq_along(bp) ){
  lines( x = c( bp[i],           bp[i] ), lwd=1.5,
         y = c( data[ i, "lo" ], data[ i, "hi" ] ) )
}

# ExWAS benchmarking set
data <- ex_pr$pr[ row.names(ex_pr$pr) != "CALDERA_m" , ]
bp <- barplot( height=data$auprc, las=2, names.arg=row.names(data),
               main="ExWAS", ylab="AUPRC (±95% CI)", ylim=c( 0, 1 ),
               col=bar_cols )
for( i in seq_along(bp) ){
  lines( x = c( bp[i],           bp[i] ), lwd=1.5,
         y = c( data[ i, "lo" ], data[ i, "hi" ] ) )
}

# Three molecular traits benchmarking set
data <- mt_pr$pr[ row.names(mt_pr$pr) != "CALDERA_m" , ]
bp <- barplot( height=data$auprc, las=2, names.arg=row.names(data),
               main="Molecular traits", ylab="AUPRC (±95% CI)", ylim=c( 0, 1 ),
               col=bar_cols )
for( i in seq_along(bp) ){
  lines( x = c( bp[i],           bp[i] ), lwd=1.5,
         y = c( data[ i, "lo" ], data[ i, "hi" ] ) )
}

# Wrap up
par( mfrow=c(1,1) )
dev.off()

# Two-sample z-test: Open Targets
two_sample_z_test( mean1 = lg_pr$pr[ "CALDERA", "auprc_logit" ], 
                   mean2 = lg_pr$pr[ "L2G",       "auprc_logit" ], 
                   se1   = lg_pr$pr[ "CALDERA", "se_logit" ], 
                   se2   = lg_pr$pr[ "L2G",       "se_logit" ] )
two_sample_z_test( mean1 = lg_pr$pr[ "CALDERA", "auprc_logit" ], 
                   mean2 = lg_pr$pr[ "FLAMES",    "auprc_logit" ], 
                   se1   = lg_pr$pr[ "CALDERA", "se_logit" ], 
                   se2   = lg_pr$pr[ "FLAMES",    "se_logit" ] )


#-------------------------------------------------------------------------------
#   Figure S4: Calibration of CALDERA and L2G in external gold standards
#-------------------------------------------------------------------------------

fig_s4_file <- file.path( fig_dir, "FigS4.tif" )
tiff( filename=fig_s4_file, width=900*4, height=300*4, res=75*4 )
grid.arrange( lg_pr$p1, ex_pr$p1, mt_pr$p1, ncol=3 )
dev.off()


#-------------------------------------------------------------------------------
#   Figure S5: distance between eQTLs/GWAS hits and their causal genes
#-------------------------------------------------------------------------------

# See file: e_distance_between_gwas_hits_and_genes_from_mostafavi2023.R


#-------------------------------------------------------------------------------
#   Figure S6: features vs. P(causal) before and after transformation
#-------------------------------------------------------------------------------

# Distance body, before
fig_s6a_file <- file.path( fig_dir, "FigS6a.tif" )
tiff( filename=fig_s6a_file, width=300*4, height=400*4, res=75*4 )
plot_logOR_relationship_smooth( data=tr, varname="distance_genebody" )
dev.off()

# Distance body, after
fig_s6b_file <- file.path( fig_dir, "FigS6b.tif" )
tiff( filename=fig_s6b_file, width=300*4, height=400*4, res=75*4 )
plot_logOR_relationship_smooth( data=tr, varname="dist_gene_glo", add_v0=FALSE )
dev.off()

# Coding sequence length, before
fig_s6c_file <- file.path( fig_dir, "FigS6c.tif" )
tiff( filename=fig_s6c_file, width=300*4, height=400*4, res=75*4 )
plot_logOR_relationship_smooth( data=tr, varname="CDS_length", add_v0=FALSE )
dev.off()

# Coding sequence length, after
fig_s6d_file <- file.path( fig_dir, "FigS6d.tif" )
tiff( filename=fig_s6d_file, width=300*4, height=400*4, res=75*4 )
plot_logOR_relationship_smooth( data=tr, varname="cds_bp_log10", add_v0=FALSE )
dev.off()

# Roadmap enhancer length, before
fig_s6e_file <- file.path( fig_dir, "FigS6e.tif" )
tiff( filename=fig_s6e_file, width=300*4, height=400*4, res=75*4 )
plot_logOR_relationship_smooth( data=tr, varname="Roadmap_length_per_type", add_v0=FALSE )
dev.off()

# Roadmap enhancer length, after
fig_s6f_file <- file.path( fig_dir, "FigS6f.tif" )
tiff( filename=fig_s6f_file, width=300*4, height=400*4, res=75*4 )
plot_logOR_relationship_smooth( data=tr, varname="roadmap_bp_log10", add_v0=FALSE )
dev.off()

# TSS
plot_logOR_relationship_smooth( data=tr, varname="distance_tss" )
plot_logOR_relationship_smooth( data=tr, varname="dist_tss_glo", add_v0=FALSE )

# Plot relationship: TWAS
plot_logOR_relationship_smooth( data=tr, varname="twas_abs_z" )

# Plot relationship: CLPP
plot_logOR_relationship_smooth( data=tr, varname="clpp_prob" )

# Plot relationship: ABC
plot_logOR_relationship_smooth( data=tr, varname="abc_score" )

# Plot relationship: SMR
plot_logOR_relationship_smooth( data=tr, varname="smr_z" )

# Plot relationship: PCHi-C Jung
plot_logOR_relationship_smooth( data=tr, varname="pchic_jung_score" )

# Plot relationship: PCHi-C Javierre
plot_logOR_relationship_smooth( data=tr, varname="pchic_javierre_score" )

# Plot relationship: E-P Liu
plot_logOR_relationship_smooth( data=tr, varname="corr_liu_score" )

# Plot relationship: E-P Andersson
plot_logOR_relationship_smooth( data=tr, varname="corr_andersson_score" )

# Plot relationship: E-P Ulirsch
plot_logOR_relationship_smooth( data=tr, varname="corr_ulirsch_score" )

# DEPICT
plot_logOR_relationship_smooth( data=tr, varname="depict_z" )

# NetWAS
plot_logOR_relationship_smooth( data=tr, varname="no_loco_nc_netwas_score" )

# NetWAS Bonferroni
plot_logOR_relationship_smooth( data=tr, varname="no_loco_nc_netwas_bon_score" )

# pLI
plot_logOR_relationship_smooth( data=tr, varname="pLI" )


#-------------------------------------------------------------------------------
#   Figure SX: Transformed and untransformed variable distributions
#-------------------------------------------------------------------------------

# Plot transformed and untransformed variables, with and without imputed values
plot_trans <- function( var_u, var_t, flip=FALSE ){
  
  # Set up
  par( mfrow = c(2,1) )
  par( mar=c(4,4,4,1) )
  col <- brewer.pal( n=6, name="Greens" )[5]
  
  # Untransformed
  hist( tr[[var_u]][ !is.na( tr[[var_u]] ) ], 
        main="Untransformed", xlab="", las=1, col=col )
  
  # Transformed
  hist( tr[[var_t]][ !is.na( tr[[var_u]] ) ], 
        main="Transformed", xlab="", las=1,col=col )
  par( mfrow = c(1,1) )
}

# Plot
plot_trans( var_u="corr_liu_score", var_t="corr_liu_glo" )
plot_trans( var_u="clpp_prob", var_t="clpp_glo" )
plot_trans( var_u="abc_score", var_t="abc_glo" )
plot_trans( var_u="smr_p", var_t="smr_glo", flip=TRUE )
plot_trans( var_u="distance_genebody", var_t="dist_gene_glo", flip=TRUE )
plot_trans( var_u="distance_tss", var_t="dist_tss_glo", flip=TRUE )

hist(tr$pops_score)
hist(tr$coding_prob)
hist(tr$corr_andersson_score)
hist(tr$corr_ulirsch_score)
hist(tr$pchic_jung_score)
hist(tr$pchic_javierre_score)
hist(tr$depict_z_glo)
hist(tr$netwas_score_glo)
hist(tr$netwas_bon_score_glo)


#-------------------------------------------------------------------------------
#   Figure SX?: CALDERA precision-recall curve versus NG+TP
#-------------------------------------------------------------------------------

# Precision and recall: POPS + nearest gene
pnn_p <- sum( tr$causal & tr$pops_and_nearest ) / sum( tr$pops_and_nearest ) #precision
pnn_r <- sum( tr$causal & tr$pops_and_nearest ) / sum( tr$causal )           #recall

# For various P(causal) thresholds, what is precision and recall for the BIL?
cal_pr0 <- list()
threshs <- c( 0, seq( 0.5, 0.9, 0.1 ) )
for( i in seq_along(threshs) ){
  tp           <- sum( bgl_las$preds$causal & 
                         bgl_las$preds$bil == 1 & 
                         bgl_las$preds$recal >= threshs[i] )
  p_denom      <- sum( bgl_las$preds$bil == 1 & 
                         bgl_las$preds$recal >= threshs[i] )
  r_denom      <- sum( bgl_las$preds$causal )
  cal_pr0[[i]] <- c( threshold = threshs[i],
                     precision = tp/p_denom, 
                     recall    = tp/r_denom )
}
cal_pr     <- as.data.frame( do.call( rbind, cal_pr0 ) )
cal_pr$col <- brewer.pal( n=NROW(cal_pr), name="Greens" )

# Figure S5
fig_s5_file <- file.path( fig_dir, "figure_s5.tif" )
tiff( filename=fig_s5_file, width=300*4, height=300*4, res=75*4 )
plot_loto_pr(bgl_las, type="recalibrated" )
points( x=cal_pr$recall, y=cal_pr$precision, pch=19, cex=1.5, col=cal_pr$col )
points( x=pnn_r,         y=pnn_p,            pch=19, cex=1.5, 
        col=brewer.pal( n=4, name="Blues" )[3] )
# better_pr( p=pnn_p, r=pnn_r, t=0.75 )
dev.off()


#-------------------------------------------------------------------------------
#   Table S2
#-------------------------------------------------------------------------------

# Read in Weeks et al. data, subset to our traits and columns
it0 <- fread("~/projects/causal_genes/weeks2023_table_s1.csv")
it  <- it0[ it0$Trait %in% unique(tr$trait) , 
            c( "Trait", "Description", "N individuals", "N cases", "N controls" ) ]

# Add number of causal/non-causal genes
it$`N causal genes` <- tapply( X=tr$causal, INDEX=tr$trait, FUN=sum )
it$`N non-causal genes` <- tapply( X=!tr$causal, INDEX=tr$trait, FUN=sum )

# Write
tab_dir <- file.path( maindir, "tables")
tab_s1_file <- file.path( tab_dir, "table_s1.tsv" )
dir.create( path=tab_dir, showWarnings=FALSE, recursive=TRUE )
fwrite( x=it, file=tab_s1_file, sep="\t" )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Statistics reported in the manuscript
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   AUPRC after removing PoPS or substituting for DEPICT or NetWAS
#-------------------------------------------------------------------------------

# Patterns
pattern1 <- "^causal$|^dist_gene_|^magma_|^coding_|^prior"
pattern2 <- "^causal$|^depict_|^dist_gene_|^magma_|^coding_|^prior"
pattern3 <- "^causal$|^netwas_|^dist_gene_|^magma_|^coding_|^prior"

# Basic feature sets
b_cols1  <- grep( pattern=pattern1, x=f_cols, value=TRUE )
b_cols2  <- grep( pattern=pattern2, x=f_cols, value=TRUE )
b_cols3  <- grep( pattern=pattern3, x=f_cols, value=TRUE )

# Basic + GLC feature sets
bg_cols1 <- c( b_cols1, bias_cols, glc_cols )
bg_cols2 <- c( b_cols2, bias_cols, glc_cols )
bg_cols3 <- c( b_cols3, bias_cols, glc_cols )

# Models
bgl_glm1 <- loto( data=tr, method="glm", feat_cols=bg_cols1, bias_cols=bias_cols )
bgl_glm2 <- loto( data=tr, method="glm", feat_cols=bg_cols2, bias_cols=bias_cols )
bgl_glm3 <- loto( data=tr, method="glm", feat_cols=bg_cols3, bias_cols=bias_cols )

# AUPRCs
pr1 <- plot_loto_pr( loto_obj=bgl_glm1,  color=viridis(11), legend=TRUE, type="scaled")
pr2 <- plot_loto_pr( loto_obj=bgl_glm2,  color=viridis(11), legend=TRUE, type="scaled")
pr3 <- plot_loto_pr( loto_obj=bgl_glm3,  color=viridis(11), legend=TRUE, type="scaled")
pr1; pr2; pr3


#-------------------------------------------------------------------------------
#   Fisher's exact test of causal gene status versus pLI < 10%
#-------------------------------------------------------------------------------

tab <- table( tr$causal, tr$pLI_lt_0.1 ) 
fisher.test(tab)


#-------------------------------------------------------------------------------
#   AUPRC of CALDERA after removing covariates
#-------------------------------------------------------------------------------

bl_glm_sca_pr <- plot_loto_pr( loto_obj=bl_glm, color=viridis(11), legend=TRUE,
                               type="scaled")
bgl_glm_sca_pr; bl_glm_sca_pr


#-------------------------------------------------------------------------------
#   Mean change in P(causal) after GLC correction (original P(causal) > 20%)
#-------------------------------------------------------------------------------

summary( ( bl_las$preds$pred - bgl_las$preds$pred )[ bl_las$preds$pred > 0.2 ] )


#-------------------------------------------------------------------------------
#   Mostafavi-esque gene feature enrichments for causal/non-causal genes
#-------------------------------------------------------------------------------

# Get ENSGIDs for genes with P(causal) > 50% and those without
c_genes <- unique( tr$ensgid[ bgl_glm$preds$scaled >  0.5 ] )
n_genes <- unique( tr$ensgid[ bgl_glm$preds$scaled <= 0.5 ] )
# c_genes <- unique( tr$ensgid[ tr$causal ] )
# n_genes <- unique( tr$ensgid[ !tr$causal ] )
length(c_genes); length(n_genes)

# Read in Mostafavi data
mo_file <- "~/projects/causal_genes/mostafavi2023_gene_annots/pc_genes.txt"
mo <- fread(mo_file)
names(mo)[ names(mo) == "GeneSymbol" ] <- "ensgid"

# Subset Mostafavi data to predicted causal/non-causal genes
mo2 <- mo[ mo$ensgid %in% c( c_genes, n_genes ) , ]
mo2$causal <- mo2$ensgid %in% c_genes
mo2$pLI_gt_0.9 <- mo2$pLI > 0.9
head(mo2)

# Compare high-pLI
x1 <- agresti_95ci_prop( X=sum( mo2$causal & mo2$pLI > 0.9 & !is.na(mo2$pLI) ),
                         n=sum( mo2$causal & !is.na(mo2$pLI) ) )
x2 <- agresti_95ci_prop( X=sum( !mo2$causal & mo2$pLI > 0.9 & !is.na(mo2$pLI) ),
                         n=sum( !mo2$causal & !is.na(mo2$pLI) ) )
rbind( x1, x2 )
mod1 <- glm( formula=pLI_gt_0.9 ~ causal, data=mo2, family="binomial" )
summary(mod1)$coef
x3 <- as.data.frame( rbind( x1, x2 ) )
x3$col <- brewer.pal( n=6, name="Greens" )[5]
# row.names(x3) <- c( "Selected genes", "Remaining genes" )
row.names(x3) <- c( "Causal genes (training)", "Non-causal genes (training)" )
forest_plot( df=x3, value.col="prop", mtext.col=NULL, colour.col="col",
             xlab="Proportion mutation intolerant (pLI > 90%)", xmax=0.4,
             margins=c(5,12,1,1) )

# Compare TFs
x1 <- agresti_95ci_prop( X=sum( mo2$causal & mo2$TF == 1 ),
                   n=sum( mo2$causal ) )
x2 <- agresti_95ci_prop( X=sum( !mo2$causal & mo2$TF == 1 ),
                   n=sum( !mo2$causal ) )
rbind( x1, x2 )
mod2 <- glm( formula=TF ~ causal, data=mo2, family="binomial" )
summary(mod2)$coef

# Compare TSS count
mean( mo2$promoter_count[ mo2$causal] )
mean( mo2$promoter_count[ !mo2$causal] )
mod3 <- glm( formula=promoter_count ~ causal, data=mo2, family="gaussian" )
summary(mod3)$coef


#-------------------------------------------------------------------------------
#   Number of genes that are causal for multiple traits
#-------------------------------------------------------------------------------

tab               <- table( tr$ensgid[tr$causal], tr$trait[tr$causal] )
n_traits_per_gene <- apply( X=tab, MARGIN=1, FUN=function(x) sum( x > 0 ) )
table( n_traits_per_gene > 1 )


#-------------------------------------------------------------------------------
#   GLM v. GLMM
#-------------------------------------------------------------------------------

bgl_glm  <- loto( data=tr, method="glm",    feat_cols=bg_cols, bias_cols=bias_cols )
bgl_glmm <- loto( data=tr, method="glmm",   feat_cols=bg_cols, bias_cols=bias_cols )
bgl_glm_pr   <- plot_loto_pr( loto_obj=bgl_glm,  color=viridis(11), legend=TRUE )
bgl_glmm_pr  <- plot_loto_pr( loto_obj=bgl_glmm, color=viridis(11), legend=TRUE )
rbind( bgl_glm_pr, bgl_glmm_pr )


#-------------------------------------------------------------------------------
#   Genebass
#-------------------------------------------------------------------------------

# Read in the Genebass data
gb <- fread( file.path( maindir, "genebass_data.txt" ) )
names(gb) <- tolower( names(gb) )
gb_cols <- c( "trait", "gene", "burden_ms", "burden_plof", "skato_ms", "skato_plof")

# Join
gb2 <- left_join( x=tr[ tr$trait %in% unique(gb$trait) ], y=gb[,..gb_cols] )

# Format
gb2$bur_mis_p <- ifelse( is.na(gb2$burden_ms),   0, -log10(gb2$burden_ms) )
gb2$bur_lof_p <- ifelse( is.na(gb2$burden_plof), 0, -log10(gb2$burden_plof) )
gb2$bur_mis_z <- ifelse( is.na(gb2$burden_ms),   0, p_to_z(gb2$burden_ms) )
gb2$bur_lof_z <- ifelse( is.na(gb2$burden_plof), 0, p_to_z(gb2$burden_plof) )
gb2$skt_mis_p <- ifelse( is.na(gb2$skato_ms),    0, -log10(gb2$skato_ms) )
gb2$skt_lof_p <- ifelse( is.na(gb2$skato_plof),  0, -log10(gb2$skato_plof) )
gb2$skt_mis_z <- ifelse( is.na(gb2$skato_ms),    0, p_to_z(gb2$skato_ms) )
gb2$skt_lof_z <- ifelse( is.na(gb2$skato_plof),  0, p_to_z(gb2$skato_plof) )

# Plot: P value (or z-score) vs. P(causal)
plot_logOR_relationship_smooth( data=gb2[ gb2$bur_mis_p != 0 , ], 
                                varname="bur_mis_p" )
plot_logOR_relationship_smooth( data=gb2[ gb2$bur_mis_z != 0 , ], 
                                varname="bur_mis_z" )
plot_logOR_relationship_smooth( data=gb2[ gb2$bur_lof_p != 0 , ], 
                                varname="bur_lof_p" )
plot_logOR_relationship_smooth( data=gb2[ gb2$bur_lof_z != 0 , ], 
                                varname="bur_lof_z" )
plot_logOR_relationship_smooth( data=gb2[ gb2$skt_mis_p != 0 , ], 
                                varname="skt_mis_p" )
plot_logOR_relationship_smooth( data=gb2[ gb2$skt_mis_z != 0 , ], 
                                varname="skt_mis_z" )
plot_logOR_relationship_smooth( data=gb2[ gb2$skt_lof_p != 0 , ], 
                                varname="skt_lof_p" )
plot_logOR_relationship_smooth( data=gb2[ gb2$skt_lof_z != 0 , ], 
                                varname="skt_lof_z" )

# What is the PPV of P < 1e-5?
tp_bur_mis    <- sum( gb2$causal & gb2$bur_mis_p >= 5 )
calls_bur_mis <- sum( gb2$bur_mis_p >= 5 )
tp_bur_lof    <- sum( gb2$causal & gb2$bur_lof_p >= 5 )
calls_bur_lof <- sum( gb2$bur_lof_p >= 5 )
tp_skt_mis    <- sum( gb2$causal & gb2$skt_mis_p >= 5 )
calls_skt_mis <- sum( gb2$skt_mis_p >= 5 )
tp_skt_lof    <- sum( gb2$causal & gb2$skt_lof_p >= 5 )
calls_skt_lof <- sum( gb2$skt_lof_p >= 5 )
agresti_95ci_prop( X=tp_bur_mis, n=calls_bur_mis )
agresti_95ci_prop( X=tp_bur_lof, n=calls_bur_lof )
agresti_95ci_prop( X=tp_skt_mis, n=calls_skt_mis )
agresti_95ci_prop( X=tp_skt_lof, n=calls_skt_lof )

# What are the CALDERA_multi scores for "non-causal" LOFs?
bad_lof <- gb2[ !gb2$causal & gb2$bur_lof_p >= 5 , ]
bad_lof_tgp <- paste( bad_lof$trait, bad_lof$ensgid, bad_lof$cs_id, sep="_" )
tr_lof_tgp  <- paste( tr$trait,      tr$ensgid,      tr$cs_id, sep="_" )
bad_lof$multi <- bgl_las$preds$pred[ match( bad_lof_tgp, tr_lof_tgp ) ]
hist( bad_lof$multi, las=1, col=brewer.pal( n=6, name="Greens" )[5],
      breaks=seq(0,1,0.05) )
mean( bad_lof$multi > 0.5 )
sum(  bad_lof$multi > 0.5 )
sort(bad_lof$multi)
bad_lof_tcp <- paste( "TRUE",    bad_lof$trait, bad_lof$region, bad_lof$cs_id, sep="_" )
tr_lof_tcp  <- paste( tr$causal, tr$trait,      tr$region,      tr$cs_id, sep="_" )
x1 <- tr[ match( bad_lof_tcp, tr_lof_tcp ) , ]
par( mfrow=c(2,1) )
hist( x1$distance_genebody, breaks=seq(0,4e5,1e4) )
hist( tr$distance_genebody[tr$causal], breaks=seq(0,4e5,1e4) )
par( mfrow=c(1,1) )
summary(x1$distance_genebody)
hist( tr$distance_genebody[tr$causal], ylim=c(0,50) )
mean( tr$distance_genebody[tr$causal] < 1e4, na.rm=TRUE )
mean( x1$distance_genebody < 1e4, na.rm=TRUE )


#-------------------------------------------------------------------------------
#   Functionally-validated genes with high CALDERA scores
#-------------------------------------------------------------------------------

# Get OT genes with CALDERA scores
# Extract study ID
ot <- lg_pr$data %>%
  select( ensgid:causal, L2G, CALDERA, )
pattern <- "(.*)-[[:digit:]]+_[[:digit:]]+_[[:alpha:]]+_[[:alpha:]]+-credset.*"
ot$study_id <- sub( pattern=pattern, replacement="\\1", x=ot$CS_ID )
head(ot,2)

# Read in Mountjoy2021 Table S8
mnt_file <- file.path( maindir, "benchmarking_datasets/mountjoy2021_table_s8_v2.txt" )
mnt <- fread(mnt_file)
names(mnt) <- tolower( names(mnt) )
names(mnt) <- gsub( pattern="\\s|-",     replacement="_", x=names(mnt) )
names(mnt) <- gsub( pattern="\\)|\\(", replacement="",  x=names(mnt) )
mnt <- mnt %>%
  # filter( gold_standard_confidence == "High" ) %>%
  # filter( grepl( pattern="^functional", gold_standard_class ) ) %>%
  rename( ensgid=gene_id ) %>%
  select( study_id, ensgid, gene_name, trait_name, gold_standard_confidence,
          gold_standard_class, gold_standard_evidence )
head(mnt,2)
mnt %>%
  select( study_id, ensgid) %>%
  duplicated() %>%
  which()

# Join
mnt_cols <- c( "study_id", "ensgid", "trait_name", "gold_standard_evidence", "gold_standard_class" )
ot2 <- inner_join( x=ot, y=mnt[,..mnt_cols], by=c( "study_id", "ensgid" ) ) %>%
  arrange( -CALDERA )
head(ot2,2)
dim(ot2)
plot( ecdf(ot2$CALDERA) )
mean( ot2$CALDERA > 0.5 ); sum( ot2$CALDERA > 0.5 )
mean( ot2$CALDERA > 0.75 ); sum( ot2$CALDERA > 0.75 )
ot2[ ot2$CALDERA > 0.75 ,  ]


#-------------------------------------------------------------------------------
#   Overlap between benchmarking causal genes
#-------------------------------------------------------------------------------

names(lg_pr$data)
names(ex_pr$data)
names(mt_pr$data)

# Shared loci: CALDERA vs. OT
shared_genes_ot <- intersect( tr$ensgid, lg_pr$data$ensgid )
ot_loci  <- lg_pr$data[ lg_pr$data$ensgid %in% shared_genes_ot , ]
cal_loci <- tr[ tr$ensgid %in% shared_genes_ot , ]
cal_loci$tcp <- paste0( cal_loci$trait, "_", cal_loci$region, "_", cal_loci$cs_id )
length( unique(cal_loci$tcp) )
length( unique(ot_loci$CS_ID) )

# Subset CALDERA training data to overlapping traits with MT
head(mt_pr$data,2)
cal_mt <- tr %>%
  filter( trait %in% tmt_traits ) %>%
  select( causal, trait, gene, ensgid, region, cs_id )
cal_mt$tgp <- paste0( cal_mt$trait, "_", cal_mt$ensgid )

# Format MT
mt <- mt_pr$data %>%
  select( ensgid:symbol, causal, CS_ID ) %>%
  filter( pheno == "IGF1" ) %>%
  rename( trait = pheno,
          gene  = symbol ) %>%
  mutate( causal = as.logical(causal) )
mt$tgp <- paste0( mt$trait, "_", mt$ensgid )

# Intersect
shared_tgp_mt <- intersect( cal_mt$tgp, mt$tgp )
cal_mt2 <- cal_mt[ match( shared_tgp_mt, cal_mt$tgp ) ]
mt2     <- mt[     match( shared_tgp_mt, mt$tgp ) ]
cal_mt2$mt_causal <- mt2$causal
cal_mt2$mt_ca <- mt2$CS_ID

# Are the causal genes the same?
table( cal_mt2$causal, mt2$causal )

#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Refine basic+GLC model
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Find most extreme decrease in P(causal) after adding GLCs
#-------------------------------------------------------------------------------

# Combine predictions with/without GLCs and pLI features
p_tol <- cbind( bl_las$preds[ , c( 1:3, 6 ) ],       
                glc=bgl_las$preds[ , c( "pred", "bil" ) ], 
                pLI_lt_0.1=tr$pLI_lt_0.1, pLI_lt_0.9=tr$pLI_lt_0.9, pLI=tr$pLI )

# Plot predicted P(causal) with and without GLCs
plot( p_tol$pred, p_tol$glc.pred, las=1, col="steelblue", lwd=1.5,
      xlab="Without covariates", ylab="With covariates"  )
abline( a=0, b=1, lty=2, lwd=2 )

# Example for manuscript: most extreme decrease in P(causal) after adding GLCs
summary(  abs( p_tol$pred - p_tol$glc.pred ) )
quantile( abs( p_tol$pred - p_tol$glc.pred ), probs=seq(0.9,1,0.01) )
diff <- abs( p_tol$pred - p_tol$glc.pred )
idx <- which( diff == max(diff) )
tr[ idx , ]
p_tol[ idx , ]

# Proportion of prioritized genes that are tolerant to mutation (pLI < 10%)
mean( p_tol$pLI_lt_0.1[ p_tol$bil     == 1 & p_tol$pred      > 0.5 ] )
mean( p_tol$pLI_lt_0.1[ p_tol$glc.bil == 1 & p_tol$glc.pred  > 0.5 ] )

# Proportion of prioritized genes that are intolerant to mutation (pLI > 90%)
mean( !p_tol$pLI_lt_0.9[ p_tol$bil     == 1 & p_tol$pred      > 0.5 ] )
mean( !p_tol$pLI_lt_0.9[ p_tol$glc.bil == 1 & p_tol$glc.pred  > 0.5 ] )

# Mean pLI of prioritized genes
mean( p_tol$pLI[        p_tol$bil     == 1 & p_tol$pred      > 0.5 ] )
mean( p_tol$pLI[        p_tol$glc.bil == 1 & p_tol$glc.pred  > 0.5 ] )


#-------------------------------------------------------------------------------
#   Do we need a covariate in the model to account for the fact that different
#   traits have different power for POPS?
#-------------------------------------------------------------------------------

# Look at correlation between mean, SD, skewness, and kurtosis across traits
tr$tcp <- paste0( tr$trait, "_", tr$region, "_", tr$cs_id )
pl <- tapply( X=tr$tcp,      INDEX=tr$trait, FUN=function(x) length( unique(x) ) )
pm <- tapply( X=tr$pops_glo, INDEX=tr$trait, FUN=mean )
pv <- tapply( X=tr$pops_glo, INDEX=tr$trait, FUN=var )
ps <- tapply( X=tr$pops_glo, INDEX=tr$trait, FUN=skewness )
pk <- tapply( X=tr$pops_glo, INDEX=tr$trait, FUN=kurtosis )
corrplot( cor( cbind( pl, pm, pv, ps, pk ) ) )

# eCDF
plot( ecdf(tr$pops_glo), col="white", las=1 )
for( i in names( sort(pl) ) ){
  sub <- tr[ tr$trait == i , ]
  plot( ecdf(sub$pops_glo), add=TRUE, 
        col=viridis( n=length(pl) )[ which( names( sort(pl) ) == i)] )
}

# Make columns for the value
tr$n_tcp     <- pl[ match( tr$trait, names(pl) ) ]
tr$pops_mean <- pm[ match( tr$trait, names(pm) ) ]
tr$pops_var  <- pv[ match( tr$trait, names(pv) ) ]
tr$pops_skew <- ps[ match( tr$trait, names(ps) ) ]
tr$pops_int  <- tr$pops_glo * tr$pops_mean

# Add trait-level covariates to the basic model
sgt_cols <- c( bg_cols, "pops_mean", "pops_int" )
sgt_las <- cv.glmnet( x=as.matrix( tr[ , ..sgt_cols ] )[,-1], 
                      foldid=as.integer( as.factor(tr$trait) ),
                      y=tr[["causal"]], family="binomial" )
sgt_coef <- as.data.frame( coef( object=sgt_las, s="lambda.min" ) )
lasso_to_forest( las=sgt_las, lambda_type="lambda.min", train=tr, rm_feats=bias_cols )

# Look at how P(causal) looks at 25th, 50th, and 75th percentile...
# ...for both pops_glo and pops_mean. Assuming other variables are at their mean.
feat_cols <- sgt_coef@Dimnames[[1]][-1]
cm <- colMeans( x=tr[ , ..feat_cols ] )
df <- as.data.frame( matrix( rep( cm, 9), nrow=9, byrow=TRUE ) )
names(df) <- names(cm)
df$pops_glo  <- rep( quantile( tr$pops_glo,  probs=c( 0.1, 0.5, 0.9 ) ), 3 )
df$pops_mean <- rep( quantile( tr$pops_mean, probs=c( 0.1, 0.5, 0.9 ) ), 
                     rep( 3, 3 ) )
df$pops_rel  <- 0
df$pops_bil  <- 1
df$dist_gene_rel <- -3
df$p_causal <- as.vector( predict( object=sgt_las, s="lambda.min", type="response",
                                   newx=as.matrix( df[ , feat_cols ] ) ) )

# Plot
plot( df$pops_glo, df$p_causal, type="n", las=1, 
      # ylim=c( 0, max(df$p_causal) ),
      xlab="pops_glo", ylab="P(causal)" )
cols <- viridis(3)
for( i in unique(df$pops_mean) ){
  sub <- df[ df$pops_mean == i , ]
  idx <- which( unique(df$pops_mean) == i )
  lines( sub$pops_glo, sub$p_causal, col=cols[idx], lwd=3 )
}
legend( "topleft", legend=paste( "Mean PoPS:", c( "10th", "50th", "90th" ) ), fill=cols )


#-------------------------------------------------------------------------------
#   Do we need a covariate in the model to account for the fact that different
#   traits have different power for MAGMA?
#-------------------------------------------------------------------------------

# Look at correlation between mean, SD, skewness, and kurtosis across traits
mm <- tapply( X=tr$magma_glo, INDEX=tr$trait, FUN=mean )
mv <- tapply( X=tr$magma_glo, INDEX=tr$trait, FUN=var )
ms <- tapply( X=tr$magma_glo, INDEX=tr$trait, FUN=skewness )
mk <- tapply( X=tr$magma_glo, INDEX=tr$trait, FUN=kurtosis )
corrplot( cor( cbind( mm, mv, ms, mk ) ) )

# eCDF
plot( ecdf(tr$magma_glo), col="white", las=1 )
for( i in names( sort(pl) ) ){
  sub <- tr[ tr$trait == i , ]
  plot( ecdf(sub$magma_glo), add=TRUE, 
        col=viridis( n=length(pl) )[ which( names( sort(pl) ) == i)] )
}

# Make columns for the value
tr$magma_mean <- pm[ match( tr$trait, names(mm) ) ]
tr$magma_var  <- pv[ match( tr$trait, names(mv) ) ]
tr$magma_skew <- ps[ match( tr$trait, names(ms) ) ]
tr$magma_int  <- tr$magma_glo * tr$magma_mean

# Add trait-level covariates to the basic model
sgt_cols <- c( bg_cols, "magma_mean", "pops_mean" )
sgt_las <- cv.glmnet( x=as.matrix( tr[ , ..sgt_cols ] )[,-1], 
                      foldid=as.integer( as.factor(tr$trait) ),
                      y=tr[["causal"]], family="binomial" )
sgt_coef <- coef( object=sgt_las, s="lambda.min" )
lasso_to_forest( las=sgt_las, lambda_type="lambda.min", train=tr, rm_feats=bias_cols )

# Look at how P(causal) looks at 25th, 50th, and 75th percentile...
# ...for both magma_glo and magma_mean. Assuming other variables are at their mean.
feat_cols <- sgt_coef@Dimnames[[1]][-1]
cm <- colMeans( x=tr[ , ..feat_cols ] )
df <- as.data.frame( matrix( rep( cm, 9), nrow=9, byrow=TRUE ) )
names(df) <- names(cm)
df$magma_glo  <- rep( quantile( tr$magma_glo,  probs=c( 0.1, 0.5, 0.9 ) ), 3 )
df$magma_mean <- rep( quantile( tr$magma_mean, probs=c( 0.1, 0.5, 0.9 ) ), 
                     rep( 3, 3 ) )
df$magma_rel  <- 0
df$magma_bil  <- 1
df$dist_gene_rel <- -3
df$p_causal <- as.vector( predict( object=sgt_las, s="lambda.min", type="response",
                                   newx=as.matrix( df[ , feat_cols ] ) ) )

# Plot
plot( df$magma_glo, df$p_causal, type="n", las=1, 
      # ylim=c( 0, max(df$p_causal) ),
      xlab="magma_glo", ylab="P(causal)" )
cols <- viridis(3)
for( i in unique(df$magma_mean) ){
  sub <- df[ df$magma_mean == i , ]
  idx <- which( unique(df$magma_mean) == i )
  lines( sub$magma_glo, sub$p_causal, col=cols[idx], lwd=3 )
}
legend( "topleft", legend=paste( "Mean MAGMA:", c( "10th", "50th", "90th" ) ), fill=cols )


#-------------------------------------------------------------------------------
#   Does adding quadratic terms improve the model? No.
#-------------------------------------------------------------------------------

# Which polynomial leads to the biggest improvement in model fit?
poly_lrt <- function(variable){
  u_form <- update( bg_glm$formula, paste0( "~ . - ", variable, " + poly(", variable, ", 2)" ) )
  u_glm  <- glm( formula=u_form, data=tr, family="binomial" )
  anova( bg_glm, u_glm, test="Chi" )$`Pr(>Chi)`[2]
}
poly_lrt( variable="pops_glo" )
poly_lrt( variable="pops_rel" )
poly_lrt( variable="dist_gene_rel" )
poly_lrt( variable="magma_rel" )
poly_lrt( variable="coding_rel" )
poly_lrt( variable="prior_n_genes_locus" )
poly_lrt( variable="ABC_count" )
poly_lrt( variable="Roadmap_count" )


#-------------------------------------------------------------------------------
#   Does adding interaction terms improve the model? Hardly.
#-------------------------------------------------------------------------------

# Add all possible interactions
i_form  <- update( bg_glm$formula, . ~ (.)^2 )

# Perform backwards stepwise model selection
# Remove features predicted to have P > 0.05/5
i_glm0 <- glm( formula=i_form, data=tr, family="binomial" )
i_glm <- stats::step( object=i_glm0, k=qchisq( 1-0.05/10, df=1 ) )

# Forest plot
glm_to_forest_p( mod=i_glm, xmax=NULL )
summary(i_glm)$coef

# Look at how P(causal) looks at 25th, 50th, and 75th percentile...
# ...for both pops_glo and Roadmap_count, assuming other variables are at their mean.
feat_cols <- all.vars(i_form)[-1]
cm <- colMeans( x=tr[ , ..feat_cols ] )
df <- as.data.frame( matrix( rep( cm, 9), nrow=9, byrow=TRUE ) )
names(df) <- names(cm)
df$pops_glo      <- rep( quantile( tr$pops_glo,      probs=c( 0.1, 0.5, 0.9 ) ), 3 )
df$Roadmap_count <- rep( quantile( tr$Roadmap_count, probs=c( 0.1, 0.5, 0.9 ) ), 
                         rep( 3, 3 ) )
df$pops_rel  <- 0
df$dist_gene_rel <- -3
pred <- predict( object=i_glm, newdata=df, se=TRUE )
df$p_causal <- logistic(pred$fit)

# Plot
plot( df$pops_glo, df$p_causal, type="n", las=1, ylim=c( 0, max(df$p_causal) ),
      xlab="pops_glo", ylab="P(causal)" )
cols <- viridis(3)
for( i in unique(df$Roadmap_count) ){
  sub <- df[ df$Roadmap_count == i , ]
  idx <- which( unique(df$Roadmap_count) == i )
  lines( sub$pops_glo, sub$p_causal, col=cols[idx], lwd=3 )
}
legend( "topleft", legend=paste( "Roadmap count:", c( "10th", "50th", "90th" ) ), fill=cols )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   POPS + local
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

# Run a naive regression using published best-in-locus
lcols <- c( "causal", "pops_and_local" )

# Run regression
l_form <- make_formula( lhs=lcols[1], rhs=lcols[-1] )
l_glm  <- glm( formula=l_form, data=tr, family="binomial" )

# Forest plot
glm_to_forest(l_glm)

# Precision and recall: training
pnl_p <- sum( tr$causal & tr$pops_and_local ) / sum( tr$pops_and_local ) #precision
pnl_r <- sum( tr$causal & tr$pops_and_local ) / sum( tr$causal )         #recall

# Precision and recall: test
sum( te$causal & te$pops_and_local ) / sum( te$pops_and_local ) #precision
sum( te$causal & te$pops_and_local ) / sum( te$causal )         #recall


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Global-only
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Feature engineering: compare definitions, find linear relationships
#-------------------------------------------------------------------------------

# Plot relationship: distance to gene body
plot_logOR_relationship_smooth( data=tr, varname="distance_genebody" )
plot_logOR_relationship_smooth( data=tr[ tr$distance_genebody > 0 ], 
                                varname="distance_genebody" )
plot_logOR_relationship_smooth( data=tr, varname="dist_gene_glo" )

# Plot relationship: distance to TSS
plot_logOR_relationship_smooth( data=tr, varname="distance_tss" )
plot_logOR_relationship_smooth( data=tr, varname="dist_tss_glo" )

# Plot relationship: POPS
plot_logOR_relationship_smooth( data=tr, varname="pops_glo" )

# Plot relationship: TWAS
mean( tr$causal[ is.na(tr$twas_z) ])
plot_logOR_relationship_smooth( data=tr[ !is.na(tr$twas_abs_z) , ], varname="twas_abs_z" )
plot_logOR_relationship_smooth( data=tr, varname="twas_glo" )

# Plot relationship: E-P Liu
mean( tr$causal[ is.na(tr$corr_liu_score) ])
plot_logOR_relationship_smooth( data=tr[ !is.na(tr$corr_liu_score) , ], 
                                varname="corr_liu_score" )
plot_logOR_relationship_smooth( data=tr, varname="corr_liu_glo" )

# Plot relationship: E-P Andersson
mean( tr$causal[ is.na(tr$corr_andersson_score) ] )
plot_logOR_relationship_smooth( data=tr[ !is.na(tr$corr_andersson_score) , ], 
                                varname="corr_andersson_score" )
plot_logOR_relationship_smooth( data=tr, varname="corr_and_glo" )

# Plot relationship: E-P Ulirsch
mean( tr$causal[ is.na(tr$corr_ulirsch_score) ] )
plot_logOR_relationship_smooth( data=tr[ !is.na(tr$corr_ulirsch_score) , ], 
                                varname="corr_ulirsch_score" )
plot_logOR_relationship_smooth( data=tr, varname="corr_uli_glo" )

# Plot relationship: PCHi-C Jung
mean( tr$causal[ is.na(tr$pchic_jung_score) ] )
plot_logOR_relationship_smooth( data=tr[ !is.na(tr$pchic_jung_score) , ], 
                                varname="pchic_jung_score" )
plot_logOR_relationship_smooth( data=tr, varname="pchic_jung_glo" )

# Plot relationship: PCHi-C Javierre
mean( tr$causal[ is.na(tr$pchic_javierre_score) ] )
plot_logOR_relationship_smooth( data=tr[ !is.na(tr$pchic_javierre_score) , ], 
                                varname="pchic_javierre_score" )
plot_logOR_relationship_smooth( data=tr, varname="pchic_jav_glo" )

# Plot relationship: CLPP
mean( tr$causal[ is.na(tr$clpp_prob) ] )
plot_logOR_relationship_smooth( data=tr[ !is.na(tr$clpp_prob) , ], 
                                varname="clpp_prob" )
plot_logOR_relationship_smooth( data=tr, varname="clpp_glo" )

# Plot relationship: ABC
mean( tr$causal[ is.na(tr$abc_score) ] )
plot_logOR_relationship_smooth( data=tr[ !is.na(tr$abc_score) , ], 
                                varname="abc_score" )
plot_logOR_relationship_smooth( data=tr, varname="abc_glo" )

# Plot relationship: MAGMA
mean( tr$causal[ is.na(tr$magma_score) ] )
plot_logOR_relationship_smooth( data=tr[ !is.na(tr$magma_score) , ], 
                                varname="magma_score" )
plot_logOR_relationship_smooth( data=tr, varname="magma_glo" )

# Plot relationship: SMR
mean( tr$causal[ is.na(tr$smr_z) ] )
plot_logOR_relationship_smooth( data=tr[ !is.na(tr$smr_z) , ], varname="smr_z" )
plot_logOR_relationship_smooth( data=tr, varname="smr_glo" )

# Plot relationship: coding
mean( tr$causal[ is.na(tr$coding_prob) ] )
plot_logOR_relationship_smooth( data=tr[ !is.na(tr$coding_prob) , ], varname="coding_prob" )
plot_logOR_relationship_smooth( data=tr, varname="coding_glo" )

# DEPICT and NetWAS
mean( tr$causal[ is.na(tr$depict_z) ] )
plot_logOR_relationship_smooth( data=tr[ !is.na(tr$depict_z) , ], 
                                varname="depict_z" )
plot_logOR_relationship_smooth( data=tr, varname="depict_z_glo" )

# NetWAS
mean( tr$causal[ is.na(tr$no_loco_nc_netwas_score) ] )
plot_logOR_relationship_smooth( data=tr[ !is.na(tr$no_loco_nc_netwas_score) , ], 
                                varname="no_loco_nc_netwas_score" )
plot_logOR_relationship_smooth( data=tr, varname="netwas_score_glo" )

# NetWAS Bonferroni
mean( tr$causal[ is.na(tr$no_loco_nc_netwas_bon_score) ] )
plot_logOR_relationship_smooth( data=tr[ !is.na(tr$no_loco_nc_netwas_bon_score) , ], 
                                varname="no_loco_nc_netwas_bon_score" )
plot_logOR_relationship_smooth( data=tr, varname="netwas_bon_score_glo" )


#-------------------------------------------------------------------------------
#   Compare competing feature definitions
#-------------------------------------------------------------------------------

# Distance to gene body
mod1 <- glm( causal ~ distance_genebody, data=tr, family="binomial" )
mod3 <- glm( causal ~ dist_gene_glo, data=tr, family="binomial" )
dev_exp(mod1)
dev_exp(mod3)
summary(mod1)$coef
summary(mod3)$coef

# Distance to TSS
mod1 <- glm( causal ~ distance_tss,     data=tr, family="binomial" )
mod3 <- glm( causal ~ dist_tss_glo, data=tr, family="binomial" )
dev_exp(mod1)
dev_exp(mod3)
summary(mod1)$coef
summary(mod3)$coef

# Number of genes in the locus
mod1 <- glm( causal ~ ngenes_nearby,       data=tr, family="binomial" )
mod2 <- glm( causal ~ prior_n_genes_locus, data=tr, family="binomial" )
dev_exp(mod1)
dev_exp(mod2)
summary(mod1)$coef
summary(mod2)$coef

# TWAS
mod1 <- glm( causal ~ twas_logp_glo,  data=tr, family="binomial" )
mod2 <- glm( causal ~ twas_glo,  data=tr, family="binomial" )
dev_exp(mod1)
dev_exp(mod2)
summary(mod1)$coef
summary(mod2)$coef

# Plot relationship: E-P Liu
mod1 <- glm( causal ~ corr_liu_raw_l2g,  data=tr, family="binomial" )
mod2 <- glm( causal ~ corr_liu_glo,  data=tr, family="binomial" )
dev_exp(mod1)
dev_exp(mod2)
summary(mod1)$coef
summary(mod2)$coef

# CLPP
mod1 <- glm( causal ~ clpp_raw_l2g,  data=tr, family="binomial" )
mod2 <- glm( causal ~ clpp_glo,  data=tr, family="binomial" )
dev_exp(mod1)
dev_exp(mod2)
summary(mod1)$coef
summary(mod2)$coef

# ABC
mod1 <- glm( causal ~ abc_raw_l2g,  data=tr, family="binomial" )
mod2 <- glm( causal ~ abc_glo,  data=tr, family="binomial" )
dev_exp(mod1)
dev_exp(mod2)
summary(mod1)$coef
summary(mod2)$coef

# DEPICT and NetWAS
mod1 <- glm( causal ~ depict_z_glo,         data=tr, family="binomial" )
mod2 <- glm( causal ~ netwas_score_glo,     data=tr, family="binomial" )
mod3 <- glm( causal ~ netwas_bon_score_glo, data=tr, family="binomial" )
dev_exp(mod1)
dev_exp(mod2)
dev_exp(mod3)


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Relative-only
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Feature engineering: compare with/without best value (captured by BIL features)
#-------------------------------------------------------------------------------

# Plot relationship: distance to gene body
plot_logOR_relationship_smooth( data=tr, varname="dist_gene_rel" )
plot_logOR_relationship_smooth( data=tr[ tr$dist_gene_rel != -3 , ], varname="dist_gene_rel" )

# Plot relationship: distance to TSS
plot_logOR_relationship_smooth( data=tr, varname="dist_tss_rel" )
plot_logOR_relationship_smooth( data=tr[ tr$dist_tss_rel != -3 , ], varname="dist_tss_rel" )

# Plot relationship: POPS
plot_logOR_relationship_smooth( data=tr, varname="pops_rel" )
plot_logOR_relationship_smooth( data=tr[ tr$pops_rel !=0 , ], varname="pops_rel" )

# Plot relationship: TWAS
plot_logOR_relationship_smooth( data=tr, varname="twas_rel" )
plot_logOR_relationship_smooth( data=tr[ tr$twas_rel !=0 , ], varname="twas_rel" )

# Plot relationship: E-P Liu
plot_logOR_relationship_smooth( data=tr[ tr$corr_liu_rel != 0 , ], varname="corr_liu_rel" )
plot_logOR_relationship_smooth( data=tr, varname="corr_liu_rel" )

# Plot relationship: E-P Andersson
plot_logOR_relationship_smooth( data=tr[ tr$corr_and_rel != 0 , ], varname="corr_and_rel" )
plot_logOR_relationship_smooth( data=tr, varname="corr_and_rel" )

# Plot relationship: E-P Ulirsch
plot_logOR_relationship_smooth( data=tr[ tr$corr_uli_rel != 0 , ], varname="corr_uli_rel" )
plot_logOR_relationship_smooth( data=tr, varname="corr_uli_rel" )

# Plot relationship: PCHi-C Jung
plot_logOR_relationship_smooth( data=tr[ tr$pchic_jung_rel != 0 , ], varname="pchic_jung_rel" )
plot_logOR_relationship_smooth( data=tr, varname="pchic_jung_rel" )

# Plot relationship: PCHi-C Javierre
plot_logOR_relationship_smooth( data=tr[ tr$pchic_jav_rel != 0 , ], varname="pchic_jav_rel" )
plot_logOR_relationship_smooth( data=tr, varname="pchic_jav_rel" )

# Plot relationship: CLPP
plot_logOR_relationship_smooth( data=tr[ tr$clpp_rel != 0 , ], varname="clpp_rel" )
plot_logOR_relationship_smooth( data=tr, varname="clpp_rel" )

# Plot relationship: ABC
plot_logOR_relationship_smooth( data=tr[ tr$abc_rel != 0 , ], varname="abc_rel" )
plot_logOR_relationship_smooth( data=tr, varname="abc_rel" )

# Plot relationship: MAGMA
plot_logOR_relationship_smooth( data=tr, varname="magma_rel" )
plot_logOR_relationship_smooth( data=tr[ tr$magma_rel != 0 , ], varname="magma_rel" )

# Plot relationship: SMR
plot_logOR_relationship_smooth( data=tr, varname="smr_rel" )
plot_logOR_relationship_smooth( data=tr[ tr$smr_rel != 0 , ], varname="smr_rel" )

# Plot relationship: coding
plot_logOR_relationship_smooth( data=tr, varname="coding_rel" )
plot_logOR_relationship_smooth( data=tr[ tr$coding_rel != 0 , ], varname="coding_rel" )

# Plot relationship: DEPICT
plot_logOR_relationship_smooth( data=tr, varname="depict_z_rel" )
plot_logOR_relationship_smooth( data=tr[ tr$depict_z_rel != 0 , ], varname="depict_z_rel" )

# Plot relationship: NetWAS
plot_logOR_relationship_smooth( data=tr, varname="netwas_score_rel" )
plot_logOR_relationship_smooth( data=tr[ tr$netwas_score_rel != 0 , ], varname="netwas_score_rel" )
plot_logOR_relationship_smooth( data=tr, varname="netwas_bon_score_rel" )
plot_logOR_relationship_smooth( data=tr[ tr$netwas_bon_score_rel != 0 , ], varname="netwas_bon_score_rel" )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   For each V2G method, compare BIL, global, and relative features
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

# Plot deviance explained by all univariate methods
all_de0 <- list()
all_cols <- unique( c( a_cols, glo_cols, rel_cols, cov_cols ) )[-1]
for( i in all_cols ){
  form <- paste0( "causal ~ ", i )
  mod  <- glm( as.formula(form), data=tr, family="binomial" )
  all_de0[[i]] <- c( ( mod$null.deviance - mod$deviance ) / mod$null.deviance )
}
all_de <- unlist(all_de0)
par( mar=c(9,4,1,1) )
barplot( sort( all_de[ all_de > 0.04 ], decreasing=TRUE ), las=2)
par( mar=c(5,5,1,1) )

# POPS
v2g_devexp( data=tr, prefix="pops" )
v2g_cor( data=tr, prefix="pops" )
corrplot( v2g_cor( data=tr, prefix="pops" ) )

# Distance: gene
v2g_devexp( data=tr, prefix="dist_gene" )
v2g_cor( data=tr, prefix="dist_gene" )
corrplot( v2g_cor( data=tr, prefix="dist_gene" ) )

# Distance: TSS
v2g_devexp( data=tr, prefix="dist_tss" )
v2g_cor( data=tr, prefix="dist_tss" )
corrplot( v2g_cor( data=tr, prefix="dist_tss" ) )

# TWAS
v2g_devexp( data=tr, prefix="twas" )
v2g_cor( data=tr, prefix="twas" )
corrplot( v2g_cor( data=tr, prefix="twas" ) )

# MAGMA
v2g_devexp( data=tr, prefix="magma" )
v2g_cor( data=tr, prefix="magma" )
corrplot( v2g_cor( data=tr, prefix="magma" ) )

# SMR
v2g_devexp( data=tr, prefix="smr" )
v2g_cor( data=tr, prefix="smr" )
corrplot( v2g_cor( data=tr, prefix="smr" ) )

# Coding
v2g_devexp( data=tr, prefix="coding" )
v2g_cor( data=tr, prefix="coding" )
corrplot( v2g_cor( data=tr, prefix="coding" ) )

# Liu
v2g_devexp( data=tr, prefix="corr_liu" )
v2g_cor( data=tr, prefix="corr_liu" )
corrplot( v2g_cor( data=tr, prefix="corr_liu" ) )

# Andersson
v2g_devexp( data=tr, prefix="corr_and" )
v2g_cor( data=tr, prefix="corr_and" )
corrplot( v2g_cor( data=tr, prefix="corr_and" ) )

# Ulirsch
v2g_devexp( data=tr, prefix="corr_uli" )
v2g_cor( data=tr, prefix="corr_uli" )
corrplot( v2g_cor( data=tr, prefix="corr_uli" ) )

# Jung
v2g_devexp( data=tr, prefix="pchic_jung" )
v2g_cor( data=tr, prefix="pchic_jung" )
corrplot( v2g_cor( data=tr, prefix="pchic_jung" ) )

# Javierre
v2g_devexp( data=tr, prefix="pchic_jav" )
v2g_cor( data=tr, prefix="pchic_jav" )
corrplot( v2g_cor( data=tr, prefix="pchic_jav" ) )

# CLPP
v2g_devexp( data=tr, prefix="clpp" )
v2g_cor( data=tr, prefix="clpp" )
corrplot( v2g_cor( data=tr, prefix="clpp" ) )

# ABC
v2g_devexp( data=tr, prefix="abc" )
v2g_cor( data=tr, prefix="abc" )
corrplot( v2g_cor( data=tr, prefix="abc" ) )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Gene-level covariates (GLCs)
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Feature engineering: number of nearby genes
#-------------------------------------------------------------------------------

# Number of genes in the locus
plot_logOR_relationship_smooth( data=tr, varname="prior_n_genes_locus" )
plot_logOR_relationship_smooth( data=tr, varname="ngenes_nearby" )


#-------------------------------------------------------------------------------
#   Feature engineering: bias features
#-------------------------------------------------------------------------------

# pLI
plot_logOR_relationship_smooth( data=tr, varname="pLI" )
plot_logOR_relationship_smooth( data=tr, varname="pLI_log10" )

# LOEUF
plot_logOR_relationship_smooth( data=tr, varname="LOEUF" )

# hs
plot_logOR_relationship_smooth( data=tr, varname="hs" )
plot_logOR_relationship_smooth( data=tr, varname="hs_log10" )

# Gene length
plot_logOR_relationship_smooth( data=tr, varname="length" )
plot_logOR_relationship_smooth( data=tr, varname="gene_bp_log10" )

# CDS length
plot_logOR_relationship_smooth( data=tr, varname="CDS_length" )
plot_logOR_relationship_smooth( data=tr, varname="cds_bp_log10" )

# ABC_count
plot_logOR_relationship_smooth( data=tr, varname="ABC_count" )

# ABC_length_per_type
plot_logOR_relationship_smooth( data=tr, varname="ABC_length_per_type" )
plot_logOR_relationship_smooth( data=tr, varname="abc_bp_log10" )

# Roadmap_count
plot_logOR_relationship_smooth( data=tr, varname="Roadmap_count" )

# Roadmap_length_per_type
plot_logOR_relationship_smooth( data=tr, varname="Roadmap_length_per_type" )
plot_logOR_relationship_smooth( data=tr, varname="roadmap_bp_log10" )


#-------------------------------------------------------------------------------
#   Feature engineering: GLCs
#-------------------------------------------------------------------------------

# promoter_count
plot_logOR_relationship_smooth( data=tr, varname="promoter_count" )
plot_logOR_relationship_smooth( data=tr, varname="promoter_count_log10" )

# connect_decile
plot_logOR_relationship_smooth( data=tr, varname="connect_decile" )

# connect_quantile
plot_logOR_relationship_smooth( data=tr, varname="connect_quantile" )

# connectedness
table( tr$causal, tr$connectedness, useNA="if" )
fisher.test( table( tr$causal, tr$connectedness ) )

# PPI_degree_decile
plot_logOR_relationship_smooth( data=tr, varname="PPI_degree_decile" )

# PPI_degree_quantile
plot_logOR_relationship_smooth( data=tr, varname="PPI_degree_quantile" )

# PPI_degree_cat
table( tr$causal, tr$PPI_degree_cat, useNA="if" )
fisher.test( table( tr$causal, tr$PPI_degree_cat ) )

# TF
table( tr$causal, tr$TF, useNA="if" )
fisher.test( table( tr$causal, tr$TF ) )


#-------------------------------------------------------------------------------
#   Compare competing feature definitions
#-------------------------------------------------------------------------------

# Number of genes in the locus
mod1 <- glm( causal ~ prior_n_genes_locus, data=tr, family="binomial" )
mod2 <- glm( causal ~ ngenes_nearby, data=tr, family="binomial" )
dev_exp(mod1)
dev_exp(mod2)
summary(mod1)$coef
summary(mod2)$coef

# Number of SNPs in the CS
mod1 <- glm( causal ~ n_cs_snps,       data=tr, family="binomial" )
mod2 <- glm( causal ~ log10_n_cs_snps, data=tr, family="binomial" )
dev_exp(mod1)
dev_exp(mod2)
summary(mod1)$coef
summary(mod2)$coef

# Width of the CS
mod1 <- glm( causal ~ cs_width,       data=tr, family="binomial" )
mod2 <- glm( causal ~ log10_cs_width, data=tr, family="binomial" )
dev_exp(mod1)
dev_exp(mod2)
summary(mod1)$coef
summary(mod2)$coef

# pLI
mod1 <- glm( causal ~ pLI,                               data=tr, family="binomial" )
mod2 <- glm( causal ~ pLI_lt_0.1,                        data=tr, family="binomial" )
mod3 <- glm( causal ~ pLI_lt_0.1 + pLI_lt_0.9,           data=tr, family="binomial" )
mod4 <- glm( causal ~ pLI_log10,                       data=tr, family="binomial" )
dev_exp(mod1)
dev_exp(mod2)
dev_exp(mod3)
dev_exp(mod4)
AIC(mod1); AIC(mod2); AIC(mod3); AIC(mod4)
summary(mod1)$coef
summary(mod2)$coef
summary(mod3)$coef
summary(mod4)$coef

# PPI
mod1 <- glm( causal ~ PPI_degree_quantile, data=tr, family="binomial" )
mod2 <- glm( causal ~ PPI_degree_decile,   data=tr, family="binomial" )
mod3 <- glm( causal ~ PPI_degree_cat,      data=tr, family="binomial" )
dev_exp(mod1)
dev_exp(mod2)
dev_exp(mod3)

# Connectedness
mod1 <- glm( causal ~ connect_quantile, data=tr, family="binomial" )
mod2 <- glm( causal ~ connect_decile,   data=tr, family="binomial" )
mod3 <- glm( causal ~ connectedness,    data=tr, family="binomial" )
dev_exp(mod1)
dev_exp(mod2)
dev_exp(mod3)

# ABC
mod1 <- glm( causal ~ ABC_count,                data=tr, family="binomial" )
mod2 <- glm( causal ~ abc_bp_log10,             data=tr, family="binomial" )
mod3 <- glm( causal ~ ABC_count + abc_bp_log10, data=tr, family="binomial" )
summary(mod1)$coef
summary(mod2)$coef
summary(mod3)$coef


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Miscellaneous
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Plot coefficient path for non-LOTO basic LASSO + GLCs
#-------------------------------------------------------------------------------

# Plot LASSO coefficient path
plot(bg_las)
plot( bg_las$glmnet.fit, xvar="lambda", label=TRUE )
abline( v=log(bg_las$lambda.min), lty=2, col="grey" )
abline( v=log(bg_las$lambda.1se), lty=2, col="grey" )

# Compare lambda.min and lambda.1se coefficients
cbind( coef( bg_las, s="lambda.min" ),
       coef( bg_las, s="lambda.1se" ),
       0:bg_las$glmnet.fit$beta@Dim[1] )


#-------------------------------------------------------------------------------
#   Compare coefficients and SEs for GLM vs. GLMM
#-------------------------------------------------------------------------------

# GLM
form1 <- make_formula( lhs=bg_cols[1], rhs=bg_cols[-1] )
bg_glm <- glm( formula=form1, data=tr, family="binomial" )

# GLMM
form2 <- update( form1, . ~ . + ( 1 | ensgid_c ) )
bg_glmm <- glmer( formula=form2, data=tr, family="binomial" )

# Compare GLM and GLMM: logOR
plot( summary(bg_glm)$coef[ -1, 1 ], 
      summary(bg_glmm)$coef[ -1, 1 ], 
      las=1, xlab="GLM", ylab="GLMM", col="steelblue", cex=1.5, lwd=2 )
abline( a=0, b=1, lty=2 )

# Compare GLM and GLMM: SE
plot( summary(bg_glm)$coef[ -1, 2 ], 
      summary(bg_glmm)$coef[ -1, 2 ], 
      las=1, xlab="GLM", ylab="GLMM", col="steelblue", cex=1.5, lwd=2 )
abline( a=0, b=1, lty=2 )


#-------------------------------------------------------------------------------
#   In-sample AUPRC of the CALDERA model
#-------------------------------------------------------------------------------

prc <- pr.curve( scores.class0 = recal_preds$recal[  recal_preds$causal ], 
                 scores.class1 = recal_preds$recal[ !recal_preds$causal ], 
                 curve = TRUE )
aucpr.conf.int.expit( estimate = prc$auc.integral, 
                      num.pos  = sum(recal_preds$causal) )


#-------------------------------------------------------------------------------
#   AUPRC separately for each trait
#-------------------------------------------------------------------------------

# Full dataset
auprc0 <- list()
for( i in bgl_las$trait ){
  prc <- pr.curve( scores.class0 = bgl_las$preds$recal[  bgl_las$preds$causal & 
                                                           bgl_las$preds$trait==i ], 
                   scores.class1 = bgl_las$preds$recal[ !bgl_las$preds$causal & 
                                                          bgl_las$preds$trait==i ], 
                   curve = TRUE )
  auprc0[[i]] <- aucpr.conf.int.expit( estimate = prc$auc.integral, 
                                       num.pos  = sum( bgl_las$preds$causal &
                                                         bgl_las$preds$trait==i ) )
}

# Independent dataset
auprc0 <- list()
for( i in bgl_las_uniq$trait ){
  prc <- pr.curve( scores.class0 = bgl_las_uniq$preds$recal[  bgl_las_uniq$preds$causal & 
                                                           bgl_las_uniq$preds$trait==i ], 
                   scores.class1 = bgl_las_uniq$preds$recal[ !bgl_las_uniq$preds$causal & 
                                                          bgl_las_uniq$preds$trait==i ], 
                   curve = TRUE )
  auprc0[[i]] <- aucpr.conf.int.expit( estimate = prc$auc.integral, 
                                       num.pos  = sum( bgl_las_uniq$preds$causal &
                                                         bgl_las_uniq$preds$trait==i ) )
}

# Plot
auprc <- as.data.frame( do.call( rbind, auprc0) )
auprc <- auprc[ order(auprc$auprc) , ]
mean_re <- rma( yi=auprc$auprc_logit, sei=auprc$se_logit )
mean_fe <- rma( yi=auprc$auprc_logit, sei=auprc$se_logit, method="FE" )
forest_plot( df=auprc, value.col="auprc", mtext.col=NULL, xlab="AUPRC",
             margins=c(5,6,1,1) )
abline( v=mean_fe$beta[1,1], lty=2, lwd=1 )


#-------------------------------------------------------------------------------
#   Plot features v. residuals
#-------------------------------------------------------------------------------

# Plot features v. residuals
plot_residuals <- function( data, model, variable ){
  
  # Get LASSO deviance residuals
  preds  <- model$preds$recal
  resids <- -2 * ( ( as.integer(model$preds$causal)  * log(preds) ) + 
                   ( as.integer(!model$preds$causal) * log(1-preds) ) )
  
  # Plot
  plot( data[[variable]], resids,
        xlim=quantile( data[[variable]], probs=c( 0.05, 0.95 ) ),
        ylim=quantile( resids, probs=c( 0.05, 0.95 ) ) )
  
  # Lines
  abline( h=0, lwd=3, col="grey" )
  lines( lowess( data[[variable]], resids ), lwd=3, col="orange2" )
}
plot_residuals( data=tr, model=bgl_las, variable="pops_glo" )
plot_residuals( data=tr, model=bgl_las, variable="pops_rel" )
plot_residuals( data=tr, model=bgl_las, variable="dist_gene_rel" )
plot_residuals( data=tr, model=bgl_las, variable="magma_glo" )
plot_residuals( data=tr, model=bgl_las, variable="magma_rel" )
plot_residuals( data=tr, model=bgl_las, variable="coding_glo" )
plot_residuals( data=tr, model=bgl_las, variable="coding_rel" )


#-------------------------------------------------------------------------------
#   In-sample deviance explained and AIC for each model
#-------------------------------------------------------------------------------

de_glm  <- c( devexp_p_glm, devexp_a_glm, devexp_f_glm, devexp_o_glm, devexp_s_glm )
aic_glm <- c( AIC(p_glm), AIC(a_glm), AIC(f_glm), AIC(o_glm), AIC(s_glm) )
names(de_glm) <- c( "Published BIL", "+ extra BIL", "+ global + relative", 
                    "- LRT-pruned", "- more" )
names(aic_glm) <- names(de_glm)
par( mar=c(7,5,1,1 ) )
cols <- viridis( length(de_glm) )
barplot( de_glm,  las=2, col=cols )
barplot( aic_glm, las=2, col=cols )
par( mar=c(5,5,1,1 ) )

# Test set AUPRC
au1 <- c( f_pr1["GLM"], o_pr1["GLM"], b_pr1["GLM"],
          f_pr1["LASSO"], o_pr1["LASSO"], b_pr1["LASSO"],
          f_pr1["XGBoost"], o_pr1["XGBoost"], b_pr1["XGBoost"] )
au2 <- c( f_pr2["GLM"],     o_pr2["GLM"],     b_pr2["GLM"],
          f_pr2["LASSO"],   o_pr2["LASSO"],   b_pr2["LASSO"],
          f_pr2["XGBoost"], o_pr2["XGBoost"], b_pr2["XGBoost"] )
names(au1) <- c( "GLM, full model",     "GLM, reduced model",     "GLM, basic model",
                 "LASSO, full model",   "LASSO, reduced model",   "LASSO, basic model",
                 "XGBoost, full model", "XGBoost, reduced model", "XGBoost, basic model" )
names(au2) <- names(au1)
par( mar=c(11,5,3,1 ) )
cols <- viridis( length(au1) )
barplot( au1, las=2, col=cols, ylim=c( 0, max(au2) ),
         main="WITHOUT protein attenuation" )
abline( h=max(au1), lty=2, lwd=2, col="grey" )
abline( h=max(au2), lty=2, lwd=2, col="grey" )
barplot( au2, las=2, col=cols, main="WITH protein attenuation" )
abline( h=max(au1), lty=2, lwd=2, col="grey" )
abline( h=max(au2), lty=2, lwd=2, col="grey" )
par( mar=c(5,5,1,1 ) )


#-------------------------------------------------------------------------------
#   P(causal) distribution (POPS+local v. not)
#-------------------------------------------------------------------------------

# Get predictions
ev <- te
pred <- predict( object=s_glm, newdata=ev, se=TRUE )

# Format fitted values
y    <- pred$fit
ci   <- pred$se.fit * qnorm( 0.025, lower.tail=FALSE )
ymin <- y - ci
ymax <- y + ci
ev$causal_p  <- logistic(y)
ev$causal_lo <- logistic(ymin)
ev$causal_hi <- logistic(ymax)
ev <- ev[ order(-ev$causal_p) , ]

# Format
bcols <- c( "causal", "causal_p", "distance_genebody", "pops_glo", "magma_glo", 
            "coding_prob", "ngenes_nearby", "dist_gene_rel", "pops_rel", 
            "magma_rel", "pops_and_local",
            "trait", "gene", "region", "cs_id" )
gcols <- c( "causal", "causal_p", "dist",              "pops_glo", "magma_glo", 
            "coding_pip",  "ngenes",        "dist_gene_rel", "pops_rel", 
            "magma_rel", "pops_and_local",
            "trait", "gene", "region", "cs_id" )
ev <- ev[,..bcols]
names(ev) <- gcols

# Plot
hist( ev$causal_p[ ev$pops_and_local ], breaks=seq(0,1,0.05), col="#70AD47", 
      main="PoPS + local", xlab="P(causal)", las=1 )
hist( ev$causal_p[ !ev$pops_and_local ], breaks=seq(0,1,0.05), col="#70AD47", 
      main="Not PoPS + local", xlab="P(causal)", las=1 )


#-------------------------------------------------------------------------------
#   Make a decision tree to see how to best improve on POPS + nearest gene
#-------------------------------------------------------------------------------

# Create a decision tree model specification
tree_spec <- decision_tree( tree_depth=3 ) %>%
  set_engine("rpart") %>%
  set_mode("classification")

# Fit the model to the training data
d_form <-  as.formula( paste( "factor(causal) ~", 
                              paste( b_cols[-1], collapse=" + " ) ) )
tree_fit <- tree_spec %>%
  fit( d_form, data = tr )

# Plot the decision tree
rpart.plot( tree_fit$fit, type = 4, under = TRUE, box.palette = "auto")


#-------------------------------------------------------------------------------
#   Tune XBGoost hyperparameters
#-------------------------------------------------------------------------------

# Extract tuned XGB hyperparameter values
model <- fl_xgb
par_names <- c( "nrounds", "eta", "max_depth", "min_child_weight", "gamma",
                "subsample", "colsample_bytree" )
hyp0 <- list()
for( i in seq_along(model$model) ){
  hyp0[[ names(model$model)[i] ]] <- unlist( model$model[[i]]$learner$par.vals[par_names] )
}
hyp <- as.data.frame( do.call( rbind, hyp0 ) )

# Plot
corrplot( cor(hyp), order="hclust" )
hist( hyp[[1]], col="steelblue", breaks=10, main=par_names[1], las=1, xlim=c( 100, 250 ) )
hist( hyp[[2]], col="steelblue", breaks=10, main=par_names[2], las=1, xlim=c( 0.01, 0.15 ) )
hist( hyp[[3]], col="steelblue", breaks=10, main=par_names[3], las=1, xlim=c( 1, 6 ) )
hist( hyp[[4]], col="steelblue", breaks=10, main=par_names[4], las=1, xlim=c( 0, 5 ) )
hist( hyp[[5]], col="steelblue", breaks=10, main=par_names[5], las=1, xlim=c( 0, 8 ) )
hist( hyp[[6]], col="steelblue", breaks=10, main=par_names[6], las=1, xlim=c( 0.5, 1 ) )
hist( hyp[[7]], col="steelblue", breaks=10, main=par_names[7], las=1, xlim=c( 0.1, 0.6 ) )
plot( hyp$nrounds, hyp$eta, pch=19, col="steelblue", cex=2, las=1,
      xlab="Number of rounds", ylab="Learning rate (eta)" )


#-------------------------------------------------------------------------------
#   Distribution of causal gene distances
#-------------------------------------------------------------------------------

# Load libraries
library(fitdistrplus)
library(logspline)

# Which distributions is this likely to fit?
adj_factor <- 0
tssc <- log10( tr$distance_tss[ tr$causal] + adj_factor )
descdist( data=tssc, discrete=FALSE )

# Find parameters that best fit the data (and compare fit AICs): eQTL TSS
fit_w <- fitdist( data=tssc, distr="weibull" )
fit_g <- fitdist( data=tssc, distr="gamma" )
fit_w$aic; fit_g$aic

# Inspect: eQTL TSS
plot(fit_w)
x <- seq( 0, max(tssc), length=100 )
plot( density(tssc) )
lines( x=x, y=dweibull( x=x, shape=fit_w$estimate["shape"], 
                        scale=fit_w$estimate["scale"] ), col="tomato2" )

# Our data
# What is the relative density of the fitted Weibull distribution?
r_dens <- dweibull( x     = log10( seq( 1e5, 5e5, 1e5 ) + adj_factor ), 
                    shape = fit_w$estimate["shape"], 
                    scale = fit_w$estimate["scale"] )
r_dens[1] / r_dens

# Fauman and Hyde 2022 analysis of eQTLGen
# What is the relative density of the fitted Weibull distribution?
r_dens <- dweibull( x     = log10( seq( 1e5, 5e5, 1e5 ) + adj_factor ), 
                    shape = 7.375, 
                    scale = 4.7 )
r_dens[1] / r_dens

# Probability that an eQTL gene is >300kb away
pweibull( q = log10( seq( 1e5, 5e5, 1e5 ) ), 
          shape = 7.375, 
          scale = 4.7 )
pweibull( q = log10( seq( 1e5, 5e5, 1e5 ) ), 
          shape = fit_w$estimate["shape"], 
          scale = fit_w$estimate["scale"] )


#-------------------------------------------------------------------------------
#   Previous Figure 1: AUPRC barplots
#-------------------------------------------------------------------------------

# Figure 1: set up
fig_dir   <- file.path( maindir, "figures" )
fig1_file <- file.path( fig_dir, "figure1.tif" )
dir.create( path=fig_dir, showWarnings=FALSE, recursive=TRUE )
tiff( filename=fig1_file, width=600*4, height=300*4, res=75*4 )
par( mfrow=c(1,2) )
par( mar=c(7,5,1,1) )

# Plot Figure 1A
f_pr <- rbind( fl_xgb_pr, fl_las_pr )
dimnames(f_pr)[[1]] <- c( "XGBoost", "LASSO" )
bp1 <- barplot( height=f_pr[,"auprc"], las=2, 
                ylab="AUPRC (±95% CI)", ylim=c( 0, 0.7 ),
                col=brewer.pal( n=6, name="Greens" )[1:2] )
for( i in seq_along(bp1) ){
  lines( x = c( bp1[i],          bp1[i] ), lwd=1.5,
         y = c( f_pr[ i, "lo" ], f_pr[ i, "hi" ] ) )
}

# Figure 1B
auprcs <- rbind( fl_las_pr, bl_las_pr, bgl_las_pr, bgl_las_recal_pr )
dimnames(auprcs)[[1]] <- c( "Full", "Basic", "+GLCs", "CALDERA" )
bp2 <- barplot( height=auprcs[,1], las=2, 
                ylab="AUPRC (±95% CI)", ylim=c( 0, 0.7 ),
                col=brewer.pal( n=6, name="Greens" )[2:5] )
for( i in seq_along(bp2) ){
  lines( x = c( bp2[i],      bp2[i] ), lwd=1.5,
         y = c( auprcs[i,2], auprcs[i,3] ) )
}

# Figure 1: wrap up
par( mfrow=c(1,1) )
dev.off()


#-------------------------------------------------------------------------------
#   Previous Figure S1: Predictions before and after GLCs
#-------------------------------------------------------------------------------

# How do predictions differ before/after adding GLCs?
fig_s1_file <- file.path( fig_dir, "figure_s1.tif" )
tiff( filename=fig_s1_file, width=300*4, height=300*4, res=75*4 )
par( mar=c(5,5,1,1) )
plot( bl_las$preds$pred, bgl_las$preds$pred, las=1, lwd=1.5,
      col=brewer.pal( n=6, name="Greens" )[5],
      xlab="Causal probability without GLCs", 
      ylab="Causal probability with GLCs" )
abline( a=0, b=1, lty=1, lwd=2 )
dev.off()

# Other
lines( lowess( bl_las$preds$pred, bgl_las$preds$pred, f=0.05 ), lty=2, lwd=2 )

# Repeat, but on the logit scale
plot( logit10(bl_las$preds$pred), logit10(bgl_las$preds$pred), las=1, lwd=1.5,
      col=brewer.pal( n=6, name="Greens" )[5],
      xlab="Causal probability without GLCs", 
      # xlim=c( logit10(0.01), logit10(0.95) ),
      # ylim=c( logit10(0.01), logit10(0.95) ),
      ylab="Causal probability with GLCs" )
abline( a=0, b=1, lty=1, lwd=2 )
lines( lowess( logit10(bl_las$preds$pred), logit10(bgl_las$preds$pred), f=0.05 ), 
       lty=2, lwd=2 )

# Compare coefficients with/without GLC correction
coef1 <- as.data.frame( as.matrix( coef( bg_las, s="lambda.min" ) ) )
coef2 <- as.data.frame( as.matrix( coef( b_las,  s="lambda.min" ) ) )
shared_feats <- intersect( row.names(coef1)[ coef1$s1 !=0 ], 
                           row.names(coef2)[ coef2$s1 !=0 ] )[-1]
df <- data.frame( feature = shared_feats,
                  without = coef2$s1[ match( shared_feats, row.names(coef1) ) ],
                  with    = coef1$s1[ match( shared_feats, row.names(coef2) ) ] )
plot( df$without, df$with, las=1, cex=1.5, col="steelblue", 
      xlab="Without", ylab="With" )
abline( a=0, b=1, lty=2, lwd=2 )


#-------------------------------------------------------------------------------
#   Simulate the effect of covariate correction: 
#   realistic feature values but no cross-feature correlation
#-------------------------------------------------------------------------------

# Generate a training dataset that contains biases
# Sample feature values for each gene randomly from real data
n <- 1e5
data_tr <- data.table( 
  pops_glo            = sample( x=tr$pops_glo,            size=n, replace=TRUE ),
  dist_gene_glo       = sample( x=tr$dist_gene_glo,       size=n, replace=TRUE ),
  coding_glo          = sample( x=tr$coding_glo,          size=n, replace=TRUE ),
  prior_n_genes_locus = sample( x=tr$prior_n_genes_locus, size=n, replace=TRUE ),
  pLI_lt_0.1          = sample( x=tr$pLI_lt_0.1,          size=n, replace=TRUE ),
  cds_bp_log10        = sample( x=tr$cds_bp_log10,        size=n, replace=TRUE ),
  roadmap_bp_log10    = sample( x=tr$roadmap_bp_log10,    size=n, replace=TRUE ) )

# Assign causal gene status based on probability
data_tr$p_causal <- predict( object=bg_glm, newdata=data_tr, type="response" )
data_tr$causal   <- rbinom( n=n, size=1, prob=data_tr$p_causal )
sum(data_tr$causal)

# Train a model that corrects for biases
# Train a model that does not correct for biases
cols_c <- data_tr %>% select( pops_glo:roadmap_bp_log10 ) %>% names()
cols_u <- data_tr %>% select( pops_glo:prior_n_genes_locus ) %>% names()
form_c <- make_formula( lhs="causal", rhs=cols_c )
form_u <- make_formula( lhs="causal", rhs=cols_u )
mod_c <- glm( formula=form_c, family="binomial", data=data_tr )
mod_u <- glm( formula=form_u, family="binomial", data=data_tr )
summary(bg_glm)$coef
summary(mod_c)$coef
summary(mod_u)$coef

# Generate a test dataset without biases
data_te <- data.table( 
  pops_glo            = sample( x=tr$pops_glo,            size=n, replace=TRUE ),
  dist_gene_glo       = sample( x=tr$dist_gene_glo,       size=n, replace=TRUE ),
  coding_glo          = sample( x=tr$coding_glo,          size=n, replace=TRUE ),
  prior_n_genes_locus = sample( x=tr$prior_n_genes_locus, size=n, replace=TRUE ),
  pLI_lt_0.1          = cov_means["pLI_lt_0.1"],
  cds_bp_log10        = cov_means["cds_bp_log10"],
  roadmap_bp_log10    = cov_means["roadmap_bp_log10"] )
data_te

# Generate predictions for the two models
data_te$p_causal_c <- predict( object=mod_c, newdata=data_te, type="response" )
data_te$p_causal_u <- predict( object=mod_u, newdata=data_te, type="response" )
data_te$p_causal_t <- predict( object=bg_glm, newdata=data_te, type="response" )
data_te$causal     <- rbinom( n=n, size=1, prob=data_te$p_causal_t )
mean( abs( data_te$p_causal_c - data_te$p_causal_t ) )
mean( abs( data_te$p_causal_u - data_te$p_causal_t ) )
plot( data_te$p_causal_t[1:1e3], data_te$p_causal_u[1:1e3] )
abline( a=0, b=1, lty=2 )
plot( data_te$p_causal_t[1:1e3], data_te$p_causal_c[1:1e3] )
abline( a=0, b=1, lty=2 )

# Compute AUPRCs for the two models
pr_c <- pr.curve( scores.class0 = data_te[["p_causal_c"]][ data_te$causal == 1 ], 
                  scores.class1 = data_te[["p_causal_c"]][ data_te$causal == 0 ], 
                  curve = TRUE )
pr_u <- pr.curve( scores.class0 = data_te[["p_causal_u"]][ data_te$causal == 1 ], 
                  scores.class1 = data_te[["p_causal_u"]][ data_te$causal == 0 ], 
                  curve = TRUE )
out_c <- aucpr.conf.int.expit( estimate = pr_c$auc.integral,
                               num.pos  = sum(data_te$causal) )
out_u <- aucpr.conf.int.expit( estimate = pr_u$auc.integral,
                               num.pos  = sum(data_te$causal) )
rbind( out_c, out_u )


#-------------------------------------------------------------------------------
#   Simulate the effect of covariate correction: 
#   preserve cross-feature correlation, but draw feature values from MVN
#-------------------------------------------------------------------------------

# Get the mean and covariance matrix for the training data
gcols <- bg_glm$model %>% select(-causal) %>% names()
mu    <- colMeans( x=tr[ , ..gcols ] )
Sigma <- cov( x=tr[ , ..gcols ] )

# Sample from a multivariate normal distribution using the means and SDs of each
# variable, and the correlations between variables
data_tr <- MASS::mvrnorm( n=n, mu=mu, Sigma=Sigma ) %>% as.data.table()

# Assign causal gene status based on probability
data_tr$p_causal <- predict( object=bg_glm, newdata=data_tr, type="response" )
data_tr$causal   <- rbinom( n=n, size=1, prob=data_tr$p_causal )
sum(data_tr$causal)

# Train a model that corrects for biases
# Train a model that does not correct for biases
cols_c <- data_tr %>% select( pops_glo:roadmap_bp_log10 ) %>% names()
cols_u <- data_tr %>% select( pops_glo:prior_n_genes_locus ) %>% names()
form_c <- make_formula( lhs="causal", rhs=cols_c )
form_u <- make_formula( lhs="causal", rhs=cols_u )
mod_c <- glm( formula=form_c, family="binomial", data=data_tr )
mod_u <- glm( formula=form_u, family="binomial", data=data_tr )
summary(bg_glm)$coef
summary(mod_c)$coef
summary(mod_u)$coef

# Generate a test dataset without biases
mu_u    <- colMeans( x=tr[ , ..cols_u ] )
Sigma_u <- cov( x=tr[ , ..cols_u ] )
data_te <- MASS::mvrnorm( n=n, mu=mu_u, Sigma=Sigma_u ) %>% as.data.table()
data_te$pLI_lt_0.1       <- cov_means["pLI_lt_0.1"]
data_te$cds_bp_log10     <- cov_means["cds_bp_log10"]
data_te$roadmap_bp_log10 <- cov_means["roadmap_bp_log10"]
data_te

# Generate predictions for the two models
data_te$p_causal_c <- predict( object=mod_c, newdata=data_te, type="response" )
data_te$p_causal_u <- predict( object=mod_u, newdata=data_te, type="response" )
data_te$p_causal_t <- predict( object=bg_glm, newdata=data_te, type="response" )
data_te$causal     <- rbinom( n=n, size=1, prob=data_te$p_causal_t )
mean( abs( data_te$p_causal_c - data_te$p_causal_t ) )
mean( abs( data_te$p_causal_u - data_te$p_causal_t ) )
plot( data_te$p_causal_t[1:1e3] |> logit(), data_te$p_causal_u[1:1e3] |> logit() )
abline( a=0, b=1, lty=2 )
plot( data_te$p_causal_t[1:1e3] |> logit(), data_te$p_causal_c[1:1e3] |> logit() )
abline( a=0, b=1, lty=2 )

# Compute AUPRCs for the two models
pr_c <- pr.curve( scores.class0 = data_te[["p_causal_c"]][ data_te$causal == 1 ], 
                  scores.class1 = data_te[["p_causal_c"]][ data_te$causal == 0 ], 
                  curve = TRUE )
pr_u <- pr.curve( scores.class0 = data_te[["p_causal_u"]][ data_te$causal == 1 ], 
                  scores.class1 = data_te[["p_causal_u"]][ data_te$causal == 0 ], 
                  curve = TRUE )
out_c <- aucpr.conf.int.expit( estimate = pr_c$auc.integral,
                               num.pos  = sum(data_te$causal) )
out_u <- aucpr.conf.int.expit( estimate = pr_u$auc.integral,
                               num.pos  = sum(data_te$causal) )
rbind( out_c, out_u )


#-------------------------------------------------------------------------------
#   Brier scores
#-------------------------------------------------------------------------------

# Function
brier_fun <- function( true, pred ){
  mean( ( pred - true )^2 )
}

# LOTO
brier_fun( true = bgl_glm$preds$causal,
           pred = bgl_glm$preds$scaled )

# Open Targets
brier_fun( true = lg_pr$data$causal,
           pred = lg_pr$data$CALDERA )
brier_fun( true = lg_pr$data$causal,
           pred = lg_pr$data$L2G )
brier_fun( true = lg_pr$data$causal,
           pred = lg_pr$data$FLAMES )
brier_fun( true = lg_pr$data$causal,
           pred = lg_pr$data$cS2G )

# ExWAS
brier_fun( true = ex_pr$data$causal,
           pred = ex_pr$data$CALDERA )
brier_fun( true = ex_pr$data$causal,
           pred = ex_pr$data$L2G )
brier_fun( true = ex_pr$data$causal,
           pred = ex_pr$data$FLAMES )
brier_fun( true = ex_pr$data$causal,
           pred = ex_pr$data$cS2G )

# MT
brier_fun( true = mt_pr$data$causal,
           pred = mt_pr$data$CALDERA )
brier_fun( true = mt_pr$data$causal,
           pred = mt_pr$data$L2G )
brier_fun( true = mt_pr$data$causal,
           pred = mt_pr$data$FLAMES )
brier_fun( true = mt_pr$data$causal,
           pred = mt_pr$data$cS2G )


#-------------------------------------------------------------------------------
#   Per-trait AUPRCs
#-------------------------------------------------------------------------------

# For this code to work, you would have to add it to the end of the 
# benchmarking_auprc() function
if( bench == "3MT" )
  pattern <- "-[[:digit:]]+_[[:digit:]]+_[ACGT]+_[ACGT]+_credset.cred1"
if( bench == "Open Targets" )
  pattern <- "-[[:digit:]]+_[[:digit:]]+_[ACGT]+_[ACGT]+-credset-intersection.cred1"
lg3$trait <- lg3$CS_ID %>%
  sub( pattern=pattern, replacement="" )
pr_pt0 <- list()
for( i in unique(lg3$trait) ){
  cal_pr <- plot_mod_pr( pred_col="CALDERA", data=lg3[ lg3$trait == i ], make_plot=FALSE )
  l2g_pr <- plot_mod_pr( pred_col="L2G",     data=lg3[ lg3$trait == i ], make_plot=FALSE )
  fla_pr <- plot_mod_pr( pred_col="FLAMES",  data=lg3[ lg3$trait == i ], make_plot=FALSE )
  s2g_pr <- plot_mod_pr( pred_col="cS2G",    data=lg3[ lg3$trait == i ], make_plot=FALSE )
  m <- rbind( cal_pr, l2g_pr, fla_pr, s2g_pr )
  m$trait <- i
  m$method <- row.names(m)
  pr_pt0[[i]] <- m
}
pr_pt <- do.call( rbind, pr_pt0 )
pr_pt

# Plot: CALDERA vs. L2G
plot( jitter( pr_pt$auprc[pr_pt$method == "L2G" ] ),
      jitter( pr_pt$auprc[pr_pt$method == "CALDERA" ] ),
      xlim=c(0,1), ylim=c(0,1), las=1, ylab="CALDERA", xlab="L2G" )
abline( a=0, b=1, lty=2 )
hist( pr_pt$auprc[pr_pt$method == "CALDERA" ] - pr_pt$auprc[pr_pt$method == "L2G" ] )
plot( density( pr_pt$auprc[pr_pt$method == "CALDERA" ] - pr_pt$auprc[pr_pt$method == "L2G" ] ) )
t.test( x = pr_pt$auprc[pr_pt$method == "CALDERA" ],
        y = pr_pt$auprc[pr_pt$method == "L2G" ],
        paired = TRUE )

# Plot: CALDERA vs. FLAMES
plot( jitter( pr_pt$auprc[pr_pt$method == "FLAMES" ] ),
      jitter( pr_pt$auprc[pr_pt$method == "CALDERA" ] ),
      xlim=c(0,1), ylim=c(0,1), las=1, ylab="CALDERA", xlab="FLAMES" )
abline( a=0, b=1, lty=2 )
hist( pr_pt$auprc[pr_pt$method == "CALDERA" ] - pr_pt$auprc[pr_pt$method == "FLAMES" ] )
plot( density( pr_pt$auprc[pr_pt$method == "CALDERA" ] - pr_pt$auprc[pr_pt$method == "FLAMES" ] ) )
t.test( x = pr_pt$auprc[pr_pt$method == "CALDERA" ],
        y = pr_pt$auprc[pr_pt$method == "FLAMES" ],
        paired = TRUE )

# Plot: CALDERA vs. cS2G
plot( jitter( pr_pt$auprc[pr_pt$method == "cS2G" ] ),
      jitter( pr_pt$auprc[pr_pt$method == "CALDERA" ] ),
      xlim=c(0,1), ylim=c(0,1), las=1, ylab="CALDERA", xlab="FLAMES" )
abline( a=0, b=1, lty=2 )
hist( pr_pt$auprc[pr_pt$method == "CALDERA" ] - pr_pt$auprc[pr_pt$method == "cS2G" ] )
plot( density( pr_pt$auprc[pr_pt$method == "CALDERA" ] - pr_pt$auprc[pr_pt$method == "cS2G" ] ) )
t.test( x = pr_pt$auprc[pr_pt$method == "CALDERA" ],
        y = pr_pt$auprc[pr_pt$method == "cS2G" ],
        paired = TRUE )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Done
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------


















