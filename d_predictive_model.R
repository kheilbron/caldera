
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
  mod3$col <- brewer.pal( n=3, name="Greens" )[3]
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
                            replacement="Coding PIP*" )
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
                        makeIntegerParam(  "nrounds",          lower = 75L,  upper = 175L ), 
                        makeNumericParam(  "eta",              lower = 0.01, upper = 0.5 ), 
                        makeIntegerParam(  "max_depth",        lower = 1L,   upper = 2L ), 
                        makeNumericParam(  "min_child_weight", lower = 0,    upper = 8 ), 
                        makeNumericParam(  "gamma",            lower = 0,    upper = 6 ), 
                        makeNumericParam(  "subsample",        lower = 0.3,  upper = 1 ), 
                        makeNumericParam(  "colsample_bytree", lower = 0.1,  upper = 1 ) )
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
  cv.glmnet( x=as.matrix( data[ , ..pred_cols ] ), 
             y=data[["causal"]], family="binomial", foldid=data$foldid )
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
  out <- list( trait    = unique(data$trait), 
               model    = list(), 
               preds    = data.table( trait      = data$trait,
                                      causal     = data$causal,
                                      pred       = NA,
                                      recal = NA,
                                      scaled     = NA ) )
  
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
    
    # Train a LASSO model to recalibrate model outputs
    cal_mod <- recalibration_model( data=preds_tr )
    
    # Apply LASSO model to get recalibrated predictions
    pred_cols <- c( "global", "relative" )
    preds_tr$modeled <- predict( object=cal_mod, s="lambda.min", type="response",
                                 newx=as.matrix( preds_tr[ , ..pred_cols ] ) )
    preds_te$modeled <- predict( object=cal_mod, s="lambda.min", type="response",
                                 newx=as.matrix( preds_te[ , ..pred_cols ] ) )
    
    
    #---------------------------------------------------------------------------
    #   Save outputs
    #---------------------------------------------------------------------------
    
    # Model
    out$model[[i]] <- l_mod
    
    # Predicted probabilities
    out$preds$pred[       out$preds$trait == i ] <- preds_te$original
    out$preds$scaled[     out$preds$trait == i ] <- preds_te$scaled
    out$preds$recal[ out$preds$trait == i ] <- preds_te$modeled
    out$preds$bil[        out$preds$trait == i ] <- preds_te$best
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
  }else if( type == "scaled" ){
    prc <- pr.curve( scores.class0 = loto_obj$preds$scaled[  loto_obj$preds$causal ], 
                       scores.class1 = loto_obj$preds$scaled[ !loto_obj$preds$causal ], 
                       curve = TRUE )
  }
  
  # Plot
  if(show_plot){
    par( mar=c(5,5,4,1) )
    plot( prc, scale.color=color, las=1, legend=legend )
  }
  
  # Return AUPRC and 95% CI
  ci <- aucpr.conf.int.expit( estimate = prc$auc.integral, 
                              num.pos  = sum(loto_obj$preds$causal) )
  names(ci) <- c( "lo", "hi" )
  out <- c( auprc=prc$auc.integral, ci )
  return(out)
}


#-------------------------------------------------------------------------------
#   plot_logOR_relationship_smooth: Smoothed P(causal) v. a variable
#-------------------------------------------------------------------------------

plot_logOR_relationship_smooth <- function( data, varname ){
  
  # Format data
  data <- as.data.frame(data)
  data <- data.frame( x=data[[varname]], y=data$causal )
  
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
  p1 <- ggplot( data, aes(x=x) ) +
    geom_histogram() +
    ggtitle(varname)
  p2 <- ggplot( loess.DF, aes( x=x, y=y) ) +
    geom_smooth( aes_auto(loess.DF), data=loess.DF, stat="identity", col="#70AD47" ) +
    scale_y_continuous( trans="logit",
                        breaks=c( 0.01, 0.05, seq(0.1,0.9,0.1), 0.95, 0.99 ) ) +
    geom_hline( yintercept=0.01 ) +
    geom_hline( yintercept=0.99 ) +
    theme_bw()
  p1 + p2 + plot_layout( nrow=2, heights=c(1,3) )
}


#-------------------------------------------------------------------------------
#   plot_logOR_relationship_model:  Modeled P(causal) v. a variable
#-------------------------------------------------------------------------------

plot_logOR_relationship_model <- function( data, model, var_prefix, 
                                           qmin=0.1, qmax=0.9, 
                                           g_trans_fun=NULL ){
  
  # Extract features from the model
  feat_cols  <- model$glmnet$beta@Dimnames[[1]]
  focal_cols <- grep( pattern=paste0( "^", var_prefix ), x=feat_cols, value=TRUE )
  g_focal    <- grep( pattern="_glo$", x=focal_cols, value=TRUE )
  r_focal    <- grep( pattern="_rel$", x=focal_cols, value=TRUE )
  b_focal    <- grep( pattern="_bil$", x=focal_cols, value=TRUE )
  
  # Format data
  xycols <- c( "causal", feat_cols )
  data <- data[ , ..xycols ]
  data <- as.data.frame( lapply( X=data[ , ..xycols ], FUN=as.numeric ) )
  
  # Create new, evenly-spaced values for the focal global feature
  n      <- 100
  xrange <- quantile( data[[g_focal]], probs=c( qmin, qmax ) )
  xseq   <- seq( from=xrange[1], to=xrange[2], length=n )
  
  # Construct a dataset to make predictions in
  # Set all global features to their mean from the training dataset
  # Set all relative and BIL features to their maximum from the training dataset
  feat_means <- apply( data[,-1], 2, mean )
  feat_maxs  <- apply( data[,-1], 2, max )
  feat_mixs  <- ifelse( grepl( pattern="_rel$|_bil$", x=names(feat_means) ),
                        feat_maxs, feat_means )
  mat <- matrix( rep( feat_mixs, n ), 
                 nrow=n, ncol=length(feat_mixs), byrow=TRUE )
  nx0 <- as.data.frame(mat)
  names(nx0) <- names(feat_means)
  
  # Except for the focal global feature, which is from the previous step
  # And the focal BIL feature, which is set to 0
  # And the focal relative feature, which is set to its mean
  nx0[[g_focal]] <- xseq
  nx0[[b_focal]] <- 0
  nx0[[r_focal]] <- feat_means[r_focal]
  nx0$Relative   <- "Mean"
  
  # Copy the previous dataset, but now set the focal relative feature to its...
  # ...best value and the focal BIL feature to 1, rbind
  nx1 <- nx0
  nx1[[r_focal]] <- feat_maxs[r_focal]
  nx1[[b_focal]] <- 1
  nx1$Relative <- "Best"
  nx <- rbind( nx0, nx1)
  
  # Get predictions
  pred <- predict( model, newx=as.matrix( nx[ , feat_cols ] ),
                   s="lambda.min", type="response" )
  nx$pred <- as.vector(pred)
  df <- data.frame( x=nx[[g_focal]], y=nx$pred, Relative=nx$Relative )
  
  # If necessary, transform the scale of the global feature
  if( !is.null(g_trans_fun) ){
    df$x <- g_trans_fun(df$x)
  }
  
  # Clean up labels
  clean_label <- g_focal
  clean_label <- sub( x=clean_label, pattern="_glo",      replacement="" )
  clean_label <- sub( x=clean_label, pattern="pops",      replacement="PoPS" )
  clean_label <- sub( x=clean_label, pattern="magma",     replacement="MAGMA z-score" )
  clean_label <- sub( x=clean_label, pattern="dist_gene", replacement="Distance (bp)" )
  clean_label <- sub( x=clean_label, pattern="coding",    replacement="Coding PIP" )
  
  # Prepare the histogram
  hdf <- data.frame( x=data[[g_focal]] )
  hdf <- hdf[ hdf$x >= xrange[1] & hdf$x <= xrange[2] , , drop=FALSE ]
  if( !is.null(g_trans_fun) ){
    hdf$x <- g_trans_fun(hdf$x)
  }
  
  # Plot
  p1 <- ggplot( hdf, aes(x=x) ) +
    geom_histogram() +
    labs( y="Count", x="" ) +
    theme_bw() + 
    ggtitle(clean_label)
  p2 <- ggplot( df, aes( x=x, y=y, color=Relative ) ) +
    geom_line( aes( y=y ), linewidth=1 ) + 
    scale_y_continuous( trans="logit",
                        breaks=c( 0.05, seq(0.1,0.7,0.1) ) ) +
    geom_hline( yintercept=0.05, col="lightgrey" ) +
    geom_hline( yintercept=0.7,  col="lightgrey" ) +
    labs( y="Predicted causal probability", x=paste0( clean_label, ", global" ) ) +
    scale_color_manual( values = brewer.pal( n=5, name="Greens" )[ c( 4, 3 ) ] ) +
    theme_bw() + 
    theme( legend.position="bottom" )
  p1 + p2 + plot_layout( nrow=2, heights=c(1,3) )
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

# Best-in-locus features
a_cols <- c( "causal", 
             "pops_bil", 
             "dist_gene_bil", "dist_tss_bil",
             "magma_bil",
             "coding_bil",
             "twas_bil", 
             "corr_liu_bil", "corr_and_bil", "corr_uli_bil",
             "pchic_jung_bil", "pchic_jav_bil", 
             "clpp_bil", 
             "smr_bil",
             "abc_bil",
             "depict_z_bil",
             "netwas_score_bil", "netwas_bon_score_bil" )

# Global features
glo_cols <- c( "causal", "pops_glo", 
               "dist_gene_glo", "dist_tss_glo",
               "magma_glo",
               "coding_glo",
               "twas_glo",
               "corr_liu_glo", "corr_and_glo", "corr_uli_glo", 
               "pchic_jung_glo", "pchic_jav_glo", 
               "clpp_glo", "smr_glo", "abc_glo",
               "depict_z_glo", "netwas_score_glo", "netwas_bon_score_glo" )

# Relative features
rel_cols <- c( "causal", "pops_rel", "dist_gene_rel", "dist_tss_rel", 
               "magma_rel", "coding_rel", "twas_rel", 
               "corr_liu_rel", "corr_and_rel", "corr_uli_rel", 
               "pchic_jung_rel", "pchic_jav_rel",  
               "clpp_rel", "smr_rel", "abc_rel",
               "depict_z_rel", "netwas_score_rel", "netwas_bon_score_rel" )

# Features that may capture bias in how causal genes were selected
bias_cols <- c( "pLI_log10",     "pLI_lt_0.1",   "pLI_lt_0.9",
                "LOEUF",         "hs_log10",
                "gene_bp_log10", "cds_bp_log10",
                "abc_bp_log10",  "roadmap_bp_log10",
                "pritchard_miss" )

# Gene-level features
# glc_cols  <- c( "ABC_count",         "Roadmap_count",
#                 "TF",                "promoter_count_log10",
#                 "connect_decile",    "connectedness",
#                 "PPI_degree_decile", "PPI_degree_cat" )
glc_cols  <- NULL

# Full feature set
f_cols <- unique( c( a_cols, glo_cols, rel_cols, "prior_n_genes_locus" ) )

# Basic feature set
pattern <- "^causal$|^pops_|^dist_gene_|^magma_|^coding_|^prior"
b_cols  <- grep( pattern=pattern, x=f_cols, value=TRUE )

# Full + GLC feature set
fg_cols <- c( f_cols, bias_cols, glc_cols )

# Basic + GLC feature set
bg_cols <- c( b_cols, bias_cols, glc_cols )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Full model: best-in-locus + global + relative
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   LOTO
#-------------------------------------------------------------------------------

# Raw feature correlations
corrplot( cor( tr[ , ..f_cols ][,-1] ), order="hclust" )

# Run LOTO
fl_las <- loto( data=tr, method="lasso2", feat_cols=f_cols )
# fl_xgb <- loto( data=tr, method="xgb",    feat_cols=f_cols, maxit=100 )
# saveRDS( object=fl_xgb, file=file.path( maindir, "loto/loto_full_xgb.rds" ) )
fl_xgb  <- readRDS( file=file.path( maindir, "loto/loto_full_xgb.rds" ) )


#-------------------------------------------------------------------------------
#   LASSO
#-------------------------------------------------------------------------------

# Make LASSO model
f_las <- cv.glmnet( x=as.matrix( tr[ , ..f_cols ] )[,-1], 
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
bl_las <- loto( data=tr, method="lasso2", feat_cols=b_cols )
# bl_xgb <- loto( data=tr, method="xgb",    feat_cols=b_cols, maxit=100 )
# saveRDS( object=bl_xgb, file=file.path( maindir, "loto/loto_basic_xgb.rds" ) )
bl_xgb  <- readRDS( file=file.path( maindir, "loto/loto_basic_xgb.rds" ) )


#-------------------------------------------------------------------------------
#   LASSO
#-------------------------------------------------------------------------------

# Make LASSO model
b_las <- cv.glmnet( x=as.matrix( tr[ , ..b_cols ] )[,-1], 
                    foldid=as.integer( as.factor(tr$trait) ),
                    y=tr[["causal"]], family="binomial" )


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
bgl_las <- loto( data=tr, method="lasso2", feat_cols=bg_cols, bias_cols=bias_cols )
bgl_xgb <- bl_xgb


#-------------------------------------------------------------------------------
#   LASSO
#-------------------------------------------------------------------------------

# Run
bg_las <- cv.glmnet( x=as.matrix( tr[ , ..bg_cols ] )[,-1], 
                     foldid=as.integer( as.factor(tr$trait) ),
                     y=tr[["causal"]], family="binomial" )

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
#   Train and save a recalibration model
#-------------------------------------------------------------------------------

# Train a recalibration model
recal_preds <- recalibrate_preds( df=tr, model=bg_las, method="lasso2", 
                                  bias_cols=bias_cols )
recal_mod   <- recalibration_model( data=recal_preds )
recal_preds$recal <- predict( object=recal_mod, s="lambda.min", type="response",
                              newx=as.matrix( 
                                recal_preds[ , c( "global", "relative" ) ] ) )

# Save the LASSO1 model and the calibration model
# saveRDS( object=bg_las,    file=file.path( maindir, "bg_las.rds" ) )
# saveRDS( object=recal_mod, file=file.path( maindir, "bg_las_recal_mod.rds" ) )

# In-sample AUPRC
prc <- pr.curve( scores.class0 = recal_preds$recal[  recal_preds$causal ], 
                 scores.class1 = recal_preds$recal[ !recal_preds$causal ], 
                 curve = TRUE )
ci <- aucpr.conf.int.expit( estimate = prc$auc.integral, 
                            num.pos  = sum(recal_preds$causal) )
names(ci) <- c( "lo", "hi" )
c( auprc=prc$auc.integral, ci )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Figures
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Figure 1: AUPRC barplots
#-------------------------------------------------------------------------------

# Extract AUPRC and 95% CI for all models
fl_las_pr   <- plot_loto_pr( loto_obj=fl_las,  color=viridis(11), legend=TRUE )
fl_xgb_pr   <- plot_loto_pr( loto_obj=fl_xgb,  color=viridis(11), legend=TRUE )
bl_las_pr   <- plot_loto_pr( loto_obj=bl_las,  color=viridis(11), legend=TRUE )
bl_xgb_pr   <- plot_loto_pr( loto_obj=bl_xgb,  color=viridis(11), legend=TRUE )
bgl_las_pr  <- plot_loto_pr( loto_obj=bgl_las, color=viridis(11), legend=TRUE )
bgl_xgb_pr  <- plot_loto_pr( loto_obj=bgl_xgb, color=viridis(11), legend=TRUE )
bgl_lasc_pr <- plot_loto_pr( loto_obj=bgl_las, color=viridis(11), legend=TRUE, 
                             type="recalibrated" )
bgl_xgbc_pr <- plot_loto_pr( loto_obj=bgl_xgb, color=viridis(11), legend=TRUE,
                             type="recalibrated" )

# Figure 1: set up
fig_dir   <- file.path( maindir, "figures" )
fig1_file <- file.path( fig_dir, "figure1.jpg" )
dir.create( path=fig_dir, showWarnings=FALSE, recursive=TRUE )
jpeg( filename=fig1_file, width=600*4, height=300*4, res=75*4 )
par( mfrow=c(1,2) )
par( mar=c(7,5,1,1) )

# Plot Figure 1A
f_pr <- rbind( fl_xgb_pr, fl_las_pr )
dimnames(f_pr)[[1]] <- c( "XGBoost", "LASSO" )
bp <- barplot( height=f_pr[,"auprc"], las=2, 
               ylab="AUPRC (±95% CI)", ylim=c( 0, 0.7 ),
               col=brewer.pal( n=5, name="Greens" )[1:2] )
for( i in seq_along(bp) ){
  lines( x = c( bp[i],           bp[i] ), lwd=1.5,
         y = c( f_pr[ i, "lo" ], f_pr[ i, "hi" ] ) )
}

# Figure 1B
auprcs <- rbind( fl_las_pr, bl_las_pr, bgl_las_pr, bgl_lasc_pr )
dimnames(auprcs)[[1]] <- c( "Full", "Basic", "+GLCs", "+Recalibrated" )
bp1 <- barplot( height=auprcs[,1], las=2, 
                ylab="AUPRC (±95% CI)", ylim=c( 0, 0.7 ),
                col=brewer.pal( n=5, name="Greens" )[2:5] )
for( i in seq_along(bp1) ){
  lines( x = c( bp1[i],      bp1[i] ), lwd=1.5,
         y = c( auprcs[i,2], auprcs[i,3] ) )
}

# Figure 1: wrap up
par( mfrow=c(1,1) )
dev.off()


#-------------------------------------------------------------------------------
#   Figure 2: Standardized coefficients
#-------------------------------------------------------------------------------

# Forest plot without bias features
fig2_file <- file.path( fig_dir, "figure2.jpg" )
jpeg( filename=fig2_file, width=375*4, height=300*4, res=75*4 )
par( mar=c(5,7,1,1) )
lasso_to_forest( las=bg_las, lambda_type="lambda.min", train=tr, rm_feats=bias_cols )
dev.off()

# Forest plot with bias features
lasso_to_forest( las=bg_las, lambda_type="lambda.min", train=tr )


#-------------------------------------------------------------------------------
#   Figure 3: P(causal) v. feature group
#-------------------------------------------------------------------------------

# Plot the effect of selected features on P(causal)
fig3a_file <- file.path( fig_dir, "figure3a.jpg" )
jpeg( filename=fig3a_file, width=300*4, height=400*4, res=75*4 )
plot_logOR_relationship_model( data=tr, model=bg_las, var_prefix="pops", 
                               qmin=0.05, qmax=0.95, g_trans_fun=NULL )
dev.off()

# Plot the effect of selected features on P(causal)
fig3b_file <- file.path( fig_dir, "figure3b.jpg" )
jpeg( filename=fig3b_file, width=300*4, height=400*4, res=75*4 )
plot_logOR_relationship_model( data=tr, model=bg_las, var_prefix="magma", 
                               qmin=0.05, qmax=0.95, g_trans_fun=NULL )
dev.off()

fig3c_file <- file.path( fig_dir, "figure3c.jpg" )
jpeg( filename=fig3c_file, width=300*4, height=400*4, res=75*4 )
plot_logOR_relationship_model( data=tr, model=bg_las, var_prefix="dist_gene", 
                               qmin=0.05, qmax=0.95, 
                               g_trans_fun=function(x) 10^-x - 1e3 )
dev.off()

fig3d_file <- file.path( fig_dir, "figure3d.jpg" )
jpeg( filename=fig3d_file, width=300*4, height=400*4, res=75*4 )
plot_logOR_relationship_model( data=tr, model=bg_las, var_prefix="coding", 
                               qmin=0, qmax=1, g_trans_fun=logistic10 )
dev.off()


#-------------------------------------------------------------------------------
#   Figure 4: Calibration plots
#-------------------------------------------------------------------------------

# Figure 4
fig4_file <- file.path( fig_dir, "figure4.jpg" )
jpeg( filename=fig4_file, width=600*4, height=300*4, res=75*4 )
p1 <- cal_plot_logistic( .data=bgl_las$preds, truth=causal, estimate=pred, 
                         conf_level=0.95 ) +
  ggtitle("Original predictions")
p2 <- cal_plot_logistic( .data=bgl_las$preds, truth=causal, estimate=recal, 
                         conf_level=0.95 ) +
  ggtitle("Recalibrated predictions")
grid.arrange( p1, p2, ncol=2 )
dev.off()

# Scaled (for comparison)
cal_plot_logistic( .data=bgl_las$preds, truth=causal, estimate=scaled, 
                   conf_level=0.95 ) +
  ggtitle("Locally-scaled predictions")

# How do predictions differ before/after calibration?
# Plot predicted P(causal) with and without GLCs
plot( bgl_las$preds$pred, bgl_las$preds$recal )
abline( a=0, b=1 )


#-------------------------------------------------------------------------------
#   Figure 5?: CALDERA precision-recall curve versus NG+TP
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

# Plot
plot_loto_pr(bgl_las, type="recalibrated" )
points( x=cal_pr$recall, y=cal_pr$precision, pch=19, cex=1.5, col=cal_pr$col )
points( x=pnn_r,         y=pnn_p,            pch=19, cex=1.5, 
        col=brewer.pal( n=4, name="Blues" )[3] )


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
plot( p_tol$pred, p_tol$glc.pred )
abline( a=0, b=1 )

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
#   Full: GLM + XGB
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   GLM
#-------------------------------------------------------------------------------

# Run regression
f_form <- make_formula( lhs=f_cols[1], rhs=f_cols[-1] )
f_glm  <- glm( formula=f_form, data=tr, family="binomial" )

# Extract deviance explained
devexp_f_glm <- dev_exp(f_glm)


#-------------------------------------------------------------------------------
#   XGBoost
#-------------------------------------------------------------------------------

# Run
t0 <- proc.time()
f_xgb <- train_xgb( data=tr, feat_cols=f_cols[-1], maxit=100 )
plot_fi(f_xgb)
t1 <- proc.time()
t1-t0


#-------------------------------------------------------------------------------
#   Assess performance of the models in the test set
#-------------------------------------------------------------------------------

# AUPRC
f_pr <- auprc_test_set( test_df  = te, 
                        rand_mod = rand_mod,
                        glm_mod  = f_glm, 
                        las_mod  = f_las, 
                        xgb_mod  = f_xgb, 
                        ymax     = 1 )

# PR curves
plot_pr( model=f_glm, test_df=te )
plot_pr( model=f_xgb, test_df=te )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Basic + GLC: GLM
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   GLM
#-------------------------------------------------------------------------------

# Run regression
bg_form <- make_formula( lhs=bg_cols[1], rhs=bg_cols[-1] )
bg_glm0 <- glm( formula=bg_form, data=tr, family="binomial" )
bg_glm  <- stats::step( object=bg_glm0, k=qchisq( 1-0.05/5, df=1 ) )

# Extract deviance explained
dev_exp( model=bg_glm, in_sample=TRUE  ) / dev_exp( model=fg_glm, in_sample=TRUE  )
dev_exp( model=bg_glm, in_sample=FALSE ) / dev_exp( model=fg_glm, in_sample=FALSE )

# Forest plot
glm_to_forest_p( mod=bg_glm, xmax=NULL )
glm_to_forest(   mod=bg_glm, suffix="", standardize=FALSE )
glm_to_forest(   mod=bg_glm, suffix="", standardize=TRUE )


#-------------------------------------------------------------------------------
#   XGBoost
#-------------------------------------------------------------------------------

bg_xgb <- train_xgb( data=tr, feat_cols=bg_cols[-1], maxit=500, 
                     bias_cols=bias_cols )
xgb_to_forest( xgb_mod=bg_xgb, suffix="", xmax=NULL )
plot_fi(bg_xgb)


#-------------------------------------------------------------------------------
#   Assess performance of the models in the test set
#-------------------------------------------------------------------------------

# AUPRC: REMOVE gene-level covariates
b_pr1 <- auprc_test_set( test_df   = te, 
                         rand_mod  = rand_mod,
                         glm_mod   = bg_glm, 
                         las_mod   = bg_las, 
                         xgb_mod   = bg_xgb1, 
                         ymax      = 1, 
                         bias_cols = bias_cols, 
                         glc_cols  = glc_cols,
                         rm_glc    = TRUE )
fg_pr1; b_pr1; b_pr1/fg_pr1
og_pr1; b_pr1; b_pr1/og_pr1

# AUPRC: KEEP gene-level covariates
b_pr2 <- auprc_test_set( test_df   = te, 
                         rand_mod  = rand_mod,
                         glm_mod   = bg_glm, 
                         las_mod   = bg_las, 
                         xgb_mod   = bg_xgb2, 
                         ymax      = 1, 
                         bias_cols = bias_cols )
fg_pr2; b_pr2; b_pr2/fg_pr2
og_pr2; b_pr2; b_pr2/og_pr2

# PR curves: REMOVE gene-level covariates
plot_pr( model=bg_glm,  test_df=te, bias_cols=bias_cols, glc_cols=glc_cols, rm_glc=TRUE )
plot_pr( model=bg_xgb1, test_df=te, bias_cols=bias_cols, glc_cols=glc_cols, rm_glc=TRUE )

# PR curves: KEEP gene-level covariates
plot_pr( model=bg_glm,  test_df=te, bias_cols=bias_cols )
plot_pr( model=bg_xgb2, test_df=te, bias_cols=bias_cols )


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
#   Published best-in-locus
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

# Run a naive regression using published best-in-locus
p_cols <- c( "causal", "pops_bil", "dist_gene_bil", "magma_bil", "twas_bil", 
             "clpp_bil", "abc_bil", "corr_any_bil", "pchic_any_bil", "smr_bil" )

# Raw feature correlations
corrplot( cor( tr[ , ..p_cols ][,-1] ), order="hclust" )

# Run regression
p_form <- make_formula( lhs=p_cols[1], rhs=p_cols[-1] )
p_glm  <- glm( formula=p_form, data=tr, family="binomial" )

# Extract deviance explained
devexp_p_glm <- dev_exp(p_glm)

# Forest plot
glm_to_forest(p_glm)


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   All best-in-locus (only)
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

# Raw feature correlations
corrplot( cor( tr[ , ..a_cols ][,-1] ), order="hclust" )

# Run regression
a_form <- make_formula( lhs=a_cols[1], rhs=a_cols[-1] )
a_glm  <- glm( formula=a_form, data=tr, family="binomial" )

# Extract deviance explained
devexp_a_glm <- dev_exp(a_glm)

# Forest plot
glm_to_forest_p( mod=a_glm, suffix="_bilTRUE" )


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
plot_logOR_relationship_smooth( data=tr, varname="dist_gene_raw_l2g" )
plot_logOR_relationship_smooth( data=tr, varname="dist_gene_glo" )

# Plot relationship: distance to TSS
plot_logOR_relationship_smooth( data=tr, varname="distance_tss" )
plot_logOR_relationship_smooth( data=tr, varname="dist_tss_raw_l2g" )
plot_logOR_relationship_smooth( data=tr, varname="dist_tss_glo" )

# Plot relationship: POPS
plot_logOR_relationship_smooth( data=tr, varname="pops_glo" )

# Plot relationship: TWAS
plot_logOR_relationship_smooth( data=tr, varname="twas_logp_glo" )
plot_logOR_relationship_smooth( data=tr, varname="twas_glo" )

# Plot relationship: E-P Liu
plot_logOR_relationship_smooth( data=tr[ tr$corr_liu_raw_l2g > 0 , ], 
                                varname="corr_liu_raw_l2g" )
plot_logOR_relationship_smooth( data=tr, varname="corr_liu_raw_l2g" )
plot_logOR_relationship_smooth( data=tr[ tr$corr_liu_glo > 0 , ], 
                                varname="corr_liu_glo" )
plot_logOR_relationship_smooth( data=tr, varname="corr_liu_glo" )

# Plot relationship: E-P Andersson
plot_logOR_relationship_smooth( data=tr[ tr$corr_and_glo > 0 , ], 
                                varname="corr_and_glo" )
plot_logOR_relationship_smooth( data=tr, varname="corr_and_glo" )

# Plot relationship: E-P Ulirsch
plot_logOR_relationship_smooth( data=tr[ tr$corr_uli_glo > 0 , ], 
                                varname="corr_uli_glo" )
plot_logOR_relationship_smooth( data=tr, varname="corr_uli_glo" )

# Plot relationship: PCHi-C Jung
plot_logOR_relationship_smooth( data=tr[ tr$pchic_jung_glo > 0 , ], 
                                varname="pchic_jung_glo" )
plot_logOR_relationship_smooth( data=tr, varname="pchic_jung_glo" )

# Plot relationship: PCHi-C Javierre
plot_logOR_relationship_smooth( data=tr[ tr$pchic_jav_glo > 0 , ], 
                                varname="pchic_jav_glo" )
plot_logOR_relationship_smooth( data=tr, varname="pchic_jav_glo" )

# Plot relationship: CLPP
plot_logOR_relationship_smooth( data=tr[ tr$clpp_raw_l2g > 0 , ], 
                                varname="clpp_raw_l2g" )
plot_logOR_relationship_smooth( data=tr, varname="clpp_raw_l2g" )
plot_logOR_relationship_smooth( data=tr[ tr$clpp_glo > log10(0.0001) , ], 
                                varname="clpp_glo" )
plot_logOR_relationship_smooth( data=tr, varname="clpp_glo" )

# Plot relationship: ABC
plot_logOR_relationship_smooth( data=tr[ tr$abc_raw_l2g > 0 , ], 
                                varname="abc_raw_l2g" )
plot_logOR_relationship_smooth( data=tr, varname="abc_raw_l2g" )
plot_logOR_relationship_smooth( data=tr[ tr$abc_raw_l2g > 0 , ], 
                                varname="abc_glo" )
plot_logOR_relationship_smooth( data=tr, varname="abc_glo" )

# Plot relationship: MAGMA
plot_logOR_relationship_smooth( data=tr, varname="magma_glo" )
plot_logOR_relationship_smooth( data=tr[ tr$magma_glo != median(tr$magma_glo) , ], 
                                varname="magma_glo" )

# Plot relationship: SMR
plot_logOR_relationship_smooth( data=tr, varname="smr_glo" )
plot_logOR_relationship_smooth( data=tr[ tr$smr_glo != median(tr$smr_glo) , ], 
                                varname="smr_glo" )

# Plot relationship: coding
plot_logOR_relationship_smooth( data=tr[ !is.na(tr$coding_prob) , ], varname="coding_prob" )
plot_logOR_relationship_smooth( data=tr, varname="coding_glo" )
plot_logOR_relationship_smooth( data=tr[ tr$coding_glo != min(tr$coding_glo) , ], 
                                varname="coding_glo" )

# DEPICT and NetWAS
plot_logOR_relationship_smooth( data=tr, varname="depict_z_glo" )
plot_logOR_relationship_smooth( data=tr[ tr$depict_z_glo > 0 , ], varname="depict_z_glo" )
plot_logOR_relationship_smooth( data=tr, varname="netwas_score_glo" )
plot_logOR_relationship_smooth( data=tr, varname="netwas_bon_score_glo" )


#-------------------------------------------------------------------------------
#   Compare competing feature definitions
#-------------------------------------------------------------------------------

# Distance to gene body
mod1 <- glm( causal ~ distance_genebody, data=tr, family="binomial" )
mod2 <- glm( causal ~ dist_gene_raw_l2g, data=tr, family="binomial" )
mod3 <- glm( causal ~ dist_gene_glo, data=tr, family="binomial" )
dev_exp(mod1)
dev_exp(mod2)
dev_exp(mod3)
summary(mod1)$coef
summary(mod2)$coef
summary(mod3)$coef

# Distance to TSS
mod1 <- glm( causal ~ distance_tss,     data=tr, family="binomial" )
mod2 <- glm( causal ~ dist_tss_raw_l2g, data=tr, family="binomial" )
mod3 <- glm( causal ~ dist_tss_glo, data=tr, family="binomial" )
dev_exp(mod1)
dev_exp(mod2)
dev_exp(mod3)
summary(mod1)$coef
summary(mod2)$coef
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
#   Plot features v. residuals
#-------------------------------------------------------------------------------

# Plot features v. residuals
plot_residuals <- function( data, model, variable ){
  plot( data[[variable]], residuals(model),
        xlim=quantile( data[[variable]], probs=c( 0.05, 0.95 ) ),
        ylim=quantile( residuals(model), probs=c( 0.05, 0.95 ) ) )
  lines( lowess( data[[variable]], residuals(model) ), lwd=3, col="orange2" )
}
plot_residuals( data=tr, model=s_glm, variable="pops_glo" )
plot_residuals( data=tr, model=s_glm, variable="pops_rel" )
plot_residuals( data=tr, model=s_glm, variable="dist_gene_rel" )
plot_residuals( data=tr, model=s_glm, variable="magma_rel" )
plot_residuals( data=tr, model=s_glm, variable="coding_glo" )
plot_residuals( data=tr, model=s_glm, variable="coding_rel" )


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
d_form <-  as.formula("factor(causal) ~ distance_genebody + pops_glo")
tree_fit <- tree_spec %>%
  fit( d_form, data = tr[ tr$pops_and_nearest , ] )

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
for( i in seq_along(model) ){
  hyp0[[ names(model)[i] ]] <- unlist( model[[i]]$learner$par.vals[par_names] )
}
hyp <- as.data.frame( do.call( rbind, hyp0 ) )

# Plot
corrplot( cor(hyp), order="hclust" )
plot( hyp$max_depth, hyp$eta, pch=19, col="steelblue", cex=2, las=1,
      xlab="Number of rounds", ylab="Learning rate (eta)" )
hist( hyp[[1]], col="steelblue", breaks=10, main=par_names[1], las=1, xlim=c( 75, 175 ) )
hist( hyp[[2]], col="steelblue", breaks=10, main=par_names[2], las=1, xlim=c( 0.01, 0.5 ) )
hist( hyp[[3]], col="steelblue", breaks=10, main=par_names[3], las=1, xlim=c( 1, 2 ) )
hist( hyp[[4]], col="steelblue", breaks=10, main=par_names[4], las=1, xlim=c( 0, 8 ) )
hist( hyp[[5]], col="steelblue", breaks=10, main=par_names[5], las=1, xlim=c( 0, 6 ) )
hist( hyp[[6]], col="steelblue", breaks=10, main=par_names[6], las=1, xlim=c( 0.3, 1 ) )
hist( hyp[[7]], col="steelblue", breaks=10, main=par_names[7], las=1, xlim=c( 0.1, 1 ) )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Done
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------


















