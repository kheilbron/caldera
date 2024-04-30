
#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Functions
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   simple_functions
#-------------------------------------------------------------------------------

# ci95_lo and ci95_hi: Return the lower/upper 95% confidence interval
ci95_lo <- function( beta, se, zcrit=qnorm(.025, lower.tail=FALSE) )  beta - zcrit * se
ci95_hi <- function( beta, se, zcrit=qnorm(.025, lower.tail=FALSE) )  beta + zcrit * se

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
  de <- c( ( mod1$null.deviance - mod1$deviance ) / mod1$null.deviance,
           ( mod2$null.deviance - mod2$deviance ) / mod2$null.deviance,
           ( mod3$null.deviance - mod3$deviance ) / mod3$null.deviance )
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

glm_to_forest <- function( mod, mtext=TRUE, xmax=NULL ){
  
  # Clean up row and column names
  mod2 <- as.data.frame( summary(mod)$coef[ -1, , drop=FALSE ] )
  names(mod2) <- c( "effect", "se", "z", "p" )
  mod3 <- mod2[ grepl( pattern="TRUE$", row.names(mod2) ) , ]
  row.names(mod3) <- sub( "TRUE$", "", row.names(mod3) )
  row.names(mod3) <- sub( "_bil$", "", row.names(mod3) )
  
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

lasso_to_forest <- function( las, xmax=NULL ){
  idx <- which( las$lambda == las$lambda.min )
  las_vec  <- exp( las$glmnet.fit$beta[ , idx ] )
  las_vec2 <- las_vec[ grepl( pattern="_bil$", names(las_vec) ) ]
  las_df  <- data.frame( row.names = sub( "_bil$", "", names(las_vec2) ),
                         or        = las_vec2,
                         col       = ifelse( las_vec2 == 1, "grey", "#70AD47" ) )
  forest_plot( df=las_df, colour.col="col", xlab="Odds ratio", 
               margins=c(5,9,1,1), vert.line.pos=1, xmax=xmax,
               value.col="or", lo.col="or", hi.col="or", mtext.col=NULL )
  par( mar=c(5,5,4,1) )
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
  # Plot
  forest_plot( df=fi3, colour.col="col", xlab="Feature importance", 
               margins=c(5,9,1,1), vert.line.pos=0, mtext.col=NULL, 
               value.col="fi", lo.col="fi", hi.col="fi", xmax=xmax )
  par( mar=c(5,5,4,1) )
}


#-------------------------------------------------------------------------------
#   auprc_test_set:  Barplot of AUPRCs from a GLM, LASSO, and XGBoost model
#-------------------------------------------------------------------------------

auprc_test_set <- function( test_df, rand_mod=NULL, glm_mod=NULL, glmm_mod=NULL,
                            las_mod=NULL, xgb_mod=NULL, ymax=NULL, avg_prot_att=TRUE ){
  
  # Fix RVIS covariates to their mean value in the test dataset
  avg_rvis <- TRUE
  if(avg_rvis){
    test_df$rvis4       <- mean(test_df$rvis4)
    test_df$rvis4_poly2 <- mean(test_df$rvis4_poly2)
  }
  
  # Fix protein attenuation covariates to their mean value in the test dataset
  if(avg_prot_att){
    test_df$prot_att_poly2 <- mean(test_df$prot_att_poly2)
    test_df$prot_att_miss  <- as.numeric(test_df$prot_att_miss)
    test_df$prot_att_miss  <- mean(test_df$prot_att_miss)
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

plot_pr <- function( model, test_df, p, r, v, avg_prot_att=TRUE ){
  
  # Fix RVIS covariates to their mean value in the test dataset
  avg_rvis <- TRUE
  if(avg_rvis){
    test_df$rvis4       <- mean(test_df$rvis4)
    test_df$rvis4_poly2 <- mean(test_df$rvis4_poly2)
  }
  
  # Fix protein attenuation covariates to their mean value in the test dataset
  if(avg_prot_att){
    test_df$prot_att_poly2 <- mean(test_df$prot_att_poly2)
    test_df$prot_att_miss  <- as.numeric(test_df$prot_att_miss)
    test_df$prot_att_miss  <- mean(test_df$prot_att_miss)
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

recalibrate_preds <- function( test_df, model, avg_prot_att=TRUE ){
  
  # Fix RVIS covariates to their mean value in the test dataset
  avg_rvis <- TRUE
  if(avg_rvis){
    test_df$rvis4       <- mean(test_df$rvis4)
    test_df$rvis4_poly2 <- mean(test_df$rvis4_poly2)
  }
  
  # Fix protein attenuation covariates to their mean value in the test dataset
  if(avg_prot_att){
    test_df$prot_att_poly2 <- mean(test_df$prot_att_poly2)
    test_df$prot_att_miss  <- as.numeric(test_df$prot_att_miss)
    test_df$prot_att_miss  <- mean(test_df$prot_att_miss)
  }
  
  # Extract model predictions
  # GLM
  if( any( class(model) == "glm" ) ){
    p_mod <- predict( object=model, newdata=test_df, type="response" )
    
  # LASSO
  }else if( class(model) == "cv.glmnet" ){
    feat_cols <- model$glmnet$beta@Dimnames[[1]]
    p_mod <- predict( model, newx=as.matrix( test_df[ , ..feat_cols ] ), 
                      s="lambda.min", type="response" )
    
  # XGBoost
  }else if( class(model) == "WrappedModel" ){
    feat_cols <- sub( pattern="TRUE$", replacement="", x=model$features[-1] )
    te_data <- data.frame( causal=as.factor( as.integer( test_df[["causal"]] ) ), 
                           model.matrix( ~.+1, data=test_df[ , ..feat_cols ] )  )
    te_task <- makeClassifTask( data=te_data, target="causal" )
    p_mod <- predict( model, te_task )$data$prob.1
  }
  test_df$pred <- p_mod
  
  # Initialize output
  out <- data.table( causal=test_df[["causal"]], original=p_mod )
  out$flames <- out$scaled <- out$scaled_up <- out$scaled_down <- out$best <- 
    out$relative <- out$global <- as.numeric(NA)
  
  # Loop through test set loci
  test_df$tcp <- paste( test_df$trait, test_df$region, test_df$cs_id, sep="_" )
  for( i in unique(test_df$tcp) ){
    
    # Subset to locus
    sub <- test_df[ test_df$tcp == i , ]
    idx <- order( -sub[["pred"]] )
    
    # Re-calibrate: scaled
    scaled <- sub[["pred"]] / sum( sub[["pred"]] )
    
    # Re-calibrate: scaled_up
    if( sum( sub[["pred"]] ) < 1 ){
      scaled_up <- sub[["pred"]] / sum( sub[["pred"]] )
    }else{
      scaled_up <- sub[["pred"]]
    }
    
    # Re-calibrate: scaled_down
    if( sum( sub[["pred"]] ) > 1 ){
      scaled_down <- sub[["pred"]] / sum( sub[["pred"]] )
    }else{
      scaled_down <- sub[["pred"]]
    }
    
    # Make variables for global, relative, and BIL scores
    glo <- log( prop_to_odds( sub[["pred"]] ) )
    rel <- glo - max(glo)
    bil <- ifelse( glo == max(glo), 1, 0 )
    
    # Re-calibrate: FLAMES
    if( NROW(sub) >=2 ){
      flames <- 2*sub[["pred"]] - sub[["pred"]][ idx[1] ] - sub[["pred"]][ idx[2] ]
    }else{
      flames <- 2*sub[["pred"]] - sub[["pred"]][ idx[1] ] - sub[["pred"]][ idx[1] ]
    }
    
    # Assign values: scaled
    set( x     = out, 
         i     = which( test_df$tcp == i ), 
         j     = "scaled", 
         value = scaled )
    
    # Assign values: scaled_up
    set( x     = out, 
         i     = which( test_df$tcp == i ), 
         j     = "scaled_up", 
         value = scaled_up )
    
    # Assign values: scaled
    set( x     = out, 
         i     = which( test_df$tcp == i ), 
         j     = "scaled_down", 
         value = scaled_down )
    
    # Assign values: FLAMES
    set( x     = out, 
         i     = which( test_df$tcp == i ), 
         j     = "flames", 
         value = flames )
    
    # Assign values: global score
    set( x     = out, 
         i     = which( test_df$tcp == i ), 
         j     = "global", 
         value = glo )
    
    # Assign values: relative score
    set( x     = out, 
         i     = which( test_df$tcp == i ), 
         j     = "relative", 
         value = rel )
    
    # Assign values: BIL
    set( x     = out, 
         i     = which( test_df$tcp == i ), 
         j     = "best", 
         value = bil )
  }
  
  # Return
  out
}


#-------------------------------------------------------------------------------
#   recalibration_model: Train a LASSO to re-calibrate model outputs
#-------------------------------------------------------------------------------

recalibration_model <- function(data){
  
  # Build the model to predict causal genes based on step 1 predictions
  pred_cols <- c( "global", "relative", "best" )
  cv.glmnet( x=as.matrix( data[ , ..pred_cols ] ), 
             y=data[["causal"]], family="binomial", standardize=FALSE )
}


#-------------------------------------------------------------------------------
#   train_xgb:       Train an XGBoost model
#-------------------------------------------------------------------------------

train_xgb <- function( data, feat_cols, rm_prot_att=TRUE, maxit=20 ){
  
  # If specified, remove RVIS features
  rm_rvis <- TRUE
  if(rm_rvis){
    feat_cols <- grep( pattern="rvis", x=feat_cols, invert=TRUE, value=TRUE )
  }
  
  # If specified, remove protein attenuation features
  if(rm_prot_att){
    feat_cols <- grep( pattern="prot_att", x=feat_cols, invert=TRUE, value=TRUE )
  }
  
  # XGBoost (basic hyperparameter tuning): raw
  tr_data <- data.frame( causal=as.factor( as.integer( data[["causal"]] ) ), 
                         model.matrix( ~.+1, data=data[ , ..feat_cols ] ) )
  tr_task <- makeClassifTask( data=tr_data, target="causal" )
  lrn <- makeLearner( cl="classif.xgboost", predict.type="prob",
                      objective="binary:logistic", verbose=1 )
  rdesc <- makeResampleDesc( "CV", stratify=TRUE, iters=5L )
  parallelStartSocket( cpus=detectCores() )
  pars <- makeParamSet( makeDiscreteParam( "booster",          values = "gbtree" ), 
                        makeIntegerParam(  "nrounds",          lower = 50L,  upper = 200L ), 
                        makeNumericParam(  "eta",              lower = 0.01, upper = 0.4 ), 
                        makeIntegerParam(  "max_depth",        lower = 3L,   upper = 10L ), 
                        makeNumericParam(  "min_child_weight", lower = 1L,   upper = 10L ), 
                        makeNumericParam(  "subsample",        lower = 0.5,  upper = 1 ), 
                        makeNumericParam(  "colsample_bytree", lower = 0.5,  upper = 1 ) )
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
library(data.table)
library(corrplot)
library(glmnet)
library(tidyverse)
library(ggrepel)
library(lme4)
library(RColorBrewer)
library(PRROC)
library(mlr)
library(xgboost)
library(parallel)
library(parallelMap) 
library(patchwork)
library(viridis)

# Read in causal/non-causal TGPs annotated with gene mapping evidence
cnc_map_file <- file.path( maindir, "causal_noncausal_trait_gene_pairs",
                           "causal_tgp_and_gene_mapping_data_300kb.tsv" )
dt <- fread(file=cnc_map_file)

# Split into training (80%) and testing (20%) sets based on trait
set.seed(1)
tr_traits <- sample( x=unique(dt$trait), size=round( length( unique(dt$trait) )*0.8 ) )
te_traits <- setdiff( unique(dt$trait), tr_traits )
tr <- dt[ dt$trait %in% tr_traits , ]
te <- dt[ dt$trait %in% te_traits , ]
NROW(tr) / NROW(dt)


#-------------------------------------------------------------------------------
#   Random predictor
#-------------------------------------------------------------------------------

# Run a regression using a random predictor
set.seed(1)
tr$rand <- rnorm( n=NROW(tr) )
te$rand <- rnorm( n=NROW(te) )
rand_mod <- glm( causal ~ rand, data=tr, family="binomial" )


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

# Extract deviance explained
devexp_l_glm <- ( l_glm$null.deviance - l_glm$deviance ) / l_glm$null.deviance

# Forest plot
glm_to_forest(l_glm)

# Precision and recall: training
sum( tr$causal & tr$pops_and_local ) / sum( tr$pops_and_local ) #precision
sum( tr$causal & tr$pops_and_local ) / sum( tr$causal )         #recall

# Precision and recall: test
sum( te$causal & te$pops_and_local ) / sum( te$pops_and_local ) #precision
sum( te$causal & te$pops_and_local ) / sum( te$causal )         #recall


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Published best-in-locus
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   GLM
#-------------------------------------------------------------------------------

# Run a naive regression using published best-in-locus
pcols <- c( "causal", "pops_bil", "dist_gene_bil", "magma_bil", "twas_bil", 
            "clpp_bil", "abc_bil", "corr_any_bil", "pchic_any_bil", "smr_bil" )

# Raw feature correlations
corrplot( cor( tr[ , ..pcols ][,-1] ), order="hclust" )

# Run regression
p_form <- make_formula( lhs=pcols[1], rhs=pcols[-1] )
p_glm  <- glm( formula=p_form, data=tr, family="binomial" )

# Extract deviance explained
devexp_p_glm <- ( p_glm$null.deviance - p_glm$deviance ) / p_glm$null.deviance

# Forest plot
glm_to_forest(p_glm)


#-------------------------------------------------------------------------------
#   LASSO
#-------------------------------------------------------------------------------

# Make LASSO model
p_las <- cv.glmnet( x=as.matrix( tr[ , ..pcols ] )[,-1], 
                    y=tr[["causal"]], family="binomial", standardize=FALSE )

# Forest plot
lasso_to_forest(p_las)

# Produce plot of deviance by lambda value
plot(p_las)


#-------------------------------------------------------------------------------
#   XGBoost
#-------------------------------------------------------------------------------

# Train XGBoost model
p_xgb <- train_xgb( data=tr, feat_cols=pcols[-1], maxit=100 )

# Forest plot
xgb_to_forest(p_xgb)


#-------------------------------------------------------------------------------
#   Assess performance of the models in the test set
#-------------------------------------------------------------------------------

auprc_test_set( test_df  = te, 
                rand_mod = rand_mod,
                glm_mod  = p_glm, 
                las_mod  = p_las, 
                xgb_mod  = p_xgb )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   All best-in-locus (only)
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

# Run a naive regression using published best-in-locus
acols <- c( "causal", 
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

# Raw feature correlations
corrplot( cor( tr[ , ..acols ][,-1] ), order="hclust" )

# Run regression
a_form <- make_formula( lhs=acols[1], rhs=acols[-1] )
a_glm  <- glm( formula=a_form, data=tr, family="binomial" )

# Extract deviance explained
devexp_a_glm <- ( a_glm$null.deviance - a_glm$deviance ) / a_glm$null.deviance

# Forest plot
glm_to_forest(a_glm)
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
plot_logOR_relationship_smooth( data=tr, varname="coding_glo" )
plot_logOR_relationship_smooth( data=tr[ tr$coding_glo != -3 , ], 
                                varname="coding_glo" )

# DEPICT and NetWAS
plot_logOR_relationship_smooth( data=tr, varname="depict_z_glo" )
plot_logOR_relationship_smooth( data=tr, varname="netwas_score_glo" )
plot_logOR_relationship_smooth( data=tr, varname="netwas_bon_score_glo" )


#-------------------------------------------------------------------------------
#   Compare competing feature definitions
#-------------------------------------------------------------------------------

# Distance to gene body
mod1 <- glm( causal ~ distance_genebody, data=tr, family="binomial" )
mod2 <- glm( causal ~ dist_gene_raw_l2g, data=tr, family="binomial" )
mod3 <- glm( causal ~ dist_gene_glo, data=tr, family="binomial" )
( mod1$null.deviance - mod1$deviance ) / mod1$null.deviance
( mod2$null.deviance - mod2$deviance ) / mod2$null.deviance
( mod3$null.deviance - mod3$deviance ) / mod3$null.deviance
summary(mod1)$coef
summary(mod2)$coef
summary(mod3)$coef

# Distance to TSS
mod1 <- glm( causal ~ distance_tss,     data=tr, family="binomial" )
mod2 <- glm( causal ~ dist_tss_raw_l2g, data=tr, family="binomial" )
mod3 <- glm( causal ~ dist_tss_glo, data=tr, family="binomial" )
( mod1$null.deviance - mod1$deviance ) / mod1$null.deviance
( mod2$null.deviance - mod2$deviance ) / mod2$null.deviance
( mod3$null.deviance - mod3$deviance ) / mod3$null.deviance
summary(mod1)$coef
summary(mod2)$coef
summary(mod3)$coef

# Number of genes in the locus
mod1 <- glm( causal ~ ngenes_nearby,       data=tr, family="binomial" )
mod2 <- glm( causal ~ prior_n_genes_locus, data=tr, family="binomial" )
( mod1$null.deviance - mod1$deviance ) / mod1$null.deviance
( mod2$null.deviance - mod2$deviance ) / mod2$null.deviance
summary(mod1)$coef
summary(mod2)$coef

# TWAS
mod1 <- glm( causal ~ twas_logp_glo,  data=tr, family="binomial" )
mod2 <- glm( causal ~ twas_glo,  data=tr, family="binomial" )
( mod1$null.deviance - mod1$deviance ) / mod1$null.deviance
( mod2$null.deviance - mod2$deviance ) / mod2$null.deviance
summary(mod1)$coef
summary(mod2)$coef

# Plot relationship: E-P Liu
mod1 <- glm( causal ~ corr_liu_raw_l2g,  data=tr, family="binomial" )
mod2 <- glm( causal ~ corr_liu_glo,  data=tr, family="binomial" )
( mod1$null.deviance - mod1$deviance ) / mod1$null.deviance
( mod2$null.deviance - mod2$deviance ) / mod2$null.deviance
summary(mod1)$coef
summary(mod2)$coef

# CLPP
mod1 <- glm( causal ~ clpp_raw_l2g,  data=tr, family="binomial" )
mod2 <- glm( causal ~ clpp_glo,  data=tr, family="binomial" )
( mod1$null.deviance - mod1$deviance ) / mod1$null.deviance
( mod2$null.deviance - mod2$deviance ) / mod2$null.deviance
summary(mod1)$coef
summary(mod2)$coef

# ABC
mod1 <- glm( causal ~ abc_raw_l2g,  data=tr, family="binomial" )
mod2 <- glm( causal ~ abc_glo,  data=tr, family="binomial" )
( mod1$null.deviance - mod1$deviance ) / mod1$null.deviance
( mod2$null.deviance - mod2$deviance ) / mod2$null.deviance
summary(mod1)$coef
summary(mod2)$coef

# DEPICT and NetWAS
mod1 <- glm( causal ~ depict_z_glo,         data=tr, family="binomial" )
mod2 <- glm( causal ~ netwas_score_glo,     data=tr, family="binomial" )
mod3 <- glm( causal ~ netwas_bon_score_glo, data=tr, family="binomial" )
( mod1$null.deviance - mod1$deviance ) / mod1$null.deviance
( mod2$null.deviance - mod2$deviance ) / mod2$null.deviance
( mod3$null.deviance - mod3$deviance ) / mod3$null.deviance


#-------------------------------------------------------------------------------
#   GLM
#-------------------------------------------------------------------------------

# Run a naive regression using published best-in-locus
glo_cols <- c( "causal", "pops_glo", 
               "dist_gene_glo", "dist_tss_glo",
               "magma_glo",
               "coding_glo",
               "twas_glo",
               "corr_liu_glo", "corr_and_glo", "corr_uli_glo", 
               "pchic_jung_glo", "pchic_jav_glo", 
               "clpp_glo", "smr_glo", "abc_glo",
               "depict_z_glo", "netwas_score_glo", "netwas_bon_score_glo" )

# Raw feature correlations
corrplot( cor( tr[ , ..glo_cols ][,-1] ), order="hclust" )

# Run regression
glo_form <- make_formula( lhs=glo_cols[1], rhs=glo_cols[-1] )
glo_glm  <- glm( formula=glo_form, data=tr, family="binomial" )

# Extract deviance explained
devexp_glo_glm <- ( glo_glm$null.deviance - glo_glm$deviance ) / glo_glm$null.deviance

# Forest plot
glm_to_forest_p( mod=glo_glm, suffix="" )


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
plot_logOR_relationship_smooth( data=tr[ tr$dist_gene_rel != 3 , ], varname="dist_gene_rel" )

# Plot relationship: distance to TSS
plot_logOR_relationship_smooth( data=tr, varname="dist_tss_rel" )
plot_logOR_relationship_smooth( data=tr[ tr$dist_tss_rel != 3 , ], varname="dist_tss_rel" )

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
#   GLM
#-------------------------------------------------------------------------------

# Run a naive regression using published best-in-locus
rel_cols <- c( "causal", "pops_rel", "dist_gene_rel", "dist_tss_rel", 
               "magma_rel", "coding_rel", "twas_rel", 
               "corr_liu_rel", "corr_and_rel", "corr_uli_rel", 
               "pchic_jung_rel", "pchic_jav_rel",  
               "clpp_rel", "smr_rel", "abc_rel",
               "depict_z_rel", "netwas_score_rel", "netwas_bon_score_rel" )

# Raw feature correlations
corrplot( cor( tr[ , ..rel_cols ][,-1] ), order="hclust" )

# Run regression
rel_form <- make_formula( lhs=rel_cols[1], rhs=rel_cols[-1] )
rel_glm  <- glm( formula=rel_form, data=tr, family="binomial" )

# Extract deviance explained
devexp_rel_glm <- ( rel_glm$null.deviance - rel_glm$deviance ) / rel_glm$null.deviance

# Forest plot
glm_to_forest_p( mod=rel_glm, suffix="_rel" )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Covariates-only
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Feature engineering 
#-------------------------------------------------------------------------------

# Protein attenuation
plot_logOR_relationship_smooth( data=tr, varname="prot_att" )
plot_logOR_relationship_smooth( data=tr[ tr$prot_att_miss == 0 , ], 
                                varname="prot_att" )
plot_logOR_relationship_smooth( data=tr, varname="prot_att_poly2" )
plot_logOR_relationship_smooth( data=tr[ tr$prot_att_miss == 0 , ], 
                                varname="prot_att_poly2" )

# Protein attenuation: What is the optimal value(s) to impute NAs to?
mod <- glm( causal ~ prot_att + prot_att_poly2, 
            data=tr[ tr$prot_att_miss == 0 , ], family="binomial" )
quad_miss_imp( data=tr, model=mod, miss_var="prot_att_miss" )

# RVIS, truncated
plot_logOR_relationship_smooth( data=tr[ tr$rvis_miss == 0 , ], varname="rvis" )
plot_logOR_relationship_smooth( data=tr[ tr$rvis_miss == 0 , ], varname="rvis4" )
plot_logOR_relationship_smooth( data=tr[ tr$rvis_miss == 0 , ], varname="rvis4_poly2" )
plot_logOR_relationship_smooth( data=tr, varname="rvis4" )

# RVIS: What is the optimal value to impute NAs to?
mod <- glm( causal ~ rvis4 + rvis4_poly2, 
            data=tr[ tr$rvis_miss == 0 , ], family="binomial" )
quad_miss_imp( data=tr, model=mod, miss_var="rvis_miss" )

# Number of genes in the locus
plot_logOR_relationship_smooth( data=tr, varname="prior_n_genes_locus" )
plot_logOR_relationship_smooth( data=tr, varname="ngenes_nearby" )

# Plot relationship: number of CSs in the region
plot_logOR_relationship_smooth( data=tr, varname="n_cs_in_region" )

# Plot relationship: number of SNPs in the CS
tr$log10_n_cs_snps <- log10(tr$n_cs_snps)
plot_logOR_relationship_smooth( data=tr, varname="n_cs_snps" )
plot_logOR_relationship_smooth( data=tr, varname="log10_n_cs_snps" )

# Plot relationship: width of the CS
tr$log10_cs_width <- log10(tr$cs_width)
plot_logOR_relationship_smooth( data=tr, varname="cs_width" )
plot_logOR_relationship_smooth( data=tr, varname="log10_cs_width" )


#-------------------------------------------------------------------------------
#   Compare competing feature definitions
#-------------------------------------------------------------------------------

# Number of genes in the locus
mod1 <- glm( causal ~ prior_n_genes_locus, data=tr, family="binomial" )
mod2 <- glm( causal ~ ngenes_nearby, data=tr, family="binomial" )
( mod1$null.deviance - mod1$deviance ) / mod1$null.deviance
( mod2$null.deviance - mod2$deviance ) / mod2$null.deviance
summary(mod1)$coef
summary(mod2)$coef

# Number of SNPs in the CS
mod1 <- glm( causal ~ n_cs_snps,       data=tr, family="binomial" )
mod2 <- glm( causal ~ log10_n_cs_snps, data=tr, family="binomial" )
( mod1$null.deviance - mod1$deviance ) / mod1$null.deviance
( mod2$null.deviance - mod2$deviance ) / mod2$null.deviance
summary(mod1)$coef
summary(mod2)$coef

# Width of the CS
mod1 <- glm( causal ~ cs_width,       data=tr, family="binomial" )
mod2 <- glm( causal ~ log10_cs_width, data=tr, family="binomial" )
( mod1$null.deviance - mod1$deviance ) / mod1$null.deviance
( mod2$null.deviance - mod2$deviance ) / mod2$null.deviance
summary(mod1)$coef
summary(mod2)$coef

# Protein attenuation
mod1 <- glm( causal ~ prot_att_class,           data=tr, family="binomial" )
mod2 <- glm( causal ~ prot_att,                 data=tr, family="binomial" )
mod3 <- glm( causal ~ prot_att + prot_att_miss, data=tr, family="binomial" )
( mod1$null.deviance - mod1$deviance ) / mod1$null.deviance
( mod2$null.deviance - mod2$deviance ) / mod2$null.deviance
( mod3$null.deviance - mod3$deviance ) / mod3$null.deviance
anova( mod2, mod3, test="Chi" )[["Pr(>Chi)"]][2]
summary(mod1)$coef
summary(mod2)$coef
summary(mod3)$coef

# RVIS
mod1 <- glm( causal ~ rvis  + rvis_miss, data=tr, family="binomial" )
mod2 <- glm( causal ~ rvis4 + rvis_miss, data=tr, family="binomial" )
( mod1$null.deviance - mod1$deviance ) / mod1$null.deviance
( mod2$null.deviance - mod2$deviance ) / mod2$null.deviance
summary(mod1)$coef
summary(mod2)$coef


#-------------------------------------------------------------------------------
#   GLM
#-------------------------------------------------------------------------------

# Run a naive regression using published best-in-locus
cov_cols <- c( "causal", "prior_n_genes_locus", 
               "prot_att", "prot_att_poly2", "prot_att_miss", 
               "rvis4", "rvis4_poly2", "rvis_miss" )

# Raw feature correlations
corrplot( cor( tr[ , ..cov_cols ][,-1] ), order="hclust" )

# Run regression
cov_form <- make_formula( lhs=cov_cols[1], rhs=cov_cols[-1] )
cov_glm  <- glm( formula=cov_form, data=tr, family="binomial" )

# Extract deviance explained
devexp_cov_glm <- ( cov_glm$null.deviance - cov_glm$deviance ) / cov_glm$null.deviance

# Forest plot
glm_to_forest_p(cov_glm, xmax=40)
summary(cov_glm)$coef


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   For each V2G method, compare BIL, global, and relative features
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

# Plot deviance explained by all univariate methods
all_de0 <- list()
all_cols <- unique( c( acols, glo_cols, rel_cols, cov_cols ) )[-1]
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
#   Full model: best-in-locus + global + relative + covariates
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   GLM
#-------------------------------------------------------------------------------

# Run a naive regression using published best-in-locus
fcols <- unique( c( acols, glo_cols, rel_cols, cov_cols ) )

# Raw feature correlations
corrplot( cor( tr[ , ..fcols ][,-1] ), order="hclust" )

# Run regression
f_form <- make_formula( lhs=fcols[1], rhs=fcols[-1] )
f_glm  <- glm( formula=f_form, data=tr, family="binomial" )

# Extract deviance explained
devexp_f_glm <- ( f_glm$null.deviance - f_glm$deviance ) / f_glm$null.deviance

# Forest plot
glm_to_forest_p( mod=f_glm, suffix="" )
glm_to_forest_p( mod=f_glm, suffix="TRUE" )
glm_to_forest_p( mod=f_glm, suffix="_glo" )
glm_to_forest_p( mod=f_glm, suffix="_rel" )


#-------------------------------------------------------------------------------
#   LASSO
#-------------------------------------------------------------------------------

# Make LASSO model
f_las <- cv.glmnet( x=as.matrix( tr[ , ..fcols ] )[,-1], 
                    y=tr[["causal"]], family="binomial", standardize=FALSE )


#-------------------------------------------------------------------------------
#   XGBoost
#-------------------------------------------------------------------------------

# WITHOUT protein attenuation
f_xgb1 <- train_xgb( data=tr, feat_cols=fcols[-1], maxit=100, rm_prot_att=TRUE )
xgb_to_forest( xgb_mod=f_xgb1 )
xgb_to_forest( xgb_mod=f_xgb1, suffix="_glo" )
xgb_to_forest( xgb_mod=f_xgb1, suffix="_rel" )

# WITH protein attenuation
f_xgb2 <- train_xgb( data=tr, feat_cols=fcols[-1], maxit=100, rm_prot_att=FALSE )
xgb_to_forest( xgb_mod=f_xgb2 )
xgb_to_forest( xgb_mod=f_xgb2, suffix="_glo" )
xgb_to_forest( xgb_mod=f_xgb2, suffix="_rel" )


#-------------------------------------------------------------------------------
#   Assess performance of the models in the test set
#-------------------------------------------------------------------------------

# AUPRC: WITHOUT protein attenuation
f_pr1 <- auprc_test_set( test_df  = te, 
                         rand_mod = rand_mod,
                         glm_mod  = f_glm, 
                         las_mod  = f_las, 
                         xgb_mod  = f_xgb1, 
                         ymax     = 1, 
                         avg_prot_att = TRUE )

# AUPRC: WITH protein attenuation
f_pr2 <- auprc_test_set( test_df  = te, 
                         rand_mod = rand_mod,
                         glm_mod  = f_glm, 
                         las_mod  = f_las, 
                         xgb_mod  = f_xgb2, 
                         ymax     = 1, 
                         avg_prot_att = FALSE )

# PR curves: WITHOUT protein attenuation
plot_pr( model=f_glm,  test_df=te, p=0.78, r=0.34, v=0.75, avg_prot_att=TRUE )
plot_pr( model=f_xgb1, test_df=te, p=0.78, r=0.34, v=0.75, avg_prot_att=TRUE )

# PR curves: WITH protein attenuation
plot_pr( model=f_glm,  test_df=te, p=0.78, r=0.34, v=0.75, avg_prot_att=FALSE )
plot_pr( model=f_xgb2, test_df=te, p=0.78, r=0.34, v=0.75, avg_prot_att=FALSE )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   "Reduced model": Maximize performance while minimizing the number of features
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Model selection
#-------------------------------------------------------------------------------

# Perform backwards stepwise model selection
# Remove features predicted to have P > 0.05/10
o_glm0 <- stats::step( object=f_glm, k=qchisq( 1-0.05/10, df=1 ) )

# Hard-coded list of features that passed model selection
# ocols0 <- c( "causal",
#              "pops_bil", "pops_glo", "pops_rel",
#              "dist_gene_rel",
#              "twas_bil", "twas_rel",
#              "magma_glo", "magma_rel",
#              "coding_glo", "coding_rel",
#              "corr_liu_rel",
#              "depict_z_glo",
#              "netwas_score_bil", "netwas_score_rel",
#              "netwas_bon_score_glo",
#              "prior_n_genes_locus",
#              "prot_att_poly2", "prot_att_miss",
#              "rvis4", "rvis4_poly2" )

# Manually make models and drop terms that are not Bonferroni-significant
ocols1 <- c( "causal",
             "pops_bil", "pops_glo", "pops_rel",
             "dist_gene_rel",
             "twas_bil", "twas_rel",
             "magma_glo", "magma_rel",
             "coding_glo", "coding_rel",
             "corr_liu_rel",
             "depict_z_glo",
             "netwas_score_bil",
             "prior_n_genes_locus",
             "prot_att_poly2", "prot_att_miss",
             "rvis4", "rvis4_poly2" )

# Run regression
o_form1 <- make_formula( lhs=ocols1[1], rhs=ocols1[-1] )
o_glm1  <- glm( formula=o_form1, data=tr, family="binomial" )

# Forest plot
glm_to_forest_p( mod=o_glm1, xmax=NULL )


#-------------------------------------------------------------------------------
#   GLM
#-------------------------------------------------------------------------------

# "Reduced" list of selected faeatures
ocols <- c( "causal",
            "pops_bil", "pops_glo", "pops_rel",
            "dist_gene_rel",
            "twas_bil", "twas_rel",
            "magma_glo", "magma_rel",
            "coding_glo", "coding_rel",
            "corr_liu_rel",
            "depict_z_glo",
            "netwas_score_bil",
            "prior_n_genes_locus",
            "prot_att_poly2", "prot_att_miss",
            "rvis4", "rvis4_poly2" )

# Raw feature correlations
corrplot( cor( tr[ , ..ocols ][,-1] ), order="hclust" )

# Run regression
o_form <- make_formula( lhs=ocols[1], rhs=ocols[-1] )
o_glm  <- glm( formula=o_form, data=tr, family="binomial" )

# Extract deviance explained
devexp_o_glm <- ( o_glm$null.deviance - o_glm$deviance ) / o_glm$null.deviance
devexp_o_glm / devexp_f_glm

# Forest plot
glm_to_forest_p( mod=o_glm, xmax=NULL )


#-------------------------------------------------------------------------------
#   LASSO
#-------------------------------------------------------------------------------

# Make LASSO model
o_las <- cv.glmnet( x=as.matrix( tr[ , ..ocols ] )[,-1], 
                    y=tr[["causal"]], family="binomial", standardize=FALSE )


#-------------------------------------------------------------------------------
#   XGBoost
#-------------------------------------------------------------------------------

# WITHOUT protein attenuation
o_xgb1 <- train_xgb( data=tr, feat_cols=ocols[-1], maxit=100, rm_prot_att=TRUE )
xgb_to_forest( xgb_mod=o_xgb1, suffix="", xmax=NULL )

# WITH protein attenuation
o_xgb2 <- train_xgb( data=tr, feat_cols=ocols[-1], maxit=100, rm_prot_att=FALSE )
xgb_to_forest( xgb_mod=o_xgb2, suffix="", xmax=NULL )


#-------------------------------------------------------------------------------
#   Assess performance of the models in the test set
#-------------------------------------------------------------------------------

# AUPRC: WITHOUT protein attenuation
o_pr1 <- auprc_test_set( test_df  = te, 
                         rand_mod = rand_mod,
                         glm_mod  = o_glm, 
                         las_mod  = o_las, 
                         xgb_mod  = o_xgb1, 
                         ymax     = 1, 
                         avg_prot_att = TRUE )
f_pr1; o_pr1; o_pr1/f_pr1

# AUPRC: WITH protein attenuation
o_pr2 <- auprc_test_set( test_df  = te, 
                         rand_mod = rand_mod,
                         glm_mod  = o_glm, 
                         las_mod  = o_las, 
                         xgb_mod  = o_xgb2, 
                         ymax     = 1, 
                         avg_prot_att = FALSE )
f_pr2; o_pr2; o_pr2/f_pr2

# PR curves: WITHOUT protein attenuation
plot_pr( model=o_glm,  test_df=te, p=0.78, r=0.34, v=0.75, avg_prot_att=TRUE )
plot_pr( model=o_xgb1, test_df=te, p=0.78, r=0.34, v=0.75, avg_prot_att=TRUE )

# PR curves: WITH protein attenuation
plot_pr( model=o_glm,  test_df=te, p=0.78, r=0.34, v=0.75, avg_prot_att=FALSE )
plot_pr( model=o_xgb2, test_df=te, p=0.78, r=0.34, v=0.75, avg_prot_att=FALSE )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Basic model: Only use features derived from 4 methods
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Model selection
#-------------------------------------------------------------------------------

# Hard-coded list of features from the final "reduced model"
# scols0 <- c( "causal",
#              "pops_bil",       "pops_glo",       "pops_rel",
#              "dist_gene_rel",
#              "magma_glo",      "magma_rel",
#              "coding_glo",     "coding_rel",
#              "prior_n_genes_locus", 
#              "prot_att_poly2", "prot_att_miss",
#              "rvis4", "rvis4_poly2" )

# Manually make models and drop terms that are not Bonferroni-significant
scols1 <- c( "causal",
             "pops_bil",       "pops_glo",       "pops_rel",
             "dist_gene_rel",
             "magma_rel",
             "coding_glo",     "coding_rel",
             "prior_n_genes_locus", 
             "prot_att_poly2", "prot_att_miss",
             "rvis4", "rvis4_poly2" )

# Run regression
s_form1 <- make_formula( lhs=scols1[1], rhs=scols1[-1] )
s_glm1  <- glm( formula=s_form1, data=tr, family="binomial" )

# Forest plot
glm_to_forest_p( mod=s_glm1, xmax=NULL )


#-------------------------------------------------------------------------------
#   GLM
#-------------------------------------------------------------------------------

# "Basic" list of selected faeatures
scols <- c( "causal",
            "pops_bil", "pops_glo", "pops_rel",
            "dist_gene_rel",
            "magma_rel",
            "coding_glo", "coding_rel",
            "prior_n_genes_locus",
            "prot_att_poly2", "prot_att_miss",
            "rvis4", "rvis4_poly2" )

# Raw feature correlations
corrplot( cor( tr[ , ..scols ][,-1] ), order="hclust" )

# Run regression
s_form <- make_formula( lhs=scols[1], rhs=scols[-1] )
s_glm  <- glm( formula=s_form, data=tr, family="binomial" )
# summary(s_glm)$coef
# saveRDS( object=s_glm, file=file.path( maindir, "scz_glm.rds" ) )

# Extract deviance explained
devexp_s_glm <- ( s_glm$null.deviance - s_glm$deviance ) / s_glm$null.deviance
devexp_s_glm / devexp_f_glm

# Forest plot
glm_to_forest_p( mod=s_glm )


#-------------------------------------------------------------------------------
#   LASSO
#-------------------------------------------------------------------------------

# Make LASSO model
s_las <- cv.glmnet( x=as.matrix( tr[ , ..scols ] )[,-1], 
                    y=tr[["causal"]], family="binomial", standardize=FALSE )


#-------------------------------------------------------------------------------
#   XGBoost
#-------------------------------------------------------------------------------

# WITHOUT protein attenuation
s_xgb1 <- train_xgb( data=tr, feat_cols=scols[-1], maxit=100, rm_prot_att=TRUE )
xgb_to_forest( xgb_mod=s_xgb1, suffix="", xmax=NULL )

# WITH protein attenuation
s_xgb2 <- train_xgb( data=tr, feat_cols=scols[-1], maxit=100, rm_prot_att=FALSE )
xgb_to_forest( xgb_mod=s_xgb2, suffix="", xmax=NULL )


#-------------------------------------------------------------------------------
#   Assess performance of the models in the test set
#-------------------------------------------------------------------------------

# AUPRC: WITHOUT protein attenuation
s_pr1 <- auprc_test_set( test_df  = te, 
                         rand_mod = rand_mod,
                         glm_mod  = s_glm, 
                         las_mod  = s_las, 
                         xgb_mod  = s_xgb1, 
                         ymax     = 1, 
                         avg_prot_att = TRUE )
f_pr1; s_pr1; s_pr1/f_pr1
o_pr1; s_pr1; s_pr1/o_pr1

# AUPRC: WITH protein attenuation
s_pr2 <- auprc_test_set( test_df  = te, 
                         rand_mod = rand_mod,
                         glm_mod  = s_glm, 
                         las_mod  = s_las, 
                         xgb_mod  = s_xgb2, 
                         ymax     = 1, 
                         avg_prot_att = FALSE )
f_pr2; s_pr2; s_pr2/f_pr2
o_pr2; s_pr2; s_pr2/o_pr2

# PR curves: WITHOUT protein attenuation
plot_pr( model=s_glm,  test_df=te, p=0.78, r=0.34, v=0.75, avg_prot_att=TRUE )
plot_pr( model=s_xgb1, test_df=te, p=0.78, r=0.34, v=0.75, avg_prot_att=TRUE )

# PR curves: WITH protein attenuation
plot_pr( model=s_glm,  test_df=te, p=0.78, r=0.34, v=0.75, avg_prot_att=FALSE )
plot_pr( model=s_xgb2, test_df=te, p=0.78, r=0.34, v=0.75, avg_prot_att=FALSE )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Calibrate model outputs
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Generate a variety of re-calibrated model predictions
#-------------------------------------------------------------------------------

# Get model predictions, plus several re-calibrations
adj_data_tr <- recalibrate_preds( test_df=tr, model=s_glm, avg_prot_att=FALSE )
adj_data_te <- recalibrate_preds( test_df=te, model=s_glm, avg_prot_att=FALSE )

# Train a LASSO model to meta-re-calibrate model outputs
adj_mod <- recalibration_model( data=adj_data_tr )

# Apply LASSO model to get meta-recalibrated predictions
pred_cols <- c( "global", "relative", "best" )
adj_data_tr$modeled <- predict( object=adj_mod, s="lambda.min", type="response",
                                newx=as.matrix( adj_data_tr[ , ..pred_cols ] ) )
adj_data_te$modeled <- predict( object=adj_mod, s="lambda.min", type="response",
                                newx=as.matrix( adj_data_te[ , ..pred_cols ] ) )


#-------------------------------------------------------------------------------
#   Plot AUPRCs and PR curves
#-------------------------------------------------------------------------------

# Compute test set AUPRC for each re-calibrated prediction
raw     <- pr.curve( scores.class0 = adj_data_te$original[  adj_data_te$causal ], 
                     scores.class1 = adj_data_te$original[ !adj_data_te$causal ], 
                     curve = TRUE )
scaled  <- pr.curve( scores.class0 = adj_data_te$scaled[  adj_data_te$causal ], 
                     scores.class1 = adj_data_te$scaled[ !adj_data_te$causal ], 
                     curve = TRUE )
scaledu  <- pr.curve( scores.class0 = adj_data_te$scaled_up[  adj_data_te$causal ], 
                      scores.class1 = adj_data_te$scaled_up[ !adj_data_te$causal ], 
                      curve = TRUE )
scaledd  <- pr.curve( scores.class0 = adj_data_te$scaled_down[  adj_data_te$causal ], 
                      scores.class1 = adj_data_te$scaled_down[ !adj_data_te$causal ], 
                      curve = TRUE )
flames  <- pr.curve( scores.class0 = adj_data_te$flames[  adj_data_te$causal ], 
                     scores.class1 = adj_data_te$flames[ !adj_data_te$causal ], 
                     curve = TRUE )
modeled <- pr.curve( scores.class0 = adj_data_te$modeled[  adj_data_te$causal ], 
                     scores.class1 = adj_data_te$modeled[ !adj_data_te$causal ], 
                     curve = TRUE )

# Barplot of AUPRCs
pr <- c( raw$auc.integral, flames$auc.integral, scaled$auc.integral,
         scaledu$auc.integral, scaledd$auc.integral, modeled$auc.integral )
names(pr) <- c( "Raw", "FLAMES", "Scaled", "ScaledUp", "ScaledDown", "Modeled" )
pr <- sort(pr)
par( mar=c(7,5,4,1) )
barplot( pr, col=viridis( length(pr) ), las=2 )

# Plot precision-recall curves
plot( raw,     las=1, scale.color=viridis(n=16, begin=0.3, end=1 ) )
plot( flames,  las=1, scale.color=viridis(n=16, begin=0.3, end=1 ) )
plot( scaled,  las=1, scale.color=viridis(n=16, begin=0.3, end=1 ) )
plot( modeled, las=1, scale.color=viridis(n=16, begin=0.3, end=1 ) )


#-------------------------------------------------------------------------------
#   Calibration plots
#-------------------------------------------------------------------------------

# Load libraries
library(tidymodels)
library(probably)
library(cowplot)

# Calibration plots: original
p1 <- cal_plot_logistic( .data=adj_data_tr, truth=causal, estimate=original ) +
  ggtitle("Training set, original predictions")
p2 <- cal_plot_logistic( .data=adj_data_te, truth=causal, estimate=original ) +
  ggtitle("Testing set, original predictions")

# Calibration plots: scaled down
p3 <- cal_plot_logistic( .data=adj_data_tr, truth=causal, estimate=scaled_down ) +
  ggtitle("Training set, locally-scaled down predictions")
p4 <- cal_plot_logistic( .data=adj_data_te, truth=causal, estimate=scaled_down ) +
  ggtitle("Testing set, locally-scaled down predictions")

# Calibration plots: modeled
p5 <- cal_plot_logistic( .data=adj_data_tr, truth=causal, estimate=modeled ) +
  ggtitle("Training set, meta-model predictions")
p6 <- cal_plot_logistic( .data=adj_data_te, truth=causal, estimate=modeled ) +
  ggtitle("Testing set, meta-model predictions")

# Original v. scaled
# Original v. modeled
plot_grid( p1, p2, p3, p4, ncol=2, nrow=2 )
plot_grid( p1, p2, p5, p6, ncol=2, nrow=2 )

# Plot original predictions v. modeled predictions, colouring by causal status
adj_data_te$col <- "black"
# adj_data_te$col[ adj_data_te$best == 1 ] <- "steelblue"
adj_data_te$col[ adj_data_te$causal ]    <- "tomato2"
plot( adj_data_te$original, adj_data_te$modeled, col=adj_data_te$col, cex=2, 
      lwd=2, xlim=c(0,1), ylim=c(0,1),
      xlab="Original P(causal)", ylab="Modeled P(causal)", las=1 )
abline( a=0, b=1, lwd=2, lty=2, col="grey" )


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Miscellaneous
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Does adding quadratic terms improve the model? No.
#-------------------------------------------------------------------------------

# Make quadratic terms for all features
tr$pops_glo_poly2      <- tr$pops_glo      * tr$pops_glo
tr$pops_rel_poly2      <- tr$pops_rel      * tr$pops_rel
tr$dist_gene_rel_poly2 <- tr$dist_gene_rel * tr$dist_gene_rel
tr$magma_rel_poly2     <- tr$magma_rel     * tr$magma_rel
tr$coding_glo_poly2    <- tr$coding_glo    * tr$coding_glo
tr$coding_rel_poly2    <- tr$coding_rel    * tr$coding_rel

# Which polynomial leads to the biggest improvement in model fit?
poly_lrt <- function(variable){
  u_form <- update( s_form, paste0( "~ . + ", variable ) )
  u_glm  <- glm( formula=u_form, data=tr, family="binomial" )
  anova( s_glm, u_glm, test="Chi" )$`Pr(>Chi)`[2]
}
poly_lrt( variable="pops_glo_poly2" )
poly_lrt( variable="pops_rel_poly2" )
poly_lrt( variable="dist_gene_rel_poly2" )
poly_lrt( variable="magma_rel_poly2" )
poly_lrt( variable="coding_glo_poly2" )
poly_lrt( variable="coding_rel_poly2" )

# "Basic" list of selected faeatures
qcols <- c( scols, "pops_glo_poly2" )

# Raw feature correlations
corrplot( cor( tr[ , ..qcols ][,-1] ), order="hclust" )

# Run regression
q_form <- make_formula( lhs=qcols[1], rhs=qcols[-1] )
q_glm  <- glm( formula=q_form, data=tr, family="binomial" )

# Extract deviance explained
devexp_q_glm <- ( q_glm$null.deviance - q_glm$deviance ) / q_glm$null.deviance
devexp_q_glm / devexp_s_glm

# Forest plot
glm_to_forest_p( mod=q_glm )


#-------------------------------------------------------------------------------
#   Does adding interaction terms improve the model? No.
#-------------------------------------------------------------------------------

# Add all possible interactions
i_form0 <- make_formula( lhs=scols[1], 
                         rhs=c( setdiff( scols[-1], cov_cols ), 
                                "prior_n_genes_locus" ) )
i_form0 <- update( i_form0, ~ (.)^2 )
i_form0 <- update( i_form0, ~ . + prot_att_poly2 + prot_att_miss + 
                     rvis4 + rvis4_poly2 )

# Perform backwards stepwise model selection
# Remove features predicted to have P > 0.05/10
i_glm0 <- glm( formula=i_form0, data=tr, family="binomial" )
i_glm1 <- stats::step( object=i_glm0, k=qchisq( 1-0.05/10, df=1 ) )

# Hard-coded list of features that passed model selection
# icols1 <- c( scols,
#             "pops_rel:prior_n_genes_locus",
#             "dist_gene_rel:coding_glo" )

# Manually make models and drop terms that are not Bonferroni-significant
icols <- c( scols, 
            "pops_rel:prior_n_genes_locus" )

# Run regression
i_form <- make_formula( lhs=icols[1], rhs=icols[-1] )
i_glm  <- glm( formula=i_form, data=tr, family="binomial" )

# Extract deviance explained
devexp_i_glm <- ( i_glm$null.deviance - i_glm$deviance ) / i_glm$null.deviance
devexp_i_glm / devexp_s_glm

# Forest plot
glm_to_forest_p( mod=i_glm, xmax=NULL )


#-------------------------------------------------------------------------------
#   Plot features v. residuals
#-------------------------------------------------------------------------------

# Plot RVIS v. residuals
# plot( tr$pubmed_log, residuals(s_glm),
#       xlim=quantile( tr$pubmed_log, probs=c( 0.05, 0.95 ) ),
#       ylim=quantile( residuals(s_glm), probs=c( 0.2, 0.8 ) ) )
# lines( lowess( tr$pubmed_log, residuals(s_glm) ), lwd=3, col="orange2" )

# Plot features v. residuals
plot_residuals <- function( data, model, variable ){
  plot( data[[variable]], residuals(model),
        xlim=quantile( data[[variable]], probs=c( 0.05, 0.95 ) ),
        ylim=quantile( residuals(model), probs=c( 0.05, 0.95 ) ) )
  lines( lowess( data[[variable]], residuals(model) ), lwd=3, col="orange2" )
}
plot_residuals( data=tr, model=s_glm, variable="rvis4" )
plot_residuals( data=tr, model=s_glm, variable="rvis4_poly2" )
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
au1 <- c( f_pr1["GLM"], o_pr1["GLM"], s_pr1["GLM"],
          f_pr1["LASSO"], o_pr1["LASSO"], s_pr1["LASSO"],
          f_pr1["XGBoost"], o_pr1["XGBoost"], s_pr1["XGBoost"] )
au2 <- c( f_pr2["GLM"],     o_pr2["GLM"],     s_pr2["GLM"],
          f_pr2["LASSO"],   o_pr2["LASSO"],   s_pr2["LASSO"],
          f_pr2["XGBoost"], o_pr2["XGBoost"], s_pr2["XGBoost"] )
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

# Libraries
library(tidymodels)
library(tidyr)

# Create a decision tree model specification
tree_spec <- decision_tree( tree_depth=3 ) %>%
  set_engine("rpart") %>%
  set_mode("classification")

# Fit the model to the training data
d_form <-  update( s_form, factor(.) ~ . - prot_att_poly2 - prot_att_miss - 
                     rvis4 - rvis4_poly2 )
tree_fit <- tree_spec %>%
  fit( d_form, data = tr[ tr$pops_and_nearest , ] )

# Plot the decision tree
library(rpart.plot)
rpart.plot( tree_fit$fit, type = 4, under = TRUE, box.palette = "auto")


#-------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////
#   Done
#///////////////////////////////////////////////////////////////////////////////
#-------------------------------------------------------------------------------


















