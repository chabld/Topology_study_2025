#adapted from Adam_G, Jun 30, 2022 at 22:32, in the thread "Understanding the process of tweaking contrasts in linear model fitting to show all levels". CC BY-SA 4.0
#https://stackoverflow.com/a/72822671

#' @description When categorical variables are set with sum contrast (contr.sum) before running gls, due to the nature of the contrast matrix coding, the last level will be "implicit" (last row made of -1) and its results will not be given by the gls() function. This script recomputes that last level's parameters and prints them along the rest of the levels.

#' @dataset The original dataset (as data.frame)
#' @model The gls() model object
#' @factor_name name of the categorical variable to print all levels from
#' @factor2_name optional: name of the second categorical variable the first categorical variable interacts with in the model

DevContrStats_gls <- function(dataset, model, factor_name, factor2_name=NULL) {
  
  #adapts whether dataset is given as matrix or data.frame
  get_column <- function(dataset, colname) {
    if (is.data.frame(dataset)) {return(dataset[[colname]])} 
    else if (is.matrix(dataset)) {return(dataset[, colname])} 
    else {stop("dataset must be a matrix or df")}}
  factor1 <- get_column(dataset, factor_name)
  
  #function to detect if factor is set as contr.sum or contr.treatment
  is_sum_contrast <- function(x) {
    contr <- contrasts(x)
    if (is.null(contr)) return(FALSE)
    # Check if the last row is all -1 and the rest are identity-like
    last_row <- tail(contr, 1)
    all(last_row == -1) &&
      all(apply(head(contr, -1), 2, function(col) sum(col == 1) == 1 && sum(col == 0) == (length(col) - 1)))
  }
  
  if (!is.null(factor2_name)) {
    message('factor2_name given. Will assume user wanted the factor1:factor2 interaction term levels')
    factor2 <- get_column(dataset, factor2_name)
    
    #will not compute matrices the same if the second variables in the interaction is not contr.sum coded but contr.treatment, but simply add the ommited last contr.sum level interacting with it
    if (is_sum_contrast(dataset[[factor_name]]) && !is_sum_contrast(dataset[[factor2_name]]))
    {   
    
    #get interaction terms from model summary (no need to change those)
    full_table <- as.data.frame(summary(model)$tTable)
    interaction_rows <- grep(paste0("^", factor_name, "\\d+:", factor2_name), rownames(full_table))
    interaction_table <- full_table[interaction_rows, , drop = FALSE]
    
    #Dynamically reconstruct the dropped interaction term
    #that level is the one ommitted as coded as all -1 
    dropped_level <- rownames(contrasts(dataset[[factor_name]]))[apply(contrasts(dataset[[factor_name]]), 1, function(x) all(x == -1))]
    # Identify the non-reference level of factor2
    non_ref_level <- colnames(contrasts(dataset[[factor2_name]]))[1]
    
    # Get the estimated interaction terms
    int_terms <- rownames(interaction_table)
    # Build contrast weights: all -1 for non omitted terms, because the dropped level is the negative sum of the others in contr.sum
    contrast_weights <- rep(-1, length(int_terms))
    names(contrast_weights) <- int_terms
    
    #build up coefficients
    # Reconstruct value, basically, the negative sum of all other coeffs under contr.sum
    estimated_coefs <- coef(model)[int_terms]
    reconstructed_value <- sum(contrast_weights * estimated_coefs)
    # Reconstruct SE
    vcov_sub <- vcov(model)[int_terms, int_terms]
    reconstructed_se <- sqrt(t(contrast_weights) %*% vcov_sub %*% contrast_weights)
    
    # t and p
    reconstructed_t <- reconstructed_value / reconstructed_se
    reconstructed_p <- 2 * pt(abs(reconstructed_t), df = model$dims$N - model$dims$p, lower.tail = FALSE)
    
    reconstructed_row <- c(reconstructed_value,
                           reconstructed_se,
                           reconstructed_t,
                           reconstructed_p)
    reconstructed_label <- paste0(factor_name, 
                                  length(levels(dataset[[factor_name]])),":", factor2_name, non_ref_level)
    
      stats.tab_b_contr <- rbind(interaction_table, reconstructed_row)
      row.names(stats.tab_b_contr)[nrow(stats.tab_b_contr)]=reconstructed_label
      
      #fix names
      # Rename rows like category1, category2, etc.
      rownames(stats.tab_b_contr) <- sapply(rownames(stats.tab_b_contr), function(term) {
        if (grepl(paste0("^",factor_name,"\\d+"), term)) {
          idx <- as.integer(sub(paste0(factor_name,"(\\d+):.*"), "\\1", term))
          gsub(pattern = paste0(factor_name,idx), levels(dataset[[factor_name]])[idx], term)
        }
      })
      
      return(stats.tab_b_contr)
      
    } else {
      
    #factor matrix
    combo <- expand.grid(
      levels(as.factor(factor1)),
      levels(as.factor(factor2)))
    rownames(combo) <- paste(combo[,1], combo[,2], sep = ":")
    colnames(combo) = c(factor_name, factor2_name)
    
    #Build a model matrix for interaction terms,
    #including the dropped last levels of each categorical variable
    fctr.mtx <- model.matrix(
      as.formula(paste0("~ ", factor_name, " * ", factor2_name)),
      data = combo,
      contrasts.arg = setNames(
        list(contr.sum(length(levels(as.factor(factor1)))),
             contr.sum(length(levels(as.factor(factor2))))),
        c(factor_name, factor2_name)
      )
    )
    
    #get index of interaction terms
    int_idx <- grep(paste0(factor_name, "\\d:", factor2_name), names(coef(model)))
    
    #Reassign names to fctr.mtx
    valid_cols <- intersect(colnames(fctr.mtx), names(coef(model)))
    fctr.mtx <- fctr.mtx[, valid_cols, drop = FALSE]
 
    # Recompute interaction estimates
    # Coefficients After Contrasts
    coef_a_contr <- coef(model)[colnames(fctr.mtx)]
    # Coefficients Before Contrasts
    ##coef_b_contr tells how each factor level deviates from the grand mean.
    coef_b_contr <- as.vector(fctr.mtx %*% coef_a_contr)
    
    # Covariance matrix After Contrasts
    var_a_contr <- vcov(model)[colnames(fctr.mtx), colnames(fctr.mtx)]
    # Covariance matrix Before Contrasts
    var_b_contr <- fctr.mtx %*% var_a_contr %*% t(fctr.mtx)
    }
    
  } else 
  {
    #factor matrix
    N <- nlevels(as.factor(factor1))
    Cmat <- contr.sum(N)
    dimnames(Cmat) <- list(levels(as.factor(factor1)), seq_len(N - 1))
    fctr.mtx <- Cmat
      
    ## coefficients After Contrasts
    coef_a_contr <- coef(model)[grep(paste0("^", factor_name, "\\d+$"), names(coef(model)))]
    
    ## coefficients Before Contrasts
    ## coef_bc tells how each factor level deviates from the grand mean.
    coef_b_contr <- (fctr.mtx %*% coef_a_contr)[, 1]
    
    ## Covariance matrix after contrasts:
    var_a_contr <- vcov(model)[grep(paste0("^", factor_name, "\\d+$"), row.names(vcov(model)))-1, grep(paste0("^", factor_name, "\\d+$"), colnames(vcov(model)))-1]
    
    ## Transform covariance matrix (after contrasts) to the one before contrasts:
    ## The diagonal of var_bc gives the estimated variance of factor-level deviations.
    var_b_contr <- fctr.mtx %*% var_a_contr %*% t(fctr.mtx)
  }
  
  
  ## standard error of point estimate `coef_bc`
  std.err_b_contr <- sqrt(diag(var_b_contr))
  
  ## t-statistics (Null Hypothesis: coef_bc = 0)
  t.stats_b_contr <- coef_b_contr / std.err_b_contr
  
  ## p-values of the t-statistics
  p.value_b_contr <- 2 * pt(abs(t.stats_b_contr),
                            df = model$dims$N - model$dims$p,
                            lower.tail = FALSE)
  
  ## construct a coefficient table that mimics `coef(summary(fit))`
  stats.tab_b_contr <- cbind(coef_b_contr,
                             std.err_b_contr,
                             t.stats_b_contr,
                             p.value_b_contr)
  colnames(stats.tab_b_contr)=colnames(colnames(summary(model)$tTable))
  
  ## extract statistics of the intercept, which = grand mean
  intercept.stats <- coef(summary(model))[1, , drop = FALSE]
  
  ## augment the coefficient table with the intercept stats
  stats.tab <- rbind(intercept.stats, stats.tab_b_contr)

  ## print stats table with significance stars
  return(as.data.frame(stats.tab))
}
