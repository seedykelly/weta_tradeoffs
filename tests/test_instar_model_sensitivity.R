#!/usr/bin/env Rscript
# Independent numerical checks of the completed sensitivity analysis.
# Run: Rscript tests/test_instar_model_sensitivity.R [--data=...] [--output=...]

script_arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
test_path <- sub('^--file=','',script_arg)
if(!file.exists(test_path)) test_path <- gsub('~+~',' ',test_path,fixed=TRUE)
root <- dirname(dirname(normalizePath(test_path,mustWork=TRUE)))
Sys.setenv(WETA_INSTAR_HELPERS_ONLY='1')
source(file.path(root,'scripts','instar_model_sensitivity.R'))
Sys.unsetenv('WETA_INSTAR_HELPERS_ONLY')
opt <- parse_options(commandArgs(trailingOnly=TRUE),root)
d <- prepare_data(opt$data)
fits <- readRDS(file.path(opt$output,'fitted_models.rds'))
results <- read.csv(file.path(opt$output,'contrast_comparison.csv'))
comparison <- read.csv(file.path(opt$output,'model_comparison_ML.csv'))
manifest <- read.csv(file.path(opt$output,'sample_manifest.csv'))
expected_n <- c(foreleg=313L,midleg=312L,hindleg=311L,ear_lin=313L,eye=313L)
independent_p <- list()
assert_close <- function(x,y,tol=1e-7) stopifnot(max(abs(x-y))<tol)

for(trait in names(fits)) {
  dt <- d[complete.cases(d[,c(trait,'logP_c')]),]
  dt$y <- log(dt[[trait]])
  stopifnot(nrow(dt)==expected_n[[trait]])
  sm <- manifest[manifest$trait==trait,]
  stopifnot(identical(sm$ID,d$ID),identical(sm$group,as.character(d$group)),
            identical(sm$ID[sm$included],dt$ID))

  # M0 must reproduce ordinary least squares, including its uncertainty.
  base <- lm(y ~ logP_c + instar_c + male,data=dt)
  original <- fits[[trait]]$M0_original$REML
  assert_close(coef(original),coef(base))
  assert_close(vcov(original),vcov(base))
  original_rows <- results[results$trait==trait & results$model=='M0_original',]
  direct <- c(average_male_step=coef(base)['instar_c'],matched_sex=coef(base)['male'],
              step_minus_sex=coef(base)['instar_c']-coef(base)['male'])
  assert_close(original_rows$estimate_log,unname(direct))

  # M3 is independently recoverable from four separate group regressions.
  independent <- lapply(split(dt,dt$group),function(z)lm(y~logP_c,data=z))
  free <- fits[[trait]]$M3_free_groups$REML
  grid <- reference_grid(dt);X <- design_matrix(free,grid)
  intercepts <- vapply(independent,function(m)unname(coef(m)[1]),numeric(1))
  variances <- vapply(independent,function(m)vcov(m)[1,1],numeric(1))
  assert_close(drop(X%*%coef(free)),intercepts)
  assert_close(X%*%vcov(free)%*%t(X),diag(variances),tol=1e-8)
  assert_close(residual_sd(free,grid),vapply(independent,sigma,numeric(1)),tol=1e-6)
  W <- contrast_weights(TRUE)
  for(i in seq_len(nrow(W))) {
    row <- results[results$trait==trait & results$model=='M3_free_groups' &
                    results$contrast==rownames(W)[i],]
    est <- sum(W[i,]*intercepts);parts <- W[i,]^2*variances;se <- sqrt(sum(parts))
    assert_close(row$estimate_log,est)
    assert_close(row$SE,se,tol=1e-6)
    # Independently calculate Welch-Satterthwaite df from separate regressions.
    df <- sum(parts)^2/sum(parts^2/vapply(independent,df.residual,numeric(1)))
    assert_close(row$df,df)
    independent_p[[length(independent_p)+1L]] <- data.frame(trait=trait,
      contrast=rownames(W)[i],p=2*pt(abs(est/se),df,lower.tail=FALSE))
  }

  for(name in names(fits[[trait]])) {
    ff <- fits[[trait]][[name]]
    stopifnot(ff$ML$method=='ML',ff$REML$method=='REML',
              nobs(ff$ML)==nrow(dt),nobs(ff$REML)==nrow(dt))
    # Verify that named group SDs reproduce nlme's normalised residuals.
    r <- (dt$y-as.numeric(predict(ff$REML,newdata=dt)))/residual_sd(ff$REML,dt)
    assert_close(r,as.numeric(residuals(ff$REML,type='normalized')))
    cmp <- comparison[comparison$trait==trait & comparison$model==name,]
    assert_close(cmp$AIC,AIC(ff$ML))
    k <- attr(logLik(ff$ML),'df')
    assert_close(cmp$AICc,AIC(ff$ML)+2*k*(k+1)/(nrow(dt)-k-1))
  }
}
independent_p <- do.call(rbind,independent_p)
independent_p$p_holm <- ave(independent_p$p,independent_p$contrast,FUN=function(x)p.adjust(x,'holm'))
for(i in seq_len(nrow(independent_p))) {
  row <- independent_p[i,]
  reported <- results[results$trait==row$trait & results$model=='M3_free_groups' &
                        results$contrast==row$contrast,]
  stopifnot((row$p_holm<.05)==reported$holm_below_05)
}
verify_original(results,opt$reference)
stopifnot(nrow(results)==75L,nrow(comparison)==20L,
          all(is.finite(results$SE)),all(results$lower_95_log<results$upper_95_log))
message('PASS: original OLS and archived results; independent group estimates and uncertainty; ',
        'independent Welch inference; sample IDs/morphs; residual SD mapping; ML comparisons.')
