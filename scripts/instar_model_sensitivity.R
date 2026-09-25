#!/usr/bin/env Rscript
# Frequentist sensitivity of the instar-matched decomposition.
# Run with source() or Rscript. Defaults resolve from this script's project.
# Rscript scripts/instar_model_sensitivity.R [--data=...] [--output=...]
#   [--reference=.../instar_matched_decomposition.csv] [--replicates=2000]
#
# Four diagnostic-motivated models are reported for each trait:
# M0: original linear male-instar gradient, common size slope and residual SD.
# M1: M0 with separate residual SDs for the four groups.
# M2: M1 with group-specific body-size slopes, retaining the linear instar
#     constraint on group means at pronotum length 7.2 mm.
# M3: M2 with unrestricted group means (ninth instar need not be the midpoint).
# Model comparisons use ML on identical rows; inference uses REML. This is a
# per-trait sensitivity, not a replacement for the joint repeated-leg model.
# Recorded morph labels are preserved; no thresholds are used to assign them.
#
# In M3 there is no single constant per-instar coefficient. The comparable
# quantity is half the adjusted eighth-to-tenth log contrast: the average male
# change per instar across those two steps. Adjacent steps are reported too.
# No model isolates causal effects of development, sex or weapon expression.
# New outputs remain separate from the maintained manuscript result tables.

resolve_script <- function() {
  paths <- Filter(Negate(is.null),lapply(sys.frames(),function(x)x$ofile))
  if(length(paths)) return(normalizePath(tail(paths,1L)[[1L]],mustWork=TRUE))
  arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
  if(length(arg)!=1L) stop('Run with Rscript or source().')
  path <- sub('^--file=','',arg)
  if(!file.exists(path)) path <- gsub('~+~',' ',path,fixed=TRUE)
  normalizePath(path,mustWork=TRUE)
}

parse_options <- function(args,root) {
  opt <- list(data=file.path(root,'data','trait_data.csv'),
    output=file.path(root,'analysis_outputs','instar_model_sensitivity'),
    reference=file.path(root,'analysis_outputs','weta_trait_analysis','tables','instar_matched_decomposition.csv'),
    replicates='2000',seed='20260924')
  for(arg in args) {
    if(!grepl('^--[^=]+=.+$',arg)) stop('Use --name=value: ',arg)
    name <- sub('^--([^=]+)=.*$','\\1',arg)
    if(!name %in% names(opt)) stop('Unknown option: ',name)
    opt[[name]] <- sub('^--[^=]+=','',arg)
  }
  for(name in c('replicates','seed')) {
    x <- suppressWarnings(as.numeric(opt[[name]]))
    if(!is.finite(x) || x!=floor(x) || x<1 || x>.Machine$integer.max-1000)
      stop('Invalid ',name)
    opt[[name]] <- as.integer(x)
  }
  if(opt$replicates<1000) stop('Use at least 1,000 diagnostic replicates.')
  opt
}

prepare_data <- function(path) {
  d <- read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
  required <- c('ID','sex','morph','pronotum','forefemur','foretibia',
                'midfemur','midtibia','hindfemur','hindtibia','ear','eye')
  if(!all(required %in% names(d))) stop('Missing columns: ',paste(setdiff(required,names(d)),collapse=', '))
  if(anyNA(d$ID) || any(!nzchar(trimws(d$ID))) || anyDuplicated(d$ID) || any(d$ID!=trimws(d$ID)))
    stop('Invalid or duplicated specimen IDs.')
  groups <- c('female','eighth','ninth','tenth')
  if(anyNA(d$morph) || !all(d$morph %in% groups) || anyNA(d$sex) || !all(d$sex %in% c('f','m')) ||
     any((d$sex=='f')!=(d$morph=='female'))) stop('Invalid morph/sex metadata.')
  for(n in setdiff(required,c('ID','sex','morph'))) {
    if(!is.numeric(d[[n]]) || any(!is.na(d[[n]]) & (!is.finite(d[[n]]) | d[[n]]<=0)))
      stop('Invalid measurements: ',n)
  }
  d$group <- factor(d$morph,levels=groups)
  d$logP_c <- log(d$pronotum/7.2)
  d$instar_c <- unname(c(female=0,eighth=-2,ninth=-1,tenth=0)[d$morph])
  d$male <- as.integer(d$sex=='m')
  d$foreleg <- d$forefemur+d$foretibia
  d$midleg <- d$midfemur+d$midtibia
  d$hindleg <- d$hindfemur+d$hindtibia
  d$ear_lin <- sqrt(d$ear)
  d
}

model_formulas <- function() {
  rhs <- c(M0_original='logP_c + instar_c + male',
    M1_group_variance='logP_c + instar_c + male',
    M2_group_slopes='instar_c + male + group:logP_c',
    M3_free_groups='group * logP_c')
  lapply(rhs,function(x)as.formula(paste('y ~',x),env=baseenv()))
}

reference_grid <- function(d) data.frame(
  group=factor(levels(d$group),levels=levels(d$group)),
  logP_c=0,instar_c=c(0,-2,-1,0),male=c(0,1,1,1))

contrast_weights <- function(extra=FALSE) {
  w <- rbind(average_male_step=c(0,-.5,0,.5),matched_sex=c(-1,0,0,1),
             step_minus_sex=c(1,-.5,0,-.5))
  if(extra) w <- rbind(w,ninth_minus_eighth=c(0,-1,1,0),
                      tenth_minus_ninth=c(0,0,-1,1),ninth_minus_midpoint=c(0,-.5,1,-.5))
  colnames(w) <- c('female','eighth','ninth','tenth')
  w
}

fit_gls <- function(formula,d,name,method) {
  # Literal call arguments let emmeans recover/refit the model without relying
  # on temporary loop variables or the caller's current working directory.
  fit <- do.call(nlme::gls,list(model=formula,data=d,
    weights=if(name=='M0_original') NULL else nlme::varIdent(form=~1|group),
    method=method,na.action=na.fail,
    control=nlme::glsControl(maxIter=200,msMaxIter=200,tolerance=1e-8,msTol=1e-8)))
  if(any(!is.finite(coef(fit))) || any(!is.finite(vcov(fit)))) stop('Nonfinite fit: ',name)
  if(name!='M0_original' && (!is.matrix(fit$apVar) || any(!is.finite(fit$apVar))))
    stop('Variance-parameter uncertainty could not be estimated: ',name)
  fit
}

design_matrix <- function(fit,grid) {
  X <- model.matrix(delete.response(terms(formula(fit))),grid,contrasts.arg=fit$contrasts)
  if(!setequal(colnames(X),names(coef(fit)))) stop('Coefficient/design mismatch.')
  X[,names(coef(fit)),drop=FALSE]
}

contrast_summary <- function(fit,d,name,trait,seed) {
  set.seed(seed)
  if(name=='M3_free_groups') {
    # This model separates into four ordinary regressions. Use their analytic
    # Welch-Satterthwaite calculation rather than a numerical Hessian for df.
    separate <- lapply(split(d,d$group),function(z)lm(y~logP_c,data=z))
    group_variance <- vapply(separate,function(m)vcov(m)[1,1],numeric(1))
    group_df <- vapply(separate,df.residual,numeric(1))
  } else {
    basis <- emmeans::ref_grid(fit,data=d,
      mode=if(name=='M0_original') 'df.error' else 'satterthwaite')
    stopifnot(isTRUE(all.equal(unname(basis@bhat),unname(coef(fit)),tolerance=1e-10)))
  }
  X <- design_matrix(fit,reference_grid(d))
  W <- contrast_weights(name=='M3_free_groups')
  do.call(rbind,lapply(seq_len(nrow(W)),function(i) {
    L <- drop(W[i,]%*%X)
    est <- sum(L*coef(fit));se <- sqrt(drop(t(L)%*%vcov(fit)%*%L))
    # emmeans' coefficient-space degrees-of-freedom function retains the GLS
    # variance-parameter uncertainty. Do not substitute the residual N-p df.
    df <- if(name=='M3_free_groups') {
      parts <- W[i,]^2*group_variance
      sum(parts)^2/sum(parts^2/group_df)
    } else basis@dffun(L,basis@dfargs)
    if(!is.finite(df) || df<=0 || !is.finite(se) || se<=0) stop('Invalid contrast uncertainty.')
    lower <- est-qt(.975,df)*se;upper <- est+qt(.975,df)*se
    data.frame(trait=trait,model=name,contrast=rownames(W)[i],n=nrow(d),
      estimate_log=est,SE=se,df=df,lower_95_log=lower,upper_95_log=upper,
      p_value=2*pt(abs(est/se),df,lower.tail=FALSE),
      percent=100*expm1(est),lower_95_percent=100*expm1(lower),upper_95_percent=100*expm1(upper),
      units=if(rownames(W)[i]=='step_minus_sex') 'percent ratio of multiplicative effects (not percentage points)' else
        if(trait=='ear_lin') 'percent linearised tympanic size' else 'percent length',
      df_method=if(name=='M0_original') 'residual (matches OLS)' else
        if(name=='M3_free_groups') 'Welch-Satterthwaite (independent group regressions)' else 'Satterthwaite (emmeans)')
  }))
}

residual_sd <- function(fit,d) {
  vs <- fit$modelStruct$varStruct
  if(is.null(vs)) return(rep(fit$sigma,nrow(d)))
  ratios <- coef(vs,unconstrained=FALSE,allCoef=TRUE)
  sd <- fit$sigma*unname(ratios[as.character(d$group)])
  if(any(!is.finite(sd) | sd<=0)) stop('Could not map residual SDs to groups.')
  sd
}

simulation_checks <- function(fit,d,name,trait,reps,seed) {
  # Conditional plug-in simulations at observed covariates. These are model-fit
  # descriptions, not posterior predictive distributions or formal hypothesis
  # tests: fitted-parameter uncertainty and diagnostic multiplicity are omitted.
  set.seed(seed)
  mu <- as.numeric(predict(fit,newdata=d));sd <- residual_sd(fit,d)
  yrep <- matrix(rnorm(nrow(d)*reps),nrow=nrow(d),ncol=reps)*sd+mu
  normres <- (d$y-mu)/sd
  do.call(rbind,lapply(levels(d$group),function(g) {
    ii <- which(d$group==g)
    means <- colMeans(yrep[ii,,drop=FALSE]);sds <- apply(yrep[ii,,drop=FALSE],2,stats::sd)
    qm <- quantile(means,c(.025,.975));qs <- quantile(sds,c(.025,.975))
    data.frame(trait=trait,model=name,group=g,n=length(ii),
      observed_mean_log=mean(d$y[ii]),fitted_mean_log=mean(mu[ii]),
      simulated_mean_lower=qm[1],simulated_mean_upper=qm[2],
      observed_sd_log=stats::sd(d$y[ii]),simulated_sd_lower=qs[1],simulated_sd_upper=qs[2],
      mean_in_simulation_band=mean(d$y[ii])>=qm[1] & mean(d$y[ii])<=qm[2],
      sd_in_simulation_band=stats::sd(d$y[ii])>=qs[1] & stats::sd(d$y[ii])<=qs[2],
      residual_sd_model=sd[ii[1]],normalised_residual_mean=mean(normres[ii]),
      normalised_residual_sd=stats::sd(normres[ii]))
  }))
}

verify_original <- function(results,path) {
  archived <- read.csv(path)
  map <- c(average_male_step='instar_component',matched_sex='weapon_component',step_minus_sex='component_difference')
  for(i in which(results$model=='M0_original')) {
    row <- results[i,];old <- archived[archived$trait==row$trait,]
    if(nrow(old)!=1L) stop('Missing original trait: ',row$trait)
    prefix <- map[[row$contrast]]
    estcol <- if(prefix=='instar_component') 'instar_component_log_per_instar' else paste0(prefix,'_log')
    stopifnot(abs(row$estimate_log-old[[estcol]])<1e-8,
              abs(row$p_value-old[[paste0(prefix,'_p')]])<1e-8,
              row$n==old$n,abs(row$df-old$df_residual)<1e-8)
  }
}

comparison_plot <- function(results) {
  d <- results[results$contrast=='step_minus_sex',]
  traits <- c('foreleg','midleg','hindleg','ear_lin','eye')
  models <- names(model_formulas())
  colours <- c('#777777','#009E73','#E69F00','#0072B2')
  par(mar=c(6.0,7,4,1),family='sans')
  plot(NA,xlim=range(c(d$lower_95_log,d$upper_95_log)),ylim=c(.4,5.6),yaxt='n',ylab='',
    xlab='',
    main='Does the instar gradient exceed the matched-sex difference?')
  axis(2,at=5:1,labels=c('Foreleg','Midleg','Hindleg','Tympanic size','Eye'),las=1)
  abline(v=0,lty=2,col='grey60')
  for(j in seq_along(models)) {
    z <- d[d$model==models[j],];y <- 6-match(z$trait,traits)+(2.5-j)*.13
    segments(z$lower_95_log,y,z$upper_95_log,y,col=colours[j],lwd=1.8)
    points(z$estimate_log,y,pch=c(1,17,15,16)[j],col=colours[j],cex=1.05)
  }
  legend('bottomright',legend=c('Original','Separate residual SDs','Also separate size slopes','Also free group means'),
    col=colours,pch=c(1,17,15,16),bty='o',box.lty=0,bg='white',cex=.85)
  mtext('Average male step minus matched-sex contrast (log scale)',side=1,line=2.7)
  mtext('Negative: matched-sex contrast larger     |     Positive: average male step larger',side=1,line=4.2,cex=.8)
}

diagnostic_plots <- function(models,data,out) {
  pdf(file.path(out,'residual_diagnostics.pdf'),width=10,height=11)
  on.exit(dev.off(),add=TRUE)
  colours <- c('#777777','#0072B2','#009E73','#D55E00')
  par(mfrow=c(3,2),mar=c(4,4,3,1),oma=c(0,0,2,0))
  for(trait in names(models)) {
    d <- data[[trait]]
    for(name in c('M0_original','M2_group_slopes','M3_free_groups')) {
      fit <- models[[trait]][[name]]$REML
      mu <- as.numeric(predict(fit,newdata=d));r <- (d$y-mu)/residual_sd(fit,d)
      plot(mu,r,col=adjustcolor(colours[as.integer(d$group)],.65),pch=16,
        xlab='Fitted log trait',ylab='Normalised residual',main=name)
      abline(h=0,lty=2)
      if(name=='M0_original') legend('bottomleft',legend=levels(d$group),col=colours,pch=16,bty='n',cex=.7)
      qqnorm(r,main=paste(name,'normal Q-Q'),pch=16,cex=.55);qqline(r,col='#D55E00')
    }
    mtext(trait,outer=TRUE,font=2)
  }
}

write_csv <- function(x,path) write.csv(x,path,row.names=FALSE,na='')

write_summary <- function(results,comparison,tests,checks,out) {
  labels <- c(foreleg='Foreleg',midleg='Midleg',hindleg='Hindleg',
              ear_lin='Tympanic linear size',eye='Eye')
  p_text <- function(p) ifelse(p<.001,formatC(p,format='e',digits=2),sprintf('%.3f',p))
  get <- function(trait,model,contrast) {
    x <- results[results$trait==trait & results$model==model & results$contrast==contrast,]
    stopifnot(nrow(x)==1L);x
  }
  lines <- c('# Flexible instar-model sensitivity: results','',
    'The analysis reproduced the original frequentist decomposition before relaxing its assumptions. The most consequential change is allowing the relationship with body size to differ among the four groups.',
    '', '## Effect estimates',
    'All comparisons below are adjusted to a pronotum length of 7.2 mm, which is inside the observed range of every group. The male step is the average eighth-to-tenth change per instar. The matched-sex contrast compares tenth-instar males with females. Percentages refer to length or linearised tympanic size.',
    '', '| Trait | Original male step | Original matched sex | M2 male step | M2 matched sex | M2 direct-difference Holm p | M3 direct-difference Holm p |',
    '|---|---:|---:|---:|---:|---:|---:|')
  for(t in names(labels)) {
    values <- c(get(t,'M0_original','average_male_step')$percent,
                get(t,'M0_original','matched_sex')$percent,
                get(t,'M2_group_slopes','average_male_step')$percent,
                get(t,'M2_group_slopes','matched_sex')$percent)
    lines <- c(lines,paste0('| ',labels[t],' | ',paste(sprintf('%.2f%%',values),collapse=' | '),' | ',
      p_text(get(t,'M2_group_slopes','step_minus_sex')$p_holm),' | ',
      p_text(get(t,'M3_free_groups','step_minus_sex')$p_holm),' |'))
  }
  lines <- c(lines,'',
    'M2 permits separate body-size slopes and residual SDs, while retaining equal male instar steps at the reference body size. M3 additionally frees the ninth-instar mean; its two adjacent steps need not be equal. Both models, and the intermediate variance-only model, are reported in full rather than selecting by hypothesis-test significance.',
    '', '## What changes for hindlegs?')
  for(m in c('M0_original','M1_group_variance','M2_group_slopes','M3_free_groups')) {
    h <- get('hindleg',m,'step_minus_sex')
    lines <- c(lines,sprintf('- %s: direct difference %.4f log units; pointwise 95%% CI %.4f to %.4f; Holm p = %s.',
      m,h$estimate_log,h$lower_95_log,h$upper_95_log,p_text(h$p_holm)))
  }
  h2 <- get('hindleg','M2_group_slopes','step_minus_sex')
  if(h2$p_holm>=.05) lines <- c(lines,'',
    'The estimated hindleg instar step remains larger, but the more flexible models do not clearly establish that it exceeds the matched-sex contrast. This reflects both a smaller estimated difference and uncertainty around it. It is not evidence that the two effects are equivalent.',
    'The original claim that the hindleg instar change exceeds the sex contrast should therefore be qualified as dependent on the common body-size slope assumption.')
  other <- setdiff(names(labels),'hindleg')
  stable <- all(vapply(other,function(t)all(vapply(c('M0_original','M1_group_variance','M2_group_slopes','M3_free_groups'),
    function(m){z<-get(t,m,'step_minus_sex');z$estimate_log<0 && z$p_holm<.05},logical(1))),logical(1)))
  if(stable) lines <- c(lines,'',
    'For forelegs, midlegs, tympanic size and eyes, the matched-sex contrast remains larger than the average male instar step under every fitted specification, including after Holm adjustment.')
  lines <- c(lines,'', '## Which changes improve fit?',
    '', '| Trait | Lowest AICc among candidates | Original minus best AICc | M2 minus best AICc | M3 minus best AICc |',
    '|---|---|---:|---:|---:|')
  for(t in names(labels)) {
    z <- comparison[comparison$trait==t,]
    lines <- c(lines,sprintf('| %s | %s | %.2f | %.2f | %.2f |',labels[t],z$model[which.min(z$AICc)],
      z$delta_AICc[z$model=='M0_original'],z$delta_AICc[z$model=='M2_group_slopes'],z$delta_AICc[z$model=='M3_free_groups']))
  }
  lines <- c(lines,'',
    'AICc values compare models fitted by maximum likelihood to the same observations within each trait. The lowest value is a relative comparison among these candidates, not proof that a model is true.',
    '', '| Trait | Add residual SDs: Holm p | Add size slopes: Holm p | Free group means: Holm p |',
    '|---|---:|---:|---:|')
  for(t in names(labels)) {
    z <- tests[tests$trait==t,]
    lines <- c(lines,paste0('| ',labels[t],' | ',paste(p_text(z$p_holm),collapse=' | '),' |'))
  }
  lines <- c(lines,'',
    'Each column tests a successive extension of the previous model. The evidence for different leg slopes is stronger than the evidence for different residual variances or for freeing the ninth-instar mean. The original specification has the lowest AICc for tympanic size and eyes among the candidates examined.',
    '', '## Diagnostics and limits',
    'The conditional simulations check whether observed group means and SDs fall within 95% bands from simulations at the observed body sizes. These use fitted parameters, omit their uncertainty, and are descriptive checks rather than formal adequacy tests.',
    '', '| Model | Group means outside bands | Group SDs outside bands | Checks per column |',
    '|---|---:|---:|---:|')
  for(m in names(model_formulas())) {
    z <- checks[checks$model==m,]
    lines <- c(lines,sprintf('| %s | %d | %d | %d |',m,sum(!z$mean_in_simulation_band),sum(!z$sd_in_simulation_band),nrow(z)))
  }
  lines <- c(lines,'',
    'The group-mean checks in M3 are satisfied by construction because its group intercepts are free. Residual plots are supplied for the original, group-slope and free-group models; passing the group summaries does not guarantee every distributional assumption.',
    'This analysis is exploratory because it was prompted by diagnostics. Confidence intervals are pointwise; Holm p-values control five traits separately for each contrast/model, not the entire process of examining several models. The contrasts describe associations and do not isolate causal effects of development, sex or weapon investment.',
    'With separate body-size slopes, the contrasts are specific to the stated 7.2 mm reference. This analysis concerns the secondary instar-sex decomposition and does not re-run or replace the main multivariate analyses or the joint repeated-leg model.',
    '', '## Re-running',
    'From R/RStudio, source the project script `scripts/instar_model_sensitivity.R`. It locates its input and output folders from the script location, reproduces the original baseline, and runs all four models for all five traits. Required R packages: nlme and emmeans. Default diagnostic simulations: 2,000 per fitted model, with a recorded random seed.',
    'Detailed estimates and all confidence intervals are in `contrast_comparison.csv`; model definitions, methods and output descriptions are in `READ_ME_FIRST.md`. Fitted models and software versions are saved alongside them.')
  # Keep a blank line after headings and between explanatory paragraphs.
  rendered <- character()
  for(line in lines) {
    if(length(rendered) && nzchar(line) && nzchar(tail(rendered,1)) &&
       (startsWith(tail(rendered,1),'#') ||
        (!grepl('^[|#-]',line) && !grepl('^[|#-]',tail(rendered,1))))) rendered <- c(rendered,'')
    rendered <- c(rendered,line)
  }
  writeLines(rendered,file.path(out,'RESULTS_SUMMARY.md'))
}

main <- function() {
  script <- resolve_script();root <- dirname(dirname(script))
  sourced <- length(Filter(Negate(is.null),lapply(sys.frames(),function(x)x$ofile)))>0
  opt <- parse_options(if(sourced)character() else commandArgs(trailingOnly=TRUE),root)
  for(p in c('nlme','emmeans')) if(!requireNamespace(p,quietly=TRUE)) stop('Required package unavailable: ',p)
  old <- options(contrasts=c('contr.treatment','contr.poly'));on.exit(options(old),add=TRUE)
  d <- prepare_data(opt$data)
  if(!file.exists(opt$reference)) stop('Original decomposition table is required to verify the baseline.')
  dir.create(opt$output,recursive=TRUE,showWarnings=FALSE);out <- normalizePath(opt$output)
  traits <- c('foreleg','midleg','hindleg','ear_lin','eye');forms <- model_formulas()
  allfits <- list();alldata <- list();contrasts <- list();comparison <- list();checks <- list()
  samples <- list();means <- list();fit_tests <- list();k <- 0L
  for(trait in traits) {
    message('Analysing ',trait)
    included <- complete.cases(d[,c(trait,'logP_c')])
    samples[[trait]] <- data.frame(ID=d$ID,trait=trait,group=d$group,included=included)
    dt <- d[included,c('ID','group','pronotum','logP_c','instar_c','male')]
    dt$y <- log(d[[trait]][included]);alldata[[trait]] <- dt
    if(any(table(dt$group)<10)) stop('Too few observations in a group.')
    for(g in levels(dt$group)) if(min(dt$pronotum[dt$group==g])>7.2 || max(dt$pronotum[dt$group==g])<7.2)
      stop('Reference pronotum length is outside a group range.')
    allfits[[trait]] <- list()
    for(name in names(forms)) {
      k <- k+1L
      fit_ml <- fit_gls(forms[[name]],dt,name,'ML')
      fit <- fit_gls(forms[[name]],dt,name,'REML')
      allfits[[trait]][[name]] <- list(ML=fit_ml,REML=fit)
      contrasts[[k]] <- contrast_summary(fit,dt,name,trait,opt$seed+k)
      ll <- logLik(fit_ml);npar <- attr(ll,'df')
      comparison[[k]] <- data.frame(trait=trait,model=name,n=nrow(dt),parameters=npar,
        logLik_ML=as.numeric(ll),AIC=AIC(fit_ml),AICc=AIC(fit_ml)+2*npar*(npar+1)/(nrow(dt)-npar-1),
        BIC=BIC(fit_ml))
      checks[[k]] <- simulation_checks(fit,dt,name,trait,opt$replicates,opt$seed+100L+k)
      grid <- reference_grid(dt);up <- grid;up$logP_c <- 1
      X <- design_matrix(fit,grid);Xup <- design_matrix(fit,up)
      means[[k]] <- data.frame(trait=trait,model=name,group=grid$group,
        adjusted_log_mean=drop(X%*%coef(fit)),adjusted_geometric_mean=exp(drop(X%*%coef(fit))),
        body_size_slope=drop((Xup-X)%*%coef(fit)))
    }
    for(j in 2:length(forms)) {
      a <- logLik(allfits[[trait]][[j-1]]$ML);b <- logLik(allfits[[trait]][[j]]$ML)
      stat <- 2*(as.numeric(b)-as.numeric(a));df <- attr(b,'df')-attr(a,'df')
      if(stat< -1e-6 || df<=0) stop('Invalid nested comparison.')
      fit_tests[[length(fit_tests)+1L]] <- data.frame(trait=trait,
        extension=paste(names(forms)[j-1],'to',names(forms)[j]),LR=max(0,stat),df=df,
        p_value=pchisq(max(0,stat),df,lower.tail=FALSE))
    }
  }
  results <- do.call(rbind,contrasts)
  results$p_holm <- ave(results$p_value,interaction(results$model,results$contrast),FUN=function(x)p.adjust(x,'holm'))
  results$holm_family <- paste('Five traits:',results$model,results$contrast)
  results$direction <- ifelse(results$estimate_log>0,'positive','negative')
  results$holm_below_05 <- results$p_holm<.05
  verify_original(results,opt$reference)
  comparison <- do.call(rbind,comparison)
  comparison$delta_AIC <- comparison$AIC-ave(comparison$AIC,comparison$trait,FUN=min)
  comparison$delta_AICc <- comparison$AICc-ave(comparison$AICc,comparison$trait,FUN=min)
  tests <- do.call(rbind,fit_tests)
  tests$p_holm <- ave(tests$p_value,tests$extension,FUN=function(x)p.adjust(x,'holm'))
  checks <- do.call(rbind,checks)
  write_summary(results,comparison,tests,checks,out)
  write_csv(results,file.path(out,'contrast_comparison.csv'))
  write_csv(comparison,file.path(out,'model_comparison_ML.csv'))
  write_csv(tests,file.path(out,'nested_model_tests_ML.csv'))
  write_csv(checks,file.path(out,'conditional_simulation_checks.csv'))
  write_csv(do.call(rbind,means),file.path(out,'adjusted_means_and_slopes.csv'))
  write_csv(do.call(rbind,samples),file.path(out,'sample_manifest.csv'))
  saveRDS(allfits,file.path(out,'fitted_models.rds'))
  saveRDS(opt,file.path(out,'run_settings.rds'))
  write_csv(data.frame(item=c('data','data_md5','script_md5','reference_md5','run_utc','replicates'),
    value=c(normalizePath(opt$data),unname(tools::md5sum(opt$data)),unname(tools::md5sum(script)),
      unname(tools::md5sum(opt$reference)),format(Sys.time(),tz='UTC',usetz=TRUE),opt$replicates)),
    file.path(out,'run_manifest.csv'))
  capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
  png(file.path(out,'contrast_comparison.png'),width=2400,height=1650,res=250)
  comparison_plot(results);invisible(dev.off())
  pdf(file.path(out,'contrast_comparison.pdf'),width=9.6,height=6.6)
  comparison_plot(results);invisible(dev.off())
  diagnostic_plots(allfits,alldata,out)
  notes <- c('# Instar-model sensitivity analysis','',
    'This is an exploratory check prompted by predictive diagnostics. Recorded morph labels and the existing manuscript outputs are preserved.',
    'The original per-trait sample sizes, estimates, p-values and residual degrees of freedom were reproduced against the archived decomposition.',
    '', '## Models',
    'M0_original: log trait ~ log pronotum + male terminal instar + sex; one residual SD.',
    'M1_group_variance: same mean model; separate residual SDs for female, eighth, ninth and tenth groups.',
    'M2_group_slopes: also permits separate body-size slopes; male means at 7.2 mm still follow equal log-scale steps.',
    'M3_free_groups: separate group means, body-size slopes and residual SDs; equivalent to independent within-group regressions.',
    '', '## Comparisons',
    'Each candidate uses the same observations within each trait. ML likelihoods/AIC/AICc and nested likelihood-ratio tests compare specifications. REML estimates and covariance matrices provide inference.',
    'For each model, average_male_step = (tenth - eighth)/2 on the log scale; matched_sex = tenth - female; step_minus_sex is their direct difference. All are evaluated at 7.2 mm pronotum length.',
    'In M0-M2 the average male step is also the fitted constant per-instar change. M3 imposes no such constant gradient; adjacent ninth-minus-eighth and tenth-minus-ninth contrasts and the ninth-minus-midpoint departure are reported separately.',
    'The sign of step_minus_sex compares two descriptive contrasts, not causal developmental and weapon contributions. Positive means the average male step is larger; negative means the matched-sex difference is larger.',
    'All ear results are for square-root-transformed area (linearised size). Exponentiating step_minus_sex yields a ratio of multiplicative effects, not a subtraction of percentage changes.',
    'Confidence intervals are pointwise 95% t intervals. M0 uses the original residual df; M1-M2 use emmeans Satterthwaite df including uncertainty in the variance parameters. For M3, the equivalent independent group regressions allow an analytic Welch-Satterthwaite calculation. Holm p-values control each five-trait contrast family separately within each model. Model-extension tests have separate five-trait Holm families.',
    'All candidates are retained; the analysis does not choose whichever model favours a biological hypothesis. This sensitivity does not re-estimate the separate joint repeated-leg model.',
    '', '## Diagnostics',
    paste('Conditional simulation checks use',opt$replicates,'replicates at observed covariates with fitted coefficients and residual SDs held fixed.'),
    'These are descriptive plug-in checks, not posterior predictive intervals, cross-validation or formal model-adequacy tests. Parameter uncertainty and multiplicity are not included. Free group intercepts reproduce group means by construction; passing those mean checks is therefore not independent validation.',
    'Inspect residual_diagnostics.pdf for residual shapes and fitted-value trends as well as conditional_simulation_checks.csv for group means/spreads.',
    '', '## Files',
    'contrast_comparison.csv contains all effects, intervals, df and Holm p-values; contrast_comparison.png/pdf displays the direct step-minus-sex contrast.',
    'model_comparison_ML.csv and nested_model_tests_ML.csv identify which assumptions affect fit. adjusted_means_and_slopes.csv records the group predictions and body-size slopes.',
    'sample_manifest.csv, run_manifest.csv, session_info.txt and fitted_models.rds retain the analysis provenance.',
    '', 'Method references: https://rvlenth.github.io/emmeans/articles/models.html ; https://stat.ethz.ch/R-manual/R-devel/library/nlme/html/anova.gls.html')
  writeLines(notes,file.path(out,'READ_ME_FIRST.md'))
  message('Verified original decomposition. Completed sensitivity analysis: ',out)
  invisible(out)
}

if(!identical(Sys.getenv('WETA_INSTAR_HELPERS_ONLY'),'1')) main()
