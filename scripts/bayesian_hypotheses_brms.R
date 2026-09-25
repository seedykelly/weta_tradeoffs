#!/usr/bin/env Rscript
# Bayesian tests of the trait-compensation hypotheses, with matched frequentist
# estimates. Standalone: deliberately NOT added to run_all.R.
#
# From the original project root:
#   Rscript scripts/bayesian_hypotheses_brms.R --mode=check
#   Rscript scripts/bayesian_hypotheses_brms.R --mode=smoke
#   Rscript scripts/bayesian_hypotheses_brms.R --mode=fit
# Optional: --models=morphology,integration,instar --prior-scales=1,2
#           --chains=4 --iter=4000 --warmup=2000 --cores=2 --seed=20260923
#           --data=/path/trait_data.csv --output=/path/new-output
#           --cv-folds=5 --slope-rope=0.1
# Values containing spaces must be shell-quoted. source() also works, with defaults.
#
# PURPOSE AND SCOPE
# * Preserve the CSV morph column exactly; do not derive or replace morph labels.
# * H1 searching: eighth > tenth for relative leg and eye dimensions. Tympanic
#   size is an additional, less direct prediction. H2 positive integration:
#   tenth > eighth in relative dimensions; positive within-group head slopes;
#   stronger head slopes in tenth-instar males are tested as separate contrasts.
# * H3 coordinated divergence: posterior uncertainty in the SAME design-weighted
#   rank-one share as the manuscript, using pooled size-residual SD scaling.
# * Report planned contrasts, anterior/posterior allocation, midpoint practical
#   equivalence, and descriptive instar/sex contrasts. These do not identify
#   causal weapon effects, functional compensation, or resource-allocation costs.
# * The five outcomes form a multivariate Gaussian response with an unstructured
#   residual covariance shared across groups. No ID random intercept is added:
#   each row already contains all five traits of one individual. This model does
#   NOT test equality of group-specific covariance matrices or exactly rank one.
# * All models use the SAME complete cases and matched mean formulas. Thus the
#   frequentist companion is a fair likelihood comparison, not a reproduction of
#   every manuscript model: the paper uses trait-specific missing-data samples,
#   common body-size slopes for selected sensory models, and repeated-leg GLS.
# * Priors regularize coefficients without forcing positive slopes or any morph
#   ordering. Morphs are fixed categories, not assumed exchangeable random groups.
# * Fit both prior scales, inspect diagnostics and predictive checks, and report
#   both. Do not select a prior or framework because it favours a hypothesis.
# * Optional specimen-level cross-validation compares prediction accuracy with
#   identical folds. It cannot decide which statistical philosophy is 'correct'.
#
# Compilation disables precompiled headers per model to avoid stale macOS SDK
# caches; the installed CmdStan toolchain is not rebuilt or changed.
#
# Outputs go to analysis_outputs/weta_brms/{check,smoke,fit}, keeping the manuscript
# outputs untouched. brms caches fits with file_refit='on_change'. Smoke estimates
# are for software checks only. No packages are installed automatically.
# Documentation: https://paulbuerkner.com/brms/articles/brms_multivariate.html
#                https://paulbuerkner.com/brms/reference/set_prior.html

resolve_script <- function() {
  sourced <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (length(sourced)) return(normalizePath(tail(sourced, 1L)[[1L]]))
  arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  if (length(arg) != 1L) stop('Run with Rscript or source().')
  path <- sub('^--file=', '', arg)
  if (!file.exists(path)) path <- gsub('~+~', ' ', path, fixed = TRUE)
  normalizePath(path, mustWork = TRUE)
}

parse_options <- function(args, root) {
  opt <- list(mode='fit', models='morphology,integration,instar',
              `prior-scales`='1,2', chains='4', iter='4000', warmup='2000',
              cores='2', seed='20260923', `cv-folds`='0', `slope-rope`='NA',
              data=file.path(root,'data','trait_data.csv'), output='', sdk='auto')
  if (length(args)) for (arg in args) {
    if (!grepl('^--[^=]+=.+$', arg)) stop('Arguments must have the form --name=value: ',arg)
    name <- sub('^--([^=]+)=.*$', '\\1', arg)
    if (!name %in% names(opt)) stop('Unknown option: ',name)
    opt[[name]] <- sub('^--[^=]+=', '', arg)
  }
  if (!opt$mode %in% c('check','smoke','fit')) stop('mode must be check, smoke or fit.')
  integer_option <- function(name, minimum) {
    x <- suppressWarnings(as.numeric(opt[[name]]))
    if (length(x)!=1L || !is.finite(x) || x!=floor(x) || x<minimum || x>.Machine$integer.max)
      stop('Invalid integer for ', name)
    as.integer(x)
  }
  for (name in c('chains','iter','warmup','cores','seed')) opt[[name]] <- integer_option(name,1)
  opt$`cv-folds` <- integer_option('cv-folds',0)
  if (opt$`cv-folds`==1L) stop('cv-folds must be zero or at least two.')
  if (opt$warmup>=opt$iter) stop('warmup must be smaller than iter.')
  opt$models <- strsplit(opt$models,',',fixed=TRUE)[[1L]]
  if (anyDuplicated(opt$models) || !all(opt$models %in% c('morphology','integration','instar')))
    stop('Invalid or duplicated model names.')
  opt$`prior-scales` <- suppressWarnings(as.numeric(strsplit(opt$`prior-scales`,',',fixed=TRUE)[[1L]]))
  if (!length(opt$`prior-scales`) || any(!is.finite(opt$`prior-scales`)) ||
      any(opt$`prior-scales`<=0) || anyDuplicated(opt$`prior-scales`)) stop('Invalid prior-scales.')
  if (opt$`slope-rope`!='NA' && is.na(suppressWarnings(as.numeric(opt$`slope-rope`)))) stop('Invalid slope-rope.')
  opt$`slope-rope` <- suppressWarnings(as.numeric(opt$`slope-rope`))
  if (length(opt$`slope-rope`)!=1L || (!is.na(opt$`slope-rope`) &&
      (!is.finite(opt$`slope-rope`) || opt$`slope-rope`<=0))) stop('Invalid slope-rope.')
  if (opt$mode=='smoke') {
    opt$chains <- 2L; opt$iter <- 600L; opt$warmup <- 300L
    if (!any(grepl('^--prior-scales=',args))) opt$`prior-scales` <- 1
    opt$`cv-folds` <- 0L
  }
  if (!nzchar(opt$output)) opt$output <- file.path(root,'analysis_outputs','weta_brms',opt$mode)
  if (opt$mode=='fit' && (opt$chains<4 || opt$iter-opt$warmup<1000))
    warning('A final analysis normally needs >=4 chains and >=1000 retained iterations per chain.')
  opt
}

prepare_data <- function(path) {
  d <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
  required <- c('ID','sex','morph','pronotum','head_length','head_width',
                'forefemur','foretibia','midfemur','midtibia','hindfemur','hindtibia','ear','eye')
  if (!all(required %in% names(d))) stop('Missing columns: ',paste(setdiff(required,names(d)),collapse=', '))
  if (anyNA(d$ID) || any(!nzchar(trimws(d$ID))) || any(d$ID!=trimws(d$ID)) || anyDuplicated(d$ID))
    stop('ID values must be unique, nonblank and without surrounding spaces.')
  groups <- c('female','eighth','ninth','tenth')
  if (anyNA(d$morph) || !all(d$morph %in% groups) || anyNA(d$sex) || !all(d$sex %in% c('f','m')))
    stop('Unexpected or missing morph/sex labels. Correct the source metadata, not this analysis.')
  if (any((d$sex=='f')!=(d$morph=='female'))) stop('Inconsistent sex and morph labels.')
  for (name in setdiff(required,c('ID','sex','morph'))) {
    if (!is.numeric(d[[name]])) stop('Measurement must be numeric: ',name)
    if (any(!is.na(d[[name]]) & (!is.finite(d[[name]]) | d[[name]]<=0)))
      stop('Nonpositive or nonfinite measurement: ',name)
  }
  included <- complete.cases(d[,required])
  manifest <- data.frame(ID=d$ID,sex=d$sex,morph=d$morph,included=included)
  d <- d[included,required]
  d$group <- factor(d$morph,levels=groups)
  if (any(table(d$group)<5)) stop('Too few complete specimens in a group.')
  d$logP <- log(d$pronotum/7.2)
  d$logH <- log(sqrt(d$head_length*d$head_width)/15)
  d$instar <- unname(c(female=10,eighth=8,ninth=9,tenth=10)[d$morph])-10
  d$male <- as.integer(d$sex=='m')
  # Fixed reference lengths (mm) centre intercept priors; no SD estimated from
  # these data is used to tune the coefficient priors. Ear is a linear dimension.
  refs <- c(fore=25,mid=25,hind=40,earlin=1.75,eye=3)
  original <- cbind(fore=d$forefemur+d$foretibia,mid=d$midfemur+d$midtibia,
                    hind=d$hindfemur+d$hindtibia,earlin=sqrt(d$ear),eye=d$eye)
  for (name in names(refs)) d[[name]] <- log(original[,name]/refs[[name]])
  list(data=d,manifest=manifest,refs=refs)
}

model_rhs <- function(name) switch(name,
  morphology='0 + group + group:logP',
  integration='0 + group + group:logP + group:logH',
  instar='1 + logP + instar + male', stop('Unknown model'))

make_spec <- function(name, data, scale=1) {
  responses <- c('fore','mid','hind','earlin','eye')
  rhs <- model_rhs(name)
  parts <- lapply(responses,function(r) brms::bf(as.formula(paste(r,'~',rhs)),family=gaussian()))
  formula <- Reduce(`+`,parts)+brms::set_rescor(TRUE)
  available <- brms::get_prior(formula,data=data)
  prior <- brms::set_prior('lkj(2)',class='rescor')
  for (r in responses) {
    prior <- c(prior,brms::set_prior('student_t(3, 0, 0.15 * prior_scale)',class='sigma',resp=r))
    if (name=='instar') prior <- c(prior,brms::set_prior('normal(0, 0.5 * prior_scale)',class='Intercept',resp=r))
    coefs <- available$coef[available$class=='b' & available$resp==r & nzchar(available$coef)]
    for (coef in unique(coefs)) {
      # Log-log slopes: N(0,1); per-instar and sex shifts: N(0,.25);
      # reference-size group means of centred log traits: N(0,.5).
      sd <- if (grepl('logP|logH',coef)) 1 else if (coef %in% c('instar','male')) .25 else .5
      prior <- c(prior,brms::set_prior(sprintf('normal(0, %s * prior_scale)',sd),
                                      class='b',coef=coef,resp=r))
    }
  }
  list(formula=formula,prior=prior,
       stanvars=brms::stanvar(x=scale,name='prior_scale'),rhs=rhs)
}

reference_grid <- function(d) {
  data.frame(group=factor(levels(d$group),levels=levels(d$group)),
             logP=0,logH=0,instar=c(0,-2,-1,0),male=c(0,1,1,1))
}

frequentist_models <- function(d,name) {
  setNames(lapply(c('fore','mid','hind','earlin','eye'),function(r)
    lm(as.formula(paste(r,'~',model_rhs(name))),data=d)),c('fore','mid','hind','earlin','eye'))
}

linear_frequentist <- function(model,grid,weights) {
  X <- model.matrix(delete.response(terms(model)),grid,contrasts.arg=model$contrasts,xlev=model$xlevels)
  a <- drop(crossprod(weights,X))
  estimate <- sum(a*coef(model));se <- sqrt(drop(t(a)%*%vcov(model)%*%a))
  df <- df.residual(model)
  c(estimate=estimate,se=se,lower=estimate-qt(.975,df)*se,
    upper=estimate+qt(.975,df)*se,p=2*pt(abs(estimate/se),df,lower.tail=FALSE),df=df)
}

summarize_draws <- function(x, exponent=NA_real_, bounds=c(NA_real_,NA_real_)) {
  q <- unname(quantile(x,c(.025,.5,.975)))
  transformed <- if (is.na(exponent)) q else 100*expm1(exponent*q)
  data.frame(median=q[2],lower_95=q[1],upper_95=q[3],
             effect_median=transformed[2],effect_lower_95=transformed[1],effect_upper_95=transformed[3],
             probability_positive=mean(x>0),probability_negative=mean(x<0),
             probability_in_region=if (anyNA(bounds)) NA_real_ else mean(x>bounds[1]&x<bounds[2]))
}

expected_draws <- function(fit,grid) {
  responses <- c('fore','mid','hind','earlin','eye')
  arrays <- lapply(responses,function(r) brms::posterior_epred(fit,newdata=grid,resp=r,re_formula=NA))
  if (!all(vapply(arrays,is.matrix,logical(1)))) stop('Unexpected prediction dimensions.')
  out <- array(NA_real_,c(nrow(arrays[[1]]),nrow(grid),length(responses)),
               dimnames=list(NULL,NULL,responses))
  for (j in seq_along(responses)) out[,,j] <- arrays[[j]]
  out
}

# The Cholesky factor supplies the SAME design weighting used by the existing
# rank-one analysis. The trait scaling is held fixed across posterior draws.
rank_geometry <- function(d) {
  Y <- as.matrix(d[,c('fore','mid','hind','earlin','eye')])
  scale <- apply(qr.resid(qr(model.matrix(~logP,d)),Y),2,sd)
  Z <- cbind(intercept=1,sweep(model.matrix(~0+group,d),1,d$logP,'*'))
  W <- model.matrix(~group,d)[,-1,drop=FALSE]
  R <- chol(crossprod(qr.resid(qr(Z),W)))
  list(R=R,scale=scale)
}
rank_fraction <- function(B,geometry) {
  singular <- svd(geometry$R %*% sweep(B,2,geometry$scale,'/'),nu=0,nv=0)$d
  if (sum(singular^2)<=0) return(NA_real_)
  singular[1]^2/sum(singular^2)
}

write_csv <- function(x,path) write.csv(x,path,row.names=FALSE,na='')

sampler_diagnostics <- function(fit,label,out) {
  draws <- posterior::as_draws_array(fit)
  draws <- posterior::subset_draws(draws,variable=grep('^(b_|sigma_|rescor_|Intercept_)',posterior::variables(draws),value=TRUE))
  stats <- posterior::summarise_draws(draws)
  write_csv(as.data.frame(stats),file.path(out,paste0(label,'_parameter_diagnostics.csv')))
  np <- brms::nuts_params(fit)
  divergences <- sum(np$Value[np$Parameter=='divergent__'])
  depth <- sum(np$Value[np$Parameter=='treedepth__']>=12)
  energy <- split(np$Value[np$Parameter=='energy__'],np$Chain[np$Parameter=='energy__'])
  ebfmi <- vapply(energy,function(e) mean(diff(e)^2)/var(e),numeric(1))
  good <- all(is.finite(stats$rhat)) && all(stats$rhat<=1.01) &&
    all(stats$ess_bulk>=400) && all(stats$ess_tail>=400) &&
    divergences==0 && depth==0 && all(is.finite(ebfmi)&ebfmi>=.3)
  data.frame(fit=label,max_rhat=max(stats$rhat),min_bulk_ess=min(stats$ess_bulk),
             min_tail_ess=min(stats$ess_tail),divergences=divergences,
             max_treedepth_hits=depth,min_ebfmi=min(ebfmi),diagnostics_pass=good)
}

predictive_checks <- function(fit,d,label,out) {
  # Group-wise location and spread checks catch problems hidden by pooled plots.
  checks <- list();index <- 0L
  pdf(file.path(out,paste0(label,'_predictive_checks.pdf')),width=9,height=8)
  on.exit(dev.off(),add=TRUE)
  for (r in c('fore','mid','hind','earlin','eye')) {
    yrep <- brms::posterior_predict(fit,resp=r,ndraws=200)
    par(mfrow=c(2,2),mar=c(4,4,3,1))
    for (g in levels(d$group)) {
      ii <- which(d$group==g);observed <- d[[r]][ii]
      means <- rowMeans(yrep[,ii,drop=FALSE])
      sds <- apply(yrep[,ii,drop=FALSE],1,sd)
      limits <- range(c(observed,quantile(yrep[,ii],c(.001,.999))))
      plot(density(observed),xlim=limits,main=paste(label,r,g),xlab='Centred log trait',lwd=2)
      for (k in seq_len(min(30,nrow(yrep)))) lines(density(yrep[k,ii]),col=adjustcolor('steelblue',.2))
      lines(density(observed),lwd=2)
      index <- index+1L
      checks[[index]] <- data.frame(fit=label,trait=r,group=g,observed_mean=mean(observed),
        replicated_mean_lower=quantile(means,.025),replicated_mean_upper=quantile(means,.975),
        observed_sd=sd(observed),replicated_sd_lower=quantile(sds,.025),replicated_sd_upper=quantile(sds,.975))
    }
  }
  write_csv(do.call(rbind,checks),file.path(out,paste0(label,'_predictive_checks.csv')))
}

summarize_hypotheses <- function(fit,d,name,scale,out,slope_rope) {
  label <- paste0(name,'_prior',scale)
  grid <- reference_grid(d);mu <- expected_draws(fit,grid)
  freq <- frequentist_models(d,name)
  traits <- dimnames(mu)[[3]];rows <- list();i <- 0L
  events <- list()
  if (name=='morphology') {
    contrasts <- list(eighth_minus_tenth=c(0,1,0,-1),
                      ninth_minus_midpoint=c(0,-.5,1,-.5),tenth_minus_female=c(-1,0,0,1))
    for (r in traits) for (cn in names(contrasts)) {
      w <- contrasts[[cn]];x <- drop(mu[,,r]%*%w);ff <- linear_frequentist(freq[[r]],grid,w)
      exponent <- if (r=='earlin') 2 else 1
      bounds <- if(cn=='ninth_minus_midpoint') log(c(.975,1.025))/exponent else c(NA_real_,NA_real_)
      i <- i+1L
      rows[[i]] <- cbind(data.frame(model=name,prior_scale=scale,trait=r,contrast=cn,
        effect_unit=if(r=='earlin') 'percent area' else 'percent length',
        region=if(cn=='ninth_minus_midpoint') 'within +/-2.5% (illustrative)' else ''),
        summarize_draws(x,exponent,bounds),data.frame(frequentist_estimate=ff['estimate'],
        frequentist_lower_95=ff['lower'],frequentist_upper_95=ff['upper'],
        frequentist_p=ff['p'],frequentist_df=ff['df']))
    }
    result <- do.call(rbind,rows)
    result$frequentist_holm_p <- ave(result$frequentist_p,result$trait,FUN=function(x)p.adjust(x,'holm'))
    write_csv(result,file.path(out,paste0(label,'_morph_contrasts.csv')))
    delta <- mu[,2,]-mu[,4,]
    events$searching_legs_and_eye_all_eighth_greater <- mean(apply(delta[,c('fore','mid','hind','eye')]>0,1,all))
    events$all_five_tenth_greater_than_eighth <- mean(apply(delta<0,1,all))
    events$all_five_increase_eighth_ninth_tenth <- mean(apply((mu[,3,]-mu[,2,]>0)&(mu[,4,]-mu[,3,]>0),1,all))
    group_w <- list(tenth_minus_eighth=c(0,-1,0,1),ninth_minus_midpoint=c(0,-.5,1,-.5),tenth_minus_female=c(-1,0,0,1))
    leg_w <- list(fore_relative_hind=c(1,0,-1),mid_relative_hind=c(0,1,-1),anterior_mean_relative_hind=c(.5,.5,-1))
    allocation <- list();k <- 0L
    # Frequentist covariance across response coefficients uses the multivariate
    # residual covariance and the common design matrix (not independent legs).
    X <- model.matrix(freq$fore);residual <- sapply(freq,residuals)
    Sigma <- crossprod(residual)/df.residual(freq$fore)
    invXX <- solve(crossprod(X));Xg <- model.matrix(delete.response(terms(freq$fore)),grid)
    for (gn in names(group_w)) for (ln in names(leg_w)) {
      gw <- group_w[[gn]];lw <- leg_w[[ln]]
      x <- drop((mu[,,1]*lw[1]+mu[,,2]*lw[2]+mu[,,3]*lw[3])%*%gw)
      a <- drop(crossprod(gw,Xg));estimate <- sum(lw*vapply(freq[1:3],function(m)sum(a*coef(m)),numeric(1)))
      se <- sqrt(drop(t(a)%*%invXX%*%a)*drop(t(lw)%*%Sigma[1:3,1:3]%*%lw));df <- df.residual(freq$fore)
      k <- k+1L
      allocation[[k]] <- cbind(data.frame(group_contrast=gn,leg_contrast=ln),summarize_draws(x,1),
        data.frame(frequentist_estimate=estimate,frequentist_lower_95=estimate-qt(.975,df)*se,
                   frequentist_upper_95=estimate+qt(.975,df)*se,frequentist_p=2*pt(abs(estimate/se),df,lower.tail=FALSE)))
    }
    allocation <- do.call(rbind,allocation);allocation$frequentist_holm_p <- p.adjust(allocation$frequentist_p,'holm')
    write_csv(allocation,file.path(out,paste0(label,'_allocation_contrasts.csv')))
    geometry <- rank_geometry(d)
    fractions <- apply(mu,1,function(m) rank_fraction(m[2:4,,drop=FALSE]-matrix(m[1,],3,5,byrow=TRUE),geometry))
    fitted_means <- sapply(freq,predict,newdata=grid)
    plugin <- rank_fraction(fitted_means[2:4,]-matrix(fitted_means[1,],3,5,byrow=TRUE),geometry)
    rank <- cbind(summarize_draws(fractions),data.frame(frequentist_point_estimate=plugin,
      probability_share_above_80_percent=mean(fractions>.8),probability_share_above_90_percent=mean(fractions>.9)))
    write_csv(rank,file.path(out,paste0(label,'_dominant_axis.csv')))
    write_csv(data.frame(draw=seq_along(fractions),rank_one_fraction=fractions),file.path(out,paste0(label,'_dominant_axis_draws.csv')))
    # Store response axes explicitly, avoiding claims that the 80/90% descriptive
    # thresholds are preregistered hypothesis tests or probabilities of exact rank.
  }
  if (name=='integration') {
    plus <- grid;plus$logH <- 1
    slopes <- expected_draws(fit,plus)-mu
    for (r in traits) for (g in seq_len(4)) {
      weights <- rep(0,8);weights[g] <- -1;weights[g+4] <- 1
      ff <- linear_frequentist(freq[[r]],rbind(grid,plus),weights)
      bounds <- if(is.na(slope_rope)) c(NA_real_,NA_real_) else c(-slope_rope,slope_rope)
      i <- i+1L;rows[[i]] <- cbind(data.frame(trait=r,group=levels(d$group)[g],
        slope_region_halfwidth=slope_rope),summarize_draws(slopes[,g,r],bounds=bounds),
        data.frame(frequentist_estimate=ff['estimate'],frequentist_lower_95=ff['lower'],
                   frequentist_upper_95=ff['upper'],frequentist_p=ff['p']))
    }
    write_csv(do.call(rbind,rows),file.path(out,paste0(label,'_head_trait_slopes.csv')))
    differences <- list();k <- 0L
    for (r in traits) for (g in c(2L,3L)) {
      w <- rep(0,8);w[c(4,g,8,g+4)] <- c(-1,1,1,-1)
      ff <- linear_frequentist(freq[[r]],rbind(grid,plus),w)
      k <- k+1L;differences[[k]] <- cbind(data.frame(trait=r,contrast=paste0('tenth_minus_',levels(d$group)[g])),
        summarize_draws(slopes[,4,r]-slopes[,g,r]),data.frame(frequentist_estimate=ff['estimate'],
          frequentist_lower_95=ff['lower'],frequentist_upper_95=ff['upper'],frequentist_p=ff['p']))
    }
    write_csv(do.call(rbind,differences),file.path(out,paste0(label,'_head_slope_differences.csv')))
    events$all_nine_male_head_leg_slopes_positive <- mean(apply(slopes[,2:4,1:3]>0,1,all))
    events$all_three_tenth_head_leg_slopes_positive <- mean(apply(slopes[,4,1:3]>0,1,all))
  }
  if (name=='instar') {
    # Use a single reference row and finite differences in this linear model.
    base <- grid[1,,drop=FALSE];instar <- base;instar$instar <- 1
    male <- base;male$male <- 1
    grid3 <- rbind(base,instar,male);ep <- expected_draws(fit,grid3)
    for (r in traits) for (cn in c('per_male_instar','matched_sex','instar_minus_sex')) {
      w <- switch(cn,per_male_instar=c(-1,1,0),matched_sex=c(-1,0,1),instar_minus_sex=c(0,1,-1))
      ff <- linear_frequentist(freq[[r]],grid3,w)
      i <- i+1L;rows[[i]] <- cbind(data.frame(trait=r,contrast=cn,
        effect_unit='percent linear dimension (ear is square root of area)'),
        summarize_draws(drop(ep[,,r]%*%w),1),data.frame(frequentist_estimate=ff['estimate'],
          frequentist_lower_95=ff['lower'],frequentist_upper_95=ff['upper'],frequentist_p=ff['p']))
    }
    result <- do.call(rbind,rows)
    result$frequentist_holm_p <- ave(result$frequentist_p,result$contrast,FUN=function(x)p.adjust(x,'holm'))
    write_csv(result,file.path(out,paste0(label,'_instar_sex.csv')))
  }
  if (length(events)) write_csv(data.frame(event=names(events),posterior_probability=unlist(events)),
                                file.path(out,paste0(label,'_joint_hypothesis_probabilities.csv')))
}

configure_macos_sdk <- function(requested='auto') {
  old <- Sys.getenv('SDKROOT',unset=NA_character_)
  if (Sys.info()[['sysname']]!='Darwin' || requested=='default') return(old)
  sdk <- requested
  if (requested=='auto') {
    compiler <- suppressWarnings(system2('/usr/bin/xcrun',c('--find','clang++'),stdout=TRUE,stderr=FALSE))
    if (length(compiler)!=1L) stop('Cannot locate the selected macOS compiler; use --sdk=default or an explicit SDK path.')
    if (grepl('/Contents/Developer/Toolchains/',compiler,fixed=TRUE)) {
      developer <- sub('/Toolchains/.*$','',compiler)
      sdk <- file.path(developer,'Platforms','MacOSX.platform','Developer','SDKs','MacOSX.sdk')
    } else {
      sdk <- suppressWarnings(system2('/usr/bin/xcrun','--show-sdk-path',stdout=TRUE,stderr=FALSE))
    }
  }
  if(length(sdk)!=1L || !dir.exists(sdk)) stop('Cannot locate a matching SDK; provide --sdk=/path/to/MacOSX.sdk.')
  # Process-local only. Keep compiler and SDK from the same developer toolchain.
  Sys.setenv(SDKROOT=normalizePath(sdk))
  message('Using macOS SDK: ',Sys.getenv('SDKROOT'))
  old
}

fit_model <- function(spec,d,settings,file,seed,sample_prior='no',compiled=NULL) {
  request <- list(chains=settings$chains,iter=settings$iter,warmup=settings$warmup,
                  seed=seed,sample_prior=sample_prior)
  cached_file <- paste0(file,'.rds')
  cached <- if(file.exists(cached_file)) readRDS(cached_file) else NULL
  refit <- if(!is.null(cached) && identical(cached$weta_sampling_request,request)) 'on_change' else 'always'
  common <- list(prior=spec$prior,stanvars=spec$stanvars,
    chains=settings$chains,cores=min(settings$cores,settings$chains),
    iter=settings$iter,warmup=settings$warmup,seed=seed,
    sample_prior=sample_prior,control=list(adapt_delta=.99,max_treedepth=12),
    file=file,file_refit=refit,refresh=500,
    stan_model_args=list(cpp_options=list(PRECOMPILED_HEADERS='false')))
  if (is.null(compiled)) {
    result <- do.call(brms::brm,c(list(formula=spec$formula,data=d,backend='cmdstanr',
                                    save_pars=brms::save_pars(all=TRUE)),common))
  } else {
    # brm(fit=...) ignores changed priors/data. update() is essential here to
    # regenerate Stan data for wider priors and for sample_prior='only'.
    # The prior expressions are identical across these fits; only the Stan data
    # variable prior_scale changes. Inherit the validated expressions, avoiding
    # brms 2.21 duplicating the transformed LKJ prior when merging specifications.
    common$prior <- NULL
    result <- do.call(stats::update,c(list(object=compiled,newdata=d),common))
  }
  used <- brms::standata(result)
  requested <- brms::make_standata(spec$formula,data=d,prior=spec$prior,
                                  stanvars=spec$stanvars,sample_prior=sample_prior)
  stopifnot(identical(used$prior_scale,requested$prior_scale),
            identical(used$prior_only,requested$prior_only),nrow(result$data)==nrow(d),
            isTRUE(all.equal(as.data.frame(result$data),d[,names(result$data),drop=FALSE],
                             check.attributes=FALSE)))
  result$weta_sampling_request <- request
  result$weta_specimen_ids <- d$ID
  saveRDS(result,cached_file)
  result
}

cross_validate <- function(d,settings,out) {
  K <- settings$`cv-folds`
  if (!K) return(invisible(NULL))
  if (K>min(table(d$group))) stop('cv-folds exceeds the smallest group size.')
  set.seed(settings$seed+7000L);fold <- integer(nrow(d))
  for (g in levels(d$group)) {
    ii <- which(d$group==g);fold[ii] <- sample(rep(seq_len(K),length.out=length(ii)))
  }
  write_csv(data.frame(ID=d$ID,group=d$group,fold=fold),file.path(out,'specimen_cv_folds.csv'))
  results <- list();diagnostics <- list();index <- 0L
  for (name in settings$models) {
   compiled <- NULL
   for (k in seq_len(K)) {
    train <- d[fold!=k,];test <- d[fold==k,]
    spec <- make_spec(name,train,1)
    fit <- fit_model(spec,train,settings,file.path(out,paste0('cv_',name,'_',k)),settings$seed+7000L+k,compiled=compiled)
    compiled <- fit
    diagnostics[[length(diagnostics)+1L]] <- sampler_diagnostics(fit,paste0('cv_',name,'_',k),out)
    ep <- expected_draws(fit,test);freq <- frequentist_models(train,name)
    for (r in dimnames(ep)[[3]]) {
      index <- index+1L
      results[[index]] <- data.frame(model=name,fold=k,ID=test$ID,trait=r,observed=test[[r]],
        bayesian=colMeans(ep[,,r]),frequentist=predict(freq[[r]],newdata=test))
    }
  }
  }
  results <- do.call(rbind,results)
  results$bayesian_squared_error <- (results$observed-results$bayesian)^2
  results$frequentist_squared_error <- (results$observed-results$frequentist)^2
  write_csv(results,file.path(out,'specimen_cv_predictions.csv'))
  grouped <- split(results,interaction(results$model,results$trait,drop=TRUE))
  summary <- do.call(rbind,lapply(grouped,function(x) data.frame(model=x$model[1],trait=x$trait[1],
    n=nrow(x),bayesian_rmse=sqrt(mean(x$bayesian_squared_error)),
    frequentist_rmse=sqrt(mean(x$frequentist_squared_error)),
    mean_paired_squared_error_difference=mean(x$bayesian_squared_error-x$frequentist_squared_error))))
  write_csv(summary,file.path(out,'specimen_cv_rmse.csv'))
  write_csv(do.call(rbind,diagnostics),file.path(out,'specimen_cv_diagnostics.csv'))
  # Positive paired differences favour frequentist prediction, negative Bayesian.
  # These are descriptive scores, not significance tests or evidence of mechanism.
}

comparison_outputs <- function(settings,out) {
  files <- list.files(out,pattern='_(morph_contrasts|head_trait_slopes|head_slope_differences|instar_sex|allocation_contrasts)[.]csv$',full.names=TRUE)
  # Use only requested fits, so stale outputs from an earlier run cannot enter
  # the comparison after the user narrows --models or --prior-scales.
  prefixes <- unlist(lapply(settings$models,function(m)paste0(m,'_prior',settings$`prior-scales`,'_')))
  files <- files[vapply(basename(files),function(f)any(startsWith(f,prefixes)),logical(1))]
  if(!length(files)) return(invisible(NULL))
  pdf(file.path(out,'bayesian_frequentist_comparison.pdf'),width=10,height=8)
  on.exit(dev.off(),add=TRUE)
  all_rows <- list()
  for (file in files) {
    d <- read.csv(file,check.names=FALSE)
    label_columns <- intersect(c('trait','group','contrast','group_contrast','leg_contrast'),names(d))
    labels <- apply(d[,label_columns,drop=FALSE],1,paste,collapse=': ')
    # The native coefficient scale is common to both methods. Percent-transformed
    # Bayesian effects are in separate columns of the detailed tables.
    limits <- range(c(d$lower_95,d$upper_95,d$frequentist_lower_95,d$frequentist_upper_95))
    y <- rev(seq_len(nrow(d)))
    par(mar=c(4,14,4,1))
    plot(d$median,y+.13,xlim=limits,ylim=c(.4,nrow(d)+.6),yaxt='n',ylab='',
         xlab='Log contrast or log-log slope (95% intervals)',pch=16,col='#0072B2',
         main=basename(file),cex=.8)
    abline(v=0,col='grey70',lty=3)
    segments(d$lower_95,y+.13,d$upper_95,y+.13,col='#0072B2',lwd=1.5)
    points(d$frequentist_estimate,y-.13,pch=17,col='#D55E00',cex=.8)
    segments(d$frequentist_lower_95,y-.13,d$frequentist_upper_95,y-.13,col='#D55E00',lwd=1.5)
    axis(2,at=y,labels=labels,las=1,cex.axis=.6)
    legend('topright',legend=c('Bayesian median / credible interval','Frequentist estimate / pointwise CI'),
           pch=c(16,17),col=c('#0072B2','#D55E00'),bty='n',cex=.65)
    all_rows[[length(all_rows)+1L]] <- data.frame(source=basename(file),quantity=labels,
      bayesian_median=d$median,bayesian_lower_95=d$lower_95,bayesian_upper_95=d$upper_95,
      frequentist_estimate=d$frequentist_estimate,frequentist_lower_95=d$frequentist_lower_95,
      frequentist_upper_95=d$frequentist_upper_95,
      posterior_probability_positive=d$probability_positive,frequentist_p=d$frequentist_p)
  }
  write_csv(do.call(rbind,all_rows),file.path(out,'bayesian_frequentist_comparison.csv'))
}

main <- function() {
  script <- resolve_script();root <- dirname(dirname(script))
  # When sourced, ignore unrelated command-line arguments of the parent session.
  sourced <- length(Filter(Negate(is.null),lapply(sys.frames(),function(x)x$ofile)))>0
  opt <- parse_options(if(sourced) character() else commandArgs(trailingOnly=TRUE),root)
  needed <- c('brms','cmdstanr','posterior')
  absent <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
  if(length(absent)) stop('Install required R packages first: ',paste(absent,collapse=', '))
  prepared <- prepare_data(opt$data);d <- prepared$data
  dir.create(opt$output,recursive=TRUE,showWarnings=FALSE)
  out <- normalizePath(opt$output)
  write_csv(prepared$manifest,file.path(out,'sample_manifest.csv'))
  write_csv(data.frame(group=levels(d$group),n=as.integer(table(d$group))),file.path(out,'sample_counts.csv'))
  write_csv(data.frame(trait=names(prepared$refs),reference_length_mm=prepared$refs),file.path(out,'response_reference_lengths.csv'))
  capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
  saveRDS(opt,file.path(out,'run_settings.rds'))
  write_csv(data.frame(item=c('input','input_md5','script_md5','run_utc','mode'),value=c(
    normalizePath(opt$data),unname(tools::md5sum(opt$data)),unname(tools::md5sum(script)),
    format(Sys.time(),tz='UTC',usetz=TRUE),opt$mode)),file.path(out,'run_manifest.csv'))
  notes <- c('# Bayesian hypothesis analysis', '',
    paste('Run mode:',opt$mode),paste('Complete specimens:',nrow(d)),
    'Morph assignments are read directly from the supplied CSV; no reclassification is performed.', '',
    'Models: morphology = group-specific body-size slopes for all five traits; integration adds group-specific raw-head slopes; instar = body size + male terminal instar + sex.',
    'Each model has a shared unstructured residual covariance. Individual slopes are fixed effects with symmetric priors, not a hierarchical model pooling male morphs.',
    'Matched frequentist models use identical complete cases and mean formulas. Their estimates can differ from the published trait-specific samples and selected common-slope sensory models.',
    'Bayesian 95% intervals are marginal equal-tailed credible intervals. Frequentist 95% intervals are pointwise t intervals; named Holm families adjust p-values only. Do not equate their coverage or posterior tail areas.',
    'Positive/negative probabilities are conditional on the likelihood and priors, not frequentist p-values. Joint probabilities are computed draw by draw and do not assume traits are independent.',
    'The +/-2.5% midpoint region is illustrative, matching the earlier sensitivity analysis; justify a biologically meaningful margin before making a substantive equivalence claim.',
    'Head-slope practical equivalence is omitted unless --slope-rope is supplied. That optional halfwidth applies to raw log-log slopes, NOT the earlier standardized +/-0.20 region.',
    'Dominant-axis shares use the manuscript design weighting and pooled size-residual SDs. Probabilities above 80% and 90% are descriptive summaries; neither is a test of exactly rank one.',
    'Group-specific covariance equality, detailed segment models and measurement-error models are outside this script.',
    'Prior scale 2 doubles all coefficient and residual-SD prior scales; LKJ(2) is held constant. Both fits must be reported rather than choosing a favourable prior.',
    'Prior predictive draws are simulated without conditioning on outcomes; inspect them alongside posterior predictive checks of group means and spreads.',
    'All final fits require satisfactory Rhat, bulk/tail ESS, divergences, tree depth and energy diagnostics. Summary tables are not automatically interpretable when diagnostics fail.',
    'Smoke-mode tables are software checks and must not be reported as scientific results.',
    'Optional cross-validation holds out entire specimens, stratifies folds by morph, and compares mean predictions on the centred log scale. It uses prior scale 1 only; predictive accuracy does not establish a causal hypothesis.',
    'No automated declaration of a winning framework or supported hypothesis is made. Review effect sizes, uncertainty, prior sensitivity, diagnostics and substantive relevance together.')
  writeLines(notes,file.path(out,'READ_ME_FIRST.md'))
  # Fast validation produces inspectable Stan code/data and matched-model estimates.
  for (name in opt$models) {
    spec <- make_spec(name,d,1)
    writeLines(brms::make_stancode(spec$formula,data=d,prior=spec$prior,stanvars=spec$stanvars),file.path(out,paste0(name,'.stan')))
    saveRDS(brms::make_standata(spec$formula,data=d,prior=spec$prior,stanvars=spec$stanvars),file.path(out,paste0(name,'_stan_data.rds')))
    write_csv(as.data.frame(spec$prior),file.path(out,paste0(name,'_priors.csv')))
    saveRDS(frequentist_models(d,name),file.path(out,paste0(name,'_matched_frequentist.rds')))
  }
  if(opt$mode=='check') {message('Input, priors and Stan generation validated. No sampling performed.');return(invisible(out))}
  previous_sdk <- configure_macos_sdk(opt$sdk)
  on.exit(if(is.na(previous_sdk)) Sys.unsetenv('SDKROOT') else Sys.setenv(SDKROOT=previous_sdk),add=TRUE)
  cmdstanr::check_cmdstan_toolchain(quiet=TRUE)
  diagnostics <- list();k <- 0L
  for (name in opt$models) {
    compiled <- NULL
    for (scale in opt$`prior-scales`) {
      label <- paste0(name,'_prior',scale);message('Fitting ',label)
      spec <- make_spec(name,d,scale)
      fit <- fit_model(spec,d,opt,file.path(out,label),opt$seed+k,compiled=compiled)
      compiled <- fit
      k <- k+1L;diagnostics[[k]] <- sampler_diagnostics(fit,label,out)
      summarize_hypotheses(fit,d,name,scale,out,opt$`slope-rope`)
      set.seed(opt$seed+2000L+k)
      predictive_checks(fit,d,label,out)
      # Prior-only fits reuse compiled Stan code. Use the same priors/scale and
      # complete-case covariates; outcomes do not enter their likelihood.
      prior_settings <- opt;prior_settings$iter <- 1000L;prior_settings$warmup <- 500L
      prior_fit <- fit_model(spec,d,prior_settings,file.path(out,paste0(label,'_prior_only')),
                             opt$seed+1000L+k,sample_prior='only',compiled=fit)
      k <- k+1L;diagnostics[[k]] <- sampler_diagnostics(prior_fit,paste0(label,'_prior_only'),out)
      set.seed(opt$seed+2000L+k)
      predictive_checks(prior_fit,d,paste0(label,'_prior_only'),out)
      write_csv(do.call(rbind,diagnostics),file.path(out,'sampler_diagnostics.csv'))
      rm(prior_fit);gc(verbose=FALSE)
    }
  }
  comparison_outputs(opt,out)
  cross_validate(d,opt,out)
  diagnostics <- do.call(rbind,diagnostics)
  if(any(!diagnostics$diagnostics_pass)) warning('Some diagnostics failed. Inspect sampler_diagnostics.csv before interpreting results.')
  writeLines(c(notes,'',paste('Sampling finished:',format(Sys.time(),tz='UTC',usetz=TRUE)),
               paste('All diagnostic screens passed:',all(diagnostics$diagnostics_pass))),file.path(out,'READ_ME_FIRST.md'))
  message('Finished. Results: ',out)
  invisible(out)
}

# Tests may source pure helpers without triggering model fitting.
if (!identical(Sys.getenv('WETA_BRMS_HELPERS_ONLY'),'1')) main()
