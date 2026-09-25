#!/usr/bin/env Rscript
# Verify reported estimates, multiplicity adjustments and precision inputs.
script_arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)
script_file <- sub('^--file=', '', script_arg[[1L]])
if (!file.exists(script_file)) script_file <- gsub('~+~', ' ', script_file, fixed = TRUE)
root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))
tables <- file.path(root, 'analysis_outputs', 'weta_trait_analysis', 'tables')
read_table <- function(name) read.csv(file.path(tables, paste0(name, '.csv')), check.names = FALSE)
close <- function(x,y,tol=1e-8) {
  stopifnot(length(x)>0, length(x)==length(y), all(is.finite(x)), all(is.finite(y)), max(abs(x-y))<tol)
}
checks <- character()
record <- function(x) {checks <<- c(checks,x);message('PASS: ',x)}
files <- c('joint_leg_males_vs_female','joint_leg_strategy_contrasts',
  'joint_anterior_posterior_difference_contrasts','head_leg_raw_strategy_slope_contrasts',
  'ear_body_strategy_contrasts','eye_body_strategy_contrasts',
  'ear_body_all_pairwise_group_contrasts','eye_body_all_pairwise_group_contrasts')
for (name in files) {
  d <- read_table(name)
  families <- if ('leg' %in% names(d)) split(d,d$leg) else list(d)
  for (f in families) {
    raw <- 2*pt(abs(f$t.ratio),f$df,lower.tail=FALSE)
    close(f$p.value,p.adjust(raw,'holm'))
    critical <- qt(1-0.05/(2*nrow(f)),f$df)
    close(f$lower.CL,f$estimate-critical*f$SE)
    close(f$upper.CL,f$estimate+critical*f$SE)
    if ('ratio' %in% names(f)) close(f$ratio,exp(f$estimate))
    if ('area_ratio' %in% names(f)) close(f$area_ratio,exp(2*f$estimate))
  }
}
record('Eight reported contrast tables reproduce Holm tests, Bonferroni intervals and response-scale ratios.')

precision <- read_table('precision_planned_leg_contrasts')
strategy <- read_table('joint_leg_strategy_contrasts')
reverse <- strategy$contrast=='Eighth instar - Tenth instar'
strategy$estimate[reverse] <- -strategy$estimate[reverse]
strategy$contrast[reverse] <- 'Tenth instar - Eighth instar'
strategy$contrast <- sub('instar,Tenth','instar, Tenth',strategy$contrast,fixed=TRUE)
allocation <- read_table('joint_anterior_posterior_difference_contrasts')
for (i in seq_len(nrow(precision))) {
  row <- precision[i,]
  if (row$contrast_type=='planned') {
    ref <- strategy[strategy$contrast==row$contrast & strategy$leg==row$leg,]
  } else ref <- allocation[allocation$contrast==row$contrast,]
  stopifnot(nrow(ref)==1L)
  close(c(row$estimate,row$SE,row$df),c(ref$estimate,ref$SE,ref$df))
}
close(precision$ci_lower,exp(precision$estimate-qt(.975,precision$df)*precision$SE))
close(precision$ci_upper,exp(precision$estimate+qt(.975,precision$df)*precision$SE))
close(precision$min_detectable_percent,100*(exp((qt(.975,precision$df)+qt(.8,precision$df))*precision$SE)-1))
record('All 12 leg-precision rows use the main contrast estimates, SEs and Satterthwaite degrees of freedom.')

slopes <- read_table('head_leg_raw_slopes')
precision_slopes <- read_table('precision_weapon_leg_slopes')
ref <- slopes[match(paste(precision_slopes$group,precision_slopes$leg),paste(slopes$group,slopes$leg)),]
close(precision_slopes$weapon_leg_slope,ref$log_head_size.trend)
close(precision_slopes$df,ref$df)
close(precision_slopes$ci_lower,ref$lower.CL)
close(precision_slopes$ci_upper,ref$upper.CL)
stopifnot(all(slopes$log_head_size.trend>0),sum(slopes$lower.CL>0)==11)
record('All 12 head-leg estimates are positive; 11 pointwise intervals exclude zero. Precision intervals match the slope table.')

# Independently reconstruct the instar/sex regressions from the analysis input.
d <- read.csv(file.path(root,'data','trait_data.csv'))
d$sex <- factor(ifelse(d$sex=='f','female','male'),levels=c('female','male'))
d$instar_c <- c(female=10,eighth=8,ninth=9,tenth=10)[d$morph]-10
d$logP_c <- log(d$pronotum/7.2)
d$foreleg <- d$forefemur+d$foretibia;d$midleg <- d$midfemur+d$midtibia
d$hindleg <- d$hindfemur+d$hindtibia;d$ear_lin <- sqrt(d$ear)
reported <- read_table('instar_matched_decomposition')
for (i in seq_len(nrow(reported))) {
 row <- reported[i,]
 fit <- lm(log(d[[row$trait]])~logP_c+instar_c+sex,data=d)
 co <- coef(summary(fit));V <- vcov(fit)
 close(row$instar_component_log_per_instar,co['instar_c','Estimate'])
 close(row$weapon_component_log,co['sexmale','Estimate'])
 close(row$instar_component_p,co['instar_c','Pr(>|t|)'])
 close(row$weapon_component_p,co['sexmale','Pr(>|t|)'])
 se_diff <- sqrt(V['instar_c','instar_c']+V['sexmale','sexmale']-2*V['instar_c','sexmale'])
 close(row$component_difference_SE,se_diff)
 difference <- co['instar_c','Estimate']-co['sexmale','Estimate']
 close(row$component_difference_log,difference)
 close(row$component_difference_p,2*pt(abs(difference/se_diff),df.residual(fit),lower.tail=FALSE))
 close(c(row$n,row$df_residual),c(nobs(fit),df.residual(fit)))
 for (term in c('instar_c','sexmale')) {
   prefix <- if (term=='instar_c') 'instar_component' else 'weapon_component'
   expected_ci <- co[term,'Estimate']+c(-1,1)*qt(.975,df.residual(fit))*co[term,'Std. Error']
   close(c(row[[paste0(prefix,'_CI_lower')]],row[[paste0(prefix,'_CI_upper')]]),expected_ci)
 }
 close(row$instar_component_pct_per_instar,100*expm1(co['instar_c','Estimate']))
 close(row$weapon_component_pct_male_vs_female,100*expm1(co['sexmale','Estimate']))
}
for (prefix in c('instar_component','weapon_component','component_difference'))
 close(reported[[paste0(prefix,'_holm_p')]],p.adjust(reported[[paste0(prefix,'_p')]],'holm'))
record('Original-model sensitivity: five instar/sex regressions reproduce coefficients, intervals, sample sizes, direct coefficient comparisons, percentage changes and three Holm families.')

# Verify the final flexible specification independently using GLS and explicit
# contrast vectors. Production uses separate group regressions instead.
final <- read_table('instar_matched_contrasts_long')
parameters <- read_table('instar_group_model_parameters')
samples <- read_table('instar_contrast_samples')
d$group <- factor(d$morph,levels=c('female','eighth','ninth','tenth'))
W <- rbind(male_step=c(0,-.5,0,.5),matched_sex=c(-1,0,0,1),contrast_difference=c(1,-.5,0,-.5))
for (trait in unique(final$trait)) {
 included <- complete.cases(d[,c(trait,'logP_c')]);dt <- d[included,];dt$y <- log(dt[[trait]])
 fit <- nlme::gls(y~group*logP_c,data=dt,weights=nlme::varIdent(form=~1|group),method='REML',
   control=nlme::glsControl(maxIter=200,msMaxIter=200,tolerance=1e-8,msTol=1e-8))
 grid <- data.frame(group=factor(levels(dt$group),levels=levels(dt$group)),logP_c=0)
 X <- model.matrix(~group*logP_c,grid)
 fitted_means <- drop(X%*%coef(fit));V <- X%*%vcov(fit)%*%t(X)
 close(fitted_means,parameters$adjusted_log_mean[parameters$trait==trait])
 sm <- samples[samples$trait==trait,]
 stopifnot(identical(sm$ID,d$ID),identical(as.logical(sm$included),included))
 group_labels <- c('Female','Eighth instar','Ninth instar','Tenth instar')
 stopifnot(identical(sm$group,group_labels[as.integer(d$group)]))
 # Closed-form regression prediction variance at x=0, independently of lm.
 pv <- numeric(4);rdf <- numeric(4)
 for (g in seq_along(levels(dt$group))) {
   z <- dt[as.integer(dt$group)==g,];x <- z$logP_c;y <- z$y
   slope <- sum((x-mean(x))*(y-mean(y)))/sum((x-mean(x))^2)
   resid <- y-(mean(y)-slope*mean(x)+slope*x)
   rdf[g] <- length(x)-2
   pv[g] <- sum(resid^2)/rdf[g]*(1/length(x)+mean(x)^2/sum((x-mean(x))^2))
 }
 for (j in seq_len(nrow(W))) {
   z <- final[final$trait==trait & final$contrast==rownames(W)[j],]
   stopifnot(nrow(z)==1L,z$n==nrow(dt))
   w <- W[j,];estimate <- sum(w*fitted_means);se <- sqrt(drop(w%*%V%*%w))
   close(z$estimate_log,estimate);close(z$SE,se,tol=1e-6)
   parts <- w^2*pv;df <- sum(parts)^2/sum(parts^2/rdf)
   close(z$df,df)
   close(z$p,2*pt(abs(z$estimate_log/z$SE),df,lower.tail=FALSE))
   close(c(z$CI_lower,z$CI_upper),z$estimate_log+c(-1,1)*qt(.975,df)*z$SE)
   close(z$percent,100*expm1(z$estimate_log))
 }
}
for (contrast in rownames(W)) {
 z <- final[final$contrast==contrast,]
 close(z$holm_p,p.adjust(z$p,'holm'))
}
wide <- read_table('instar_matched_contrasts')
for (contrast in rownames(W)) {
 z <- final[final$contrast==contrast,];z <- z[match(wide$trait,z$trait),]
 for (column in c('estimate_log','SE','df','CI_lower','CI_upper','percent','p','holm_p'))
   close(wide[[paste0(contrast,'_',column)]],z[[column]])
}
hind <- final[final$trait=='hindleg' & final$contrast=='contrast_difference',]
stopifnot(hind$CI_lower<0,hind$CI_upper>0,hind$holm_p>.05)
record('Final flexible contrasts independently reproduce GLS estimates, analytic Welch-Satterthwaite df, intervals, Holm tests and recorded sample assignments; wide and long tables agree.')

sensitivity <- read_table('group_divergence_classification_sensitivity')
stopifnot(identical(as.integer(sensitivity$n),c(311L,301L)),sensitivity$excluded[2]==10L,
  sensitivity$permutations[2]==9999L,sensitivity$p_value[2]<.001)
close(sensitivity$rank_one_fraction[1],read_table('female_reference_rank_one_tests')$rank_one_fraction[1])
axes <- read_table('group_divergence_weighted_axes')
close(axes$axis_1_weighted_percent,rep(100*sensitivity$rank_one_fraction[1],nrow(axes)))
record('Classification sensitivity and the new figure use the same design-weighted group-effect definition as the main analysis.')

# Cross-analysis audit labels and p-values.
audit <- read_table('focal_multiplicity_sensitivity')
close(audit$p_Holm,p.adjust(audit$p_value,'holm'));close(audit$p_BH,p.adjust(audit$p_value,'BH'))
stopifnot(audit$tier[audit$analysis=='Five-trait correlation-matrix equality']=='Exploratory')
record('The ten-test audit reproduces Holm and Benjamini–Hochberg adjustments; the matrix test is labelled exploratory.')

report <- c('# Statistical reporting verification', '',
 'The checks below passed against the maintained analysis input and result tables.', '',
 paste0('- ',checks), '',
 sprintf('Excluding ten category-disagreement males changed the dominant-axis share from %.1f%% to %.1f%%; the residual rank-one test remained p = %.4f (%s randomizations).',
 100*sensitivity$rank_one_fraction[1],100*sensitivity$rank_one_fraction[2],sensitivity$p_value[2],sensitivity$permutations[2]),'',
 'Holm-adjusted p-values and Bonferroni-adjusted simultaneous intervals are explicitly distinguished in the documents. Precision intervals and minimum detectable effects refer to single comparisons without multiplicity adjustment.', '',
 'Leg-precision estimates now reuse the main contrast/slope tables rather than substituting residual observation degrees of freedom. Primary-model estimates and Holm p-values are retained. The final secondary instar comparison allows group-specific means, slopes and residual variances and treats the hindleg contrast ordering as uncertain.', '',
 'Morph categories are the externally derived categories supplied with the dataset; no mixture model was fitted anew. The exclusion analysis is a sensitivity to historical category definitions, not a reclassification of the main sample.', '',
 'Scope: targeted verification of reported contrasts, transformations, precision calculations, original and final instar comparisons and the classification sensitivity. This is not an independent rerun of every model or a verification of raw measurements.')
writeLines(report,file.path(root,'analysis_outputs','weta_trait_analysis','statistical_reporting_audit.md'))
