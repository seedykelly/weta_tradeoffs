#!/usr/bin/env Rscript
arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
path <- sub('^--file=','',arg)
if(!file.exists(path)) path <- gsub('~+~',' ',path,fixed=TRUE)
root <- dirname(dirname(normalizePath(path)))
Sys.setenv(WETA_BRMS_HELPERS_ONLY='1')
source(file.path(root,'scripts','bayesian_hypotheses_brms.R'))
Sys.unsetenv('WETA_BRMS_HELPERS_ONLY')
input <- commandArgs(trailingOnly=TRUE)
data_file <- if(length(input)) input[[1]] else file.path(root,'data','trait_data.csv')
d <- prepare_data(data_file)$data
stopifnot(nrow(d)==311L,identical(as.character(d$group),d$morph))
close <- function(x,y,tol=1e-9) stopifnot(length(x)==length(y),all(is.finite(x)),all(is.finite(y)),max(abs(x-y))<tol)
# The comparison uses the original morph assignments, even if an alternative
# threshold would put a specimen in a different group.
raw <- read.csv(data_file);close(sum(d$morph=='tenth'),sum(raw$morph=='tenth'))
stopifnot(d$morph[d$ID=='GDWH']=='tenth')
# Bayesian rank summary must match the published definition, not ordinary PCA.
freq <- frequentist_models(d,'morphology');grid <- reference_grid(d)
m <- sapply(freq,predict,newdata=grid);B <- m[2:4,]-matrix(m[1,],3,5,byrow=TRUE)
geometry <- rank_geometry(d)
close(rank_fraction(B,geometry),.88842213451758)
rank_one <- outer(c(1,2,3),c(1,2,-1,.5,3));close(rank_fraction(rank_one,geometry),1)
close(rank_fraction(2*B,geometry),rank_fraction(B,geometry))
# Contrast algebra: exact directions, geometric midpoint, and log-area units.
w <- c(0,1,0,-1)
f <- linear_frequentist(freq$fore,grid,w)
close(f['estimate'],m[2,'fore']-m[4,'fore'])
f <- linear_frequentist(freq$fore,grid,c(0,-.5,1,-.5))
close(f['estimate'],m[3,'fore']-.5*(m[2,'fore']+m[4,'fore']))
s <- summarize_draws(rep(log(1.1),50),2)
close(s$effect_median,21);close(s$probability_positive,1)
s <- summarize_draws(c(-.02,0,.02,.04),bounds=c(-.025,.025));close(s$probability_in_region,.75)
# A perfect synthetic model makes the coefficient-difference contrast known.
synthetic <- d;synthetic$fore <- .3+1.2*d$logP+.07*d$instar+.15*d$male+sin(seq_len(nrow(d)))*.001
fit <- frequentist_models(synthetic,'instar')$fore
base <- grid[1,,drop=FALSE];instar <- base;instar$instar <- 1;male <- base;male$male <- 1
contrast <- linear_frequentist(fit,rbind(base,instar,male),c(0,1,-1))
close(contrast['estimate'],coef(fit)['instar']-coef(fit)['male'])
# Input validation rejects malformed arguments and duplicated IDs.
expect_error <- function(expr) stopifnot(inherits(tryCatch({force(expr);NULL},error=identity),'error'))
expect_error(parse_options('--chains=1.5',root));expect_error(parse_options('--cv-folds=1',root))
expect_error(parse_options('--prior-scales=0',root));expect_error(parse_options('--warmup=5000',root))
expect_error(parse_options('--slope-rope=-0.1',root))
smoke <- parse_options(c('--mode=smoke','--prior-scales=1,2'),root)
close(smoke$`prior-scales`,c(1,2));close(smoke$chains,2)
tmp <- tempfile(fileext='.csv');write.csv(rbind(raw,raw[1,]),tmp,row.names=FALSE)
expect_error(prepare_data(tmp));unlink(tmp)
# No likelihood change from centring head length by a fixed reference.
fhead <- frequentist_models(d,'integration')$fore
rawhead <- lm(fore~0+group+group:logP+group:I(logH+log(15)),d)
close(fitted(fhead),fitted(rawhead),1e-8)
cat('Bayesian helper checks passed: labels, weighted rank, contrasts, area transformation, equivalence region and validation.\n')
