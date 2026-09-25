#!/usr/bin/env Rscript
# Classification sensitivity and a display of the design-weighted group effect.
# Reuses the rank-one functions from the maintained multivariate analysis.
file_arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)
script_file <- sub('^--file=', '', file_arg[[1L]])
if (!file.exists(script_file)) script_file <- gsub('~+~', ' ', script_file, fixed = TRUE)
script_dir <- dirname(normalizePath(script_file, mustWork = TRUE))
project_root <- dirname(script_dir)
source(file.path(script_dir, 'workflow_helpers.R'))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3L) stop('Expected data, output directory and permutations.')
data_file <- if (length(args) >= 1L) args[[1L]] else file.path(project_root, 'data', 'trait_data.csv')
output_dir <- if (length(args) >= 2L) args[[2L]] else file.path(project_root, 'analysis_outputs', 'weta_trait_analysis')
permutations <- parse_replicates(if (length(args) >= 3L) args[[3L]] else 9999L, 'Residual randomizations', 999L)
tables_dir <- file.path(output_dir, 'tables')
figures_dir <- file.path(output_dir, 'figures')
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Load only pure function definitions; do not rerun the full upstream script.
wanted <- c('matrix_fit', 'attributes_from_vectors', 'rank_one_fit', 'rank_one_randomization')
expressions <- as.list(parse(file.path(script_dir, 'female_reference_multivariate_reanalysis.R')))
definitions <- Filter(function(x) is.call(x) && identical(x[[1]], as.name('<-')) &&
  is.symbol(x[[2]]) && as.character(x[[2]]) %in% wanted, expressions)
stopifnot(length(definitions) == length(wanted))
for (x in definitions) eval(x)

d <- read.csv(data_file)
levels_group <- c('female', 'eighth', 'ninth', 'tenth')
d$group <- factor(d$morph, levels = levels_group)
d$logP_c <- log(d$pronotum / 7.2)
d$foreleg <- d$forefemur + d$foretibia
d$midleg <- d$midfemur + d$midtibia
d$hindleg <- d$hindfemur + d$hindtibia
d$ear_lin <- sqrt(d$ear)
traits <- c('foreleg', 'midleg', 'hindleg', 'ear_lin', 'eye')
legacy <- ifelse(d$head_length < 19.04, 'eighth', ifelse(d$head_length < 24.20, 'ninth', 'tenth'))
d$classification_sensitive <- d$sex == 'm' & as.character(d$group) != legacy
complete <- complete.cases(d[, c(traits, 'logP_c', 'group', 'head_length')])
d <- d[complete, ]

prepare <- function(data) {
  Y <- log(as.matrix(data[, traits]))
  scale_sd <- apply(qr.resid(qr(model.matrix(~ logP_c, data)), Y), 2, sd)
  Y <- sweep(Y, 2, scale_sd, '/')
  Z <- cbind(intercept = 1, sweep(model.matrix(~ 0 + group, data), 1, data$logP_c, '*'))
  W <- model.matrix(~ group, data)[, -1, drop = FALSE]
  colnames(W) <- c('D8', 'D9', 'D10')
  list(Y = Y, Z = Z, W = W)
}
full <- prepare(d)
observed <- rank_one_fit(full$Y, full$Z, full$W)
archived <- read.csv(file.path(tables_dir, 'female_reference_rank_one_tests.csv'))
base <- archived[archived$scenario == 'standardized_head_excluded', ]
stopifnot(nrow(base) == 1L, abs(observed$rank_one_fraction - base$rank_one_fraction) < 1e-10)
filtered <- prepare(d[!d$classification_sensitive, ])
sensitivity <- rank_one_randomization(filtered$Y, filtered$Z, filtered$W, permutations, 5301L)
result <- data.frame(
  scenario = c('Full sample', 'Exclude disagreements with 2010 categories'),
  n = c(nrow(d), sum(!d$classification_sensitive)),
  excluded = c(0L, sum(d$classification_sensitive)),
  rank_one_fraction = c(observed$rank_one_fraction, sensitivity$observed$rank_one_fraction),
  proportion_beyond_rank_one = c(observed$proportion_beyond_rank_one, sensitivity$observed$proportion_beyond_rank_one),
  null_95 = c(base$null_95, unname(sensitivity$null_quantiles['95%'])),
  permutations = c(base$permutations, permutations),
  p_value = c(base$p_value, sensitivity$p_value)
)
write.csv(result, file.path(tables_dir, 'group_divergence_classification_sensitivity.csv'), row.names = FALSE)

# Right singular axes of the same design-weighted matrix used by rank_one_fit.
residual_W <- qr.resid(qr(full$Z), full$W)
R <- chol(crossprod(residual_W))
Q <- residual_W %*% solve(R)
A <- crossprod(Q, qr.resid(qr(full$Z), full$Y))
decomposition <- svd(A)
B <- observed$B_unrestricted
coordinates <- B %*% decomposition$v[, 1:2, drop = FALSE]
if (coordinates['D10', 1] < 0) coordinates[, 1] <- -coordinates[, 1]
if (coordinates['D8', 2] < 0) coordinates[, 2] <- -coordinates[, 2]
shares <- 100 * decomposition$d^2 / sum(decomposition$d^2)
stopifnot(abs(shares[1] / 100 - base$rank_one_fraction) < 1e-10)
write.csv(data.frame(vector = rownames(B), axis_1 = coordinates[,1], axis_2 = coordinates[,2],
  axis_1_weighted_percent = shares[1], axis_2_weighted_percent = shares[2]),
  file.path(tables_dir, 'group_divergence_weighted_axes.csv'), row.names = FALSE)
colours <- c('#0072B2', '#009E73', '#D55E00')
labels <- c('Eighth instar', 'Ninth instar', 'Tenth instar')
png(file.path(figures_dir, 'group_divergence.png'), width = 2700, height = 1350, res = 300, bg = 'white')
par(mfrow = c(1,2), mar = c(5.1,4.8,2.8,1), oma = c(3,0,0,0), family = 'serif', las = 1)
plot(coordinates, type = 'n', xlim = range(c(0,coordinates[,1])) + c(-0.6,0.8),
  ylim = range(c(0,coordinates[,2])) + c(-0.65,0.65), bty = 'l',
  xlab = sprintf('Axis 1 (%.1f%%)', shares[1]),
  ylab = sprintf('Axis 2 (%.1f%%)', shares[2]), main = 'A  Adjusted group differences', cex.main = 1)
abline(h = 0,v = 0,lty = 3,col = 'grey75')
for (i in 1:3) {
  arrows(0,0,coordinates[i,1],coordinates[i,2],length = .09,col = colours[i],lwd = 2)
  points(coordinates[i,1],coordinates[i,2],pch = c(16,17,15)[i],col = colours[i],cex = 1.2)
  text(coordinates[i,1],coordinates[i,2], labels = c('Eighth','Ninth','Tenth')[i],pos = 3,offset = .6,cex = .85)
}
points(0,0,pch = 21,bg = 'white',cex = 1.1)
text(0.25,0,'Female reference',pos = 1,offset = .7,cex = .8)
matplot(t(B),type = 'b',pch = c(16,17,15),lty = 1,col = colours,lwd = 2,
  xaxt = 'n',bty = 'l',xlab = '',ylab = 'Male–female difference (size-adjusted SD)',
  ylim = range(c(0,B)) + c(-.3,.5),main = 'B  Contributions of individual traits',cex.main = 1)
axis(1,at = 1:5,labels = c('Foreleg','Midleg','Hindleg','Tympanic\nsize','Eye'),cex.axis = .8)
abline(h = 0,lty = 3,col = 'grey65')
par(fig = c(0,1,0,1),new = TRUE,mar = c(0,0,0,0),oma = c(0,0,0,0))
plot.new();legend('bottom',legend = labels,col = colours,pch = c(16,17,15),lty = 1,
  horiz = TRUE,bty = 'n',cex = 1,inset = c(0,.015))
invisible(dev.off())
print(result)
