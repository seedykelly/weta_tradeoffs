# ============================================================
# WELLINGTON TREE WETA MORPHOLOGICAL ALLOCATION
# Females and eighth-, ninth-, and tenth-instar males
# ============================================================

source_dir <- "/Users/uqam/Documents/Research_Admin/research projects/Trait compensation/data/clean"

# 1. PACKAGES ------------------------------------------------

library(tidyverse)
library(emmeans)
library(broom)
library(patchwork)

# Optional package for model diagnostics
# install.packages("performance")
library(performance)


# 2. IMPORT DATA ---------------------------------------------

dat_raw <- read_csv(
  "data/clean/trait_data.csv",
  na = c("", "NA"),
  show_col_types = FALSE
)

glimpse(dat_raw)


# 3. BASIC DATA CHECKS ---------------------------------------

# Dimensions
dim(dat_raw)

# Number of unique individuals
n_distinct(dat_raw$ID)

# Duplicate IDs
dat_raw %>%
  count(ID) %>%
  filter(n > 1)

# Group sample sizes
dat_raw %>%
  count(sex, morph)

# Missing values by variable
dat_raw %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "n_missing"
  )

# Inspect ranges of all measured traits
dat_raw %>%
  select(
    pronotum,
    head_length,
    head_width,
    forefemur,
    foretibia,
    midfemur,
    midtibia,
    hindfemur,
    hindtibia,
    ear,
    eye
  ) %>%
  summarise(
    across(
      everything(),
      list(
        minimum = ~ min(.x, na.rm = TRUE),
        maximum = ~ max(.x, na.rm = TRUE)
      )
    )
  )


# 4. CREATE ANALYSIS VARIABLES -------------------------------

# Reference pronotum length for adjusted comparisons.
# This lies within the observed body-size range of all groups.
reference_pronotum <- 7.2

dat <- dat_raw %>%
  mutate(
    # Four biological groups
    group = factor(
      morph,
      levels = c("female", "eighth", "ninth", "tenth")
    ),
    
    # Ordered version for descriptive purposes only
    male_morph_ordered = case_when(
      group == "eighth" ~ 8,
      group == "ninth"  ~ 9,
      group == "tenth"  ~ 10,
      TRUE              ~ NA_real_
    ),
    
    # Total femur + tibia length for each leg pair
    foreleg = forefemur + foretibia,
    midleg  = midfemur + midtibia,
    hindleg = hindfemur + hindtibia,
    
    # Composite linear head dimension
    head_size = sqrt(head_length * head_width),
    
    # Ear is an area, so square root places it on a linear scale
    ear_linear = sqrt(ear),
    
    # Centred log body size
    # A value of zero represents a pronotum length of 7.2 mm
    logP_c = log(pronotum) - log(reference_pronotum),
    
    # Log-transformed traits
    log_foreleg       = log(foreleg),
    log_midleg        = log(midleg),
    log_hindleg       = log(hindleg),
    
    log_forefemur     = log(forefemur),
    log_foretibia     = log(foretibia),
    log_midfemur      = log(midfemur),
    log_midtibia      = log(midtibia),
    log_hindfemur     = log(hindfemur),
    log_hindtibia     = log(hindtibia),
    
    log_head_size     = log(head_size),
    log_head_length   = log(head_length),
    log_head_width    = log(head_width),
    
    log_eye           = log(eye),
    log_ear_linear    = log(ear_linear)
  )


# Verify that all values to be logged are positive
dat %>%
  summarise(
    across(
      c(
        pronotum,
        foreleg,
        midleg,
        hindleg,
        head_size,
        ear_linear,
        eye
      ),
      ~ sum(.x <= 0, na.rm = TRUE)
    )
  )


# 5. DESCRIPTIVE STATISTICS ----------------------------------

group_summary <- dat %>%
  group_by(group) %>%
  summarise(
    n = n(),
    
    pronotum_mean = mean(pronotum, na.rm = TRUE),
    pronotum_sd   = sd(pronotum, na.rm = TRUE),
    pronotum_min  = min(pronotum, na.rm = TRUE),
    pronotum_max  = max(pronotum, na.rm = TRUE),
    
    foreleg_mean = mean(foreleg, na.rm = TRUE),
    midleg_mean  = mean(midleg, na.rm = TRUE),
    hindleg_mean = mean(hindleg, na.rm = TRUE),
    
    ear_mean  = mean(ear, na.rm = TRUE),
    eye_mean  = mean(eye, na.rm = TRUE),
    
    .groups = "drop"
  )

group_summary


# Determine common body-size overlap among all four groups
body_size_ranges <- dat %>%
  group_by(group) %>%
  summarise(
    minimum = min(pronotum, na.rm = TRUE),
    maximum = max(pronotum, na.rm = TRUE),
    .groups = "drop"
  )

body_size_ranges

common_lower <- max(body_size_ranges$minimum)
common_upper <- min(body_size_ranges$maximum)

c(
  common_lower = common_lower,
  common_upper = common_upper,
  reference_pronotum = reference_pronotum
)


# 6. BODY-SIZE DISTRIBUTIONS ---------------------------------

p_body <- ggplot(
  dat,
  aes(x = pronotum, fill = group)
) +
  geom_density(alpha = 0.35) +
  geom_vline(
    xintercept = reference_pronotum,
    linetype = 2
  ) +
  labs(
    x = "Pronotum length",
    y = "Density",
    fill = "Group"
  ) +
  theme_classic()

p_body


# 7. LEG DATA IN LONG FORMAT ---------------------------------

legs_long <- dat %>%
  select(
    ID,
    group,
    pronotum,
    logP_c,
    foreleg,
    midleg,
    hindleg
  ) %>%
  pivot_longer(
    cols = c(foreleg, midleg, hindleg),
    names_to = "leg",
    values_to = "leg_length"
  ) %>%
  mutate(
    leg = factor(
      leg,
      levels = c("foreleg", "midleg", "hindleg"),
      labels = c("Foreleg", "Midleg", "Hindleg")
    ),
    log_leg_length = log(leg_length),
    log_pronotum = log(pronotum)
  )


# 8. ALLOMETRIC PLOT -----------------------------------------

p_leg_allometry <- ggplot(
  legs_long,
  aes(
    x = log_pronotum,
    y = log_leg_length,
    colour = group
  )
) +
  geom_point(alpha = 0.50) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE
  ) +
  facet_wrap(
    ~ leg,
    scales = "free_y"
  ) +
  labs(
    x = "Log pronotum length",
    y = "Log leg length",
    colour = "Group"
  ) +
  theme_classic()

p_leg_allometry

# 9. FORELEG ALLOMETRY ---------------------------------------

# Different slopes among groups
mod_fore_interaction <- lm(
  log_foreleg ~ logP_c * group,
  data = dat
)

# Common slope among groups
mod_fore_common <- lm(
  log_foreleg ~ logP_c + group,
  data = dat
)

# Partial F-test for whether group-specific slopes improve fit
anova(
  mod_fore_common,
  mod_fore_interaction
)

# Full interaction-model coefficients
summary(mod_fore_interaction)


# 10. MIDLEG ALLOMETRY ---------------------------------------

mod_mid_interaction <- lm(
  log_midleg ~ logP_c * group,
  data = dat
)

mod_mid_common <- lm(
  log_midleg ~ logP_c + group,
  data = dat
)

anova(
  mod_mid_common,
  mod_mid_interaction
)

summary(mod_mid_interaction)


# 11. HINDLEG ALLOMETRY --------------------------------------

mod_hind_interaction <- lm(
  log_hindleg ~ logP_c * group,
  data = dat
)

mod_hind_common <- lm(
  log_hindleg ~ logP_c + group,
  data = dat
)

anova(
  mod_hind_common,
  mod_hind_interaction
)

summary(mod_hind_interaction)

# 12. ESTIMATED ALLOMETRIC SLOPES ----------------------------

fore_slopes <- emtrends(
  mod_fore_interaction,
  ~ group,
  var = "logP_c"
)

mid_slopes <- emtrends(
  mod_mid_interaction,
  ~ group,
  var = "logP_c"
)

hind_slopes <- emtrends(
  mod_hind_interaction,
  ~ group,
  var = "logP_c"
)

fore_slopes
mid_slopes
hind_slopes

# Pairwise slope comparisons
pairs(fore_slopes, adjust = "holm")
pairs(mid_slopes, adjust = "holm")
pairs(hind_slopes, adjust = "holm")

# 13. ADJUSTED GROUP MEANS -----------------------------------

emm_fore <- emmeans(
  mod_fore_common,
  ~ group,
  at = list(logP_c = 0)
)

emm_mid <- emmeans(
  mod_mid_common,
  ~ group,
  at = list(logP_c = 0)
)

emm_hind <- emmeans(
  mod_hind_common,
  ~ group,
  at = list(logP_c = 0)
)

emm_fore
emm_mid
emm_hind


# Pairwise group comparisons
contrast_fore <- pairs(
  emm_fore,
  adjust = "holm"
)

contrast_mid <- pairs(
  emm_mid,
  adjust = "holm"
)

contrast_hind <- pairs(
  emm_hind,
  adjust = "holm"
)

contrast_fore
contrast_mid
contrast_hind

# 14. BACK-TRANSFORM ADJUSTED MEANS --------------------------

emm_fore_original <- as.data.frame(emm_fore) %>%
  mutate(
    predicted_length = exp(emmean),
    lower_CL = exp(lower.CL),
    upper_CL = exp(upper.CL),
    leg = "Foreleg"
  )

emm_mid_original <- as.data.frame(emm_mid) %>%
  mutate(
    predicted_length = exp(emmean),
    lower_CL = exp(lower.CL),
    upper_CL = exp(upper.CL),
    leg = "Midleg"
  )

emm_hind_original <- as.data.frame(emm_hind) %>%
  mutate(
    predicted_length = exp(emmean),
    lower_CL = exp(lower.CL),
    upper_CL = exp(upper.CL),
    leg = "Hindleg"
  )

adjusted_leg_means <- bind_rows(
  emm_fore_original,
  emm_mid_original,
  emm_hind_original
) %>%
  mutate(
    leg = factor(
      leg,
      levels = c("Foreleg", "Midleg", "Hindleg")
    )
  )

adjusted_leg_means

# 15. EXPRESS CONTRASTS AS RATIOS AND PERCENT DIFFERENCES ----

ratio_table <- function(contrast_object, leg_name) {
  
  as.data.frame(contrast_object) %>%
    mutate(
      leg = leg_name,
      
      # contrast = log(mean1) - log(mean2)
      ratio = exp(estimate),
      
      percent_difference = 100 * (ratio - 1),
      
      ratio_lower = exp(estimate - 1.96 * SE),
      ratio_upper = exp(estimate + 1.96 * SE)
    ) %>%
    select(
      leg,
      contrast,
      estimate,
      SE,
      df,
      t.ratio,
      p.value,
      ratio,
      ratio_lower,
      ratio_upper,
      percent_difference
    )
}

leg_contrasts <- bind_rows(
  ratio_table(contrast_fore, "Foreleg"),
  ratio_table(contrast_mid,  "Midleg"),
  ratio_table(contrast_hind, "Hindleg")
)

leg_contrasts

# 16. PLOT ADJUSTED MEANS ------------------------------------

p_adjusted_legs <- ggplot(
  adjusted_leg_means,
  aes(
    x = group,
    y = predicted_length
  )
) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(
      ymin = lower_CL,
      ymax = upper_CL
    ),
    width = 0.12
  ) +
  facet_wrap(
    ~ leg,
    scales = "free_y"
  ) +
  labs(
    x = NULL,
    y = paste0(
      "Predicted leg length at pronotum = ",
      reference_pronotum
    )
  ) +
  theme_classic()

p_adjusted_legs

# 17. REPEATED-TRAIT MODEL -----------------------------------

library(lme4)
library(lmerTest)

mod_leg_allocation <- lmer(
  log_leg_length ~ logP_c * group * leg +
    (1 | ID),
  data = legs_long,
  REML = FALSE
)

anova(
  mod_leg_allocation,
  type = 3,
  ddf = "Satterthwaite"
)

summary(mod_leg_allocation)

# Adjusted group means within each leg pair
emm_allocation <- emmeans(
  mod_leg_allocation,
  ~ group | leg,
  at = list(logP_c = 0)
)

emm_allocation

# Compare groups separately for each leg
pairs(
  emm_allocation,
  by = "leg",
  adjust = "holm"
)

# 18. SEGMENT-LEVEL DATA -------------------------------------

segments_long <- dat %>%
  select(
    ID,
    group,
    logP_c,
    forefemur,
    foretibia,
    midfemur,
    midtibia,
    hindfemur,
    hindtibia
  ) %>%
  pivot_longer(
    cols = c(
      forefemur,
      foretibia,
      midfemur,
      midtibia,
      hindfemur,
      hindtibia
    ),
    names_to = "trait",
    values_to = "segment_length"
  ) %>%
  mutate(
    leg = case_when(
      str_starts(trait, "fore") ~ "Foreleg",
      str_starts(trait, "mid")  ~ "Midleg",
      str_starts(trait, "hind") ~ "Hindleg"
    ),
    
    segment = case_when(
      str_ends(trait, "femur") ~ "Femur",
      str_ends(trait, "tibia") ~ "Tibia"
    ),
    
    leg = factor(
      leg,
      levels = c("Foreleg", "Midleg", "Hindleg")
    ),
    
    segment = factor(
      segment,
      levels = c("Femur", "Tibia")
    ),
    
    log_segment_length = log(segment_length)
  )


# Plot segment allometries
p_segments <- ggplot(
  segments_long,
  aes(
    x = logP_c,
    y = log_segment_length,
    colour = group
  )
) +
  geom_point(alpha = 0.45) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE
  ) +
  facet_grid(
    segment ~ leg,
    scales = "free_y"
  ) +
  labs(
    x = paste0(
      "Log pronotum length, centred at ",
      reference_pronotum,
      " mm"
    ),
    y = "Log segment length",
    colour = "Group"
  ) +
  theme_classic()

p_segments

# Forefemur
mod_forefemur <- lm(
  log_forefemur ~ logP_c * group,
  data = dat
)

# Foretibia
mod_foretibia <- lm(
  log_foretibia ~ logP_c * group,
  data = dat
)

anova(mod_forefemur)
anova(mod_foretibia)

# 19. EAR–FORETIBIA RELATIONSHIP -----------------------------

p_ear <- ggplot(
  dat,
  aes(
    x = log_foretibia,
    y = log_ear_linear,
    colour = group
  )
) +
  geom_point(alpha = 0.55) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE
  ) +
  labs(
    x = "Log foretibia length",
    y = "Log linearized ear size",
    colour = "Group"
  ) +
  theme_classic()

p_ear

mod_ear_full <- lm(
  log_ear_linear ~
    logP_c * group +
    log_foretibia * group,
  data = dat
)

anova(mod_ear_full)
summary(mod_ear_full)

mod_ear_common <- lm(
  log_ear_linear ~
    logP_c * group +
    log_foretibia,
  data = dat
)

anova(
  mod_ear_common,
  mod_ear_full
)

ear_tibia_slopes <- emtrends(
  mod_ear_full,
  ~ group,
  var = "log_foretibia"
)

ear_tibia_slopes

pairs(
  ear_tibia_slopes,
  adjust = "holm"
)

# 20. EYE–HEAD RELATIONSHIP ----------------------------------

p_eye <- ggplot(
  dat,
  aes(
    x = log_head_size,
    y = log_eye,
    colour = group
  )
) +
  geom_point(alpha = 0.55) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE
  ) +
  labs(
    x = "Log head size",
    y = "Log eye length",
    colour = "Group"
  ) +
  theme_classic()

p_eye

mod_eye_full <- lm(
  log_eye ~
    logP_c * group +
    log_head_size * group,
  data = dat
)

anova(mod_eye_full)
summary(mod_eye_full)

mod_eye_common <- lm(
  log_eye ~
    logP_c * group +
    log_head_size,
  data = dat
)

anova(
  mod_eye_common,
  mod_eye_full
)

eye_head_slopes <- emtrends(
  mod_eye_full,
  ~ group,
  var = "log_head_size"
)

eye_head_slopes

pairs(
  eye_head_slopes,
  adjust = "holm"
)

# 21. DIAGNOSTICS --------------------------------------------

performance::check_model(mod_fore_common)
check_model(mod_mid_common)
check_model(mod_hind_common)

check_model(mod_ear_full)
check_model(mod_eye_full)


# 22. SAVE OUTPUTS -------------------------------------------

write_csv(
  group_summary,
  "weta_group_descriptive_statistics.csv"
)

write_csv(
  adjusted_leg_means,
  "weta_adjusted_leg_means.csv"
)

write_csv(
  leg_contrasts,
  "weta_leg_pairwise_contrasts.csv"
)

ggsave(
  "weta_leg_allometries.png",
  p_leg_allometry,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  "weta_adjusted_leg_means.png",
  p_adjusted_legs,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  "weta_ear_foretibia.png",
  p_ear,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  "weta_eye_head.png",
  p_eye,
  width = 7,
  height = 5,
  dpi = 300
)
