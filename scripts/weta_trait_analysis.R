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
