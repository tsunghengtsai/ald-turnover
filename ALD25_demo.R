# Source code to demonstrate a nonlinear regression-based approach for the 
# analysis of kinetic data in terms of inference and testing about turnover rate

# Required packages -------------------------------------------------------

library(tidyverse)
library(broom)

# Functions ---------------------------------------------------------------

fit_sep <- function(df) {
    df$peptide <- factor(df$peptide)
    init <- df |> 
        group_by(peptide) |> 
        summarise(min = min(labeling), amp = max(labeling) - min(labeling))
    
    nls(
        labeling ~ label_min[peptide] + label_amp[peptide] * (1 - exp(-k[condition] * time)), 
        data = df, 
        start = list(label_min = init$min, label_amp = init$amp, k = rep(0.02, 2))
    )
}

fit_shared <- function(df) {
    df$peptide <- factor(df$peptide)
    init <- df |> 
        group_by(peptide) |> 
        summarise(min = min(labeling), amp = max(labeling) - min(labeling))
    
    nls(
        labeling ~ label_min[peptide] + label_amp[peptide] * (1 - exp(-k * time)), 
        data = df, 
        start = list(label_min = init$min, label_amp = init$amp, k = 0.02)
    )
}

# Sample size calculation with simulations --------------------------------

# Time in days
t_grid <- (0:504) / 24

# Consider 9 time points
t_obs <- c(0, 1, 3, 8, 24, 72, 168, 288, 504) / 24

# Baseline and plateau labeling
sim_min <- 0.55
sim_max <- 0.85

# Turnover rate
sim_k <- 0.15

# Standard deviation of random error
sim_s <- 0.015

# Underlying curve based on the one-compartment model
sim_curve <- sim_min + (sim_max - sim_min) * (1 - exp(-sim_k * t_grid))

ggplot() + 
    geom_line(aes(x = t_grid, y = sim_curve)) + 
    labs(x = "Time", y = "Total labeling")

# Simulate kinetic data for 3 peptides corresponding to the same protein,
# with the same turnover rate and baseline and plateau labeling.
# (The labeling condition is unnecessary, as the labeling parameters for 
# different peptides will be estimated separately.)

set.seed(1)
n_peptide <- 3
l_peptide <- vector("list", n_peptide)
for (i in 1:n_peptide) {
    l_peptide[[i]] <- data.frame(time = t_grid, labeling = sim_curve) |> 
        filter(time %in% t_obs) |> 
        mutate(labeling = labeling + rnorm(length(t_obs), mean = 0, sd = sim_s)) |> 
        mutate(peptide = str_c("P", i))
}
df <- bind_rows(l_peptide)

# Plot the simulated kinetic data (colored points)
ggplot(df, aes(time, labeling)) + 
    geom_point(aes(color = peptide)) + 
    geom_line(data = data.frame(time = t_grid, labeling = sim_curve))

# Simulate kinetic data with an increased turnover rate (sim_k + 0.05), for another group
sim_curve2 <- sim_min + (sim_max - sim_min) * (1 - exp(-(sim_k + 0.05) * t_grid))

l_peptide2 <- vector("list", n_peptide)
for (i in 1:n_peptide) {
    l_peptide2[[i]] <- data.frame(time = t_grid, labeling = sim_curve2) |> 
        filter(time %in% t_obs) |> 
        mutate(labeling = labeling + rnorm(length(t_obs), mean = 0, sd = sim_s)) |> 
        mutate(peptide = str_c("P", i))
}
df2 <- bind_rows(l_peptide2)

# Show data in both groups
combined <- bind_rows(
    df |> mutate(condition = "C1"),
    df2 |> mutate(condition = "C2")
) |> 
    mutate(condition = factor(condition))

ggplot(combined, aes(time, labeling)) + 
    geom_point(aes(color = peptide)) + 
    facet_wrap(~ condition)

# Fit 2 models: one with separate turnover rates, the other with shared turnover rate
mod_full <- fit_sep(combined)
mod_reduced <- fit_shared(combined)

# Display the parameter estimates in the full model
#  - baseline labeling for each peptide
#  - difference between plateau and baseline labeling for each peptide
#  - turnover rates in both conditions (k1 & k2)
tidy(mod_full)

# Display the parameter estimates in the reduced model
tidy(mod_reduced)

# Compare the model fits with the reduced model being the null model
# (if null model is rejected -> the notion of shared turnover rate is rejected
# -> the difference in turnover rate is significant)
anova(mod_reduced, mod_full)

