###############################################################################
# Bayesian model-based meta-analysis of targeted biologics in generalised
# myasthenia gravis (gMG): hierarchical nonlinear Emax models (brms/Stan)
#
# Run from the repository root:  source("R/gMG_MBMA_analysis.R")
# Input:   data/gMG_MBMA_extraction.xlsx
# Output:  outputs/
###############################################################################

# SETUP
data_path <- normalizePath(file.path("data", "gMG_MBMA_extraction.xlsx"), mustWork = TRUE)
dir.create("outputs", showWarnings = FALSE)
output_dir <- normalizePath("outputs", mustWork = TRUE)
setwd(output_dir)

library(readxl)
library(dplyr)
library(brms)
library(tidybayes)
library(ggplot2)
library(tidyr)
library(stringr)
library(patchwork)
library(gridExtra)
library(grid)
library(pracma)
library(forcats)

set.seed(2024)

###############################################################################
# SECTION 1: MG-ADL DATA LOADING & PREPARATION
###############################################################################

data_raw <- read_excel(data_path, sheet = "MGADL_change")

data <- data_raw %>%
  mutate(
    study = as.factor(study),
    Trial = Trial_acronym,
    Drug = as.factor(ifelse(Drug == "PLACEBO", "Placebo", as.character(Drug))),
    Mechanism = as.factor(ifelse(is.na(Mechanism) | Mechanism == "", "Placebo", Mechanism)),
    Arm = as.factor(Arm),
    Regimen_type = as.factor(Regimen_type),
    time = as.numeric(time_weeks),
    change = mgadl_change,
    se = as.numeric(se_mgadl),
    dose_eq_perweek = as.numeric(Dose_eq_perweek),
    study_arm = as.factor(paste(study, Arm, sep = "_")),
    weight = 1 / se^2
  ) %>%
  filter(time > 0, !is.na(change), !is.na(se), se > 0) %>%
  arrange(study_arm, time)

# Treatment arms
trt <- data %>% filter(Arm == "Treatment")

# Split regimens
trt_continuous <- trt %>% filter(!str_detect(Regimen_type, "Cyclical|Weekly"))
trt_cyclical <- trt %>% filter(str_detect(Regimen_type, "Cyclical|Weekly"))

# Continuous: Cumulative dose exposure
trt_continuous <- trt_continuous %>%
  group_by(study_arm) %>%
  mutate(
    interval = time - lag(time, default = 0),
    cum_dose = cumsum(dose_eq_perweek * interval)
  ) %>%
  ungroup()

# Cyclical: On-period average dose exposure
trt_cycl_onperiod <- trt_cyclical %>%
  group_by(study_arm) %>%
  mutate(
    dose_eq_on = mean(dose_eq_perweek[dose_eq_perweek > 0], na.rm = TRUE),
    dose_eq = dose_eq_on
  ) %>%
  ungroup()

# Placebo data
plac <- data %>% filter(Arm == "Control", time > 0)
stopifnot(all(plac$se > 0))


###############################################################################
# SECTION 2: MG-ADL EXPLORATORY PLOTS
###############################################################################

# 1. Drug vs Placebo by Trial
p_drug_vs_placebo <- ggplot(data, aes(x = time, y = change, color = Drug, linetype = Arm)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_line(size = 1.2) +
  geom_errorbar(aes(ymin = change - 1.96*se, ymax = change + 1.96*se), width = 0.6, alpha = 0.8) +
  facet_wrap(~ Trial, scales = "free_y", ncol = 3) +
  scale_linetype_manual(values = c("Treatment" = "solid", "Control" = "dashed")) +
  scale_color_brewer(palette = "Paired") +
  labs(
    title = "MG-ADL Trajectories: Each Drug vs Trial-Specific Placebo",
    subtitle = "Solid line = Treatment | Dashed line = Placebo | Points with 95% CI",
    x = "Weeks", y = "MG-ADL Change from Baseline",
    color = "Drug", linetype = "Arm"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12))

print(p_drug_vs_placebo)
ggsave("mgadl_drug_vs_placebo_by_trial.png", p_drug_vs_placebo, width = 18, height = 14, dpi = 400, bg = "white")

# 2. Mechanism-Specific
p_mechanism <- ggplot(data %>% filter(Arm == "Treatment"), aes(x = time, y = change, color = Drug, group = study_arm)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_line(size = 1.2) +
  geom_errorbar(aes(ymin = change - 1.96*se, ymax = change + 1.96*se), width = 0.6, alpha = 0.8) +
  facet_wrap(~ Mechanism, scales = "free_y", ncol = 2) +
  scale_color_brewer(palette = "Paired") +
  labs(
    title = "MG-ADL Trajectories by Mechanism (Treatment Arms Only)",
    subtitle = "Points with 95% CI | Lines connect study arms",
    x = "Weeks", y = "MG-ADL Change from Baseline", color = "Drug"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12))

print(p_mechanism)
ggsave("mgadl_trajectories_by_mechanism.png", p_mechanism, width = 14, height = 10, dpi = 400, bg = "white")

# 3. Regimen-Specific
p_regimen <- ggplot(data %>% filter(Arm == "Treatment"), aes(x = time, y = change, color = Drug, group = study_arm)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_line(size = 1.2) +
  geom_errorbar(aes(ymin = change - 1.96*se, ymax = change + 1.96*se), width = 0.6, alpha = 0.8) +
  facet_wrap(~ Regimen_type, scales = "free_y", ncol = 2) +
  scale_color_brewer(palette = "Paired") +
  labs(
    title = "MG-ADL Trajectories by Regimen Type (Treatment Arms Only)",
    subtitle = "Points with 95% CI | Lines connect study arms",
    x = "Weeks", y = "MG-ADL Change from Baseline", color = "Drug"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12))

print(p_regimen)
ggsave("mgadl_trajectories_by_regimen.png", p_regimen, width = 14, height = 10, dpi = 400, bg = "white")

# Combined collage
combined_plot <- (p_drug_vs_placebo +
                    theme(plot.margin = margin(10, 10, 10, 10),
                          panel.border = element_rect(color = "grey70", fill = NA, size = 1)) +
                    ggtitle("Drug vs Trial-Specific Placebo")) +
  (p_mechanism +
     theme(plot.margin = margin(10, 10, 10, 10),
           panel.border = element_rect(color = "grey70", fill = NA, size = 1)) +
     ggtitle("Mechanism-Specific Trajectories")) +
  (p_regimen +
     theme(plot.margin = margin(10, 10, 10, 10),
           panel.border = element_rect(color = "grey70", fill = NA, size = 1)) +
     ggtitle("Regimen-Specific Trajectories")) +
  plot_layout(ncol = 1) +
  plot_annotation(
    title = "MG-ADL Trajectory Comparisons",
    subtitle = "Post-baseline data | Points with 95% CI | Grey borders separate panels",
    caption = "Top: Drug vs Placebo by Trial | Middle: By Mechanism | Bottom: By Regimen Type"
  ) &
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        text = element_text(size = 12))

print(combined_plot)
ggsave("mgadl_trajectory_collage.png", combined_plot, width = 18, height = 24, dpi = 400, bg = "white")


###############################################################################
# SECTION 3: MG-ADL MODEL FITTING
###############################################################################

cat("\n=== Fitting MG-ADL Continuous Model ===\n")
fit_continuous <- brm(
  bf(
    change | se(se, sigma = TRUE) ~ -emax * (cum_dose / (ed50 + cum_dose)) * (1 - exp(-k * time)),
    emax ~ 1 + (1 | ID | Drug),
    ed50 ~ 1,
    k ~ 1,
    nl = TRUE
  ),
  data = trt_continuous,
  prior = c(
    prior(normal(5, 3), nlpar = "emax", lb = 0),
    prior(exponential(2), class = "sd", nlpar = "emax"),
    prior(normal(5000, 5000), nlpar = "ed50", lb = 0),
    prior(normal(0.1, 0.1), nlpar = "k", lb = 0)
  ),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4, iter = 12000, warmup = 4000,
  seed = 123
)

cat("\n=== Fitting MG-ADL Cyclical Model ===\n")
fit_cyclical <- brm(
  bf(
    change | se(se, sigma = TRUE) ~ -emax * (dose_eq / (ed50 + dose_eq)) * (1 - exp(-k * time)),
    emax ~ 1 + (1 | ID | Drug),
    ed50 ~ 1,
    k ~ 1,
    nl = TRUE
  ),
  data = trt_cycl_onperiod,
  prior = c(
    prior(normal(6, 3), nlpar = "emax", lb = 0),
    prior(exponential(2), class = "sd", nlpar = "emax"),
    prior(normal(200, 400), nlpar = "ed50", lb = 0),
    prior(normal(0.4, 0.2), nlpar = "k", lb = 0)
  ),
  control = list(adapt_delta = 0.995, max_treedepth = 15),
  chains = 4, iter = 15000, warmup = 5000,
  seed = 123
)

cat("\n=== Fitting MG-ADL Placebo Model ===\n")
fit_plac_exp <- brm(
  bf(
    change | se(se, sigma = TRUE) ~ plateau * (1 - exp(-k * time)),
    plateau ~ 1 + (1 | ID | study),
    k ~ 1 + (1 | ID | study),
    nl = TRUE
  ),
  data = plac,
  prior = c(
    prior(normal(-3, 1), nlpar = "plateau", lb = -10),
    prior(normal(0.1, 0.1), nlpar = "k", lb = 0),
    prior(exponential(2), class = "sd", nlpar = "plateau"),
    prior(exponential(2), class = "sd", nlpar = "k")
  ),
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  chains = 4, iter = 20000, warmup = 10000,
  seed = 123
)


###############################################################################
# SECTION 4: MG-ADL MODEL DIAGNOSTICS
###############################################################################

diag_cont <- as.data.frame(summary(fit_continuous)$fixed) %>%
  tibble::rownames_to_column("Parameter") %>%
  mutate(`95% CrI` = paste(round(`l-95% CI`, 2), "–", round(`u-95% CI`, 2))) %>%
  select(Parameter, Estimate, `Est.Error`, `95% CrI`, Rhat, Bulk_ESS, Tail_ESS)

diag_cycl <- as.data.frame(summary(fit_cyclical)$fixed) %>%
  tibble::rownames_to_column("Parameter") %>%
  mutate(`95% CrI` = paste(round(`l-95% CI`, 2), "–", round(`u-95% CI`, 2))) %>%
  select(Parameter, Estimate, `Est.Error`, `95% CrI`, Rhat, Bulk_ESS, Tail_ESS)

print("Continuous Model Diagnostics")
print(diag_cont)
write.csv(diag_cont, "continuous_model_diagnostics.csv", row.names = FALSE)

print("Cyclical Model Diagnostics")
print(diag_cycl)
write.csv(diag_cycl, "cyclical_model_diagnostics.csv", row.names = FALSE)

# Diagnostics tables as PDF
table_cont <- tableGrob(diag_cont, rows = NULL, theme = ttheme_default(base_size = 10))
title_cont <- textGrob("Continuous Regimens Model Diagnostics", gp = gpar(fontface = "bold", fontsize = 14))
table_cycl <- tableGrob(diag_cycl, rows = NULL, theme = ttheme_default(base_size = 10))
title_cycl <- textGrob("Cyclical/Weekly Regimens Model Diagnostics", gp = gpar(fontface = "bold", fontsize = 14))

combined_diag_layout <- arrangeGrob(
  title_cont, table_cont, title_cycl, table_cycl,
  ncol = 1, heights = c(0.5, 2, 0.5, 2)
)
ggsave("model_diagnostics_combined.pdf", combined_diag_layout, width = 12, height = 16, dpi = 400, device = "pdf")

# PP checks
pp_check(fit_continuous, ndraws = 50) + ggtitle("Continuous Model PP Check")
ggsave("pp_check_continuous.png", width = 10, height = 8, dpi = 400)

pp_check(fit_cyclical, ndraws = 50) + ggtitle("Cyclical Model PP Check")
ggsave("pp_check_cyclical.png", width = 10, height = 8, dpi = 400)

# Placebo diagnostics
cat("\n=== Placebo Model Summary ===\n")
print(summary(fit_plac_exp))


###############################################################################
# SECTION 5: MG-ADL VPCs
###############################################################################

trt_pred_cont <- trt_continuous %>%
  add_predicted_draws(fit_continuous, ndraws = 100, allow_new_levels = TRUE)

trt_pred_summary_cont <- trt_pred_cont %>%
  group_by(Drug, time) %>%
  median_hdci(.prediction)

p_vpc_cont <- ggplot() +
  geom_errorbar(data = trt_continuous, aes(x = time, ymin = change - 1.96*se, ymax = change + 1.96*se, color = Drug), width = 0.5, alpha = 0.8) +
  geom_point(data = trt_continuous, aes(x = time, y = change, color = Drug), size = 3) +
  geom_line(data = trt_pred_summary_cont, aes(x = time, y = .prediction, color = Drug), size = 1.2) +
  geom_ribbon(data = trt_pred_summary_cont, aes(x = time, y = .prediction, ymin = .lower, ymax = .upper, fill = Drug), alpha = 0.3) +
  facet_wrap(~ Drug, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") +
  labs(title = "VPC: Continuous Regimens (Cumulative Dose Model)",
       subtitle = "Observed points with 95% CI | Predicted median line + 95% prediction interval ribbon",
       x = "Weeks", y = "MG-ADL Change from Baseline") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 12, color = "darkblue"),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12), panel.grid.minor = element_blank())

print(p_vpc_cont)
ggsave("vpc_continuous_regimens.png", p_vpc_cont, width = 16, height = 10, dpi = 400, bg = "white")

trt_pred_cycl <- trt_cycl_onperiod %>%
  add_predicted_draws(fit_cyclical, ndraws = 100, allow_new_levels = TRUE)

trt_pred_summary_cycl <- trt_pred_cycl %>%
  group_by(Drug, time) %>%
  median_hdci(.prediction)

p_vpc_cycl <- ggplot() +
  geom_errorbar(data = trt_cycl_onperiod, aes(x = time, ymin = change - 1.96*se, ymax = change + 1.96*se, color = Drug), width = 0.5, alpha = 0.8) +
  geom_point(data = trt_cycl_onperiod, aes(x = time, y = change, color = Drug), size = 3) +
  geom_line(data = trt_pred_summary_cycl, aes(x = time, y = .prediction, color = Drug), size = 1.2) +
  geom_ribbon(data = trt_pred_summary_cycl, aes(x = time, y = .prediction, ymin = .lower, ymax = .upper, fill = Drug), alpha = 0.3) +
  facet_wrap(~ Drug, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") +
  labs(title = "VPC: Cyclical/Weekly Regimens (On-Period Average Dose Model)",
       subtitle = "Observed points with 95% CI | Predicted median line + 95% prediction interval ribbon",
       x = "Weeks", y = "MG-ADL Change from Baseline") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 12, color = "darkblue"),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12), panel.grid.minor = element_blank())

print(p_vpc_cycl)
ggsave("vpc_cyclical_regimens.png", p_vpc_cycl, width = 16, height = 8, dpi = 400, bg = "white")

p_combined_vpc <- p_vpc_cont / p_vpc_cycl +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Visual Predictive Checks – MG-ADL Models",
    subtitle = "Continuous Regimens (Cumulative Dose) vs Cyclical/Weekly Regimens (On-Period Average Dose)",
    caption = "Observed points with 95% CI | Predicted median line + 95% prediction interval ribbon",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 14, hjust = 0.5),
                  plot.caption = element_text(size = 11, hjust = 1, colour = "gray50")))

print(p_combined_vpc)
ggsave("mgadl_vpc_combined.png", p_combined_vpc, width = 16, height = 12, dpi = 400, bg = "white")


###############################################################################
# SECTION 6: MG-ADL FITTED TRAJECTORIES
###############################################################################

trt_fitted_cont <- trt_continuous %>%
  add_epred_draws(fit_continuous, ndraws = 100, allow_new_levels = TRUE)
trt_fitted_summary_cont <- trt_fitted_cont %>% mean_hdci(.epred)

p_fitted_cont <- ggplot() +
  geom_errorbar(data = trt_continuous, aes(x = time, ymin = change - 1.96*se, ymax = change + 1.96*se, color = Drug), width = 0.5, alpha = 0.8) +
  geom_point(data = trt_continuous, aes(x = time, y = change, color = Drug), size = 3) +
  geom_line(data = trt_fitted_summary_cont, aes(x = time, y = .epred, color = Drug), size = 1.2) +
  geom_ribbon(data = trt_fitted_summary_cont, aes(x = time, y = .epred, ymin = .lower, ymax = .upper, fill = Drug), alpha = 0.3) +
  facet_wrap(~ Drug, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") +
  labs(title = "Fitted Trajectories: Continuous Regimens (Cumulative Dose Model)",
       subtitle = "Observed points with 95% CI | Predicted mean line + 95% HDCI ribbon",
       x = "Weeks", y = "MG-ADL Change from Baseline") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 12, color = "darkblue"),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12), panel.grid.minor = element_blank())

print(p_fitted_cont)
ggsave("fitted_trajectories_continuous.png", p_fitted_cont, width = 16, height = 10, dpi = 400, bg = "white")

trt_fitted_cycl <- trt_cycl_onperiod %>%
  add_epred_draws(fit_cyclical, ndraws = 100, allow_new_levels = TRUE)
trt_fitted_summary_cycl <- trt_fitted_cycl %>% mean_hdci(.epred)

p_fitted_cycl <- ggplot() +
  geom_errorbar(data = trt_cycl_onperiod, aes(x = time, ymin = change - 1.96*se, ymax = change + 1.96*se, color = Drug), width = 0.5, alpha = 0.8) +
  geom_point(data = trt_cycl_onperiod, aes(x = time, y = change, color = Drug), size = 3) +
  geom_line(data = trt_fitted_summary_cycl, aes(x = time, y = .epred, color = Drug), size = 1.2) +
  geom_ribbon(data = trt_fitted_summary_cycl, aes(x = time, y = .epred, ymin = .lower, ymax = .upper, fill = Drug), alpha = 0.3) +
  facet_wrap(~ Drug, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") +
  labs(title = "Fitted Trajectories: Cyclical/Weekly Regimens (On-Period Average Dose Model)",
       subtitle = "Observed points with 95% CI | Predicted mean line + 95% HDCI ribbon",
       x = "Weeks", y = "MG-ADL Change from Baseline") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 12, color = "darkblue"),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12), panel.grid.minor = element_blank())

print(p_fitted_cycl)
ggsave("fitted_trajectories_cyclical.png", p_fitted_cycl, width = 16, height = 8, dpi = 400, bg = "white")


###############################################################################
# SECTION 7: MG-ADL EMAX FOREST PLOTS
###############################################################################

# Emax from Continuous model
emax_cont <- as_draws_df(fit_continuous) %>%
  select(.draw, starts_with("r_Drug"), b_emax_Intercept) %>%
  pivot_longer(cols = starts_with("r_Drug"), names_to = "temp", values_to = "offset") %>%
  mutate(
    Drug = str_extract(temp, "(?<=\\[)[^,]+"),
    emax = b_emax_Intercept + offset,
    regimen = "Continuous"
  ) %>% filter(!is.na(Drug))

# Emax from Cyclical model
emax_cycl <- as_draws_df(fit_cyclical) %>%
  select(.draw, starts_with("r_Drug"), b_emax_Intercept) %>%
  pivot_longer(cols = starts_with("r_Drug"), names_to = "temp", values_to = "offset") %>%
  mutate(
    Drug = str_extract(temp, "(?<=\\[)[^,]+"),
    emax = b_emax_Intercept + offset,
    regimen = "Cyclical/Weekly"
  ) %>% filter(!is.na(Drug))

emax_combined <- bind_rows(emax_cont, emax_cycl)

emax_summary_combined <- emax_combined %>%
  group_by(Drug, regimen) %>%
  mean_qi(emax, .width = 0.95) %>%
  mutate(label = sprintf("%.1f [%.1f, %.1f]", emax, .lower, .upper)) %>%
  arrange(desc(emax))

label_x_pos <- max(emax_summary_combined$.upper) + 0.5

p_emax_combined <- ggplot(emax_combined, aes(x = emax, color = regimen)) +
  stat_halfeye(aes(y = reorder(Drug, emax, mean)), .width = c(0.66, 0.95), slab_alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray50") +
  geom_text(data = emax_summary_combined, aes(y = reorder(Drug, emax, mean), label = label, x = label_x_pos),
            hjust = "left", size = 4.5, fontface = "bold") +
  labs(title = "Drug-Specific Posterior Emax for MG-ADL (Distinct Models by Regimen)",
       subtitle = "Higher = greater long-term efficacy | Continuous vs Cyclical/Weekly",
       x = "Emax (points improvement)", y = NULL, color = "Regimen Type") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.margin = margin(10, 140, 10, 10))

print(p_emax_combined)
ggsave("emax_per_drug_distinct_models_forest.png", p_emax_combined,
       width = 14, height = max(6, 0.8 * n_distinct(emax_combined$Drug)), dpi = 400, bg = "white")


###############################################################################
# SECTION 8: MG-ADL NET ADDED BENEFIT OVER PLACEBO
###############################################################################

placebo_plateau_draws <- as_draws_df(fit_plac_exp) %>%
  select(.draw, b_plateau_Intercept) %>%
  mutate(placebo_mag = abs(b_plateau_Intercept))

placebo_mag_mean <- mean(placebo_plateau_draws$placebo_mag)

emax_cont_net <- emax_cont %>%
  mutate(net_emax = emax - placebo_mag_mean)

emax_cycl_net <- emax_cycl %>%
  mutate(net_emax = emax - placebo_mag_mean)

emax_net_all <- bind_rows(emax_cont_net, emax_cycl_net)

emax_net_summary <- emax_net_all %>%
  group_by(Drug, regimen) %>%
  mean_qi(net_emax, .width = 0.95) %>%
  mutate(label = sprintf("%.1f [%.1f, %.1f]", net_emax, .lower, .upper)) %>%
  arrange(desc(net_emax))

label_x_pos <- max(emax_net_summary$.upper) + 0.5

p_net_emax <- ggplot(emax_net_all, aes(x = net_emax, color = regimen)) +
  stat_halfeye(aes(y = reorder(Drug, net_emax, mean)), .width = c(0.66, 0.95), slab_alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray50") +
  geom_text(data = emax_net_summary, aes(y = reorder(Drug, net_emax, mean), label = label, x = label_x_pos),
            hjust = "left", size = 4.5, fontface = "bold") +
  labs(title = "Drug-Specific Net Added Improvement vs Placebo for MG-ADL",
       subtitle = "Net = Drug Emax magnitude - Placebo plateau magnitude | Higher = greater added benefit",
       x = "Net Added Improvement over Placebo (points)", y = NULL, color = "Regimen Type") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.margin = margin(10, 140, 10, 10))

print(p_net_emax)
ggsave("net_added_emax_vs_placebo_forest.png", p_net_emax, width = 14, height = 10, dpi = 400, bg = "white")


###############################################################################
# SECTION 9: MG-ADL WEEK-26 PLACEBO-ADJUSTED + SUPERIORITY
###############################################################################

placebo_week26 <- plac %>% mutate(time = 26)
placebo_pred <- posterior_predict(fit_plac_exp, newdata = placebo_week26, ndraws = 1000)
placebo_mean <- mean(placebo_pred)

# Continuous drugs
pred_week26_cont <- trt_continuous %>%
  group_by(Drug, Mechanism) %>%
  summarise(dose_eq_mean = mean(dose_eq_perweek, na.rm = TRUE), .groups = "drop") %>%
  mutate(time = 26, cum_dose = dose_eq_mean * 26, dose_eq = dose_eq_mean, se = 0)

pred_cont_matrix <- posterior_predict(fit_continuous, newdata = pred_week26_cont, ndraws = 1000, allow_new_levels = TRUE)

pred_cont_post <- pred_cont_matrix %>%
  as.data.frame() %>% mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "obs_str", values_to = "pred_change") %>%
  mutate(obs = as.integer(str_remove(obs_str, "^V"))) %>% select(-obs_str) %>%
  left_join(pred_week26_cont %>% mutate(obs = row_number()), by = "obs") %>%
  mutate(regimen = "Continuous")

# Cyclical drugs
pred_week26_cycl <- trt_cycl_onperiod %>%
  group_by(Drug, Mechanism) %>%
  summarise(dose_eq = mean(dose_eq, na.rm = TRUE), .groups = "drop") %>%
  mutate(time = 26, se = 0)

pred_cycl_matrix <- posterior_predict(fit_cyclical, newdata = pred_week26_cycl, ndraws = 1000, allow_new_levels = TRUE)

pred_cycl_post <- pred_cycl_matrix %>%
  as.data.frame() %>% mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "obs_str", values_to = "pred_change") %>%
  mutate(obs = as.integer(str_remove(obs_str, "^V"))) %>% select(-obs_str) %>%
  left_join(pred_week26_cycl %>% mutate(obs = row_number()), by = "obs") %>%
  mutate(regimen = "Cyclical/Weekly")

pred_26_all_post <- bind_rows(pred_cont_post, pred_cycl_post)
pred_26_net_post <- pred_26_all_post %>% mutate(net_change = pred_change - placebo_mean)

pred_26_net_summary <- pred_26_net_post %>%
  group_by(Drug, Mechanism, regimen) %>%
  mean_qi(net_change, .width = 0.95) %>%
  mutate(label = sprintf("%.1f [%.1f, %.1f]", net_change, .lower, .upper)) %>%
  arrange(net_change)

p_net_week26 <- ggplot(pred_26_net_post, aes(x = net_change)) +
  stat_halfeye(aes(y = reorder(Drug, net_change, mean), fill = Mechanism),
               .width = c(0.66, 0.95), slab_alpha = 0.85, point_interval = mean_qi) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray40", size = 1) +
  geom_text(data = pred_26_net_summary,
            aes(y = reorder(Drug, net_change, mean), label = label, x = min(.lower) - 0.3),
            hjust = "right", size = 4.5, fontface = "bold", colour = "black") +
  scale_fill_manual(values = c("Complement inhibitor" = "#1f77b4", "FcRn inhibitor" = "#ff7f0e")) +
  labs(title = "Placebo-Adjusted Predicted MG-ADL Change at Week 26 by Drug",
       subtitle = "Net added benefit over placebo (posterior mean [95% credible interval]) | More negative = greater improvement",
       caption = "On-period average dose for cyclical regimens | Ranked by greatest net benefit",
       x = "Net Change in MG-ADL Score at Week 26 (points)", y = NULL, fill = "Mechanism") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, colour = "gray30"),
        plot.caption = element_text(size = 10, hjust = 1, colour = "gray50"),
        axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12),
        legend.position = "bottom", legend.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(20, 140, 20, 20))

print(p_net_week26)
ggsave("placebo_adjusted_week26_forest_publication.png", p_net_week26, width = 14, height = 10, dpi = 400, bg = "white")

# Superiority probabilities
placebo_draws <- posterior_predict(fit_plac_exp, newdata = plac %>% mutate(time = 26), ndraws = 1000) %>% as.vector()

pred_week26_cont_sup <- trt_continuous %>%
  group_by(Drug) %>%
  summarise(dose_eq_mean = mean(dose_eq_perweek, na.rm = TRUE), .groups = "drop") %>%
  mutate(time = 26, cum_dose = dose_eq_mean * 26, dose_eq = dose_eq_mean, se = 0)

pred_cont_matrix_sup <- posterior_predict(fit_continuous, newdata = pred_week26_cont_sup, ndraws = 1000, allow_new_levels = TRUE)

pred_cont_post_sup <- pred_cont_matrix_sup %>%
  as.data.frame() %>% mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "obs_str", values_to = "drug_change") %>%
  mutate(obs = as.integer(str_remove(obs_str, "^V"))) %>% select(-obs_str) %>%
  left_join(pred_week26_cont_sup %>% mutate(obs = row_number()), by = "obs") %>%
  mutate(regimen = "Continuous")

pred_week26_cycl_sup <- trt_cycl_onperiod %>%
  group_by(Drug) %>%
  summarise(dose_eq = mean(dose_eq, na.rm = TRUE), .groups = "drop") %>%
  mutate(time = 26, se = 0)

pred_cycl_matrix_sup <- posterior_predict(fit_cyclical, newdata = pred_week26_cycl_sup, ndraws = 1000, allow_new_levels = TRUE)

pred_cycl_post_sup <- pred_cycl_matrix_sup %>%
  as.data.frame() %>% mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "obs_str", values_to = "drug_change") %>%
  mutate(obs = as.integer(str_remove(obs_str, "^V"))) %>% select(-obs_str) %>%
  left_join(pred_week26_cycl_sup %>% mutate(obs = row_number()), by = "obs") %>%
  mutate(regimen = "Cyclical/Weekly")

pred_26_drug_post <- bind_rows(pred_cont_post_sup, pred_cycl_post_sup)
pred_26_net_sup_post <- pred_26_drug_post %>%
  mutate(net_change = drug_change - placebo_draws[draw])

superior_prob <- pred_26_net_sup_post %>%
  group_by(Drug, regimen) %>%
  summarise(prob_superior = mean(net_change < 0) * 100, .groups = "drop")

net_summary <- pred_26_net_sup_post %>%
  group_by(Drug, regimen) %>%
  mean_qi(net_change, .width = 0.95) %>%
  mutate(label = sprintf("%.1f [%.1f, %.1f] (P=%.1f%%)", net_change, .lower, .upper,
                         superior_prob$prob_superior[match(Drug, superior_prob$Drug)])) %>%
  arrange(net_change)

p_net_superior <- ggplot(pred_26_net_sup_post, aes(x = net_change, color = regimen)) +
  stat_halfeye(aes(y = reorder(Drug, net_change, mean)), .width = c(0.66, 0.95), slab_alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray5") +
  geom_text(data = net_summary, aes(y = reorder(Drug, net_change, mean), label = label, x = min(.lower) - 0.5),
            hjust = "right", size = 4.5, fontface = "bold") +
  labs(title = "Net Added Benefit vs Placebo at Week 26 with Superiority Probability",
       x = "Net ΔMG-ADL at Week 26 (points)", y = NULL, color = "Regimen Type") +
  scale_color_manual(values = c("Continuous" = "#1B9E77", "Cyclical/Weekly" = "#D95F02")) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 13, hjust = 0.5),
        legend.position = "bottom", legend.title = element_text(face = "bold"),
        axis.text = element_text(size = 8),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        plot.margin = margin(20, 30, 20, 20))

print(p_net_superior)
ggsave("net_benefit_with_superiority_week26.png", p_net_superior, width = 14, height = 10, dpi = 400)

###############################################################################
# SECTION 10: MG-ADL TIME TO MCID
###############################################################################

MCID <- -2
time_grid <- seq(1, 26, by = 0.5)

# Continuous
pred_time_grid_cont <- trt_continuous %>%
  group_by(Drug) %>%
  summarise(dose_eq_mean = mean(dose_eq_perweek, na.rm = TRUE), .groups = "drop") %>%
  crossing(time = time_grid) %>%
  mutate(cum_dose = dose_eq_mean * time, dose_eq = dose_eq_mean, se = 0)

pred_time_cont <- posterior_predict(fit_continuous, newdata = pred_time_grid_cont, ndraws = 500, allow_new_levels = TRUE)

pred_time_cont_post <- pred_time_cont %>%
  as.data.frame() %>% mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "obs_str", values_to = "pred_change") %>%
  mutate(obs = as.integer(str_remove(obs_str, "^V"))) %>% select(-obs_str) %>%
  left_join(pred_time_grid_cont %>% mutate(obs = row_number()), by = "obs") %>%
  mutate(regimen = "Continuous")

# Cyclical
pred_time_grid_cycl <- trt_cycl_onperiod %>%
  group_by(Drug) %>%
  summarise(dose_eq = mean(dose_eq, na.rm = TRUE), .groups = "drop") %>%
  crossing(time = time_grid) %>%
  mutate(se = 0)

pred_time_cycl <- posterior_predict(fit_cyclical, newdata = pred_time_grid_cycl, ndraws = 500, allow_new_levels = TRUE)

pred_time_cycl_post <- pred_time_cycl %>%
  as.data.frame() %>% mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "obs_str", values_to = "pred_change") %>%
  mutate(obs = as.integer(str_remove(obs_str, "^V"))) %>% select(-obs_str) %>%
  left_join(pred_time_grid_cycl %>% mutate(obs = row_number()), by = "obs") %>%
  mutate(regimen = "Cyclical/Weekly")

pred_time_all_post <- bind_rows(pred_time_cont_post, pred_time_cycl_post)

time_to_mcid <- pred_time_all_post %>%
  group_by(draw, Drug, regimen) %>%
  arrange(time) %>%
  filter(pred_change <= MCID) %>%
  slice(1) %>% ungroup()

time_mcid_summary <- time_to_mcid %>%
  group_by(Drug, regimen) %>%
  mean_qi(time, .width = 0.95) %>%
  mutate(label = sprintf("%.1f [%.1f, %.1f] weeks", time, .lower, .upper)) %>%
  arrange(time)

p_time_mcid <- ggplot(time_to_mcid, aes(x = time, color = regimen)) +
  stat_halfeye(aes(y = reorder(Drug, time, mean)), .width = c(0.66, 0.95), slab_alpha = 0.8) +
  geom_text(data = time_mcid_summary,
            aes(y = reorder(Drug, time, mean), label = label, x = max(.upper) + 2),
            hjust = "right", size = 4.5, fontface = "bold") +
  labs(title = "Drug-Specific Time to Achieve MCID (≥2-Point MG-ADL Improvement)",
       subtitle = "Posterior median [95% CrI] weeks | Shorter = faster meaningful benefit",
       x = "Time to MCID (weeks)", y = NULL, color = "Regimen Type") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "bottom", plot.margin = margin(10, 140, 10, 10))

print(p_time_mcid)
ggsave("time_to_mcid_per_drug_distinct_models.png", p_time_mcid, width = 16, height = 10, dpi = 400, bg = "white")

###############################################################################
# SECTION 11: MG-ADL AUEC
###############################################################################

time_grid_auec <- seq(0, 26, by = 1)

# Continuous
pred_auec_cont <- trt_continuous %>%
  group_by(Drug) %>%
  summarise(dose_eq_mean = mean(dose_eq_perweek, na.rm = TRUE), .groups = "drop") %>%
  crossing(time = time_grid_auec) %>%
  mutate(cum_dose = dose_eq_mean * time, dose_eq = dose_eq_mean, se = 0)

pred_auec_cont_matrix <- posterior_predict(fit_continuous, newdata = pred_auec_cont, ndraws = 500, allow_new_levels = TRUE)

pred_auec_cont_post <- pred_auec_cont_matrix %>%
  as.data.frame() %>% mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "obs_str", values_to = "pred_change") %>%
  mutate(obs = as.integer(str_remove(obs_str, "^V"))) %>% select(-obs_str) %>%
  left_join(pred_auec_cont %>% mutate(obs = row_number()), by = "obs") %>%
  mutate(regimen = "Continuous")

# Cyclical
pred_auec_cycl <- trt_cycl_onperiod %>%
  group_by(Drug) %>%
  summarise(dose_eq = mean(dose_eq, na.rm = TRUE), .groups = "drop") %>%
  crossing(time = time_grid_auec) %>%
  mutate(se = 0)

pred_auec_cycl_matrix <- posterior_predict(fit_cyclical, newdata = pred_auec_cycl, ndraws = 500, allow_new_levels = TRUE)

pred_auec_cycl_post <- pred_auec_cycl_matrix %>%
  as.data.frame() %>% mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "obs_str", values_to = "pred_change") %>%
  mutate(obs = as.integer(str_remove(obs_str, "^V"))) %>% select(-obs_str) %>%
  left_join(pred_auec_cycl %>% mutate(obs = row_number()), by = "obs") %>%
  mutate(regimen = "Cyclical/Weekly")

pred_auec_all_post <- bind_rows(pred_auec_cont_post, pred_auec_cycl_post)

auec_all <- pred_auec_all_post %>%
  group_by(draw, Drug, regimen) %>%
  arrange(time) %>%
  summarise(auec = trapz(time, abs(pred_change)), .groups = "drop")

auec_summary <- auec_all %>%
  group_by(Drug, regimen) %>%
  mean_qi(auec, .width = 0.95) %>%
  mutate(label = sprintf("%.0f [%.0f, %.0f]", auec, .lower, .upper)) %>%
  arrange(desc(auec))

label_x_pos <- max(auec_summary$.upper) + 50

p_auec <- ggplot(auec_all, aes(x = auec, color = regimen)) +
  stat_halfeye(aes(y = reorder(Drug, auec, mean)), .width = c(0.66, 0.95), slab_alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray50") +
  geom_text(data = auec_summary, aes(y = reorder(Drug, auec, mean), label = label, x = label_x_pos),
            hjust = "right", size = 4.5, fontface = "bold") +
  labs(title = "Drug-Specific Area Under the Effect Curve (AUEC, Week 0-26)",
       subtitle = "Total integrated improvement | Higher = greater cumulative benefit",
       x = "AUEC (point·weeks)", y = NULL, color = "Regimen Type") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.margin = margin(10, 140, 10, 10))

print(p_auec)
ggsave("auec_per_drug_distinct_models.png", p_auec, width = 14, height = 10, dpi = 400, bg = "white")

# Combined dashboard
p_combined_mgadl_dashboard <- p_time_mcid / p_auec +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(
    title = "Clinical Benefit Dashboard – MG-ADL Outcome",
    subtitle = "Time to meaningful improvement + Total cumulative benefit over 26 weeks",
    caption = "Posterior median [95% CrI] | Colour by regimen type",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 14, hjust = 0.5),
                  plot.caption = element_text(size = 11, hjust = 1, colour = "gray50")))

print(p_combined_mgadl_dashboard)
ggsave("mgadl_clinical_benefit_dashboard_combined.png", p_combined_mgadl_dashboard,
       width = 14, height = 12, dpi = 400, bg = "white")

###############################################################################
# SECTION 12: MG-ADL POOLED ED50 AND k
###############################################################################

ed50_k_continuous <- as_draws_df(fit_continuous) %>%
  select(.draw, b_ed50_Intercept, b_k_Intercept) %>%
  rename(ed50 = b_ed50_Intercept, k = b_k_Intercept) %>%
  mutate(regimen = "Continuous")

ed50_k_cyclical <- as_draws_df(fit_cyclical) %>%
  select(.draw, b_ed50_Intercept, b_k_Intercept) %>%
  rename(ed50 = b_ed50_Intercept, k = b_k_Intercept) %>%
  mutate(regimen = "Cyclical/Weekly")

ed50_k_combined <- bind_rows(ed50_k_continuous, ed50_k_cyclical)

ed50_k_table <- ed50_k_combined %>%
  group_by(regimen) %>%
  mean_qi(ed50, k, .width = 0.95) %>%
  mutate(
    ED50 = sprintf("%.0f [%.0f, %.0f]", ed50, ed50.lower, ed50.upper),
    `Onset Rate k` = sprintf("%.2f [%.2f, %.2f]", k, k.lower, k.upper)
  ) %>%
  select(regimen, ED50, `Onset Rate k`) %>% arrange(regimen)

print("Pooled ED50 and k by Regimen")
print(ed50_k_table)
write.csv(ed50_k_table, "pooled_ed50_k_by_regimen.csv", row.names = FALSE)

p_ed50 <- ggplot(ed50_k_combined, aes(x = ed50, fill = regimen)) +
  stat_halfeye(.width = c(0.66, 0.95), slab_alpha = 0.8) +
  labs(title = "Pooled ED50", x = "ED50 (mg/week equivalent)", y = "Density", fill = "Regimen") +
  theme_minimal() + theme(legend.position = "none")

p_k <- ggplot(ed50_k_combined, aes(x = k, fill = regimen)) +
  stat_halfeye(.width = c(0.66, 0.95), slab_alpha = 0.8) +
  labs(title = "Pooled Onset Rate k", x = "k (per week)", y = "Density", fill = "Regimen") +
  theme_minimal() + theme(legend.position = "bottom")

plots_combined <- p_ed50 + p_k + plot_layout(guides = "collect") & theme(legend.position = "bottom")

table_grob <- tableGrob(ed50_k_table, rows = NULL,
                         theme = ttheme_default(base_size = 12, core = list(fg_params = list(fontface = "bold"))))

final_figure <- plots_combined / table_grob +
  plot_annotation(
    title = "Pooled Pharmacodynamic Parameters by Regimen Type",
    subtitle = "ED50 (lower = higher potency) | k (higher = faster onset)",
    caption = "Posterior mean [95% CrI]"
  ) &
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 20, hjust = 0.5))

print(final_figure)
ggsave("pooled_ed50_k_by_regimen_combined.png", final_figure, width = 16, height = 12, dpi = 400, bg = "white")


###############################################################################
# SECTION 13: MG-ADL FULL MODEL DIAGNOSTICS PDF
###############################################################################

extract_summary <- function(fit) {
  fixed <- as.data.frame(summary(fit)$fixed) %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(`95% CrI` = paste(round(`l-95% CI`, 2), "–", round(`u-95% CI`, 2))) %>%
    select(Parameter, Estimate, `Est.Error`, `95% CrI`, Rhat, Bulk_ESS, Tail_ESS)

  random_param <- if ("Drug" %in% names(summary(fit)$random)) "Drug" else "study"
  random <- as.data.frame(summary(fit)$random[[random_param]]) %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(`95% CrI` = paste(round(`l-95% CI`, 2), "–", round(`u-95% CI`, 2))) %>%
    select(Parameter, Estimate, `Est.Error`, `95% CrI`)

  sigma <- as.data.frame(summary(fit)$spec_pars) %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(`95% CrI` = paste(round(`l-95% CI`, 2), "–", round(`u-95% CI`, 2))) %>%
    select(Parameter, Estimate, `Est.Error`, `95% CrI`, Rhat, Bulk_ESS, Tail_ESS)

  list(fixed = fixed, random = random, sigma = sigma)
}

plac_s <- extract_summary(fit_plac_exp)
cont_s <- extract_summary(fit_continuous)
cycl_s <- extract_summary(fit_cyclical)

create_table_grob <- function(df) tableGrob(df, rows = NULL, theme = ttheme_default(base_size = 10))

plac_section <- arrangeGrob(
  textGrob("Placebo Model (Exponential Asymptotic)", gp = gpar(fontface = "bold", fontsize = 14)),
  textGrob("Population-Level Effects", gp = gpar(fontface = "bold", fontsize = 12)),
  create_table_grob(plac_s$fixed),
  textGrob("Group-Level Effects", gp = gpar(fontface = "bold", fontsize = 12)),
  create_table_grob(plac_s$random),
  textGrob("Residual Variance", gp = gpar(fontface = "bold", fontsize = 12)),
  create_table_grob(plac_s$sigma),
  ncol = 1, heights = c(0.5, 0.3, 2, 0.3, 1.5, 0.3, 1))

cont_section <- arrangeGrob(
  textGrob("Continuous Regimens Model (Cumulative Dose)", gp = gpar(fontface = "bold", fontsize = 14)),
  textGrob("Population-Level Effects", gp = gpar(fontface = "bold", fontsize = 12)),
  create_table_grob(cont_s$fixed),
  textGrob("Group-Level Effects", gp = gpar(fontface = "bold", fontsize = 12)),
  create_table_grob(cont_s$random),
  textGrob("Residual Variance", gp = gpar(fontface = "bold", fontsize = 12)),
  create_table_grob(cont_s$sigma),
  ncol = 1, heights = c(0.5, 0.3, 2, 0.3, 1.5, 0.3, 1))

cycl_section <- arrangeGrob(
  textGrob("Cyclical/Weekly Regimens Model (On-Period Average Dose)", gp = gpar(fontface = "bold", fontsize = 14)),
  textGrob("Population-Level Effects", gp = gpar(fontface = "bold", fontsize = 12)),
  create_table_grob(cycl_s$fixed),
  textGrob("Group-Level Effects", gp = gpar(fontface = "bold", fontsize = 12)),
  create_table_grob(cycl_s$random),
  textGrob("Residual Variance", gp = gpar(fontface = "bold", fontsize = 12)),
  create_table_grob(cycl_s$sigma),
  ncol = 1, heights = c(0.5, 0.3, 2, 0.3, 1.5, 0.3, 1))

combined_mgadl_diag <- arrangeGrob(
  textGrob("Model Diagnostics Summary (MG-ADL)", gp = gpar(fontface = "bold", fontsize = 16)),
  plac_section, cont_section, cycl_section,
  ncol = 1, heights = c(0.5, 6, 6, 6))

pdf("mgadl_model_diagnostics_combined.pdf", width = 14, height = 28)
grid.draw(combined_mgadl_diag)
dev.off()


###############################################################################
# SECTION 14: QMG DATA LOADING & PREPARATION
###############################################################################

qmg_raw <- read_excel(data_path, sheet = "QMG_change")

qmg <- qmg_raw %>%
  mutate(
    study = as.factor(study),
    Trial = Trial_acronym,
    Drug = as.factor(ifelse(Drug == "PLACEBO", "Placebo", as.character(Drug))),
    Mechanism = as.factor(ifelse(is.na(Mechanism) | Mechanism == "", "Placebo", Mechanism)),
    Arm = as.factor(Arm),
    Regimen_type = as.factor(Regimen_type),
    time = as.numeric(time_weeks),
    change = qmg_change,
    se = as.numeric(se_qmg),
    dose_eq_perweek = as.numeric(Dose_eq_perweek),
    study_arm = as.factor(paste(study, Arm, sep = "_")),
    weight = 1 / se^2
  ) %>%
  filter(time > 0, !is.na(change), change != 0, !is.na(se), se > 0) %>%
  arrange(study_arm, time) %>%
  group_by(study_arm) %>%
  mutate(
    interval = time - lag(time, default = 0),
    cum_dose = cumsum(dose_eq_perweek * interval)
  ) %>% ungroup()

trt_qmg <- qmg %>% filter(Arm == "Treatment", Drug != "INEBILIZUMAB")

trt_qmg_continuous <- trt_qmg %>% filter(!str_detect(Regimen_type, "Cyclical|Weekly"))
trt_qmg_cyclical   <- trt_qmg %>% filter(str_detect(Regimen_type, "Cyclical|Weekly"))

trt_qmg_cycl_onperiod <- trt_qmg_cyclical %>%
  group_by(study_arm) %>%
  mutate(dose_eq_on = mean(dose_eq_perweek[dose_eq_perweek > 0], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(dose_eq = dose_eq_on)

plac_qmg <- qmg %>% filter(Arm == "Control", Trial != "MINT")
stopifnot(all(plac_qmg$se > 0))


###############################################################################
# SECTION 15: QMG EXPLORATORY PLOTS
###############################################################################

p_qmg_drug_vs_placebo <- ggplot(qmg, aes(x = time, y = change, color = Drug, linetype = Arm)) +
  geom_point(size = 3, alpha = 0.9) + geom_line(size = 1.2) +
  geom_errorbar(aes(ymin = change - 1.96*se, ymax = change + 1.96*se), width = 0.6, alpha = 0.8) +
  facet_wrap(~ Trial, scales = "free_y", ncol = 3) +
  scale_linetype_manual(values = c("Treatment" = "solid", "Control" = "dashed")) +
  scale_color_brewer(palette = "Paired") +
  labs(title = "QMG Trajectories: Each Drug vs Trial-Specific Placebo",
       subtitle = "Solid line = Treatment | Dashed line = Placebo | Points with 95% CI",
       x = "Weeks", y = "QMG Change from Baseline", color = "Drug", linetype = "Arm") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 12))

print(p_qmg_drug_vs_placebo)
ggsave("qmg_drug_vs_placebo_by_trial.png", p_qmg_drug_vs_placebo, width = 18, height = 14, dpi = 400, bg = "white")

p_qmg_mechanism <- ggplot(qmg %>% filter(Arm == "Treatment"), aes(x = time, y = change, color = Drug, group = study_arm)) +
  geom_point(size = 3, alpha = 0.9) + geom_line(size = 1.2) +
  geom_errorbar(aes(ymin = change - 1.96*se, ymax = change + 1.96*se), width = 0.6, alpha = 0.8) +
  facet_wrap(~ Mechanism, scales = "free_y", ncol = 2) +
  scale_color_brewer(palette = "Paired") +
  labs(title = "QMG Trajectories by Mechanism (Treatment Arms Only)",
       subtitle = "Points with 95% CI | Lines connect study arms",
       x = "Weeks", y = "QMG Change from Baseline", color = "Drug") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12))

print(p_qmg_mechanism)
ggsave("qmg_trajectories_by_mechanism.png", p_qmg_mechanism, width = 14, height = 10, dpi = 400, bg = "white")

p_qmg_regimen <- ggplot(qmg %>% filter(Arm == "Treatment"), aes(x = time, y = change, color = Drug, group = study_arm)) +
  geom_point(size = 3, alpha = 0.9) + geom_line(size = 1.2) +
  geom_errorbar(aes(ymin = change - 1.96*se, ymax = change + 1.96*se), width = 0.6, alpha = 0.8) +
  facet_wrap(~ Regimen_type, scales = "free_y", ncol = 2) +
  scale_color_brewer(palette = "Paired") +
  labs(title = "QMG Trajectories by Regimen Type (Treatment Arms Only)",
       subtitle = "Points with 95% CI | Lines connect study arms",
       x = "Weeks", y = "QMG Change from Baseline", color = "Drug") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12))

print(p_qmg_regimen)
ggsave("qmg_trajectories_by_regimen.png", p_qmg_regimen, width = 14, height = 10, dpi = 400, bg = "white")

combined_qmg_plot <- (p_qmg_drug_vs_placebo +
                        theme(plot.margin = margin(10,10,10,10), panel.border = element_rect(color = "grey70", fill = NA, size = 1)) +
                        ggtitle("Drug vs Trial-Specific Placebo")) +
  (p_qmg_mechanism +
     theme(plot.margin = margin(10,10,10,10), panel.border = element_rect(color = "grey70", fill = NA, size = 1)) +
     ggtitle("Mechanism-Specific Trajectories")) +
  (p_qmg_regimen +
     theme(plot.margin = margin(10,10,10,10), panel.border = element_rect(color = "grey70", fill = NA, size = 1)) +
     ggtitle("Regimen-Specific Trajectories")) +
  plot_layout(ncol = 1) +
  plot_annotation(title = "QMG Trajectory Comparisons",
                  subtitle = "Post-baseline data | Points with 95% CI | Grey borders separate panels",
                  caption = "Top: Drug vs Placebo by Trial | Middle: By Mechanism | Bottom: By Regimen Type") &
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5), text = element_text(size = 12))

print(combined_qmg_plot)
ggsave("qmg_trajectory_collage.png", combined_qmg_plot, width = 18, height = 24, dpi = 400, bg = "white")


###############################################################################
# SECTION 16: QMG MODEL FITTING
###############################################################################

cat("\n=== Fitting QMG Continuous Model ===\n")
fit_qmg_continuous <- brm(
  bf(
    change | se(se, sigma = TRUE) ~ -emax * (cum_dose / (ed50 + cum_dose)) * (1 - exp(-k * time)),
    emax ~ 1 + (1 | ID | Drug),
    ed50 ~ 1,
    k ~ 1,
    nl = TRUE
  ),
  data = trt_qmg_continuous,
  prior = c(
    prior(normal(8, 4), nlpar = "emax", lb = 0),
    prior(exponential(2), class = "sd", nlpar = "emax"),
    prior(normal(5000, 5000), nlpar = "ed50", lb = 0),
    prior(normal(0.1, 0.1), nlpar = "k", lb = 0)
  ),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4, iter = 12000, warmup = 4000,
  seed = 123
)

cat("\n=== Fitting QMG Cyclical Model ===\n")
fit_qmg_cyclical <- brm(
  bf(
    change | se(se, sigma = TRUE) ~ -emax * (dose_eq / (ed50 + dose_eq)) * (1 - exp(-k * time)),
    emax ~ 1 + (1 | ID | Drug),
    ed50 ~ 1,
    k ~ 1,
    nl = TRUE
  ),
  data = trt_qmg_cycl_onperiod,
  prior = c(
    prior(normal(10, 5), nlpar = "emax", lb = 0),
    prior(exponential(2), class = "sd", nlpar = "emax"),
    prior(normal(200, 400), nlpar = "ed50", lb = 0),
    prior(normal(0.4, 0.2), nlpar = "k", lb = 0)
  ),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4, iter = 12000, warmup = 4000,
  seed = 123
)

cat("\n=== Fitting QMG Placebo Model ===\n")
fit_plac_qmg <- brm(
  bf(
    change | se(se, sigma = TRUE) ~ plateau * (1 - exp(-k * time)),
    plateau ~ 1 + (1 | ID | study),
    k ~ 1 + (1 | ID | study),
    nl = TRUE
  ),
  data = plac_qmg,
  prior = c(
    prior(normal(-3, 1), nlpar = "plateau", lb = -10),
    prior(normal(0.1, 0.1), nlpar = "k", lb = 0),
    prior(exponential(2), class = "sd", nlpar = "plateau"),
    prior(exponential(2), class = "sd", nlpar = "k")
  ),
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  chains = 4, iter = 20000, warmup = 5000,
  seed = 123
)


###############################################################################
# SECTION 17: QMG VPCs
###############################################################################

trt_pred_qmg_cont <- trt_qmg_continuous %>%
  add_predicted_draws(fit_qmg_continuous, ndraws = 100, allow_new_levels = TRUE)
trt_pred_summary_qmg_cont <- trt_pred_qmg_cont %>% group_by(Drug, time) %>% median_hdci(.prediction)

p_vpc_qmg_cont <- ggplot() +
  geom_errorbar(data = trt_qmg_continuous, aes(x = time, ymin = change - 1.96*se, ymax = change + 1.96*se, color = Drug), width = 0.5, alpha = 0.8) +
  geom_point(data = trt_qmg_continuous, aes(x = time, y = change, color = Drug), size = 3) +
  geom_line(data = trt_pred_summary_qmg_cont, aes(x = time, y = .prediction, color = Drug), size = 1.2) +
  geom_ribbon(data = trt_pred_summary_qmg_cont, aes(x = time, y = .prediction, ymin = .lower, ymax = .upper, fill = Drug), alpha = 0.3) +
  facet_wrap(~ Drug, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") +
  labs(title = "VPC: QMG Continuous Regimens (Cumulative Dose Model)",
       subtitle = "Observed points with 95% CI | Predicted median line + 95% prediction interval ribbon",
       x = "Weeks", y = "QMG Change from Baseline") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 12, color = "darkblue"),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12), panel.grid.minor = element_blank())

print(p_vpc_qmg_cont)
ggsave("vpc_qmg_continuous_regimens.png", p_vpc_qmg_cont, width = 16, height = 10, dpi = 400, bg = "white")

trt_pred_qmg_cycl <- trt_qmg_cycl_onperiod %>%
  add_predicted_draws(fit_qmg_cyclical, ndraws = 100, allow_new_levels = TRUE)
trt_pred_summary_qmg_cycl <- trt_pred_qmg_cycl %>% group_by(Drug, time) %>% median_hdci(.prediction)

p_vpc_qmg_cycl <- ggplot() +
  geom_errorbar(data = trt_qmg_cycl_onperiod, aes(x = time, ymin = change - 1.96*se, ymax = change + 1.96*se, color = Drug), width = 0.5, alpha = 0.8) +
  geom_point(data = trt_qmg_cycl_onperiod, aes(x = time, y = change, color = Drug), size = 3) +
  geom_line(data = trt_pred_summary_qmg_cycl, aes(x = time, y = .prediction, color = Drug), size = 1.2) +
  geom_ribbon(data = trt_pred_summary_qmg_cycl, aes(x = time, y = .prediction, ymin = .lower, ymax = .upper, fill = Drug), alpha = 0.3) +
  facet_wrap(~ Drug, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") +
  labs(title = "VPC: QMG Cyclical/Weekly Regimens (On-Period Average Dose Model)",
       subtitle = "Observed points with 95% CI | Predicted median line + 95% prediction interval ribbon",
       x = "Weeks", y = "QMG Change from Baseline") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 12, color = "darkblue"),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12), panel.grid.minor = element_blank())

print(p_vpc_qmg_cycl)
ggsave("vpc_qmg_cyclical_regimens.png", p_vpc_qmg_cycl, width = 16, height = 8, dpi = 400, bg = "white")

p_combined_vpc_qmg <- p_vpc_qmg_cont / p_vpc_qmg_cycl +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(title = "Visual Predictive Checks – QMG Models",
                  subtitle = "Continuous Regimens (Cumulative Dose) vs Cyclical/Weekly Regimens (On-Period Average Dose)",
                  caption = "Observed points with 95% CI | Predicted median line + 95% prediction interval ribbon",
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
                                plot.subtitle = element_text(size = 14, hjust = 0.5),
                                plot.caption = element_text(size = 11, hjust = 1, colour = "gray50")))

print(p_combined_vpc_qmg)
ggsave("qmg_vpc_combined.png", p_combined_vpc_qmg, width = 16, height = 12, dpi = 400, bg = "white")


###############################################################################
# SECTION 18: QMG FITTED TRAJECTORIES
###############################################################################

trt_fitted_qmg_cont <- trt_qmg_continuous %>%
  add_epred_draws(fit_qmg_continuous, ndraws = 100, allow_new_levels = TRUE)
trt_fitted_summary_qmg_cont <- trt_fitted_qmg_cont %>% mean_hdci(.epred)

p_fitted_qmg_cont <- ggplot() +
  geom_errorbar(data = trt_qmg_continuous, aes(x = time, ymin = change - 1.96*se, ymax = change + 1.96*se, color = Drug), width = 0.5, alpha = 0.8) +
  geom_point(data = trt_qmg_continuous, aes(x = time, y = change, color = Drug), size = 3) +
  geom_line(data = trt_fitted_summary_qmg_cont, aes(x = time, y = .epred, color = Drug), size = 1.2) +
  geom_ribbon(data = trt_fitted_summary_qmg_cont, aes(x = time, y = .epred, ymin = .lower, ymax = .upper, fill = Drug), alpha = 0.3) +
  facet_wrap(~ Drug, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") +
  labs(title = "Fitted Trajectories: QMG Continuous Regimens (Cumulative Dose Model)",
       subtitle = "Observed points with 95% CI | Predicted mean line + 95% HDCI ribbon",
       x = "Weeks", y = "QMG Change from Baseline") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 12, color = "darkblue"),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12), panel.grid.minor = element_blank())

print(p_fitted_qmg_cont)
ggsave("fitted_trajectories_qmg_continuous.png", p_fitted_qmg_cont, width = 16, height = 10, dpi = 400, bg = "white")

trt_fitted_qmg_cycl <- trt_qmg_cycl_onperiod %>%
  add_epred_draws(fit_qmg_cyclical, ndraws = 100, allow_new_levels = TRUE)
trt_fitted_summary_qmg_cycl <- trt_fitted_qmg_cycl %>% mean_hdci(.epred)

p_fitted_qmg_cycl <- ggplot() +
  geom_errorbar(data = trt_qmg_cycl_onperiod, aes(x = time, ymin = change - 1.96*se, ymax = change + 1.96*se, color = Drug), width = 0.5, alpha = 0.8) +
  geom_point(data = trt_qmg_cycl_onperiod, aes(x = time, y = change, color = Drug), size = 3) +
  geom_line(data = trt_fitted_summary_qmg_cycl, aes(x = time, y = .epred, color = Drug), size = 1.2) +
  geom_ribbon(data = trt_fitted_summary_qmg_cycl, aes(x = time, y = .epred, ymin = .lower, ymax = .upper, fill = Drug), alpha = 0.3) +
  facet_wrap(~ Drug, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") +
  labs(title = "Fitted Trajectories: QMG Cyclical/Weekly Regimens (On-Period Average Dose Model)",
       subtitle = "Observed points with 95% CI | Predicted mean line + 95% HDCI ribbon",
       x = "Weeks", y = "QMG Change from Baseline") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 12, color = "darkblue"),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12), panel.grid.minor = element_blank())

print(p_fitted_qmg_cycl)
ggsave("fitted_trajectories_qmg_cyclical.png", p_fitted_qmg_cycl, width = 16, height = 8, dpi = 400, bg = "white")


###############################################################################
# SECTION 19: QMG EMAX + NET BENEFIT + WEEK-26 + TIME TO MCID + AUEC
###############################################################################

# QMG Emax
emax_qmg_cont <- as_draws_df(fit_qmg_continuous) %>%
  select(.draw, starts_with("r_Drug"), b_emax_Intercept) %>%
  pivot_longer(cols = starts_with("r_Drug"), names_to = "temp", values_to = "offset") %>%
  mutate(Drug = str_extract(temp, "(?<=\\[)[^,]+"), emax = b_emax_Intercept + offset, regimen = "Continuous") %>%
  filter(!is.na(Drug))

emax_qmg_cycl <- as_draws_df(fit_qmg_cyclical) %>%
  select(.draw, starts_with("r_Drug"), b_emax_Intercept) %>%
  pivot_longer(cols = starts_with("r_Drug"), names_to = "temp", values_to = "offset") %>%
  mutate(Drug = str_extract(temp, "(?<=\\[)[^,]+"), emax = b_emax_Intercept + offset, regimen = "Cyclical/Weekly") %>%
  filter(!is.na(Drug))

emax_qmg_combined <- bind_rows(emax_qmg_cont, emax_qmg_cycl)

emax_qmg_summary <- emax_qmg_combined %>%
  group_by(Drug, regimen) %>%
  mean_qi(emax, .width = 0.95) %>%
  mutate(label = sprintf("%.1f [%.1f, %.1f]", emax, .lower, .upper)) %>%
  arrange(desc(emax))

label_x_pos_qmg <- max(emax_qmg_summary$.upper) + 1

p_emax_qmg_combined <- ggplot(emax_qmg_combined, aes(x = emax, color = regimen)) +
  stat_halfeye(aes(y = reorder(Drug, emax, mean)), .width = c(0.66, 0.95), slab_alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray50") +
  geom_text(data = emax_qmg_summary, aes(y = reorder(Drug, emax, mean), label = label, x = label_x_pos_qmg),
            hjust = "left", size = 4.5, fontface = "bold") +
  labs(title = "Drug-Specific Posterior Emax for QMG (Distinct Models by Regimen)",
       subtitle = "Higher = greater long-term efficacy | Continuous vs Cyclical/Weekly",
       x = "Emax (points improvement)", y = NULL, color = "Regimen Type") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.margin = margin(10, 140, 10, 10))

print(p_emax_qmg_combined)
ggsave("emax_qmg_per_drug_distinct_models_forest.png", p_emax_qmg_combined,
       width = 14, height = max(6, 0.8 * n_distinct(emax_qmg_combined$Drug)), dpi = 400, bg = "white")

# QMG Net added benefit
placebo_plateau_draws_qmg <- as_draws_df(fit_plac_qmg) %>%
  select(.draw, b_plateau_Intercept) %>%
  mutate(placebo_mag = abs(b_plateau_Intercept))
placebo_mag_mean_qmg <- mean(placebo_plateau_draws_qmg$placebo_mag)

emax_qmg_net_all <- bind_rows(
  emax_qmg_cont %>% mutate(net_emax = emax - placebo_mag_mean_qmg),
  emax_qmg_cycl %>% mutate(net_emax = emax - placebo_mag_mean_qmg)
)

emax_qmg_net_summary <- emax_qmg_net_all %>%
  group_by(Drug, regimen) %>% mean_qi(net_emax, .width = 0.95) %>%
  mutate(label = sprintf("%.1f [%.1f, %.1f]", net_emax, .lower, .upper)) %>%
  arrange(desc(net_emax))

label_x_pos_qmg_net <- max(emax_qmg_net_summary$.upper) + 1

p_net_emax_qmg <- ggplot(emax_qmg_net_all, aes(x = net_emax, color = regimen)) +
  stat_halfeye(aes(y = reorder(Drug, net_emax, mean)), .width = c(0.66, 0.95), slab_alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray50") +
  geom_text(data = emax_qmg_net_summary, aes(y = reorder(Drug, net_emax, mean), label = label, x = label_x_pos_qmg_net),
            hjust = "left", size = 4.5, fontface = "bold") +
  labs(title = "Drug-Specific Net Added Improvement vs Placebo for QMG",
       subtitle = "Net = Drug Emax magnitude - Placebo plateau magnitude | Higher = greater added benefit",
       x = "Net Added Improvement over Placebo (points)", y = NULL, color = "Regimen Type") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.margin = margin(10, 140, 10, 10))

print(p_net_emax_qmg)
ggsave("net_added_emax_qmg_vs_placebo_forest.png", p_net_emax_qmg, width = 14, height = 10, dpi = 400, bg = "white")

# QMG Week-26 with superiority
placebo_draws_qmg <- posterior_predict(fit_plac_qmg, newdata = plac_qmg %>% mutate(time = 26), ndraws = 1000, allow_new_levels = TRUE) %>% as.vector()

pred_week26_qmg_cont <- trt_qmg_continuous %>%
  group_by(Drug, Mechanism) %>%
  summarise(dose_eq_mean = mean(dose_eq_perweek, na.rm = TRUE), .groups = "drop") %>%
  mutate(time = 26, cum_dose = dose_eq_mean * 26, dose_eq = dose_eq_mean, se = 0)

pred_qmg_cont_matrix <- posterior_predict(fit_qmg_continuous, newdata = pred_week26_qmg_cont, ndraws = 1000, allow_new_levels = TRUE)

pred_qmg_cont_post <- pred_qmg_cont_matrix %>%
  as.data.frame() %>% mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "obs_str", values_to = "drug_change") %>%
  mutate(obs = as.integer(str_remove(obs_str, "^V"))) %>% select(-obs_str) %>%
  left_join(pred_week26_qmg_cont %>% mutate(obs = row_number()), by = "obs") %>%
  mutate(regimen = "Continuous")

pred_week26_qmg_cycl <- trt_qmg_cycl_onperiod %>%
  group_by(Drug, Mechanism) %>%
  summarise(dose_eq = mean(dose_eq, na.rm = TRUE), .groups = "drop") %>%
  mutate(time = 26, se = 0)

pred_qmg_cycl_matrix <- posterior_predict(fit_qmg_cyclical, newdata = pred_week26_qmg_cycl, ndraws = 1000, allow_new_levels = TRUE)

pred_qmg_cycl_post <- pred_qmg_cycl_matrix %>%
  as.data.frame() %>% mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "obs_str", values_to = "drug_change") %>%
  mutate(obs = as.integer(str_remove(obs_str, "^V"))) %>% select(-obs_str) %>%
  left_join(pred_week26_qmg_cycl %>% mutate(obs = row_number()), by = "obs") %>%
  mutate(regimen = "Cyclical/Weekly")

pred_26_qmg_post <- bind_rows(pred_qmg_cont_post, pred_qmg_cycl_post)
pred_26_qmg_net_post <- pred_26_qmg_post %>% mutate(net_change = drug_change - placebo_draws_qmg[draw])

superior_prob_qmg <- pred_26_qmg_net_post %>%
  group_by(Drug, Mechanism, regimen) %>%
  summarise(prob_superior = mean(net_change < 0) * 100, .groups = "drop")

pred_26_qmg_net_summary <- pred_26_qmg_net_post %>%
  group_by(Drug, Mechanism, regimen) %>%
  mean_qi(net_change, .width = 0.95) %>%
  left_join(superior_prob_qmg, by = c("Drug", "Mechanism", "regimen")) %>%
  mutate(label = sprintf("%.1f [%.1f, %.1f] (P=%.1f%%)", net_change, .lower, .upper, prob_superior)) %>%
  arrange(net_change)

p_net_superior_qmg <- ggplot(pred_26_qmg_net_post, aes(x = net_change, color = regimen)) +
  stat_halfeye(aes(y = reorder(Drug, net_change, mean)), .width = c(0.66, 0.95), slab_alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray5") +
  geom_text(data = pred_26_qmg_net_summary,
            aes(y = reorder(Drug, net_change, mean), label = label, x = min(.lower) - 0.5),
            hjust = "right", size = 4.5, fontface = "bold") +
  labs(title = "Net Added Benefit vs Placebo at Week 26 with Superiority Probability",
       x = "Net ΔQMG at Week 26 (points)", y = NULL, color = "Regimen Type") +
  scale_color_manual(values = c("Continuous" = "#1B9E77", "Cyclical/Weekly" = "#D95F02")) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 13, hjust = 0.5),
        legend.position = "bottom", legend.title = element_text(face = "bold"),
        axis.text = element_text(size = 8),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        plot.margin = margin(20, 30, 20, 20))

print(p_net_superior_qmg)
ggsave("qmg_net_benefit_week26.png", p_net_superior_qmg, width = 14, height = 10, dpi = 400)

# QMG Time to MCID
MCID_QMG <- -3
time_grid_qmg <- seq(1, 26, by = 0.5)

pred_time_grid_qmg_cont <- trt_qmg_continuous %>%
  group_by(Drug) %>% summarise(dose_eq_mean = mean(dose_eq_perweek, na.rm = TRUE), .groups = "drop") %>%
  crossing(time = time_grid_qmg) %>%
  mutate(cum_dose = dose_eq_mean * time, dose_eq = dose_eq_mean, se = 0)

pred_time_qmg_cont <- posterior_predict(fit_qmg_continuous, newdata = pred_time_grid_qmg_cont, ndraws = 500, allow_new_levels = TRUE)

pred_time_qmg_cont_post <- pred_time_qmg_cont %>%
  as.data.frame() %>% mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "obs_str", values_to = "pred_change") %>%
  mutate(obs = as.integer(str_remove(obs_str, "^V"))) %>% select(-obs_str) %>%
  left_join(pred_time_grid_qmg_cont %>% mutate(obs = row_number()), by = "obs") %>%
  mutate(regimen = "Continuous")

pred_time_grid_qmg_cycl <- trt_qmg_cycl_onperiod %>%
  group_by(Drug) %>% summarise(dose_eq = mean(dose_eq, na.rm = TRUE), .groups = "drop") %>%
  crossing(time = time_grid_qmg) %>%
  mutate(se = 0)

pred_time_qmg_cycl <- posterior_predict(fit_qmg_cyclical, newdata = pred_time_grid_qmg_cycl, ndraws = 500, allow_new_levels = TRUE)

pred_time_qmg_cycl_post <- pred_time_qmg_cycl %>%
  as.data.frame() %>% mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "obs_str", values_to = "pred_change") %>%
  mutate(obs = as.integer(str_remove(obs_str, "^V"))) %>% select(-obs_str) %>%
  left_join(pred_time_grid_qmg_cycl %>% mutate(obs = row_number()), by = "obs") %>%
  mutate(regimen = "Cyclical/Weekly")

pred_time_qmg_all_post <- bind_rows(pred_time_qmg_cont_post, pred_time_qmg_cycl_post)

time_to_mcid_qmg <- pred_time_qmg_all_post %>%
  group_by(draw, Drug, regimen) %>% arrange(time) %>%
  filter(pred_change <= MCID_QMG) %>% slice(1) %>% ungroup()

time_mcid_qmg_summary <- time_to_mcid_qmg %>%
  group_by(Drug, regimen) %>% mean_qi(time, .width = 0.95) %>%
  mutate(label = sprintf("%.1f [%.1f, %.1f] weeks", time, .lower, .upper)) %>% arrange(time)

p_time_mcid_qmg <- ggplot(time_to_mcid_qmg, aes(x = time, color = regimen)) +
  stat_halfeye(aes(y = reorder(Drug, time, mean)), .width = c(0.66, 0.95), slab_alpha = 0.8) +
  geom_text(data = time_mcid_qmg_summary,
            aes(y = reorder(Drug, time, mean), label = label),
            x = Inf, hjust = -0.05,
            size = 4.5, fontface = "bold") +
  coord_cartesian(clip = "off") +
  labs(title = "Drug-Specific Time to Achieve MCID (≥3-Point QMG Improvement)",
       subtitle = "Posterior median [95% CrI] weeks | Shorter = faster meaningful benefit",
       x = "Time to MCID (weeks)", y = NULL, color = "Regimen Type") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "bottom", plot.margin = margin(10, 160, 10, 10))

print(p_time_mcid_qmg)
ggsave("time_to_mcid_qmg.png", p_time_mcid_qmg, width = 16, height = 10, dpi = 400, bg = "white")

# QMG AUEC
time_grid_auec_qmg <- seq(0, 26, by = 1)

pred_auec_qmg_cont <- trt_qmg_continuous %>%
  group_by(Drug) %>% summarise(dose_eq_mean = mean(dose_eq_perweek, na.rm = TRUE), .groups = "drop") %>%
  crossing(time = time_grid_auec_qmg) %>%
  mutate(cum_dose = dose_eq_mean * time, dose_eq = dose_eq_mean, se = 0)

pred_auec_qmg_cont_matrix <- posterior_predict(fit_qmg_continuous, newdata = pred_auec_qmg_cont, ndraws = 500, allow_new_levels = TRUE)

pred_auec_qmg_cont_post <- pred_auec_qmg_cont_matrix %>%
  as.data.frame() %>% mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "obs_str", values_to = "pred_change") %>%
  mutate(obs = as.integer(str_remove(obs_str, "^V"))) %>% select(-obs_str) %>%
  left_join(pred_auec_qmg_cont %>% mutate(obs = row_number()), by = "obs") %>%
  mutate(regimen = "Continuous")

pred_auec_qmg_cycl <- trt_qmg_cycl_onperiod %>%
  group_by(Drug) %>% summarise(dose_eq = mean(dose_eq, na.rm = TRUE), .groups = "drop") %>%
  crossing(time = time_grid_auec_qmg) %>%
  mutate(se = 0)

pred_auec_qmg_cycl_matrix <- posterior_predict(fit_qmg_cyclical, newdata = pred_auec_qmg_cycl, ndraws = 500, allow_new_levels = TRUE)

pred_auec_qmg_cycl_post <- pred_auec_qmg_cycl_matrix %>%
  as.data.frame() %>% mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "obs_str", values_to = "pred_change") %>%
  mutate(obs = as.integer(str_remove(obs_str, "^V"))) %>% select(-obs_str) %>%
  left_join(pred_auec_qmg_cycl %>% mutate(obs = row_number()), by = "obs") %>%
  mutate(regimen = "Cyclical/Weekly")

pred_auec_qmg_all_post <- bind_rows(pred_auec_qmg_cont_post, pred_auec_qmg_cycl_post)

auec_qmg_all <- pred_auec_qmg_all_post %>%
  group_by(draw, Drug, regimen) %>% arrange(time) %>%
  summarise(auec = trapz(time, abs(pred_change)), .groups = "drop")

auec_qmg_summary <- auec_qmg_all %>%
  group_by(Drug, regimen) %>% mean_qi(auec, .width = 0.95) %>%
  mutate(label = sprintf("%.0f [%.0f, %.0f]", auec, .lower, .upper)) %>% arrange(desc(auec))

label_x_pos_auec_qmg <- max(auec_qmg_summary$.upper) + 50

p_auec_qmg <- ggplot(auec_qmg_all, aes(x = auec, color = regimen)) +
  stat_halfeye(aes(y = reorder(Drug, auec, mean)), .width = c(0.66, 0.95), slab_alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray50") +
  geom_text(data = auec_qmg_summary, aes(y = reorder(Drug, auec, mean), label = label, x = label_x_pos_auec_qmg),
            hjust = "right", size = 4.5, fontface = "bold") +
  labs(title = "Drug-Specific Area Under the Effect Curve (AUEC, Week 0-26) for QMG",
       subtitle = "Total integrated improvement | Higher = greater cumulative benefit",
       x = "AUEC (point·weeks)", y = NULL, color = "Regimen Type") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.margin = margin(10, 140, 10, 10))

print(p_auec_qmg)
ggsave("auec_qmg_per_drug_distinct_models.png", p_auec_qmg, width = 14, height = 10, dpi = 400, bg = "white")

# QMG dashboard
p_combined_qmg_dashboard <- p_time_mcid_qmg / p_auec_qmg +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(
    title = "Clinical Benefit Dashboard – QMG Outcome",
    subtitle = "Time to meaningful improvement + Total cumulative benefit over 26 weeks",
    caption = "Posterior median [95% CrI] | Colour by regimen type",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 14, hjust = 0.5),
                  plot.caption = element_text(size = 11, hjust = 1, colour = "gray50")))

print(p_combined_qmg_dashboard)
ggsave("qmg_clinical_benefit_dashboard_combined.png", p_combined_qmg_dashboard,
       width = 14, height = 12, dpi = 400, bg = "white")

###############################################################################
# SECTION 20: QMG POOLED ED50 AND k
###############################################################################

ed50_k_qmg_continuous <- as_draws_df(fit_qmg_continuous) %>%
  select(.draw, b_ed50_Intercept, b_k_Intercept) %>%
  rename(ed50 = b_ed50_Intercept, k = b_k_Intercept) %>%
  mutate(regimen = "Continuous")

ed50_k_qmg_cyclical <- as_draws_df(fit_qmg_cyclical) %>%
  select(.draw, b_ed50_Intercept, b_k_Intercept) %>%
  rename(ed50 = b_ed50_Intercept, k = b_k_Intercept) %>%
  mutate(regimen = "Cyclical/Weekly")

ed50_k_qmg_combined <- bind_rows(ed50_k_qmg_continuous, ed50_k_qmg_cyclical)

ed50_k_qmg_summary <- ed50_k_qmg_combined %>%
  group_by(regimen) %>% mean_qi(ed50, k, .width = 0.95) %>%
  mutate(ED50 = sprintf("%.0f [%.0f, %.0f]", ed50, ed50.lower, ed50.upper),
         `Onset Rate k` = sprintf("%.2f [%.2f, %.2f]", k, k.lower, k.upper)) %>%
  select(regimen, ED50, `Onset Rate k`) %>% arrange(regimen)

print("Pooled ED50 and k by Regimen (QMG)")
print(ed50_k_qmg_summary)
write.csv(ed50_k_qmg_summary, "pooled_ed50_k_qmg_by_regimen.csv", row.names = FALSE)

p_ed50_qmg <- ggplot(ed50_k_qmg_combined, aes(x = ed50, fill = regimen)) +
  stat_halfeye(.width = c(0.66, 0.95), slab_alpha = 0.8) +
  labs(title = "Pooled ED50 (QMG)", x = "ED50 (mg/week equivalent)", y = "Density", fill = "Regimen") +
  theme_minimal() + theme(legend.position = "none")

p_k_qmg <- ggplot(ed50_k_qmg_combined, aes(x = k, fill = regimen)) +
  stat_halfeye(.width = c(0.66, 0.95), slab_alpha = 0.8) +
  labs(title = "Pooled Onset Rate k (QMG)", x = "k (per week)", y = "Density", fill = "Regimen") +
  theme_minimal() + theme(legend.position = "bottom")

plots_combined_qmg <- p_ed50_qmg + p_k_qmg + plot_layout(guides = "collect") & theme(legend.position = "bottom")

table_grob_qmg <- tableGrob(ed50_k_qmg_summary, rows = NULL,
                              theme = ttheme_default(base_size = 12, core = list(fg_params = list(fontface = "bold"))))

final_figure_qmg <- plots_combined_qmg / table_grob_qmg +
  plot_annotation(title = "Pooled Pharmacodynamic Parameters by Regimen Type (QMG)",
                  subtitle = "ED50 (lower = higher potency) | k (higher = faster onset)",
                  caption = "Posterior mean [95% CrI]") &
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5))

print(final_figure_qmg)
ggsave("pooled_ed50_k_qmg_by_regimen_combined.png", final_figure_qmg,
       width = 16, height = 12, dpi = 400, bg = "white")


###############################################################################
# SECTION 21: QMG MODEL DIAGNOSTICS PDF
###############################################################################

plac_qmg_s <- extract_summary(fit_plac_qmg)
cont_qmg_s <- extract_summary(fit_qmg_continuous)
cycl_qmg_s <- extract_summary(fit_qmg_cyclical)

create_section_grob <- function(summary_list, title_text) {
  arrangeGrob(
    textGrob(title_text, gp = gpar(fontface = "bold", fontsize = 14)),
    textGrob("Population-Level Effects", gp = gpar(fontface = "bold", fontsize = 12)),
    create_table_grob(summary_list$fixed),
    textGrob("Group-Level Effects", gp = gpar(fontface = "bold", fontsize = 12)),
    create_table_grob(summary_list$random),
    textGrob("Residual Variance (sigma)", gp = gpar(fontface = "bold", fontsize = 12)),
    create_table_grob(summary_list$sigma),
    ncol = 1, heights = c(0.5, 0.3, 2.2, 0.3, 1.8, 0.3, 1.2))
}

combined_qmg_pdf <- arrangeGrob(
  textGrob("QMG Model Diagnostics Summary", gp = gpar(fontface = "bold", fontsize = 16)),
  create_section_grob(plac_qmg_s, "Placebo Model (QMG)"),
  create_section_grob(cont_qmg_s, "Continuous Regimens Model (QMG)"),
  create_section_grob(cycl_qmg_s, "Cyclical/Weekly Regimens Model (QMG)"),
  ncol = 1, heights = c(0.6, 6, 6, 6))

ggsave("qmg_model_diagnostics_combined.pdf", combined_qmg_pdf,
       width = 14, height = 28, dpi = 400, device = "pdf")


###############################################################################
# SECTION 22: SAVE ALL OBJECTS BEFORE MIDD MODULES
###############################################################################

cat("\n=== Saving all fitted model objects ===\n")
save(fit_continuous, fit_cyclical, fit_plac_exp,
     fit_qmg_continuous, fit_qmg_cyclical, fit_plac_qmg,
     trt_continuous, trt_cycl_onperiod,
     trt_qmg_continuous, trt_qmg_cycl_onperiod,
     plac, plac_qmg, data, qmg,
     time_to_mcid, time_to_mcid_qmg,
     auec_all, auec_qmg_all,
     emax_cont, emax_cycl, emax_combined,
     emax_qmg_cont, emax_qmg_cycl, emax_qmg_combined,
     file = file.path(output_dir, "all_fitted_models.RData"))
cat("Models saved to:", file.path(output_dir, "all_fitted_models.RData"), "\n")


###############################################################################
# SECTION 23: MIDD MODULE 1 – CLINICAL TRIAL SIMULATION
###############################################################################

cat("\n=== MIDD Module 1: Endpoint Sensitivity Simulation ===\n")

n_sim_trials    <- 500
n_per_arm       <- 60
time_grid_cts   <- seq(0, 26, by = 1)
cyclical_drugs  <- c("EFGARTIGIMOD", "BATOCLIMAB", "ROZANOLIXIZUMAB")
assessment_weeks <- c(12, 16, 20, 24, 26)

extract_drug_emax_cts <- function(fit, drug_name) {
  draws <- as_draws_df(fit)
  intercept <- draws$b_emax_Intercept
  drug_col <- grep(paste0("r_Drug__emax\\[", drug_name), names(draws), value = TRUE)
  if (length(drug_col) == 0) return(intercept)
  return(intercept + draws[[drug_col[1]]])
}

extract_pop_params <- function(fit) {
  draws <- as_draws_df(fit)
  list(ed50 = draws$b_ed50_Intercept, k = draws$b_k_Intercept, sigma = draws$sigma)
}

extract_placebo_params <- function(fit_plac) {
  draws <- as_draws_df(fit_plac)
  list(plateau = draws$b_plateau_Intercept, k = draws$b_k_Intercept,
       sigma = if ("sigma" %in% names(draws)) draws$sigma else rep(1, nrow(draws)))
}

params_mgadl_cycl <- extract_pop_params(fit_cyclical)
params_mgadl_plac <- extract_placebo_params(fit_plac_exp)
params_qmg_cycl   <- extract_pop_params(fit_qmg_cyclical)
params_qmg_plac   <- extract_placebo_params(fit_plac_qmg)

simulate_trial <- function(drug_name, fit_drug, fit_plac, params_drug, params_plac,
                            n_per_arm, time_grid_sim, scale_name = "MG-ADL") {
  idx <- sample(1:length(params_drug$ed50), 1)
  emax_true  <- extract_drug_emax_cts(fit_drug, drug_name)[idx]
  ed50_true  <- params_drug$ed50[idx]
  k_true     <- params_drug$k[idx]
  sigma_drug <- params_drug$sigma[idx]
  plac_plateau <- params_plac$plateau[idx]
  plac_k       <- params_plac$k[idx]
  sigma_plac   <- params_plac$sigma[min(idx, length(params_plac$sigma))]

  dose_eq <- mean(trt_cycl_onperiod$dose_eq[trt_cycl_onperiod$Drug == drug_name], na.rm = TRUE)
  if (is.na(dose_eq) || dose_eq == 0)
    dose_eq <- mean(trt_qmg_cycl_onperiod$dose_eq[trt_qmg_cycl_onperiod$Drug == drug_name], na.rm = TRUE)

  drug_effect <- -emax_true * (dose_eq / (ed50_true + dose_eq)) * (1 - exp(-k_true * time_grid_sim))
  plac_effect <- plac_plateau * (1 - exp(-plac_k * time_grid_sim))

  trt_data <- matrix(NA, nrow = n_per_arm, ncol = length(time_grid_sim))
  plac_data <- matrix(NA, nrow = n_per_arm, ncol = length(time_grid_sim))
  for (i in 1:n_per_arm) {
    trt_data[i, ] <- drug_effect + rnorm(length(time_grid_sim), 0, abs(sigma_drug))
    plac_data[i, ] <- plac_effect + rnorm(length(time_grid_sim), 0, abs(sigma_plac))
  }

  endpoint_results <- list()
  for (wk in assessment_weeks) {
    wk_idx <- which(time_grid_sim == wk)
    if (length(wk_idx) == 0) next
    tt <- t.test(trt_data[, wk_idx], plac_data[, wk_idx], alternative = "less")
    endpoint_results[[paste0("week", wk)]] <- data.frame(
      drug = drug_name, scale = scale_name,
      endpoint = paste0("Week ", wk, " Change"),
      p_value = tt$p.value,
      effect_size = mean(trt_data[, wk_idx]) - mean(plac_data[, wk_idx]),
      stringsAsFactors = FALSE)
  }

  trt_auec <- apply(trt_data, 1, function(row) trapz(time_grid_sim, abs(row)))
  plac_auec <- apply(plac_data, 1, function(row) trapz(time_grid_sim, abs(row)))
  tt_auec <- t.test(trt_auec, plac_auec, alternative = "greater")
  endpoint_results[["auec"]] <- data.frame(
    drug = drug_name, scale = scale_name,
    endpoint = "AUEC (Week 0-26)", p_value = tt_auec$p.value,
    effect_size = mean(trt_auec) - mean(plac_auec), stringsAsFactors = FALSE)

  trt_mean_change <- rowMeans(trt_data)
  plac_mean_change <- rowMeans(plac_data)
  tt_avg <- t.test(trt_mean_change, plac_mean_change, alternative = "less")
  endpoint_results[["time_avg"]] <- data.frame(
    drug = drug_name, scale = scale_name,
    endpoint = "Time-Averaged Change", p_value = tt_avg$p.value,
    effect_size = mean(trt_mean_change) - mean(plac_mean_change), stringsAsFactors = FALSE)

  bind_rows(endpoint_results)
}

all_sim_results <- list()
counter <- 0
for (drug in cyclical_drugs) {
  cat("  Simulating:", drug, "\n")
  for (sim in 1:n_sim_trials) {
    tryCatch({
      res_mgadl <- simulate_trial(drug, fit_cyclical, fit_plac_exp, params_mgadl_cycl, params_mgadl_plac, n_per_arm, time_grid_cts, "MG-ADL")
      counter <- counter + 1; all_sim_results[[counter]] <- res_mgadl
    }, error = function(e) NULL)
    tryCatch({
      res_qmg <- simulate_trial(drug, fit_qmg_cyclical, fit_plac_qmg, params_qmg_cycl, params_qmg_plac, n_per_arm, time_grid_cts, "QMG")
      counter <- counter + 1; all_sim_results[[counter]] <- res_qmg
    }, error = function(e) NULL)
  }
}

sim_df <- bind_rows(all_sim_results)

alpha <- 0.05
power_summary <- sim_df %>%
  group_by(drug, scale, endpoint) %>%
  summarise(n_sims = n(), power = mean(p_value < alpha, na.rm = TRUE) * 100,
            mean_effect = mean(effect_size, na.rm = TRUE), .groups = "drop") %>%
  arrange(drug, scale, desc(power))

cat("\n=== STATISTICAL POWER BY ENDPOINT ===\n")
print(power_summary, n = Inf)
write.csv(power_summary, "midd_endpoint_sensitivity_power.csv", row.names = FALSE)

# Power heatmap
key_endpoints <- c("Week 12 Change", "Week 16 Change", "Week 20 Change",
                   "Week 24 Change", "Week 26 Change", "AUEC (Week 0-26)", "Time-Averaged Change")

power_plot_data <- power_summary %>%
  filter(endpoint %in% key_endpoints) %>%
  mutate(endpoint = factor(endpoint, levels = key_endpoints),
         drug = factor(drug, levels = cyclical_drugs))

p_power_heatmap <- ggplot(power_plot_data, aes(x = endpoint, y = drug, fill = power)) +
  geom_tile(colour = "white", size = 1.2) +
  geom_text(aes(label = sprintf("%.0f%%", power)), size = 4.5, fontface = "bold") +
  facet_wrap(~ scale, ncol = 1) +
  scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850", midpoint = 80, limits = c(0, 100), name = "Power (%)") +
  labs(title = "Statistical Power by Endpoint Definition for Cyclically Dosed Agents",
       subtitle = paste0("Simulated RCTs (n=", n_per_arm, "/arm, ", n_sim_trials, " trials) | α = 0.05 one-sided"),
       x = "Endpoint", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
        axis.text.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(face = "bold", size = 14),
        legend.position = "right", panel.grid = element_blank())

print(p_power_heatmap)
ggsave("midd_endpoint_sensitivity_heatmap.png", p_power_heatmap, width = 16, height = 10, dpi = 400, bg = "white")

# Power curves
snapshot_power <- power_summary %>%
  filter(grepl("^Week", endpoint)) %>%
  mutate(week = as.numeric(str_extract(endpoint, "\\d+")))

p_power_curves <- ggplot(snapshot_power, aes(x = week, y = power, colour = drug, linetype = scale)) +
  geom_line(size = 1.3, alpha = 0.9) + geom_point(size = 3.5) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "gray40", size = 0.8) +
  scale_colour_brewer(palette = "Set1") +
  scale_linetype_manual(values = c("MG-ADL" = "solid", "QMG" = "dashed")) +
  labs(title = "Power of Snapshot Endpoints by Assessment Week",
       x = "Assessment Week", y = "Statistical Power (%)", colour = "Drug", linetype = "Scale") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        legend.position = "bottom")

print(p_power_curves)
ggsave("midd_power_curves_by_week.png", p_power_curves, width = 12, height = 8, dpi = 400, bg = "white")

# Snapshot vs integrated comparison
comparison_data <- power_summary %>%
  filter(endpoint %in% c("Week 26 Change", "AUEC (Week 0-26)", "Time-Averaged Change")) %>%
  mutate(endpoint_type = case_when(
    endpoint == "Week 26 Change" ~ "Snapshot\n(Week 26)",
    endpoint == "AUEC (Week 0-26)" ~ "Integrated\n(AUEC)",
    endpoint == "Time-Averaged Change" ~ "Integrated\n(Time-Averaged)"
  ),
  endpoint_type = factor(endpoint_type, levels = c("Snapshot\n(Week 26)", "Integrated\n(AUEC)", "Integrated\n(Time-Averaged)")))

p_power_bar <- ggplot(comparison_data, aes(x = endpoint_type, y = power, fill = drug)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.0f%%", power)), position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3.8, fontface = "bold") +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "gray40") +
  facet_wrap(~ scale) + scale_fill_brewer(palette = "Set2") +
  labs(title = "Snapshot vs Integrated Endpoint Power for Cyclically Dosed Agents",
       x = "Endpoint Type", y = "Statistical Power (%)", fill = "Drug") +
  coord_cartesian(ylim = c(0, 105)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        legend.position = "bottom", strip.text = element_text(face = "bold", size = 14))

print(p_power_bar)
ggsave("midd_snapshot_vs_integrated_power.png", p_power_bar, width = 14, height = 8, dpi = 400, bg = "white")


###############################################################################
# SECTION 24: MIDD MODULE 2 – TREATMENT SELECTION FRAMEWORK
###############################################################################

cat("\n=== MIDD Module 2: Treatment Selection Framework ===\n")

# Helper for Emax extraction
get_emax_draws_m2 <- function(fit, scale_label) {
  draws <- as_draws_df(fit)
  intercept <- draws$b_emax_Intercept
  drug_cols <- grep("^r_Drug__emax\\[", names(draws), value = TRUE)
  results <- list()
  for (col in drug_cols) {
    drug_name <- gsub("r_Drug__emax\\[([^,]+),.*", "\\1", col)
    emax_vals <- intercept + draws[[col]]
    results[[drug_name]] <- data.frame(draw = seq_along(emax_vals), Drug = drug_name,
                                        emax = emax_vals, scale = scale_label, stringsAsFactors = FALSE)
  }
  bind_rows(results)
}

emax_mgadl_m2 <- bind_rows(get_emax_draws_m2(fit_continuous, "MG-ADL"), get_emax_draws_m2(fit_cyclical, "MG-ADL"))
emax_qmg_m2   <- bind_rows(get_emax_draws_m2(fit_qmg_continuous, "QMG"), get_emax_draws_m2(fit_qmg_cyclical, "QMG"))

# Pairwise superiority
compute_pairwise_superiority <- function(emax_df, scale_label) {
  drugs <- unique(emax_df$Drug)
  n_drugs <- length(drugs)
  min_draws <- emax_df %>% group_by(Drug) %>% summarise(n = n()) %>% pull(n) %>% min()
  emax_aligned <- emax_df %>% group_by(Drug) %>% slice_head(n = min_draws) %>% mutate(draw = row_number()) %>% ungroup()
  mat <- matrix(NA, nrow = n_drugs, ncol = n_drugs, dimnames = list(drugs, drugs))
  for (i in 1:n_drugs) for (j in 1:n_drugs) {
    if (i == j) { mat[i,j] <- 50 } else {
      di <- emax_aligned %>% filter(Drug == drugs[i]) %>% pull(emax)
      dj <- emax_aligned %>% filter(Drug == drugs[j]) %>% pull(emax)
      n_c <- min(length(di), length(dj))
      mat[i,j] <- mean(di[1:n_c] > dj[1:n_c]) * 100
    }
  }
  as.data.frame(mat) %>% tibble::rownames_to_column("Drug_row") %>%
    pivot_longer(-Drug_row, names_to = "Drug_col", values_to = "P_superior") %>%
    mutate(scale = scale_label)
}

pairwise_qmg   <- compute_pairwise_superiority(emax_qmg_m2, "QMG")
pairwise_mgadl <- compute_pairwise_superiority(emax_mgadl_m2, "MG-ADL")

# SUCRA
compute_rank_probs <- function(emax_df, scale_label) {
  drugs <- unique(emax_df$Drug)
  min_draws <- emax_df %>% group_by(Drug) %>% summarise(n = n()) %>% pull(n) %>% min()
  emax_wide <- emax_df %>% group_by(Drug) %>% slice_head(n = min_draws) %>%
    mutate(draw = row_number()) %>% ungroup() %>%
    select(draw, Drug, emax) %>% pivot_wider(names_from = Drug, values_from = emax)
  n_draws <- nrow(emax_wide); n_drugs <- length(drugs)
  rank_matrix <- matrix(0, nrow = n_drugs, ncol = n_drugs, dimnames = list(drugs, paste0("Rank_", 1:n_drugs)))
  for (d in 1:n_draws) {
    vals <- as.numeric(emax_wide[d, -1]); names(vals) <- drugs
    ranks <- rank(-vals, ties.method = "random")
    for (drug in drugs) rank_matrix[drug, ranks[drug]] <- rank_matrix[drug, ranks[drug]] + 1
  }
  rank_probs <- rank_matrix / n_draws * 100
  sucra <- apply(rank_probs, 1, function(row) { n <- length(row); sum(cumsum(row)[1:(n-1)]) / (n-1) })
  as.data.frame(rank_probs) %>% tibble::rownames_to_column("Drug") %>%
    mutate(SUCRA = sucra, scale = scale_label) %>%
    pivot_longer(cols = starts_with("Rank_"), names_to = "Rank", values_to = "Probability") %>%
    mutate(Rank = as.integer(str_extract(Rank, "\\d+")))
}

rank_qmg   <- compute_rank_probs(emax_qmg_m2, "QMG")
rank_mgadl <- compute_rank_probs(emax_mgadl_m2, "MG-ADL")
rank_all   <- bind_rows(rank_qmg, rank_mgadl)

sucra_summary <- rank_all %>% select(Drug, SUCRA, scale) %>% distinct() %>% arrange(scale, desc(SUCRA))
cat("\n=== SUCRA Rankings ===\n")
print(sucra_summary, n = Inf)
write.csv(sucra_summary, "midd_sucra_rankings.csv", row.names = FALSE)

# Pairwise heatmaps
p_pairwise_qmg <- ggplot(pairwise_qmg, aes(x = Drug_col, y = Drug_row, fill = P_superior)) +
  geom_tile(colour = "white", size = 0.8) +
  geom_text(aes(label = sprintf("%.0f%%", P_superior)), size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 50, limits = c(0, 100)) +
  labs(title = "QMG: Pairwise Probability of Superiority (Emax)", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold", hjust = 0.5), panel.grid = element_blank())

p_pairwise_mgadl <- ggplot(pairwise_mgadl, aes(x = Drug_col, y = Drug_row, fill = P_superior)) +
  geom_tile(colour = "white", size = 0.8) +
  geom_text(aes(label = sprintf("%.0f%%", P_superior)), size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 50, limits = c(0, 100)) +
  labs(title = "MG-ADL: Pairwise Probability of Superiority (Emax)", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold", hjust = 0.5), panel.grid = element_blank())

p_pairwise_combined <- p_pairwise_qmg + p_pairwise_mgadl +
  plot_annotation(title = "Pairwise Posterior Probability of Superiority",
                  theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)))

print(p_pairwise_combined)
ggsave("midd_pairwise_superiority_heatmaps.png", p_pairwise_combined, width = 20, height = 10, dpi = 400, bg = "white")

# SUCRA rankogram
p_rankogram <- ggplot(rank_all, aes(x = Rank, y = Probability, fill = Drug)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.8) +
  facet_wrap(~ scale, ncol = 1) + scale_fill_brewer(palette = "Paired") +
  labs(title = "Rankogram: Probability of Each Efficacy Rank by Drug",
       x = "Rank", y = "Probability (%)", fill = "Drug") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        strip.text = element_text(face = "bold", size = 14), legend.position = "bottom")

print(p_rankogram)
ggsave("midd_sucra_rankogram.png", p_rankogram, width = 16, height = 12, dpi = 400, bg = "white")

# SUCRA bar
sucra_plot_data <- sucra_summary %>% mutate(Drug = fct_reorder(Drug, SUCRA))

p_sucra_bar <- ggplot(sucra_plot_data, aes(x = Drug, y = SUCRA, fill = scale)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.0f%%", SUCRA)), position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3.8, fontface = "bold") +
  coord_flip(ylim = c(0, 105)) +
  scale_fill_manual(values = c("QMG" = "#1f77b4", "MG-ADL" = "#ff7f0e")) +
  labs(title = "SUCRA by Drug and Scale", x = NULL, y = "SUCRA (%)", fill = "Scale") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5), legend.position = "bottom")

print(p_sucra_bar)
ggsave("midd_sucra_bar_chart.png", p_sucra_bar, width = 14, height = 10, dpi = 400, bg = "white")


###############################################################################
# SECTION 25: MIDD MODULE 2 ENHANCED – BENEFIT-RISK-CONVENIENCE
###############################################################################

cat("\n=== MIDD Module 2 Enhanced: Benefit-Risk-Convenience ===\n")

drug_chars <- read_excel(data_path, sheet = "Tolerability_safety_profiles")
drug_chars <- drug_chars %>%
  mutate(Drug = toupper(Drug)) %>%
  mutate(Drug = case_when(Drug == "NIPOCLIMAB" ~ "NIPOCALIMAB", TRUE ~ Drug))

numeric_cols <- c("any_ae_rate_percent", "serious_ae_rate_percent", "infection_rate_percent",
                  "severe_infection_rate_percent", "headache_rate_percent",
                  "infusion_reaction_rate_percent", "discontinuation_ae_rate",
                  "Completion_rate_percentage", "dropout_any_reason_percentage",
                  "admin_per_year", "doses_per_cycle", "cycle_length_weeks")

for (col in numeric_cols) {
  if (col %in% names(drug_chars))
    drug_chars[[col]] <- suppressWarnings(as.numeric(as.character(drug_chars[[col]])))
}

# Emax summaries
emax_qmg_sum_m2 <- emax_qmg_m2 %>% group_by(Drug) %>% summarise(emax_qmg = median(emax), .groups = "drop")
emax_mgadl_sum_m2 <- emax_mgadl_m2 %>% group_by(Drug) %>% summarise(emax_mgadl = median(emax), .groups = "drop")

# Onset and AUEC summaries
onset_qmg_sum <- time_to_mcid_qmg %>% group_by(Drug) %>% summarise(onset_qmg = median(time), .groups = "drop")
onset_mgadl_sum <- time_to_mcid %>% group_by(Drug) %>% summarise(onset_mgadl = median(time), .groups = "drop")
auec_qmg_sum <- auec_qmg_all %>% group_by(Drug) %>% summarise(auec_qmg = median(auec), .groups = "drop")
auec_mgadl_sum <- auec_all %>% group_by(Drug) %>% summarise(auec_mgadl = median(auec), .groups = "drop")

normalise <- function(x, invert = FALSE) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0 || max(x_clean) == min(x_clean)) return(rep(0.5, length(x)))
  r <- (x - min(x_clean)) / (max(x_clean) - min(x_clean))
  if (invert) r <- 1 - r
  return(pmax(0, pmin(1, r)))
}

drug_metrics <- drug_chars %>%
  select(Drug, Mechanism, Route, self_admin, admin_per_year, Weight_based, Vaccination_req,
         Completion_rate_percentage, infusion_duration_min,
         serious_ae_rate_percent, infection_rate_percent,
         headache_rate_percent, infusion_reaction_rate_percent,
         discontinuation_ae_rate, Menningococcal_risk) %>%
  left_join(emax_qmg_sum_m2, by = "Drug") %>%
  left_join(emax_mgadl_sum_m2, by = "Drug") %>%
  left_join(onset_qmg_sum, by = "Drug") %>%
  left_join(onset_mgadl_sum, by = "Drug") %>%
  left_join(auec_qmg_sum, by = "Drug") %>%
  left_join(auec_mgadl_sum, by = "Drug") %>%
  mutate(
    infusion_dur_numeric = suppressWarnings(as.numeric(as.character(infusion_duration_min))),
    infusion_dur_numeric = ifelse(is.na(infusion_dur_numeric) & Route == "s.c", 5, infusion_dur_numeric),
    infusion_dur_numeric = ifelse(is.na(infusion_dur_numeric), 60, infusion_dur_numeric)
  )

drug_scores <- drug_metrics %>%
  mutate(
    d1_qmg_efficacy = normalise(emax_qmg),
    d2_mgadl_efficacy = normalise(emax_mgadl),
    d3_onset_speed = normalise(rowMeans(cbind(
      ifelse(is.na(onset_qmg), NA, onset_qmg),
      ifelse(is.na(onset_mgadl), NA, onset_mgadl)), na.rm = TRUE), invert = TRUE),
    d4_cumulative_benefit = normalise(rowMeans(cbind(
      ifelse(is.na(auec_qmg), NA, auec_qmg),
      ifelse(is.na(auec_mgadl), NA, auec_mgadl)), na.rm = TRUE)),
    safety_composite = rowMeans(cbind(
      ifelse(is.na(serious_ae_rate_percent), 10, serious_ae_rate_percent),
      ifelse(is.na(infection_rate_percent), 20, infection_rate_percent),
      ifelse(is.na(discontinuation_ae_rate), 3, discontinuation_ae_rate)), na.rm = TRUE),
    d5_safety = normalise(safety_composite, invert = TRUE),
    tolerability_composite = Completion_rate_percentage -
      0.3 * ifelse(is.na(headache_rate_percent), 15, headache_rate_percent) -
      0.2 * ifelse(is.na(infusion_reaction_rate_percent), 5, infusion_reaction_rate_percent),
    d6_tolerability = normalise(tolerability_composite),
    convenience_composite = 25 * self_admin +
      normalise(admin_per_year, invert = TRUE) * 25 +
      (1 - Weight_based) * 15 + (1 - Vaccination_req) * 15 +
      normalise(infusion_dur_numeric, invert = TRUE) * 20,
    d7_convenience = normalise(convenience_composite),
    d8_mening_free = 1 - Menningococcal_risk
  )

dim_summary <- drug_scores %>%
  select(Drug, d1_qmg_efficacy, d2_mgadl_efficacy, d3_onset_speed,
         d4_cumulative_benefit, d5_safety, d6_tolerability, d7_convenience, d8_mening_free)

cat("\n=== Normalised Dimension Scores ===\n")
print(dim_summary, n = Inf, width = Inf)
write.csv(dim_summary, "midd_dimension_scores.csv", row.names = FALSE)

# Patient profiles
patient_profiles <- tribble(
  ~profile, ~w1, ~w2, ~w3, ~w4, ~w5, ~w6, ~w7, ~w8,
  "QMG-Dominant", 0.35, 0.05, 0.10, 0.10, 0.15, 0.10, 0.05, 0.10,
  "MG-ADL-Dominant", 0.05, 0.35, 0.10, 0.10, 0.15, 0.10, 0.05, 0.10,
  "Balanced", 0.15, 0.15, 0.10, 0.10, 0.15, 0.10, 0.15, 0.10,
  "Rapid Onset Priority", 0.10, 0.10, 0.35, 0.15, 0.10, 0.05, 0.05, 0.10,
  "Cumulative Benefit Priority", 0.10, 0.10, 0.05, 0.35, 0.15, 0.10, 0.05, 0.10,
  "Safety-Conscious", 0.10, 0.10, 0.05, 0.05, 0.35, 0.15, 0.05, 0.15,
  "Convenience-First", 0.05, 0.05, 0.05, 0.05, 0.10, 0.15, 0.40, 0.15,
  "Resource-Constrained Setting", 0.10, 0.15, 0.05, 0.10, 0.10, 0.10, 0.25, 0.15
)

composite_results <- list()
for (i in 1:nrow(patient_profiles)) {
  prof <- patient_profiles[i, ]
  scores <- drug_scores %>%
    mutate(composite = prof$w1 * d1_qmg_efficacy + prof$w2 * d2_mgadl_efficacy +
             prof$w3 * d3_onset_speed + prof$w4 * d4_cumulative_benefit +
             prof$w5 * d5_safety + prof$w6 * d6_tolerability +
             prof$w7 * d7_convenience + prof$w8 * d8_mening_free,
           profile = prof$profile) %>%
    select(Drug, Mechanism, composite, profile)
  composite_results[[i]] <- scores
}

composite_df <- bind_rows(composite_results)
write.csv(composite_df, "midd_enhanced_composite_scores.csv", row.names = FALSE)

composite_plot_df <- composite_df %>% filter(!is.na(composite))
drug_order <- composite_plot_df %>% group_by(Drug) %>%
  summarise(mean_score = mean(composite, na.rm = TRUE)) %>% arrange(mean_score) %>% pull(Drug)
composite_plot_df <- composite_plot_df %>% mutate(Drug = factor(Drug, levels = drug_order))

p_heatmap <- ggplot(composite_plot_df,
                     aes(x = factor(profile, levels = patient_profiles$profile), y = Drug, fill = composite)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", composite)), size = 3.8, fontface = "bold") +
  scale_fill_gradient2(low = "#d73027", mid = "#ffffbf", high = "#1a9850", midpoint = 0.5, limits = c(0, 1), name = "Composite\nScore") +
  labs(title = "Enhanced Treatment Selection Framework: Benefit-Risk-Convenience",
       x = "Patient Priority Profile", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
        axis.text.y = element_text(size = 11, face = "bold"),
        legend.position = "right", panel.grid = element_blank())

print(p_heatmap)
ggsave("midd_enhanced_treatment_heatmap.png", p_heatmap, width = 16, height = 10, dpi = 400, bg = "white")

# Radar plots
radar_data <- dim_summary %>%
  filter(!if_any(c(d1_qmg_efficacy, d2_mgadl_efficacy), is.na)) %>%
  pivot_longer(cols = c(d1_qmg_efficacy, d2_mgadl_efficacy, d3_onset_speed,
                         d4_cumulative_benefit, d5_safety, d6_tolerability,
                         d7_convenience, d8_mening_free),
               names_to = "dimension", values_to = "score") %>%
  mutate(dimension_label = factor(case_when(
    dimension == "d1_qmg_efficacy" ~ "QMG\nEfficacy",
    dimension == "d2_mgadl_efficacy" ~ "MG-ADL\nEfficacy",
    dimension == "d3_onset_speed" ~ "Onset\nSpeed",
    dimension == "d4_cumulative_benefit" ~ "Cumulative\nBenefit",
    dimension == "d5_safety" ~ "Safety",
    dimension == "d6_tolerability" ~ "Tolerability",
    dimension == "d7_convenience" ~ "Convenience",
    dimension == "d8_mening_free" ~ "No Mening.\nRisk"
  ), levels = c("QMG\nEfficacy", "MG-ADL\nEfficacy", "Onset\nSpeed", "Cumulative\nBenefit",
                "Safety", "Tolerability", "Convenience", "No Mening.\nRisk")))

p_radar <- ggplot(radar_data, aes(x = dimension_label, y = score, group = Drug)) +
  geom_polygon(aes(fill = Drug), alpha = 0.15, colour = NA) +
  geom_line(aes(colour = Drug), linewidth = 0.8) +
  geom_point(aes(colour = Drug), size = 2) +
  coord_polar(start = 0) + facet_wrap(~ Drug, ncol = 3) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(title = "Multi-Dimensional Drug Profile: Radar Plots") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "none", panel.grid.major = element_line(colour = "gray85", linewidth = 0.3))

print(p_radar)
ggsave("midd_radar_plots.png", p_radar, width = 14, height = 16, dpi = 400, bg = "white")

# Tornado sensitivity
balanced_weights <- as.numeric(patient_profiles[patient_profiles$profile == "Balanced", 2:9])
dim_names <- c("QMG Efficacy", "MG-ADL Efficacy", "Onset Speed", "Cumulative Benefit",
               "Safety", "Tolerability", "Convenience", "Mening. Risk-Free")

tornado_results <- list()
for (d in 1:8) {
  for (direction in c("high", "low")) {
    w_mod <- balanced_weights
    if (direction == "high") w_mod[d] <- w_mod[d] * 1.5 else w_mod[d] <- w_mod[d] * 0.5
    w_mod <- w_mod / sum(w_mod)
    scores <- drug_scores %>%
      mutate(composite = w_mod[1]*d1_qmg_efficacy + w_mod[2]*d2_mgadl_efficacy +
               w_mod[3]*d3_onset_speed + w_mod[4]*d4_cumulative_benefit +
               w_mod[5]*d5_safety + w_mod[6]*d6_tolerability +
               w_mod[7]*d7_convenience + w_mod[8]*d8_mening_free)
    top_drug <- scores %>% slice_max(composite, n = 1)
    tornado_results[[length(tornado_results) + 1]] <- data.frame(
      dimension = dim_names[d], direction = direction,
      top_drug = top_drug$Drug[1], top_score = top_drug$composite[1], stringsAsFactors = FALSE)
  }
}

tornado_df <- bind_rows(tornado_results) %>%
  mutate(dimension = factor(dimension, levels = rev(dim_names)))

baseline_scores_tornado <- drug_scores %>%
  mutate(composite = balanced_weights[1]*d1_qmg_efficacy + balanced_weights[2]*d2_mgadl_efficacy +
           balanced_weights[3]*d3_onset_speed + balanced_weights[4]*d4_cumulative_benefit +
           balanced_weights[5]*d5_safety + balanced_weights[6]*d6_tolerability +
           balanced_weights[7]*d7_convenience + balanced_weights[8]*d8_mening_free)
baseline_top <- baseline_scores_tornado %>% slice_max(composite, n = 1)

tornado_wide <- tornado_df %>% select(dimension, direction, top_score) %>%
  pivot_wider(names_from = direction, values_from = top_score)

p_tornado <- ggplot(tornado_wide, aes(y = dimension)) +
  geom_segment(aes(x = low, xend = high, yend = dimension), linewidth = 6, colour = "#4292c6", alpha = 0.7) +
  geom_point(aes(x = low), size = 4, colour = "#08519c") +
  geom_point(aes(x = high), size = 4, colour = "#08519c") +
  geom_vline(xintercept = baseline_top$composite[1], linetype = "dashed", colour = "red", linewidth = 0.8) +
  labs(title = "Tornado Sensitivity Analysis: Balanced Profile",
       x = "Composite Score of Top-Ranked Drug", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        axis.text.y = element_text(size = 12, face = "bold"), panel.grid.minor = element_blank())

print(p_tornado)
ggsave("midd_tornado_sensitivity.png", p_tornado, width = 12, height = 8, dpi = 400, bg = "white")


###############################################################################
# SECTION 26: FORMAL MODEL COMPARISON
###############################################################################

cat("=== Formal Model Comparison (Pooled vs Split) ===\n")

trt_all_pooled <- data %>%
  filter(Arm == "Treatment", time > 0, !is.na(change), !is.na(se), se > 0) %>%
  group_by(study_arm) %>%
  mutate(interval = time - lag(time, default = 0),
         cum_dose = cumsum(dose_eq_perweek * interval)) %>%
  ungroup()

fit_pooled_mgadl <- brm(
  bf(change | se(se, sigma = TRUE) ~ -emax * (cum_dose / (ed50 + cum_dose)) * (1 - exp(-k * time)),
     emax ~ 1 + (1 | ID | Drug), ed50 ~ 1, k ~ 1, nl = TRUE),
  data = trt_all_pooled,
  prior = c(prior(normal(5, 3), nlpar = "emax", lb = 0),
            prior(exponential(2), class = "sd", nlpar = "emax"),
            prior(normal(5000, 5000), nlpar = "ed50", lb = 0),
            prior(normal(0.1, 0.1), nlpar = "k", lb = 0)),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4, iter = 12000, warmup = 4000, seed = 123)

loo_pooled  <- loo(fit_pooled_mgadl)
loo_cont    <- loo(fit_continuous)
loo_cycl    <- loo(fit_cyclical)

elpd_pooled <- loo_pooled$estimates["elpd_loo", "Estimate"]
elpd_split  <- loo_cont$estimates["elpd_loo", "Estimate"] + loo_cycl$estimates["elpd_loo", "Estimate"]

cat("Pooled ELPD_LOO:", round(elpd_pooled, 2), "\n")
cat("Split ELPD_LOO:", round(elpd_split, 2), "\n")
cat("Difference (split - pooled):", round(elpd_split - elpd_pooled, 2), "\n")

model_comparison <- data.frame(
  Model = c("Pooled (single cumulative dose)", "Split (continuous + cyclical)"),
  ELPD_LOO = c(round(elpd_pooled, 2), round(elpd_split, 2)))
write.csv(model_comparison, "model_comparison_loo.csv", row.names = FALSE)
print(summary(fit_pooled_mgadl))


###############################################################################
# SECTION 27: PRIOR SENSITIVITY
###############################################################################

cat("\n=== Prior Sensitivity Analysis ===\n")

fit_cont_wide <- brm(
  bf(change | se(se, sigma = TRUE) ~ -emax * (cum_dose / (ed50 + cum_dose)) * (1 - exp(-k * time)),
     emax ~ 1 + (1 | ID | Drug), ed50 ~ 1, k ~ 1, nl = TRUE),
  data = trt_continuous,
  prior = c(prior(normal(5, 10), nlpar = "emax", lb = 0),
            prior(exponential(1), class = "sd", nlpar = "emax"),
            prior(normal(5000, 10000), nlpar = "ed50", lb = 0),
            prior(normal(0.1, 0.2), nlpar = "k", lb = 0)),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4, iter = 12000, warmup = 4000, seed = 123)

fit_cont_narrow <- brm(
  bf(change | se(se, sigma = TRUE) ~ -emax * (cum_dose / (ed50 + cum_dose)) * (1 - exp(-k * time)),
     emax ~ 1 + (1 | ID | Drug), ed50 ~ 1, k ~ 1, nl = TRUE),
  data = trt_continuous,
  prior = c(prior(normal(5, 1.5), nlpar = "emax", lb = 0),
            prior(exponential(3), class = "sd", nlpar = "emax"),
            prior(normal(5000, 2000), nlpar = "ed50", lb = 0),
            prior(normal(0.1, 0.05), nlpar = "k", lb = 0)),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4, iter = 12000, warmup = 4000, seed = 123)

extract_emax_ranking <- function(fit, label) {
  draws <- as_draws_df(fit)
  intercept <- draws$b_emax_Intercept
  drug_cols <- grep("^r_Drug__emax\\[", names(draws), value = TRUE)
  results <- list()
  for (col in drug_cols) {
    drug_name <- gsub("r_Drug__emax\\[([^,]+),.*", "\\1", col)
    emax_vals <- intercept + draws[[col]]
    results[[drug_name]] <- data.frame(Drug = drug_name, median_emax = median(emax_vals),
                                        lower = quantile(emax_vals, 0.025),
                                        upper = quantile(emax_vals, 0.975), stringsAsFactors = FALSE)
  }
  bind_rows(results) %>% arrange(desc(median_emax)) %>% mutate(rank = row_number(), prior_set = label)
}

rank_primary <- extract_emax_ranking(fit_continuous, "Primary")
rank_wide    <- extract_emax_ranking(fit_cont_wide, "Wide")
rank_narrow  <- extract_emax_ranking(fit_cont_narrow, "Narrow")

prior_sensitivity <- bind_rows(rank_primary, rank_wide, rank_narrow)
cat("\n--- Prior Sensitivity: MG-ADL Continuous Emax Rankings ---\n")
print(as_tibble(prior_sensitivity %>% select(Drug, prior_set, rank, median_emax) %>% arrange(Drug, prior_set)), n = Inf)
write.csv(prior_sensitivity, "prior_sensitivity_emax_rankings.csv", row.names = FALSE)

rank_stability <- prior_sensitivity %>% select(Drug, prior_set, rank) %>%
  pivot_wider(names_from = prior_set, values_from = rank)
cat("\n--- Rank Stability ---\n")
print(rank_stability)
write.csv(rank_stability, "prior_sensitivity_rank_stability.csv", row.names = FALSE)


###############################################################################
# SECTION 28: LEAVE-ONE-STUDY-OUT
###############################################################################

cat("\n=== Leave-One-Study-Out Sensitivity ===\n")

studies_cont <- unique(trt_continuous$study)
loo_study_results <- list()

for (s in studies_cont) {
  cat("  Excluding study:", as.character(s), "\n")
  data_loo <- trt_continuous %>% filter(study != s)
  if (n_distinct(data_loo$Drug) < 2) { cat("    Skipping\n"); next }
  tryCatch({
    fit_loo <- update(fit_continuous, newdata = data_loo, chains = 4, iter = 8000, warmup = 3000, seed = 123)
    loo_study_results[[as.character(s)]] <- extract_emax_ranking(fit_loo, paste0("Excl_", s))
  }, error = function(e) cat("    Error:", conditionMessage(e), "\n"))
}

if (length(loo_study_results) > 0) {
  loo_sensitivity_df <- bind_rows(loo_study_results)
  loo_rank_summary <- loo_sensitivity_df %>% select(Drug, prior_set, rank) %>%
    pivot_wider(names_from = prior_set, values_from = rank)
  cat("\n--- LOO Rank Stability ---\n")
  print(loo_rank_summary, width = Inf)
  write.csv(loo_rank_summary, "loo_study_rank_stability.csv", row.names = FALSE)
  write.csv(loo_sensitivity_df, "loo_study_emax_results.csv", row.names = FALSE)
}
###############################################################################
# SECTION 27b: PRIOR SENSITIVITY – QMG CONTINUOUS
###############################################################################

cat("\n=== Prior Sensitivity Analysis (QMG Continuous) ===\n")

fit_qmg_cont_wide <- brm(
  bf(change | se(se, sigma = TRUE) ~ -emax * (cum_dose / (ed50 + cum_dose)) * (1 - exp(-k * time)),
     emax ~ 1 + (1 | ID | Drug), ed50 ~ 1, k ~ 1, nl = TRUE),
  data = trt_qmg_continuous,
  prior = c(prior(normal(8, 16), nlpar = "emax", lb = 0),
            prior(exponential(1), class = "sd", nlpar = "emax"),
            prior(normal(5000, 10000), nlpar = "ed50", lb = 0),
            prior(normal(0.1, 0.2), nlpar = "k", lb = 0)),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4, iter = 12000, warmup = 4000, seed = 123)

fit_qmg_cont_narrow <- brm(
  bf(change | se(se, sigma = TRUE) ~ -emax * (cum_dose / (ed50 + cum_dose)) * (1 - exp(-k * time)),
     emax ~ 1 + (1 | ID | Drug), ed50 ~ 1, k ~ 1, nl = TRUE),
  data = trt_qmg_continuous,
  prior = c(prior(normal(8, 2), nlpar = "emax", lb = 0),
            prior(exponential(3), class = "sd", nlpar = "emax"),
            prior(normal(5000, 2000), nlpar = "ed50", lb = 0),
            prior(normal(0.1, 0.05), nlpar = "k", lb = 0)),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4, iter = 12000, warmup = 4000, seed = 123)

rank_qmg_primary <- extract_emax_ranking(fit_qmg_continuous, "Primary")
rank_qmg_wide    <- extract_emax_ranking(fit_qmg_cont_wide, "Wide")
rank_qmg_narrow  <- extract_emax_ranking(fit_qmg_cont_narrow, "Narrow")

prior_sensitivity_qmg <- bind_rows(rank_qmg_primary, rank_qmg_wide, rank_qmg_narrow)
cat("\n--- Prior Sensitivity: QMG Continuous Emax Rankings ---\n")
print(as_tibble(prior_sensitivity_qmg %>% select(Drug, prior_set, rank, median_emax) %>% arrange(Drug, prior_set)), n = Inf)
write.csv(prior_sensitivity_qmg, "prior_sensitivity_qmg_emax_rankings.csv", row.names = FALSE)

rank_stability_qmg <- prior_sensitivity_qmg %>% select(Drug, prior_set, rank) %>%
  pivot_wider(names_from = prior_set, values_from = rank)
cat("\n--- QMG Rank Stability Across Priors ---\n")
print(as_tibble(rank_stability_qmg), n = Inf)
write.csv(rank_stability_qmg, "prior_sensitivity_qmg_rank_stability.csv", row.names = FALSE)

###############################################################################
# SECTION 28b: LEAVE-ONE-STUDY-OUT – QMG CONTINUOUS
###############################################################################

cat("\n=== Leave-One-Study-Out Sensitivity (QMG Continuous) ===\n")

studies_qmg_cont <- unique(trt_qmg_continuous$study)
loo_qmg_study_results <- list()

for (s in studies_qmg_cont) {
  cat("  Excluding study:", as.character(s), "\n")
  data_loo_qmg <- trt_qmg_continuous %>% filter(study != s)
  if (n_distinct(data_loo_qmg$Drug) < 2) { cat("    Skipping (< 2 drugs remaining)\n"); next }
  tryCatch({
    fit_loo_qmg <- update(fit_qmg_continuous, newdata = data_loo_qmg,
                          chains = 4, iter = 8000, warmup = 3000, seed = 123)
    loo_qmg_study_results[[as.character(s)]] <- extract_emax_ranking(fit_loo_qmg, paste0("Excl_", s))
  }, error = function(e) cat("    Error:", conditionMessage(e), "\n"))
}

if (length(loo_qmg_study_results) > 0) {
  loo_qmg_sensitivity_df <- bind_rows(loo_qmg_study_results)
  loo_qmg_rank_summary <- loo_qmg_sensitivity_df %>% select(Drug, prior_set, rank) %>%
    pivot_wider(names_from = prior_set, values_from = rank)
  cat("\n--- QMG LOO Rank Stability ---\n")
  print(as_tibble(loo_qmg_rank_summary), n = Inf, width = Inf)
  write.csv(loo_qmg_rank_summary, "loo_qmg_study_rank_stability.csv", row.names = FALSE)
  write.csv(loo_qmg_sensitivity_df, "loo_qmg_study_emax_results.csv", row.names = FALSE)
}

###############################################################################
# SECTION 29: CROSS-SCALE CONCORDANCE
###############################################################################

cat("\n=== Cross-Scale Concordance Analysis ===\n")

drugs_both <- intersect(unique(emax_mgadl_m2$Drug), unique(emax_qmg_m2$Drug))
cat("Drugs in both scales:", paste(drugs_both, collapse = ", "), "\n")

set.seed(2024)
n_boot <- 2000
concordance_rho <- numeric(n_boot)

for (d in 1:n_boot) {
  mgadl_s <- emax_mgadl_m2 %>% filter(Drug %in% drugs_both) %>% group_by(Drug) %>% slice_sample(n = 1) %>% ungroup()
  qmg_s   <- emax_qmg_m2 %>% filter(Drug %in% drugs_both) %>% group_by(Drug) %>% slice_sample(n = 1) %>% ungroup()
  merged <- inner_join(mgadl_s %>% select(Drug, emax_mgadl = emax),
                       qmg_s %>% select(Drug, emax_qmg = emax), by = "Drug")
  concordance_rho[d] <- if (nrow(merged) >= 3) cor(rank(-merged$emax_mgadl), rank(-merged$emax_qmg), method = "spearman") else NA
}

concordance_rho <- concordance_rho[!is.na(concordance_rho)]

cat("Posterior median rho:", round(median(concordance_rho), 3), "\n")
cat("95% CrI:", round(quantile(concordance_rho, 0.025), 3), "to", round(quantile(concordance_rho, 0.975), 3), "\n")
cat("P(rho > 0):", round(mean(concordance_rho > 0) * 100, 1), "%\n")

emax_mgadl_med <- emax_mgadl_m2 %>% filter(Drug %in% drugs_both) %>% group_by(Drug) %>% summarise(emax_mgadl = median(emax))
emax_qmg_med   <- emax_qmg_m2 %>% filter(Drug %in% drugs_both) %>% group_by(Drug) %>% summarise(emax_qmg = median(emax))

concordance_point <- inner_join(emax_mgadl_med, emax_qmg_med, by = "Drug") %>%
  mutate(rank_mgadl = rank(-emax_mgadl), rank_qmg = rank(-emax_qmg), rank_diff = abs(rank_mgadl - rank_qmg)) %>%
  arrange(rank_mgadl)

cat("\n--- Drug-Level Rank Comparison ---\n")
print(concordance_point, n = Inf)
write.csv(concordance_point, "cross_scale_concordance.csv", row.names = FALSE)

p_concordance <- ggplot(concordance_point, aes(x = emax_mgadl, y = emax_qmg, label = Drug)) +
  geom_point(size = 4, colour = "#2166ac") +
  geom_text(hjust = -0.15, vjust = -0.5, size = 3.5, fontface = "bold") +
  geom_smooth(method = "lm", se = TRUE, colour = "gray50", linetype = "dashed") +
  labs(title = "Cross-Scale Concordance: MG-ADL vs QMG Emax",
       subtitle = paste0("Spearman rho = ", round(median(concordance_rho), 2),
                         " [", round(quantile(concordance_rho, 0.025), 2), ", ",
                         round(quantile(concordance_rho, 0.975), 2), "]"),
       x = "MG-ADL Emax", y = "QMG Emax") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))

print(p_concordance)
ggsave("cross_scale_concordance_scatter.png", p_concordance, width = 10, height = 8, dpi = 400, bg = "white")


###############################################################################
# SECTION 30: MECHANISM-LEVEL CLASS CONTRAST
###############################################################################

cat("\n=== Mechanism-Level Class Contrast ===\n")

mechanism_map <- data %>% filter(Arm == "Treatment") %>% select(Drug, Mechanism) %>%
  distinct() %>% mutate(Drug = as.character(Drug))

# MG-ADL
emax_mgadl_mech <- emax_mgadl_m2 %>% left_join(mechanism_map, by = "Drug") %>%
  filter(Mechanism %in% c("Complement inhibitor", "FcRn inhibitor"))

class_mgadl <- emax_mgadl_mech %>% group_by(draw, Mechanism) %>%
  summarise(class_emax = mean(emax), .groups = "drop") %>%
  pivot_wider(names_from = Mechanism, values_from = class_emax) %>% drop_na()

p_comp_mgadl <- mean(class_mgadl$`Complement inhibitor` > class_mgadl$`FcRn inhibitor`) * 100

cat("MG-ADL: P(Complement > FcRn) =", round(p_comp_mgadl, 1), "%\n")

# QMG
emax_qmg_mech <- emax_qmg_m2 %>% left_join(mechanism_map, by = "Drug") %>%
  filter(Mechanism %in% c("Complement inhibitor", "FcRn inhibitor"))

class_qmg <- emax_qmg_mech %>% group_by(draw, Mechanism) %>%
  summarise(class_emax = mean(emax), .groups = "drop") %>%
  pivot_wider(names_from = Mechanism, values_from = class_emax) %>% drop_na()

p_fcrn_qmg_class <- mean(class_qmg$`FcRn inhibitor` > class_qmg$`Complement inhibitor`) * 100

cat("QMG: P(FcRn > Complement) =", round(p_fcrn_qmg_class, 1), "%\n")

class_summary_table <- data.frame(
  Scale = c("MG-ADL", "MG-ADL", "QMG", "QMG"),
  Class = c("Complement inhibitor", "FcRn inhibitor", "Complement inhibitor", "FcRn inhibitor"),
  Mean_Emax = round(c(mean(class_mgadl$`Complement inhibitor`), mean(class_mgadl$`FcRn inhibitor`),
                        mean(class_qmg$`Complement inhibitor`), mean(class_qmg$`FcRn inhibitor`)), 2))

print(class_summary_table)
write.csv(class_summary_table, "mechanism_class_contrast.csv", row.names = FALSE)

p_class_mgadl <- ggplot(emax_mgadl_mech, aes(x = emax, fill = Mechanism)) +
  stat_halfeye(.width = c(0.66, 0.95), slab_alpha = 0.7) +
  labs(title = "MG-ADL: Mechanism-Level Emax",
       subtitle = paste0("P(Complement > FcRn) = ", round(p_comp_mgadl, 1), "%"),
       x = "Class-Level Emax", fill = "Mechanism") +
  scale_fill_manual(values = c("Complement inhibitor" = "#1f77b4", "FcRn inhibitor" = "#ff7f0e")) +
  theme_minimal(base_size = 14) + theme(legend.position = "bottom")

p_class_qmg <- ggplot(emax_qmg_mech, aes(x = emax, fill = Mechanism)) +
  stat_halfeye(.width = c(0.66, 0.95), slab_alpha = 0.7) +
  labs(title = "QMG: Mechanism-Level Emax",
       subtitle = paste0("P(FcRn > Complement) = ", round(p_fcrn_qmg_class, 1), "%"),
       x = "Class-Level Emax", fill = "Mechanism") +
  scale_fill_manual(values = c("Complement inhibitor" = "#1f77b4", "FcRn inhibitor" = "#ff7f0e")) +
  theme_minimal(base_size = 14) + theme(legend.position = "bottom")

p_class_combined <- p_class_mgadl + p_class_qmg +
  plot_annotation(title = "Mechanism-Level Class Contrast",
                  theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)))

print(p_class_combined)
ggsave("mechanism_class_contrast_density.png", p_class_combined, width = 16, height = 8, dpi = 400, bg = "white")


###############################################################################
# SECTION 31: HETEROGENEITY + ED50 IDENTIFIABILITY
###############################################################################

cat("\n=== Heterogeneity Quantification ===\n")

extract_tau <- function(fit, label) {
  draws <- as_draws_df(fit)
  tau_col <- grep("^sd_Drug__emax_Intercept", names(draws), value = TRUE)
  if (length(tau_col) == 0) tau_col <- grep("^sd_Drug", names(draws), value = TRUE)[1]
  if (length(tau_col) == 0 || is.na(tau_col)) return(NULL)
  tau_vals <- draws[[tau_col]]
  data.frame(Model = label, tau_median = round(median(tau_vals), 3),
             tau_lower = round(quantile(tau_vals, 0.025), 3),
             tau_upper = round(quantile(tau_vals, 0.975), 3))
}

tau_results <- bind_rows(
  extract_tau(fit_continuous, "MG-ADL Continuous"),
  extract_tau(fit_cyclical, "MG-ADL Cyclical"),
  extract_tau(fit_qmg_continuous, "QMG Continuous"),
  extract_tau(fit_qmg_cyclical, "QMG Cyclical"))

if (nrow(tau_results) > 0) {
  tau_results$label <- paste0(tau_results$tau_median, " [", tau_results$tau_lower, ", ", tau_results$tau_upper, "]")
  cat("--- Between-Drug Heterogeneity (tau) ---\n")
  print(tau_results %>% select(Model, `tau [95% CrI]` = label))
  write.csv(tau_results, "heterogeneity_tau.csv", row.names = FALSE)
}

cat("\n=== ED50 Identifiability Check ===\n")

ed50_diagnostics <- list()
for (info in list(list(fit = fit_continuous, label = "MG-ADL Continuous"),
                  list(fit = fit_cyclical, label = "MG-ADL Cyclical"),
                  list(fit = fit_qmg_continuous, label = "QMG Continuous"),
                  list(fit = fit_qmg_cyclical, label = "QMG Cyclical"))) {
  ed50 <- as_draws_df(info$fit)$b_ed50_Intercept
  ed50_diagnostics[[info$label]] <- data.frame(
    Model = info$label, ED50_median = round(median(ed50), 1),
    ED50_lower = round(quantile(ed50, 0.025), 1),
    ED50_upper = round(quantile(ed50, 0.975), 1),
    CV_percent = round(sd(ed50) / mean(ed50) * 100, 1))
}

ed50_table <- bind_rows(ed50_diagnostics)
ed50_table$label <- paste0(ed50_table$ED50_median, " [", ed50_table$ED50_lower, ", ", ed50_table$ED50_upper, "] CV=", ed50_table$CV_percent, "%")
cat("--- ED50 Posterior Summary ---\n")
print(ed50_table %>% select(Model, `ED50 [95% CrI]` = label))
write.csv(ed50_table, "ed50_identifiability.csv", row.names = FALSE)


###############################################################################
# SECTION 32: SENSITIVITY – MCID THRESHOLD
###############################################################################

cat("\n=== MCID Threshold Sensitivity Analysis ===\n")

# ---- MG-ADL: Alternative MCID = -3 (vs primary MCID = -2) ----
MCID_alt_mgadl <- -3

time_to_mcid_alt_mgadl <- pred_time_all_post %>%
  group_by(draw, Drug, regimen) %>%
  arrange(time) %>%
  filter(pred_change <= MCID_alt_mgadl) %>%
  slice(1) %>% ungroup()

time_mcid_alt_mgadl_summary <- time_to_mcid_alt_mgadl %>%
  group_by(Drug, regimen) %>%
  mean_qi(time, .width = 0.95) %>%
  mutate(label = sprintf("%.1f [%.1f, %.1f] weeks", time, .lower, .upper)) %>%
  arrange(time)

p_time_mcid_alt_mgadl <- ggplot(time_to_mcid_alt_mgadl, aes(x = time, color = regimen)) +
  stat_halfeye(aes(y = reorder(Drug, time, mean)), .width = c(0.66, 0.95), slab_alpha = 0.8) +
  geom_text(data = time_mcid_alt_mgadl_summary,
            aes(y = reorder(Drug, time, mean), label = label),
            x = Inf, hjust = -0.05, size = 4.5, fontface = "bold") +
  coord_cartesian(clip = "off") +
  labs(title = "Sensitivity: Time to MCID (≥3-Point MG-ADL Improvement)",
       subtitle = "Alternative threshold | Posterior median [95% CrI] weeks",
       x = "Time to MCID (weeks)", y = NULL, color = "Regimen Type") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "bottom", plot.margin = margin(10, 160, 10, 10))

print(p_time_mcid_alt_mgadl)
ggsave("sensitivity_mcid3_mgadl_time_to_mcid.png", p_time_mcid_alt_mgadl,
       width = 16, height = 10, dpi = 400, bg = "white")

# Compare rank orderings
rank_primary_mgadl <- time_mcid_summary %>%
  arrange(time) %>% mutate(rank_primary = row_number()) %>%
  select(Drug, regimen, time_primary = time, rank_primary)

rank_alt_mgadl <- time_mcid_alt_mgadl_summary %>%
  arrange(time) %>% mutate(rank_alt = row_number()) %>%
  select(Drug, regimen, time_alt = time, rank_alt)

mcid_rank_comparison_mgadl <- left_join(rank_primary_mgadl, rank_alt_mgadl,
                                         by = c("Drug", "regimen")) %>%
  mutate(rank_shift = abs(rank_primary - rank_alt))

cat("\n--- MG-ADL MCID Rank Comparison (MCID -2 vs -3) ---\n")
print(as_tibble(mcid_rank_comparison_mgadl), n = Inf)
write.csv(mcid_rank_comparison_mgadl, "sensitivity_mcid_mgadl_rank_comparison.csv", row.names = FALSE)

# ---- QMG: Alternative MCID = -5 (vs primary MCID = -3) ----
MCID_alt_qmg <- -5

time_to_mcid_alt_qmg <- pred_time_qmg_all_post %>%
  group_by(draw, Drug, regimen) %>%
  arrange(time) %>%
  filter(pred_change <= MCID_alt_qmg) %>%
  slice(1) %>% ungroup()

time_mcid_alt_qmg_summary <- time_to_mcid_alt_qmg %>%
  group_by(Drug, regimen) %>%
  mean_qi(time, .width = 0.95) %>%
  mutate(label = sprintf("%.1f [%.1f, %.1f] weeks", time, .lower, .upper)) %>%
  arrange(time)

p_time_mcid_alt_qmg <- ggplot(time_to_mcid_alt_qmg, aes(x = time, color = regimen)) +
  stat_halfeye(aes(y = reorder(Drug, time, mean)), .width = c(0.66, 0.95), slab_alpha = 0.8) +
  geom_text(data = time_mcid_alt_qmg_summary,
            aes(y = reorder(Drug, time, mean), label = label),
            x = Inf, hjust = -0.05, size = 4.5, fontface = "bold") +
  coord_cartesian(clip = "off") +
  labs(title = "Sensitivity: Time to MCID (≥5-Point QMG Improvement)",
       subtitle = "Alternative threshold | Posterior median [95% CrI] weeks",
       x = "Time to MCID (weeks)", y = NULL, color = "Regimen Type") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "bottom", plot.margin = margin(10, 160, 10, 10))

print(p_time_mcid_alt_qmg)
ggsave("sensitivity_mcid5_qmg_time_to_mcid.png", p_time_mcid_alt_qmg,
       width = 16, height = 10, dpi = 400, bg = "white")

rank_primary_qmg <- time_mcid_qmg_summary %>%
  arrange(time) %>% mutate(rank_primary = row_number()) %>%
  select(Drug, regimen, time_primary = time, rank_primary)

rank_alt_qmg <- time_mcid_alt_qmg_summary %>%
  arrange(time) %>% mutate(rank_alt = row_number()) %>%
  select(Drug, regimen, time_alt = time, rank_alt)

mcid_rank_comparison_qmg <- left_join(rank_primary_qmg, rank_alt_qmg,
                                       by = c("Drug", "regimen")) %>%
  mutate(rank_shift = abs(rank_primary - rank_alt))

cat("\n--- QMG MCID Rank Comparison (MCID -3 vs -5) ---\n")
print(as_tibble(mcid_rank_comparison_qmg), n = Inf)
write.csv(mcid_rank_comparison_qmg, "sensitivity_mcid_qmg_rank_comparison.csv", row.names = FALSE)

cat("\nMG-ADL max rank shift:", max(mcid_rank_comparison_mgadl$rank_shift, na.rm = TRUE), "\n")
cat("QMG max rank shift:", max(mcid_rank_comparison_qmg$rank_shift, na.rm = TRUE), "\n")


###############################################################################
# SECTION 33: SENSITIVITY – MINT PLACEBO EXCLUSION (MG-ADL)
###############################################################################

cat("\n=== MINT Placebo Exclusion Sensitivity (MG-ADL) ===\n")

plac_excl_mint <- data %>% filter(Arm == "Control", time > 0, Trial != "MINT")
stopifnot(all(plac_excl_mint$se > 0))

cat("Primary placebo N arms:", n_distinct(plac$study_arm),
    "| Excl-MINT N arms:", n_distinct(plac_excl_mint$study_arm), "\n")

fit_plac_excl_mint <- brm(
  bf(
    change | se(se, sigma = TRUE) ~ plateau * (1 - exp(-k * time)),
    plateau ~ 1 + (1 | ID | study),
    k ~ 1 + (1 | ID | study),
    nl = TRUE
  ),
  data = plac_excl_mint,
  prior = c(
    prior(normal(-3, 1), nlpar = "plateau", lb = -10),
    prior(normal(0.1, 0.1), nlpar = "k", lb = 0),
    prior(exponential(2), class = "sd", nlpar = "plateau"),
    prior(exponential(2), class = "sd", nlpar = "k")
  ),
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  chains = 4, iter = 20000, warmup = 10000,
  seed = 123
)

# Compare placebo parameter estimates
plac_primary_summary <- as.data.frame(summary(fit_plac_exp)$fixed) %>%
  tibble::rownames_to_column("Parameter") %>%
  mutate(model = "Primary (all trials)")

plac_excl_summary <- as.data.frame(summary(fit_plac_excl_mint)$fixed) %>%
  tibble::rownames_to_column("Parameter") %>%
  mutate(model = "Excl. MINT")

plac_comparison <- bind_rows(plac_primary_summary, plac_excl_summary) %>%
  select(model, Parameter, Estimate, `Est.Error`, `l-95% CI`, `u-95% CI`, Rhat)

cat("\n--- Placebo Model Parameter Comparison ---\n")
print(as_tibble(plac_comparison), n = Inf, width = Inf)
write.csv(plac_comparison, "sensitivity_mint_placebo_comparison.csv", row.names = FALSE)

# Re-derive net benefit with excluded-MINT placebo
# Predict population-level placebo at week 26 only (scalar)
plac_wk26_grid <- data.frame(
  time = 26, se = 0,
  study = as.character(unique(plac_excl_mint$study)[1])
)
plac_wk26_pred <- posterior_predict(fit_plac_excl_mint, newdata = plac_wk26_grid,
                                    ndraws = 500, allow_new_levels = TRUE, re_formula = NA)
plac_wk26_sens <- as.numeric(median(plac_wk26_pred))

time_grid_sens <- seq(1, 26, by = 1)

# Continuous treatment predictions
cont_pred_grid_sens <- trt_continuous %>%
  group_by(Drug) %>%
  summarise(dose_eq_mean = mean(dose_eq_perweek, na.rm = TRUE), .groups = "drop") %>%
  crossing(time = time_grid_sens) %>%
  mutate(cum_dose = dose_eq_mean * time, dose_eq = dose_eq_mean, se = 0)

cont_pred_sens <- posterior_predict(fit_continuous, newdata = cont_pred_grid_sens,
                                    ndraws = 500, allow_new_levels = TRUE)

cont_pred_df_sens <- cont_pred_sens %>%
  as.data.frame() %>% mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "obs_str", values_to = "pred_change") %>%
  mutate(obs = as.integer(str_remove(obs_str, "^V"))) %>% select(-obs_str) %>%
  left_join(cont_pred_grid_sens %>% mutate(obs = row_number()), by = "obs")

# Week-26 net benefit (excl-MINT placebo)
net_benefit_sens_cont <- cont_pred_df_sens %>%
  filter(time == 26) %>%
  group_by(Drug) %>%
  summarise(
    trt_median = median(pred_change),
    net_benefit = median(pred_change) - plac_wk26_sens,
    net_lower = quantile(pred_change, 0.025) - plac_wk26_sens,
    net_upper = quantile(pred_change, 0.975) - plac_wk26_sens,
    .groups = "drop"
  ) %>% mutate(model = "Excl. MINT placebo", regimen = "Continuous")

# Cyclical treatment predictions
cycl_pred_grid_sens <- trt_cycl_onperiod %>%
  group_by(Drug) %>%
  summarise(dose_eq = mean(dose_eq, na.rm = TRUE), .groups = "drop") %>%
  crossing(time = time_grid_sens) %>%
  mutate(se = 0)

cycl_pred_sens <- posterior_predict(fit_cyclical, newdata = cycl_pred_grid_sens,
                                    ndraws = 500, allow_new_levels = TRUE)

cycl_pred_df_sens <- cycl_pred_sens %>%
  as.data.frame() %>% mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "obs_str", values_to = "pred_change") %>%
  mutate(obs = as.integer(str_remove(obs_str, "^V"))) %>% select(-obs_str) %>%
  left_join(cycl_pred_grid_sens %>% mutate(obs = row_number()), by = "obs")

net_benefit_sens_cycl <- cycl_pred_df_sens %>%
  filter(time == 26) %>%
  group_by(Drug) %>%
  summarise(
    trt_median = median(pred_change),
    net_benefit = median(pred_change) - plac_wk26_sens,
    net_lower = quantile(pred_change, 0.025) - plac_wk26_sens,
    net_upper = quantile(pred_change, 0.975) - plac_wk26_sens,
    .groups = "drop"
  ) %>% mutate(model = "Excl. MINT placebo", regimen = "Cyclical")

net_benefit_sens <- bind_rows(net_benefit_sens_cont, net_benefit_sens_cycl) %>%
  arrange(regimen, desc(abs(net_benefit)))

cat("\n--- MG-ADL Week-26 Net Benefit (Excl. MINT Placebo) ---\n")
print(as_tibble(net_benefit_sens), n = Inf, width = Inf)
write.csv(net_benefit_sens, "sensitivity_mint_net_benefit_mgadl.csv", row.names = FALSE)

# Rank comparison
rank_sens <- net_benefit_sens %>%
  group_by(regimen) %>%
  arrange(desc(abs(net_benefit))) %>%
  mutate(rank_sens = row_number()) %>%
  select(Drug, regimen, net_benefit_sens = net_benefit, rank_sens)

cat("\n--- MINT Exclusion Sensitivity: Rank Comparison ---\n")
print(as_tibble(rank_sens), n = Inf)
write.csv(rank_sens, "sensitivity_mint_rank_comparison.csv", row.names = FALSE)


###############################################################################
# SECTION 34: CONSOLIDATED MODEL SUMMARIES FOR MANUSCRIPT
###############################################################################

cat("\nMG-ADL continuous\n"); print(summary(fit_continuous))
cat("\nMG-ADL cyclical\n"); print(summary(fit_cyclical))
cat("\nMG-ADL placebo\n"); print(summary(fit_plac_exp))
cat("\nQMG continuous\n"); print(summary(fit_qmg_continuous))
cat("\nQMG cyclical\n"); print(summary(fit_qmg_cyclical))
cat("\nQMG placebo\n"); print(summary(fit_plac_qmg))


###############################################################################
# SECTION 35: FINAL SAVE
###############################################################################

save.image(file = file.path(output_dir, "gMG_MBMA_complete_workspace.RData"))

cat("Analysis complete. Outputs in:", output_dir, "\n")

