if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)
pacman::p_load(ggplot2, 
               dplyr,
               readr,
               tidyr,
               here,
               rio
               )
options(scipen = 999)

source(here("Codes", "SCCS_main_functions.R"))
source(here("Codes", "SCCS_sim_workflow_functions.R"))
source(here("Codes", "SCCS_helper_functions.R"))
source(here("Codes", "nestedloop.R"))

################################# 
# 1 - Table of all scenarios #### 
################################# 

scenarios <- scenario_table(p_C_vals = c(5e-4, 1e-3, 5e-3), 
                            C_on_E_vals = c(0.25, 1, 2, 5),
                            IRR_E_vals = c(1, 2), 
                            chosen_risk_window_C = c(10, 30, 60),
                            U_on_Y_vals = c(0.25, 1, 2, 15))

##############################
# 2 - List of methods #### 
##############################

# Need to be among 'ignore_C', 'exclude_C', 'censor_C', 'adjust_C', 'adjust_C_obs_at_exp', and 'adjust_calendar'
methods <- c("ignore_C",              # ignore C in model
             "exclude_C",             # exclude subjects with C
             "censor_C",              # censor observation time at first C occurrence
             "adjust_C")              # adjust for C in SCCS model

##############################
# 3 - Quantify bias and summarise results ####
##############################
bias_results <- summarise_simulation_results(method_scen = method_scen(scenario_table = scenarios),
                                             results_dir = here("Results"),
                                             summary_dir = file.path(here("Results"), "Summary"),
                                             summary_file_name = "Summary_all_scens_2026_05_01")

##############################
# 4 - Exploratory analysis ####
##############################

# Extract data --------------
data_exploratory <- function(method_scen_ = method_scen(scenario_table = scenarios),
                             results_dir = here("Results")){
  scenarios <- unique(method_scen_$scen_name)
  methods <- unique(method_scen_$methods)
  
  all_data <- list()
  
  for (i in seq_len(nrow(method_scen_))) {
    method_i <- method_scen_$methods[i]
    scen_file <- paste0(method_scen_$scen_name[i], ".csv")
    
    file_path <- file.path(results_dir, method_i, scen_file)
    
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }
    
    df <- tryCatch(read_csv(file_path, show_col_types = FALSE),
                   error = function(e) NULL)
    
    if (is.null(df) || !"est_E" %in% names(df)) {
      warning("gamma_1_hat missing in: ", file_path)
      next
    }
    
    all_data[[length(all_data) + 1]] <- data.frame(
      scen_name = method_scen_$scen_name[i],
      method = method_i,
      gamma_1_hat = df$est_E,
      se_gamma_1 = df$se_E,
      n_event = df$n_event
    )
  }
  
  data_all <- bind_rows(all_data)
  return(data_all)
}

data_explor <- data_exploratory()
non_convergence <- data_explor %>% filter(gamma_1_hat <= -10)
data_explor_converge <- data_explor %>% filter(!is.na(gamma_1_hat) & gamma_1_hat > -10)

# Histogram of gamma_1, n_event and se_gamma_1 ---------------
plot_exploratory <- function(
    data,
    output_dir = file.path(here(),"Plots"),
    width = 7,
    height = 6,
    dpi = 300
){
  scenarios <- unique(data$scen_name)
  
  output_dir1 = file.path(output_dir,"Hist_gamma_1_hat")
  output_dir2 = file.path(output_dir,"Hist_se_gamma_1")
  output_dir3 = file.path(output_dir,"Hist_n_event")
  create_directory(output_dir1)
  create_directory(output_dir2)
  create_directory(output_dir3)
  
  
  for (sc in scenarios) {
    df_sc <- data %>% filter(scen_name == sc)
    
    if (nrow(df_sc) == 0) {
      warning("No data for scenario: ", sc)
      next
    }
    
    p <- ggplot(df_sc, aes(x = method, y = gamma_1_hat, fill = method)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.9) +
      labs(
        title = paste("Scenario:", sc),
        x = "Method",
        y = "gamma_1_hat"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "none"
      )
    
    ggsave(
      filename = file.path(output_dir1, paste0("plot_", sc, ".png")),
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
    
    p2 <- ggplot(df_sc, aes(x = method, y = se_gamma_1, fill = method)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.9) +
      labs(
        title = paste("Scenario:", sc),
        x = "Method",
        y = "se_gamma_1"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "none"
      )
    
    ggsave(
      filename = file.path(output_dir2, paste0("plot_", sc, ".png")),
      plot = p2,
      width = width,
      height = height,
      dpi = dpi
    )
    
    p3 <- ggplot(df_sc, aes(x = method, y = n_event, fill = method)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.9) +
      labs(
        title = paste("Scenario:", sc),
        x = "Method",
        y = "n_event"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "none"
      )
    
    ggsave(
      filename = file.path(output_dir3, paste0("plot_", sc, ".png")),
      plot = p3,
      width = width,
      height = height,
      dpi = dpi
    )
  }
  
  message("All plots saved to: ", output_dir)
}

plot_exploratory(data = data_explor_converge)

# Agreement between empirical se_gamma_1 and model-based se_gamma_1

se_agreement <- ggplot(data = bias_results, aes(x = mod_se_gamma_1, y = se_gamma_1)) +
  geom_point() + geom_abline(slope = 1, intercept = 0)
se_agreement

##############################
# 5 - Visualising results
##############################

#bias_results <- import(here("Results", "Summary", "Summary_all_scens_2025_12_02.xlsx"))

bias_results$methods <- factor(bias_results$methods, levels = c("ignore_C", "exclude_C", "censor_C", "adjust_C"))

# For plotting IRR
bias_results <- bias_results %>% 
  mutate(log_IRR_E_hat = log(IRR_E_hat),
         log_IRR_E_low_CI = log(IRR_E_hat_low_CI),
         log_IRR_E_up_CI = log(IRR_E_hat_up_CI),
         exp_gamma_1_hat = exp(gamma_1_hat),
         exp_gamma_1_hat_low_CI = exp(gamma_1_hat_low_CI),
         exp_gamma_1_hat_up_CI = exp(gamma_1_hat_up_CI),
         bias_exp_gamma_1_hat = exp(gamma_1_hat) - IRR_E,
         )

bias_results1 <- bias_results %>% mutate(across(c(p_C,ConE,IRR_E,risk_win_C, UonY), as.factor))

bias_results1 <- bias_results1 %>%
  mutate(
    p_C_num = as.numeric(p_C),
    method_offset = as.numeric(methods) * 0.15 - 0.30,   # Manual dodge: create offsets for each method
    risk_win_C_labs = factor(risk_win_C, 
                             levels = c("10", "30", "60"), 
                             labels = c("Risk window C: 10 days", "Risk window C: 30 days", "Risk window C: 60 days")),
    y_dodged = p_C_num + method_offset
  )

# 5.1. No Confounding effect of C and U ----------------------------------------------------------------

# Expectation: downward bias for ignore_C due to mediating effect of C.
# No bias for other methods

results_no_UC1 <- bias_results1 %>% filter(UonY == 1 & ConE == 1 & IRR_E == 1)
results_no_UC2 <- bias_results1 %>% filter(UonY == 1 & ConE == 1 & IRR_E == 2)

## Lollipop plot for bias/ estimated IRR -----------------------------

lollipop_plot <- function(data,
                          aes_x, aes_x_low_ci, aes_x_up_ci,
                          mode = c("gamma", "irr"),   # "gamma" = linear scale, "irr" = log scale
                          refline = 0,               # 0 for bias, 1 or 2 for IRR
                          irr_breaks = NULL,         # vector of IRR ticks, e.g. c(0.5,1,2,4)
                          irr_limits = NULL,         # two-element vector, IRR scale
                          x_break = NULL,            # for gamma mode
                          x_limits = NULL,           # for gamma mode
                          xlabel = "",
                          plot_name) {
  
  mode <- match.arg(mode)
  
  png(here("Plots", paste0(plot_name, ".png")), 
      width = 18, height = 10, units = "cm", res = 300)
  
  p <- ggplot(
    data,
    aes(x = .data[[aes_x]], y = y_dodged, color = methods)
  ) +
    geom_segment(aes(
      x = if (mode == "gamma") refline else log(refline),
      xend = .data[[aes_x]],
      y = y_dodged,
      yend = y_dodged
    ),
    linewidth = 0.6) +
    
    geom_point(size = 1.5) +
    
    geom_errorbarh(
      aes(xmin = .data[[aes_x_low_ci]], xmax = .data[[aes_x_up_ci]]),
      height = 0.1, alpha = 0.6
    ) +
    
    scale_y_continuous(
      breaks = unique(data$p_C_num),
      labels = levels(data$p_C)
    ) +
    
    scale_color_manual(
      name = "Method",
      labels = c("Ignore C", "Exclude C", "Censor C", "Adjust for C"),
      values = c("#bf3729","#2f357c","#e69b00", "#999999")
    ) +
    
    facet_wrap(~ risk_win_C_labs) +
    theme_bw() +
    theme(
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 12),
      legend.position = "bottom",
      strip.text   = element_text(size = 12)
    ) +
    labs(y = "Daily probability of C", x = xlabel) +
    coord_flip()
  
  # ---- SCALE LOGIC ----------------------------------------------------------
  
  if (mode == "gamma") {
    # linear scale
    p <- p +
      scale_x_continuous(
        breaks = x_break,
        limits = x_limits
      ) +
      geom_vline(xintercept = refline, linewidth = 1.2)
    
  } else if (mode == "irr") {
    # log scale, but axis shows IRR values
    p <- p +
      scale_x_continuous(
        breaks = log(irr_breaks),
        labels = irr_breaks,
        limits = log(irr_limits)
      ) +
      geom_vline(xintercept = log(refline), linewidth = 1.2)
  }
  
  print(p)
  dev.off()
}

lollipop_plot_noerror <- function(data,
                          aes_x,
                          mode = c("gamma", "irr"),   # "gamma" = linear scale, "irr" = log scale
                          refline = 0,               # 0 for bias, 1 or 2 for IRR
                          irr_breaks = NULL,         # vector of IRR ticks, e.g. c(0.5,1,2,4)
                          irr_limits = NULL,         # two-element vector, IRR scale
                          x_break = NULL,            # for gamma mode
                          x_limits = NULL,           # for gamma mode
                          xlabel = "",
                          plot_name) {
  
  mode <- match.arg(mode)
  
  png(here("Plots", paste0(plot_name, ".png")), 
      width = 20, height = 9, units = "cm", res = 300)
  
  p <- ggplot(
    data,
    aes(x = .data[[aes_x]], y = y_dodged, color = methods)
  ) +
    geom_segment(aes(
      x = if (mode == "gamma") refline else log(refline),
      xend = .data[[aes_x]],
      y = y_dodged,
      yend = y_dodged
    ),
    linewidth = 0.6) +
    
    geom_point(size = 1.5) +
    
    scale_y_continuous(
      breaks = unique(data$p_C_num),
      labels = levels(data$p_C)
    ) +
    
    scale_color_manual(
      name = "Method",
      labels = c("Ignore C", "Exclude C", "Censor C", "Adjust for C"),
      values = c("#bf3729","#2f357c","#e69b00", "#999999")
    ) +
    
    facet_wrap(~ risk_win_C_labs) +
    theme_bw() +
    theme(
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 12),
      strip.text   = element_text(size = 12)
    ) +
    labs(y = "Daily probability of C", x = xlabel) +
    coord_flip()
  
  # ---- SCALE LOGIC ----------------------------------------------------------
  
  if (mode == "gamma") {
    # linear scale
    p <- p +
      scale_x_continuous(
        breaks = x_break,
        limits = x_limits
      ) +
      geom_vline(xintercept = refline, linewidth = 1.2)
    
  } else if (mode == "irr") {
    # log scale, but axis shows IRR values
    p <- p +
      scale_x_continuous(
        breaks = log(irr_breaks),
        labels = irr_breaks,
        limits = log(irr_limits)
      ) +
      geom_vline(xintercept = log(refline), linewidth = 1.2)
  }
  
  print(p)
  dev.off()
}


### Plot for bias gamma_1 -------------
lollipop_plot(data = results_no_UC1, aes_x ="bias_gamma_1", 
              aes_x_low_ci ="bias_gamma_1_low_CI", aes_x_up_ci = "bias_gamma_1_up_CI",
              mode = "gamma",
              refline = 0,
              x_break = seq(from = -0.2, to = 0.12, by = 0.02),
              x_limits = c(-0.2, 0.12),
              xlabel = expression("Bias of " * gamma[1] * " (" * log(IRR[E]) * ")"),
              plot_name = "bias_gamma1_no_UC_IRRE1")
lollipop_plot(data = results_no_UC2, aes_x ="bias_gamma_1", 
              aes_x_low_ci ="bias_gamma_1_low_CI", aes_x_up_ci = "bias_gamma_1_up_CI",
              mode = "gamma",
              refline = 0,
              x_break = seq(from = -0.2, to = 0.12, by = 0.02),
              x_limits = c(-0.2, 0.12),
              xlabel = expression("Bias of " * gamma[1] * " (" * log(IRR[E]) * ")"),
              plot_name = "bias_gamma1_no_UC_IRRE2")

### Plot for bias IRR_E -----
lollipop_plot(data = results_no_UC1, aes_x ="bias_IRR_E", 
              aes_x_low_ci ="bias_IRR_E_low_CI", aes_x_up_ci = "bias_IRR_E_up_CI",
              mode = "gamma",
              refline = 0,
              x_break = seq(from = -0.2, to = 0.12, by = 0.02),
              x_limits = c(-0.2, 0.12),
              xlabel  = expression("Bias of " * IRR[E]),
              plot_name = "bias_IRRE_no_UC_IRRE1")
lollipop_plot(data = results_no_UC2, aes_x ="bias_IRR_E", 
              aes_x_low_ci ="bias_IRR_E_low_CI", aes_x_up_ci = "bias_IRR_E_up_CI",
              mode = "gamma",
              refline = 0,
              x_break = seq(from = -0.2, to = 0.12, by = 0.02),
              x_limits = c(-0.2, 0.12),
              xlabel <- expression("Bias of " * IRR[E]),
              plot_name = "bias_IRRE_no_UC_IRRE2")

### Plot for bias exp(gamma_1_hat)
lollipop_plot_noerror(data = results_no_UC1, aes_x ="bias_exp_gamma_1_hat", 
              mode = "gamma",
              refline = 0,
              x_break = seq(from = -0.14, to = 0, by = 0.02),
              x_limits = c(-0.14, 0),
              xlabel = expression("Bias of " * hat(IRR[E])),
              plot_name = "bias_IRRE_expg1_no_UC_IRRE1")

lollipop_plot_noerror(data = results_no_UC2, aes_x ="bias_exp_gamma_1_hat", 
              mode = "gamma",
              refline = 0,
              x_break = seq(from = -0.2, to = 0, by = 0.02),
              x_limits = c(-0.2, 0),
              xlabel = expression("Bias of " * hat(IRR[E])),
              plot_name = "bias_IRRE_expg1_no_UC_IRRE2")

### Relative bias ----
lollipop_plot(data = results_no_UC2, aes_x ="relative_bias_g1", 
              aes_x_low_ci ="relative_bias_g1_low_CI", aes_x_up_ci = "relative_bias_g1_up_CI",
              mode = "gamma",
              refline = 0,
              x_break = seq(from = -0.18, to = 0.1, by = 0.02),
              x_limits = c(-0.19, 0.1),
              xlabel = expression("Relative bias of " * gamma[1]),
              plot_name = "relative_bias_g1_no_UC_IRRE2")

### Plot for IRR -----------------------------

lollipop_plot(
  data       = results_no_UC1,
  aes_x      = "log_IRR_E_hat",
  aes_x_low_ci = "log_IRR_E_low_CI",
  aes_x_up_ci  = "log_IRR_E_up_CI",
  mode       = "irr",
  refline    = 1,                     
  irr_breaks = seq(0.8, 1.12, by = 0.04),
  irr_limits = c(0.8, 1.12),
  xlabel     = expression("Estimated " * IRR[E]),
  plot_name  = "mean_IRRE_no_UC_IRRE1"
)

lollipop_plot(
  data       = results_no_UC2,
  aes_x      = "log_IRR_E_hat",
  aes_x_low_ci = "log_IRR_E_low_CI",
  aes_x_up_ci  = "log_IRR_E_up_CI",
  mode       = "irr",
  refline    = 2,                     
  irr_breaks = seq(1.8, 2.12, by = 0.04),
  irr_limits = c(1.8, 2.12),
  xlabel     = expression("Estimated " * IRR[E]),
  plot_name  = "mean_IRRE_no_UC_IRRE2"
)

# IRR corresponding to the average gamma_1
lollipop_plot(
  data       = results_no_UC1,
  aes_x      = "gamma_1_hat",
  aes_x_low_ci = "gamma_1_hat_low_CI",
  aes_x_up_ci  = "gamma_1_hat_up_CI",
  mode       = "irr",
  refline    = 1,                     
  irr_breaks = seq(0.8, 1.12, by = 0.04),
  irr_limits = c(0.8, 1.12),
  xlabel     = expression("Estimated " * IRR[E] *" corresponding to " * bar(gamma[1])),
  plot_name  = "mean_IRRE_expg1_no_UC_IRRE1"
)

lollipop_plot(
  data       = results_no_UC2,
  aes_x      = "gamma_1_hat",
  aes_x_low_ci = "gamma_1_hat_low_CI",
  aes_x_up_ci  = "gamma_1_hat_up_CI",
  mode       = "irr",
  refline    = 2,                     
  irr_breaks = seq(1.76, 2.12, by = 0.04),
  irr_limits = c(1.76, 2.12),
  xlabel     = expression("Estimated " * IRR[E] *" corresponding to " * bar(gamma[1])),
  plot_name  = "mean_IRRE_expg1_no_UC_IRRE2"
)


## Plot for variance ------------------------

lollipop_plot_var <- function(data, aes_x, aes_x_low_ci, aes_x_up_ci, plot_name){
  
  png(here("Plots", paste0(plot_name, ".png")), width = 18, height = 10, units = "cm", res = 300)
  
  p <- ggplot(data,
              aes(x = .data[[aes_x]],
                  y = y_dodged,
                  color = methods)) +
    
    geom_point(size = 1.5) +
    
    geom_errorbarh(aes(xmin = .data[[aes_x_low_ci]],
                       xmax = .data[[aes_x_up_ci]]),
                   width = 0.1,
                   alpha = 0.4) +
    
    scale_y_continuous(
      breaks = unique(data$p_C_num),
      labels = levels(data$p_C)
    ) +
    
    scale_x_continuous(
      limits = c(0, 0.8)) +
    
    scale_color_manual(
      name = "Method",
      labels = c("Ignore C", "Exclude C", "Censor C", "Adjust for C"),
      values = c("#bf3729","#2f357c","#e69b00", "#ada43b")
    ) +
    
    facet_wrap(~ risk_win_C_labs) +
    theme_bw() +
    theme(
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 12),
      legend.position = "bottom",
      strip.text   = element_text(size = 12)
    ) +
    labs(y = "Daily probability of C", x = expression("Empirical standard error of " * hat(gamma)[1])) +
    coord_flip()
  
  print(p)
  dev.off()
}

lollipop_plot_var(data = results_no_UC1, aes_x ="se_gamma_1", 
              aes_x_low_ci ="se_gamma_1_low_CI", aes_x_up_ci = "se_gamma_1_up_CI",
              plot_name = "empirical_se_gamma_1_no_UC_IRRE1")
lollipop_plot_var(data = results_no_UC2, aes_x ="se_gamma_1", 
                  aes_x_low_ci ="se_gamma_1_low_CI", aes_x_up_ci = "se_gamma_1_up_CI",
                  plot_name = "empirical_se_gamma_1_no_UC_IRRE2")

## Plot for mean number of events ----------------------------

nr_event_plot <- function(data, aes_x, plot_name){
  
  png(here("Plots", paste0(plot_name, ".png")), width = 18, height = 10, units = "cm", res = 300)
  
  p <- ggplot(data,
              aes(y = .data[[aes_x]],
                  x = y_dodged,
                  fill = methods)) +
    
    geom_bar(stat = "identity") +
    
    scale_x_continuous(
      breaks = unique(data$p_C_num),
      labels = levels(data$p_C)
    ) +
    
    scale_y_continuous(breaks = seq(from = 0, to = 2000, by = 200)) +
    
    scale_fill_manual(
      name = "Method",
      labels = c("Ignore C", "Exclude C", "Censor C", "Adjust for C"),
      values = c("#bf3729","#2f357c","#e69b00", "#ada43b")
    ) +
    
    facet_wrap(~ risk_win_C_labs) +
    theme_bw() +
    theme(
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 12),
      legend.position = "bottom",
      strip.text   = element_text(size = 12)
    ) +
    labs(x = "Daily probability of C", y = "Mean number of events") #+
#    coord_flip()
  
  print(p)
  dev.off()
}

nr_event_plot(data = results_no_UC2, aes_x ="mean_nr_event", 
                  plot_name = "nr_event_no_UC_IRRE2")

### Coverage ----



# 5.2. No Confounding effect of U, confounding effect of C presents----------------------------------------

# Expectation: Bias when ignoring C: combination of confounding and mediating path E-C-Y
# No bias for other methods

# Reshape data to wide format
method_cols <- setdiff(
  names(bias_results),
  c("methods", "scen_name", "p_C", "ConE", "IRR_E", "risk_win_C",
    "UonY", "scen_id"))
bias_results_wide <- pivot_wider(data = bias_results, 
                                 names_from = methods, values_from = all_of(method_cols), 
                                 names_vary = "fastest")

# Filter data of scenarios U->Y = 1
results_noU1 <- bias_results_wide %>% filter(UonY == 1 & IRR_E == 1)
results_noU2 <- bias_results_wide %>% filter(UonY == 1 & IRR_E == 2)

## IRR_E = 1 -------------------------------
# Reorder the parameters for nested loop plot
results_noU1_reorder <- nestedloop(results_noU1, 
                                   varnames = c("ConE", "p_C", "risk_win_C"), #order of the loop
                                   varlabels = c("Effect of C on E (OR)", 
                                                 "Daily probability of C",  
                                                 "Risk window of C (days)"),
                                   sign = c(1, 1, 1)) # Order of the levels of each parameter, 
# 1 = increasing, -1 = decreasing

# Create nested loop plot
parOrig <- par()
nested_loop_plot <- function(data,
                             yvars,
                             ylab = expression("Bias"),
                             outfile = NULL,
                             refline = 0,
                             ylim,
                             ymax_ref,
                             y_log = FALSE,
                             ybreaks,
                             x_lab = "4 × 3 × 3 = 36 ordered scenarios",
                             legend_labels = c("Ignore C", "Exclude C", "Censor C", "Adjust for C"),
                             colours = c("#bf3729", "#2f357c", "#e69b00", "#999999"),
                             legend_pos = c(1, 0.55),
                             plot_width = 20,
                             plot_height = 12) {
  
  # Colors
  cols <- colours
  if (length(cols) < length(yvars)) {
    cols <- rep(cols, length.out = length(yvars))
  }
  
  # Open PNG if requested
  if (!is.null(outfile)) {
    png(outfile, width = plot_width, height = plot_height, units = "cm", res = 300)
  }
  
  par(pty = "m", bg=NA)
  ## Skeleton plot
  
  if (y_log == FALSE) {
    plot(data$p_C,
         type = "n",
         ylim = ylim,
         bty = "n",
         xlab = x_lab,
         ylab = ylab,
         las = 1,
         xaxt = "n",
         yaxt = "n")
  }
  else {
    plot(data$IRR_E,
         log=c("y"),
         type = "n",
         ylim = ylim,
         bty = "n",
         xlab = x_lab,
         ylab = ylab,
         las = 1,
         xaxt = "n",
         yaxt = "n")
  }

  
  axis(side = 2, at = round(ybreaks, 1))
  
  ## Background nested-loop grid
  lines.nestedloop(data, which = "v",
                   col = c("#454545", "#BBBBBB", "#DDDDDD", "#EEEEEE"))
  
  lines.nestedloop(data, which = "r",
                   ymin.refline = ylim[1],
                   ymax.refline = ymax_ref,
                   log = y_log)
  
  ## Horizontal reference line at 0
  abline(h = refline, col = "black", lwd = 2)
  
  ## Plot the bias for different methods
  for (i in seq_along(yvars)) {
    varname <- yvars[i]
    lines(data[[varname]], col = cols[i], lwd = 1.5, type = "s")
  }
  
  ## Legend
  legend(legend_pos[1], legend_pos[2],
         lwd = rep(1.5, length(yvars)),
         col = cols[seq_along(yvars)],
         cex = 0.9,
         bty = "n",
         legend = legend_labels)
  
  if (!is.null(outfile)) dev.off()
}
### Bias gamma
nested_loop_plot(
  data = results_noU1_reorder,
  yvars = c("bias_gamma_1_ignore_C",
            "bias_gamma_1_exclude_C",
            "bias_gamma_1_censor_C",
            "bias_gamma_1_adjust_C"),
  ylab = expression("Bias of " * gamma[1]),
  refline = 0,
  ylim = c(-0.7, 0.5),
  ymax_ref = -0.3,
  y_log = FALSE,
  ybreaks = seq(-0.6, 0.5, by = 0.05),
  legend_pos = c(1, 0.5),
  outfile = here("Plots", "bias_gamma1_no_U_IRRE1_loop.png"),
  plot_height = 14
)

### IRR
nested_loop_plot(
  data = results_noU1_reorder,
  yvars = c("IRR_E_hat_ignore_C",
            "IRR_E_hat_exclude_C",
            "IRR_E_hat_censor_C",
            "IRR_E_hat_adjust_C"),
  ylab = expression(hat(IRR[E])),
  refline = 1,
  ylim = c(0.4, 2),
  ymax_ref = 0.7,
  y_log = TRUE,
  ybreaks = seq(0.4, 2, by = 0.1),
  legend_pos = c(1, 2.2),
  outfile = here("Plots", "mean_IRRE_no_U_IRRE1_loop.png"),
  plot_height = 14
)

### IRR corresponding to average gamma_1
nested_loop_plot(
  data = results_noU1_reorder,
  yvars = c("exp_gamma_1_hat_ignore_C",
            "exp_gamma_1_hat_exclude_C",
            "exp_gamma_1_hat_censor_C",
            "exp_gamma_1_hat_adjust_C"),
  ylab = expression("Estimated " * IRR[E] *" corresponding to " * bar(gamma[1])),
  refline = 1,
  ylim = c(0.4, 2),
  ymax_ref = 0.7,
  y_log = TRUE,
  ybreaks = seq(0.4, 2, by = 0.1),
  legend_pos = c(1, 2.2),
  outfile = here("Plots", "mean_IRRE_expg1_no_U_IRRE1_loop.png"),
  plot_height = 14
)

### Variance 
nested_loop_plot(
  data = results_noU1_reorder,
  yvars = c("se_gamma_1_ignore_C",
            "se_gamma_1_exclude_C",
            "se_gamma_1_censor_C",
            "se_gamma_1_adjust_C"),
  ylab = expression("Empirical SE of "*hat(gamma[1])),
  refline = 0,
  ylim = c(-0.6, 0.8),
  ymax_ref = -0.1,
  y_log = FALSE,
  ybreaks = seq(0, 0.6, by = 0.1),
  legend_pos = c(1, 0.9),
  outfile = here("Plots", "empirical_se_gamma_1_no_U_IRRE1_loop.png"),
  plot_height = 14
)

## IRR_E = 2 ------------------------
results_noU2_reorder <- nestedloop(results_noU2, 
                                   varnames = c("ConE", "p_C", "risk_win_C"), #order of the loop
                                   varlabels = c("Effect of C on E (OR)", 
                                                 "Daily probability of C",  
                                                 "Risk window of C (days)"),
                                   sign = c(1, 1, 1)) # Order of the levels of each parameter, 
### Bias gamma
nested_loop_plot(
  data = results_noU2_reorder,
  yvars = c("bias_gamma_1_ignore_C",
            "bias_gamma_1_exclude_C",
            "bias_gamma_1_censor_C",
            "bias_gamma_1_adjust_C"),
  ylab = expression("Bias of " * gamma[1]),
  refline = 0,
  ylim = c(-0.7, 0.5),
  ymax_ref = -0.3,
  y_log = FALSE,
  ybreaks = seq(-0.6, 0.5, by = 0.05),
  legend_pos = c(1, 0.5),
  outfile = here("Plots", "bias_gamma1_no_U_IRRE2_loop.png"),
  plot_height = 14
)

### Relative bias
nested_loop_plot(
  data = results_noU2_reorder,
  yvars = c("relative_bias_g1_ignore_C",
            "relative_bias_g1_exclude_C",
            "relative_bias_g1_censor_C",
            "relative_bias_g1_adjust_C"),
  ylab = expression("Relative bias of " * gamma[1]),
  refline = 0,
  ylim = c(-0.9, 0.6),
  ymax_ref = -0.5,
  y_log = FALSE,
  ybreaks = seq(-0.6, 0.6, by = 0.05),
  legend_pos = c(1, 0.6),
  outfile = here("Plots", "relative_bias_gamma1_no_U_IRRE2_loop.png"),
  plot_height = 14
)

### IRR
nested_loop_plot(
  data = results_noU2_reorder,
  yvars = c("IRR_E_hat_ignore_C",
            "IRR_E_hat_exclude_C",
            "IRR_E_hat_censor_C",
            "IRR_E_hat_adjust_C"),
  ylab = expression(hat(IRR[E])),
  refline = 2,
  ylim = c(1.2, 2.6),
  ymax_ref = 1.5,
  y_log = TRUE,
  ybreaks = seq(1.4, 2.6, by = 0.1),
  legend_pos = c(1, 2.6),
  outfile = here("Plots", "mean_IRRE_no_U_IRRE2_loop.png"),
  plot_height = 14
)

### IRR corresponding to average gamma_1
nested_loop_plot(
  data = results_noU2_reorder,
  yvars = c("exp_gamma_1_hat_ignore_C",
            "exp_gamma_1_hat_exclude_C",
            "exp_gamma_1_hat_censor_C",
            "exp_gamma_1_hat_adjust_C"),
  ylab = expression("Estimated " * IRR[E] *" corresponding to " * bar(gamma[1])),
  refline = 2,
  ylim = c(1.2, 3),
  ymax_ref = 1.5,
  y_log = TRUE,
  ybreaks = seq(1.4, 3, by = 0.1),
  legend_pos = c(1, 3),
  outfile = here("Plots", "mean_IRRE_expg1_no_U_IRRE2_loop.png"),
  plot_height = 14
)

### Variance 
nested_loop_plot(
  data = results_noU2_reorder,
  yvars = c("se_gamma_1_ignore_C",
            "se_gamma_1_exclude_C",
            "se_gamma_1_censor_C",
            "se_gamma_1_adjust_C"),
  ylab = expression("Empirical SE of "*hat(gamma[1])),
  refline = 0,
  ylim = c(-0.6, 0.8),
  ymax_ref = -0.1,
  y_log = FALSE,
  ybreaks = seq(0, 0.6, by = 0.1),
  legend_pos = c(1, 0.9),
  outfile = here("Plots", "empirical_se_gamma_1_no_U_IRRE2_loop.png"),
  plot_height = 14
)

# 5.3. All scenarios ----------------------------------------------------------------------------

## Filter scenarios IRR_E = 1 and 2 
bias_results_wide_IRRE1 <-  bias_results_wide %>% filter(IRR_E == 1)
bias_results_wide_IRRE2 <-  bias_results_wide %>% filter(IRR_E == 2)

## Reorder the parameters for nested loop plot

bias_results_wide_IRRE1_reorder <- nestedloop(bias_results_wide_IRRE1, 
                                   varnames = c("UonY", "ConE", "p_C", "risk_win_C"), #order of the loop
                                   varlabels = c("Effect of U on Y (OR)",
                                                 "Effect of C on E (OR)", 
                                                 "Daily probability of C",  
                                                 "Risk window of C (days)"),
                                   sign = c(1, 1, 1, 1)) 

bias_results_wide_IRRE2_reorder <- nestedloop(bias_results_wide_IRRE2, 
                                              varnames = c("UonY", "ConE", "p_C", "risk_win_C"), #order of the loop
                                              varlabels = c("Effect of U on Y (OR)",
                                                            "Effect of C on E (OR)", 
                                                            "Daily probability of C",  
                                                            "Risk window of C (days)"),
                                              sign = c(1, 1, 1, 1)) 
## Create nested loop plot

### True IRR_E = 1 --------------------------
#### Bias gamma
nested_loop_plot(
  data = bias_results_wide_IRRE1_reorder,
  yvars = c("bias_gamma_1_ignore_C",
            "bias_gamma_1_exclude_C",
            "bias_gamma_1_censor_C",
            "bias_gamma_1_adjust_C"),
  ylab = expression("Bias of " * gamma[1]),
  refline = 0,
  ylim = c(-1, 0.5),
  ymax_ref = -0.4,
  y_log = FALSE,
  ybreaks = seq(-1, 0.5, by = 0.05),
  x_lab = "4 x 4 x 3 x 3 = 144 ordered scenarios",
  legend_pos = c(1, 0.6),
  outfile = here("Plots", "bias_gamma1_all_scens_IRRE1_loop.png"),
  plot_height = 13
)

#### IRR
nested_loop_plot(
  data = bias_results_wide_IRRE1_reorder,
  yvars = c("IRR_E_hat_ignore_C",
            "IRR_E_hat_exclude_C",
            "IRR_E_hat_censor_C",
            "IRR_E_hat_adjust_C"),
  ylab = expression(hat(IRR[E])),
  refline = 1,
  ylim = c(0.4, 2),
  ymax_ref = 0.7,
  y_log = TRUE,
  ybreaks = seq(0.4, 2, by = 0.1),
  x_lab = "4 x 4 x 3 x 3 = 144 ordered scenarios",
  legend_pos = c(1, 2.2),
  outfile = here("Plots", "mean_IRRE_all_scens_IRRE1_loop.png"),
  plot_height = 12
)

### IRR corresponding to average gamma_1
nested_loop_plot(
  data = bias_results_wide_IRRE1_reorder,
  yvars = c("exp_gamma_1_hat_ignore_C",
            "exp_gamma_1_hat_exclude_C",
            "exp_gamma_1_hat_censor_C",
            "exp_gamma_1_hat_adjust_C"),
  ylab = expression("Estimated " * IRR[E] *" corresponding to " * bar(gamma[1])),
  refline = 1,
  ylim = c(0.4, 2),
  ymax_ref = 0.7,
  y_log = TRUE,
  ybreaks = seq(0.4, 2, by = 0.1),
  x_lab = "4 x 4 x 3 x 3 = 144 ordered scenarios",
  legend_pos = c(1, 2.2),
  outfile = here("Plots", "mean_IRRE_expg1_all_scens_IRRE1_loop.png"),
  plot_height = 14
)

#### Mean squared error
nested_loop_plot(
  data = bias_results_wide_IRRE1_reorder,
  yvars = c("mse_g1_ignore_C",
            "mse_g1_exclude_C",
            "mse_g1_censor_C",
            "mse_g1_adjust_C"),
  ylab = expression("Mean squared error of "*gamma[1]),
  refline = 0,
  ylim = c(-0.6, 0.7),
  ymax_ref = -0.1,
  y_log = FALSE,
  ybreaks = seq(0, 0.4, by = 0.01),
  x_lab = "4 x 4 x 3 x 3 = 144 ordered scenarios",
  legend_pos = c(1, 0.8),
  outfile = here("Plots", "mse_all_scens_IRRE1_loop.png"),
  plot_height = 14
)


#### Coverage

nested_loop_plot(
  data = bias_results_wide_IRRE1_reorder,
  yvars = c("coverage_irr_E_ignore_C",
            "coverage_irr_E_exclude_C",
            "coverage_irr_E_censor_C",
            "coverage_irr_E_adjust_C"),
  ylab = expression("Coverage of 95% CI of " * hat(IRR[E])),
  refline = 0.95,
  ylim = c(0, 1.2),
  ymax_ref = 0.5,
  y_log = FALSE,
  ybreaks = seq(0, 1, by = 0.1),
  x_lab = "4 x 4 x 3 x 3 = 144 ordered scenarios",
  legend_pos = c(1, 1.3),
  outfile = here("Plots", "coverage_all_scens_IRRE1_loop.png"),
  plot_height = 13
)

 

### True IRR_E = 2 --------------------------
#### Bias gamma ----
nested_loop_plot(
  data = bias_results_wide_IRRE2_reorder,
  yvars = c("bias_gamma_1_ignore_C",
            "bias_gamma_1_exclude_C",
            "bias_gamma_1_censor_C",
            "bias_gamma_1_adjust_C"),
  ylab = expression("Bias of " * gamma[1]),
  refline = 0,
  ylim = c(-1, 0.5),
  ymax_ref = -0.4,
  y_log = FALSE,
  ybreaks = seq(-1, 0.5, by = 0.05),
  x_lab = "4 x 4 x 3 x 3 = 144 ordered scenarios",
  legend_pos = c(1, 0.6),
  outfile = here("Plots", "bias_gamma1_all_scens_IRRE2_loop.png"),
  plot_height = 13
)

#### Relative bias ----
nested_loop_plot(
  data = bias_results_wide_IRRE2_reorder,
  yvars = c("relative_bias_g1_ignore_C",
            "relative_bias_g1_exclude_C",
            "relative_bias_g1_censor_C",
            "relative_bias_g1_adjust_C"),
  ylab = expression("Relative bias of " * gamma[1]),
  refline = 0,
  ylim = c(-1.5, 0.7),
  ymax_ref = -0.6,
  y_log = FALSE,
  ybreaks = seq(-1, 0.7, by = 0.05),
  x_lab = "4 x 4 x 3 x 3 = 144 ordered scenarios",
  legend_pos = c(1, 0.8),
  outfile = here("Plots", "relative_bias_gamma1_all_scens_IRRE2_loop.png"),
  plot_height = 14
)


#### IRR ----
nested_loop_plot(
  data = bias_results_wide_IRRE2_reorder,
  yvars = c("IRR_E_hat_ignore_C",
            "IRR_E_hat_exclude_C",
            "IRR_E_hat_censor_C",
            "IRR_E_hat_adjust_C"),
  ylab = expression(hat(IRR[E])),
  refline = 2,
  ylim = c(1, 3.2),
  ymax_ref = 1.5,
  y_log = TRUE,
  ybreaks = seq(1, 3.2, by = 0.1),
  x_lab = "4 x 4 x 3 x 3 = 144 ordered scenarios",
  legend_pos = c(1, 3.4),
  outfile = here("Plots", "mean_IRRE_all_scens_IRRE2_loop.png"),
  plot_height = 14
)

### IRR corresponding to average gamma_1 ----
nested_loop_plot(
  data = bias_results_wide_IRRE2_reorder,
  yvars = c("exp_gamma_1_hat_ignore_C",
            "exp_gamma_1_hat_exclude_C",
            "exp_gamma_1_hat_censor_C",
            "exp_gamma_1_hat_adjust_C"),
  ylab = expression("Estimated " * IRR[E] *" corresponding to " * bar(gamma[1])),
  refline = 2,
  ylim = c(1, 3.2),
  ymax_ref = 1.5,
  y_log = TRUE,
  ybreaks = seq(1, 3.2, by = 0.1),
  x_lab = "4 x 4 x 3 x 3 = 144 ordered scenarios",
  legend_pos = c(1, 3.4),
  outfile = here("Plots", "mean_IRRE_expg1_all_scens_IRRE2_loop.png"),
  plot_height = 14
)

#### Mean squared error
nested_loop_plot(
  data = bias_results_wide_IRRE2_reorder,
  yvars = c("mse_g1_ignore_C",
            "mse_g1_exclude_C",
            "mse_g1_censor_C",
            "mse_g1_adjust_C"),
  ylab = expression("Mean squared error of "*gamma[1]),
  refline = 0,
  ylim = c(-0.6, 0.7),
  ymax_ref = -0.1,
  y_log = FALSE,
  ybreaks = seq(0, 0.4, by = 0.01),
  x_lab = "4 x 4 x 3 x 3 = 144 ordered scenarios",
  legend_pos = c(1, 0.8),
  outfile = here("Plots", "mse_all_scens_IRRE2_loop.png"),
  plot_height = 14
)

#### Coverage ----

nested_loop_plot(
  data = bias_results_wide_IRRE2_reorder,
  yvars = c("coverage_irr_E_ignore_C",
            "coverage_irr_E_exclude_C",
            "coverage_irr_E_censor_C",
            "coverage_irr_E_adjust_C"),
  ylab = expression("Coverage of 95% CI of " * hat(IRR[E])),
  refline = 0.95,
  ylim = c(0, 1.2),
  ymax_ref = 0.5,
  y_log = FALSE,
  ybreaks = seq(0, 1, by = 0.1),
  x_lab = "4 x 4 x 3 x 3 = 144 ordered scenarios",
  legend_pos = c(1, 1.3),
  outfile = here("Plots", "coverage_all_scens_IRRE2_loop.png"),
  plot_height = 13
)


### Number of events ----------------------------
nested_loop_plot(
  data = bias_results_wide_IRRE1_reorder,
  yvars = c("mean_nr_event_ignore_C",
            "mean_nr_event_exclude_C",
            "mean_nr_event_censor_C",
            "mean_nr_event_adjust_C"),
  ylab = "Mean number of events",
  refline = NULL,
  ylim = c(-2000, 5600),
  ymax_ref = -100,
  y_log = FALSE,
  ybreaks = seq(0, 5600, by = 500),
  x_lab = "4 x 4 x 3 x 3 = 144 ordered scenarios",
  legend_pos = c(1, 5600),
  outfile = here("Plots", "nr_event_all_scens_IRRE1_loop.png"),
  plot_height = 15,
  plot_width = 25
)

nested_loop_plot(
  data = bias_results_wide_IRRE2_reorder,
  yvars = c("mean_nr_event_ignore_C",
            "mean_nr_event_exclude_C",
            "mean_nr_event_censor_C",
            "mean_nr_event_adjust_C"),
  ylab = "Mean number of events",
  refline = NULL,
  ylim = c(-2000, 5600),
  ymax_ref = -100,
  y_log = FALSE,
  ybreaks = seq(0, 5600, by = 200),
  x_lab = "4 x 4 x 3 x 3 = 144 ordered scenarios",
  legend_pos = c(1, 5600),
  outfile = here("Plots", "nr_event_all_scens_IRRE2_loop.png"),
  plot_height = 15,
  plot_width = 25
)

# 5.4. Zoom in UonY = 15 --------------------------------------------------
bias_results_wide_UY1_15_IRRE1_reorder <- bias_results_wide_IRRE1_reorder %>% 
  filter (UonY %in% c(1, 15))

bias_results_wide_UY1_15_IRRE2_reorder <- bias_results_wide_IRRE2_reorder %>% 
  filter (UonY %in% c(1, 15))

## IRR
nested_loop_plot(
  data = bias_results_wide_UY1_15_IRRE1_reorder,
  yvars = c("IRR_E_hat_exclude_C",
            "IRR_E_hat_censor_C",
            "IRR_E_hat_adjust_C"),
  ylab = expression(hat(IRR[E])),
  refline = 1,
  ylim = c(0.4, 1.8),
  ymax_ref = 0.8,
  y_log = TRUE,
  ybreaks = seq(0.4, 1.8, by = 0.1),
  x_lab = "2x 4 x 3 x 3 = 36 ordered scenarios",
  legend_labels = c("Exclude C", "Censor C", "Adjust for C"),
  colours = c("#2f357c", "#e69b00", "#999999"),
  legend_pos = c(1, 2),
  outfile = here("Plots", "mean_IRRE_UY1_15_IRRE1_loop.png"),
  plot_height = 14
)

nested_loop_plot(
  data = bias_results_wide_UY1_15_IRRE2_reorder,
  yvars = c("IRR_E_hat_exclude_C",
            "IRR_E_hat_censor_C",
            "IRR_E_hat_adjust_C"),
  ylab = expression(hat(IRR[E])),
  refline = 2,
  ylim = c(1.2, 2.8),
  ymax_ref = 1.8,
  y_log = TRUE,
  ybreaks = seq(1.2, 2.8, by = 0.1),
  x_lab = "2x 4 x 3 x 3 = 36 ordered scenarios",
  legend_labels = c("Exclude C", "Censor C", "Adjust for C"),
  colours = c("#2f357c", "#e69b00", "#999999"),
  legend_pos = c(1, 2.9),
  outfile = here("Plots", "mean_IRRE_UY1_15_IRRE2_loop.png"),
  plot_height = 14
)

## Coverage
nested_loop_plot(
  data = bias_results_wide_UY1_15_IRRE1_reorder,
  yvars = c("coverage_irr_E_exclude_C",
            "coverage_irr_E_censor_C",
            "coverage_irr_E_adjust_C"),
  ylab = expression("Coverage of 95% CI of " * hat(IRR[E])),
  refline = 0.95,
  ylim = c(0.4, 1.2),
  ymax_ref = 0.8,
  y_log = FALSE,
  ybreaks = seq(0.4, 1, by = 0.1),
  x_lab = "2x 4 x 3 x 3 = 36 ordered scenarios",
  legend_labels = c("Exclude C", "Censor C", "Adjust for C"),
  colours = c("#2f357c", "#e69b00", "#999999"),  
  legend_pos = c(1, 1.2),
  outfile = here("Plots", "coverage_UY1_15_IRRE1_loop.png"),
  plot_height = 13
)

nested_loop_plot(
  data = bias_results_wide_UY1_15_IRRE2_reorder,
  yvars = c("coverage_irr_E_exclude_C",
            "coverage_irr_E_censor_C",
            "coverage_irr_E_adjust_C"),
  ylab = expression("Coverage of 95% CI of " * hat(IRR[E])),
  refline = 0.95,
  ylim = c(0.4, 1.2),
  ymax_ref = 0.8,
  y_log = FALSE,
  ybreaks = seq(0.4, 1, by = 0.1),
  x_lab = "2x 4 x 3 x 3 = 36 ordered scenarios",
  legend_labels = c("Exclude C", "Censor C", "Adjust for C"),
  colours = c("#2f357c", "#e69b00", "#999999"), 
  legend_pos = c(1, 1.2),
  outfile = here("Plots", "coverage_UY1_15_IRRE2_loop.png"),
  plot_height = 13
)
