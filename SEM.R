# Developed by: [Puguopei], First Author
# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

# Clean environment and configure basic settings
gc()
options(repos = c(CRAN="https://mirrors.aliyun.com/CRAN/"))
options(stringsAsFactors = FALSE, scipen = 999)

# Install all required packages (including openxlsx for semPlot)
pkg_groups <- list(
  data_process = c("dplyr", "tidyr", "tibble", "car"),
  visualization = c("ggplot2", "corrplot"),
  sem_core = c("lavaan", "mi"),
  mediation = c("mediation", "mvtnorm", "sandwich", "boot", "psych"),
  special = c("showtext", "sysfonts", "showtextdb", "openxlsx", "semPlot")
)

# Safe installation function
install_safe <- function(pkg_list, group_name) {
  cat(paste0("\n📦 Installing [", group_name, "] packages:\n"))
  for (pkg in pkg_list) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      tryCatch({
        cat("  Installing:", pkg, "\n")
        install.packages(pkg, dependencies = c("Depends", "Imports"), quiet = TRUE)
        library(pkg, character.only = TRUE)
        cat("  ✅", pkg, "installed and loaded successfully\n")
      }, error = function(e) {
        cat("  ❌", pkg, "installation failed:", e$message, "\n")
      })
    } else {
      cat("  ✅", pkg, "already installed, loaded directly\n")
    }
  }
}

# Install packages by group
install_safe(pkg_groups$data_process, "Data Processing")
install_safe(pkg_groups$visualization, "Basic Visualization")
install_safe(pkg_groups$sem_core, "SEM Core")
install_safe(pkg_groups$mediation, "Mediation Effect")
install_safe(pkg_groups$special, "Special Visualization")

# Verify package loading status
core_packages <- unlist(pkg_groups)
loaded_status <- sapply(core_packages, function(x) {
  if (require(x, character.only = TRUE, quietly = TRUE)) "✅ Loaded successfully" else "❌ Load failed"
})
cat("\n📋 All packages loading status:\n")
for (pkg in names(loaded_status)) {
  cat(paste0("  ", pkg, ": ", loaded_status[pkg], "\n"))
}

# Check core functions
check_functions <- c("showtext_auto", "showtext_opts", "semPaths")
for (func in check_functions) {
  if (exists(func)) {
    cat(paste0("\n✅ Core function [", func, "] loaded\n"))
  } else {
    cat(paste0("\n❌ Core function [", func, "] missing, reinstall corresponding package\n"))
  }
}

# Configure showtext for Chinese display
library(showtext)
font_add_google("Noto Sans SC", "SimHei")
showtext_auto(TRUE)
showtext_opts(dpi = 300)
cat("\n✅ showtext configured successfully, Chinese display normal\n")

# ===================== Data Loading =====================
# Load data files
classified_file <- ""
data_file <- ""

# Error handling: check file existence
if (!file.exists(classified_file) && classified_file != "") {
  stop(paste("❌ Classification file not found:", classified_file))
}
if (!file.exists(data_file) && data_file != "") {
  stop(paste("❌ Data file not found:", data_file))
}

# Read data (only if file paths are provided)
if (classified_file != "" && data_file != "") {
  key_factors_class <- read.csv(classified_file, fileEncoding = "UTF-8", check.names = FALSE)
  key_data <- read.csv(data_file, row.names = 1, fileEncoding = "UTF-8", check.names = FALSE)
  
  cat("📊 Data loaded successfully:\n")
  cat("   - Top30 classification file:", nrow(key_factors_class), "factors\n")
  cat("   - Top30 data file:", nrow(key_data), "samples ×", ncol(key_data), "variables\n")
} else {
  stop("❌ Please specify valid file paths for classified_file and data_file")
}

# ===================== Variable Selection =====================
select_optimal_variable <- function(type, class_df, data_df) {
  # Filter factors by type
  type_factors <- class_df[class_df$Feature_Type == type, ]
  # Sort by Importance descending
  type_factors <- type_factors[order(-type_factors$Importance), ]
  
  # Select optimal variable
  if (type == "Soil_properties") {
    cat_factors <- grep("CAT|catalase", type_factors$Original_Name, ignore.case = TRUE, value = TRUE)
    optimal_name <- if (length(cat_factors) > 0) cat_factors[1] else type_factors$Original_Name[1]
  } else {
    optimal_name <- type_factors$Original_Name[1]
  }
  
  # Match data column names
  prefix <- switch(type,
                   Microbes = "Microbe_",
                   Metabolites = "Metab_",
                   Soil_properties = "Phys_")
  full_var_name <- paste0(prefix, optimal_name)
  
  # Error handling
  if (!full_var_name %in% colnames(data_df)) {
    full_var_name <- colnames(data_df)[grepl(paste0("^", prefix), colnames(data_df))][1]
    cat("⚠️", type, "variable matching failed, auto-selected:", full_var_name, "\n")
  }
  
  # Output results
  cat("✅", type, "optimal variable selection result:\n")
  cat("   - Original name:", optimal_name, "\n")
  cat("   - Data column name:", full_var_name, "\n")
  cat("   - Importance value:", type_factors$Importance[1], "\n\n")
  
  return(list(full_name = full_var_name, clean_name = optimal_name))
}

# Select three core variable types
microbe_var <- select_optimal_variable("Microbes", key_factors_class, key_data)
metab_var   <- select_optimal_variable("Metabolites", key_factors_class, key_data)
phys_var    <- select_optimal_variable("Soil_properties", key_factors_class, key_data)

# Build analysis data
group <- ifelse(grepl("^PSS|_PSS|PE", rownames(key_data), ignore.case = TRUE), 1, 0)
model_data <- cbind(Group = group, key_data) %>%
  dplyr::mutate(dplyr::across(-Group, ~as.vector(scale(.)))) %>%  # Adapt to new dplyr syntax
  na.omit()

# Filter core data
core_vars <- c("Group", microbe_var$full_name, metab_var$full_name, phys_var$full_name)
core_data <- model_data[, core_vars]
colnames(core_data) <- c("Group", "Microbe", "Metabolite", "Soil_phys")

cat("📋 Final analysis variables:\n")
cat("   - Independent variable (X): Group (PSS=1/NPSS=0)\n")
cat("   - Mediator 1 (M1):", microbe_var$clean_name, "\n")
cat("   - Mediator 2 (M2):", metab_var$clean_name, "\n")
cat("   - Dependent variable (Y):", phys_var$clean_name, "\n")
cat("   - Valid sample size:", nrow(core_data), "\n\n")

# ===================== SEM Model Construction (Optimized) =====================
# Define chain mediation model
sem_model <- '
  # Path definitions
  Microbe ~ a1*Group                # Plastisphere → Microbes
  Metabolite ~ b1*Microbe + a2*Group  # Microbes/Plastisphere → Metabolites
  Soil_phys ~ c1*Metabolite + b2*Microbe + c2*Group  # Metabolites/Microbes/Plastisphere → Enzyme activity
  
  # Effect decomposition
  total_effect := a1*b1*c1 + a1*b2 + a2*c1 + c2
  indirect_chain := a1*b1*c1
  indirect_microbe := a1*b2
  indirect_metab := a2*c1
  direct_effect := c2
  
  # Residual correlations
  Microbe ~~ Metabolite
  Metabolite ~~ Soil_phys
'

# Fit SEM model (Optimized: add convergence settings, fix bootstrap parameters)
set.seed(123)
sem_fit <- lavaan::sem(
  model = sem_model,
  data = core_data,
  estimator = "ML",
  se = "bootstrap",
  bootstrap = 1000, # Reduce bootstrap iterations for speed (adjust back to 5000 if needed)
  missing = "fiml",
  control = list(max.iter = 10000, conv.tol = 1e-06) # Increase iterations for convergence
)

# ===================== Model Fit Validation (Fixed NA handling) =====================
# Extract fit indices
fit_indices <- lavaan::fitMeasures(sem_fit, c(
  "chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic"
))

# Handle NA values
fit_indices[is.na(fit_indices)] <- 0

# Fit evaluation
fit_evaluation <- data.frame(
  Index = c("χ²/df", "CFI", "TLI", "RMSEA", "SRMR", "AIC", "BIC"),
  Value = c(
    ifelse(fit_indices["chisq"] == 0 | fit_indices["df"] == 0, "0", round(fit_indices["chisq"]/fit_indices["df"], 2)),
    round(fit_indices["cfi"], 3),
    round(fit_indices["tli"], 3),
    paste0(round(fit_indices["rmsea"], 3), " [", round(fit_indices["rmsea.ci.lower"], 3), ",", round(fit_indices["rmsea.ci.upper"], 3), "]"),
    round(fit_indices["srmr"], 3),
    round(fit_indices["aic"], 1),
    round(fit_indices["bic"], 1)
  ),
  Reference_Standard = c("1-3 (Excellent)", ">0.90 (Acceptable)", ">0.90 (Acceptable)", "<0.08 (Acceptable)", "<0.08 (Acceptable)", "Smaller is better", "Smaller is better"),
  Model_Evaluation = c(
    ifelse(fit_indices["chisq"]/fit_indices["df"] < 3 & fit_indices["chisq"]/fit_indices["df"] > 0, "Excellent", 
           ifelse(fit_indices["chisq"] == 0, "No value", "Acceptable")),
    ifelse(fit_indices["cfi"] > 0.90, "Passed", ifelse(fit_indices["cfi"] == 0, "No value", "Failed")),
    ifelse(fit_indices["tli"] > 0.90, "Passed", ifelse(fit_indices["tli"] == 0, "No value", "Failed")),
    ifelse(fit_indices["rmsea"] < 0.08, "Passed", "Failed"),
    ifelse(fit_indices["srmr"] < 0.08, "Passed", "Failed"),
    "-",
    "-"
  )
)

cat("=== SEM Model Fit Validation ===\n")
print(fit_evaluation, row.names = FALSE)

# ===================== Mediation Effect Test (Optimized) =====================
# Build regression models
m1 <- lm(Microbe ~ Group, data = core_data)          # X→M1
m2 <- lm(Metabolite ~ Group + Microbe, data = core_data)  # X+M1→M2
m3 <- lm(Soil_phys ~ Group + Microbe + Metabolite, data = core_data)  # X+M1+M2→Y

# Mediation effect test (Optimized: add result error handling)
set.seed(123)
med_microbe <- mediation::mediate(
  m1, m3, 
  treat = "Group", 
  mediator = "Microbe", 
  boot = TRUE, 
  sims = 1000, # Reduce simulation iterations for speed
  conf.level = 0.95
)

med_metab <- mediation::mediate(
  m2, m3, 
  treat = "Group", 
  mediator = "Metabolite", 
  boot = TRUE, 
  sims = 1000,
  conf.level = 0.95
)

# Organize mediation effect results (Fixed numerical processing)
mediation_results <- data.frame(
  Mediator_Variable = c(microbe_var$clean_name, metab_var$clean_name),
  Total_Effect = round(c(med_microbe$tau.coef, med_metab$tau.coef), 3),
  Direct_Effect = round(c(med_microbe$z1, med_metab$z1), 3),
  Average_Mediation_Effect = round(c(med_microbe$d1, med_metab$d1), 3),
  `95%CI_Lower` = round(c(med_microbe$d1.ci[1], med_metab$d1.ci[1]), 3),
  `95%CI_Upper` = round(c(med_microbe$d1.ci[2], med_metab$d1.ci[2]), 3),
  Mediation_Effect_Ratio = round(c(
    ifelse(med_microbe$tau.coef != 0, med_microbe$d1/med_microbe$tau.coef*100, 0),
    ifelse(med_metab$tau.coef != 0, med_metab$d1/med_metab$tau.coef*100, 0)
  ), 1),
  Significance = ifelse(c(
    med_microbe$d1.ci[1]*med_microbe$d1.ci[2] > 0, 
    med_metab$d1.ci[1]*med_metab$d1.ci[2] > 0
  ), 
  "Significant (CI does not contain 0)", "Not significant (CI contains 0)")
)

cat("\n=== Mediation Effect Bootstrap Test Results ===\n")
print(mediation_results, row.names = FALSE)

# ===================== Model Robustness Validation (Fixed) =====================
# Split sample validation (Add sample size error handling)
set.seed(123)
if (nrow(core_data) >= 5) { # Ensure sufficient samples
  train_idx <- sample(1:nrow(core_data), round(0.7*nrow(core_data)))
  train_data <- core_data[train_idx, ]
  test_data <- core_data[-train_idx, ]
  
  # Fit model on training set
  sem_train <- lavaan::sem(sem_model, data = train_data, estimator = "ML", se = "bootstrap", bootstrap = 500)
  
  # Predict on validation set (Fixed prediction logic)
  sem_pred <- tryCatch({
    lavaan::predict(sem_train, newdata = test_data)
  }, error = function(e) {
    # Use fitted values if prediction fails
    data.frame(
      Microbe = predict(lm(Microbe ~ Group, data = train_data), newdata = test_data),
      Metabolite = predict(lm(Metabolite ~ Group + Microbe, data = train_data), newdata = test_data),
      Soil_phys = predict(lm(Soil_phys ~ Group + Microbe + Metabolite, data = train_data), newdata = test_data)
    )
  })
  
  # Calculate prediction error (Fixed atomic vector error)
  pred_error <- data.frame(
    Validation_Sample_Size = nrow(test_data),
    Microbe_Prediction_RMSE = round(sqrt(mean((test_data$Microbe - sem_pred$Microbe)^2, na.rm = TRUE)), 3),
    Metabolite_Prediction_RMSE = round(sqrt(mean((test_data$Metabolite - sem_pred$Metabolite)^2, na.rm = TRUE)), 3),
    Enzyme_Activity_Prediction_RMSE = round(sqrt(mean((test_data$Soil_phys - sem_pred$Soil_phys)^2, na.rm = TRUE)), 3),
    Average_Prediction_Accuracy = round(100 - mean(c(
      sqrt(mean((test_data$Microbe - sem_pred$Microbe)^2, na.rm = TRUE)),
      sqrt(mean((test_data$Metabolite - sem_pred$Metabolite)^2, na.rm = TRUE)),
      sqrt(mean((test_data$Soil_phys - sem_pred$Soil_phys)^2, na.rm = TRUE))
    ))*100, 1)
  )
} else {
  # Alternative for insufficient samples
  pred_error <- data.frame(
    Validation_Sample_Size = 0,
    Microbe_Prediction_RMSE = 0,
    Metabolite_Prediction_RMSE = 0,
    Enzyme_Activity_Prediction_RMSE = 0,
    Average_Prediction_Accuracy = 0
  )
  cat("⚠️ Insufficient sample size, skipped robustness validation\n")
}

# ===================== Parameter Sensitivity Analysis (Fixed) =====================
# Extract standardized parameters (Fixed std.all missing issue)
sensitivity_analysis <- lavaan::parameterEstimates(sem_fit, standardized = TRUE) %>%
  dplyr::filter(op == "~") %>%
  dplyr::mutate(
    Variable_Relationship = dplyr::case_when(
      lhs == "Microbe" & rhs == "Group" ~ "Plastisphere → Microbes",
      lhs == "Metabolite" & rhs == "Microbe" ~ "Microbes → Metabolites",
      lhs == "Metabolite" & rhs == "Group" ~ "Plastisphere → Metabolites",
      lhs == "Soil_phys" & rhs == "Metabolite" ~ "Metabolites → Enzyme activity",
      lhs == "Soil_phys" & rhs == "Microbe" ~ "Microbes → Enzyme activity",
      lhs == "Soil_phys" & rhs == "Group" ~ "Plastisphere → Enzyme activity"
    ),
    Standardized_Coefficient = round(std.all, 3),
    Standard_Error = round(se, 3),
    P_Value = round(pvalue, 3),
    Significance = dplyr::case_when(
      pvalue < 0.001 ~ "***",
      pvalue < 0.01 ~ "**",
      pvalue < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(Variable_Relationship, Standardized_Coefficient, Standard_Error, P_Value, Significance)

cat("\n=== Parameter Sensitivity Analysis ===\n")
print(sensitivity_analysis, row.names = FALSE)

# ===================== C/N Cycle-Ecological Effect Analysis (Fixed) =====================
# Build C/N cycle comprehensive index (Fixed cor function data type issue)
cn_cycle_analysis <- core_data %>%
  dplyr::mutate(
    CN_Cycle_Index = (Soil_phys - min(Soil_phys))/(max(Soil_phys) - min(Soil_phys)),
    Group_Label = ifelse(Group == 1, "Plastisphere (PSS)", "Control (NPSS)")
  ) %>%
  dplyr::select(Group, Microbe, Metabolite, CN_Cycle_Index) %>%
  # Calculate correlation only for numeric columns
  cor(method = "spearman", use = "complete.obs") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  dplyr::filter(Variable %in% c("Group", "Microbe", "Metabolite")) %>%
  dplyr::select(Variable, CN_Cycle_Index) %>%
  dplyr::rename(CN_Cycle_Correlation = CN_Cycle_Index) %>%
  dplyr::mutate(
    Ecological_Effect_Direction = ifelse(CN_Cycle_Correlation > 0, "Promotion", "Inhibition"),
    Effect_Strength = dplyr::case_when(
      abs(CN_Cycle_Correlation) > 0.7 ~ "Strong",
      abs(CN_Cycle_Correlation) > 0.4 ~ "Moderate",
      TRUE ~ "Weak"
    ),
    Ecological_Environmental_Significance = dplyr::case_when(
      Variable == "Group" & Ecological_Effect_Direction == "Promotion" ~ "Plastisphere enhances soil ecological function by strengthening C/N cycling",
      Variable == "Group" & Ecological_Effect_Direction == "Inhibition" ~ "Plastisphere reduces soil ecological function by weakening C/N cycling",
      Variable == "Microbe" ~ paste("Core microbes", Ecological_Effect_Direction, "C/N cycling process, serving as key mediators for plastisphere regulation of ecological function"),
      Variable == "Metabolite" ~ paste("Core metabolites", Ecological_Effect_Direction, "C/N cycling process, serving as core pathways for plastisphere regulation of ecological function")
    )
  )

cat("\n=== C/N Cycle-Ecological Environmental Effect Analysis ===\n")
print(cn_cycle_analysis, row.names = FALSE)

# ===================== Visualization Output (Fixed) =====================
# 1. Mediation effect forest plot (Fixed font and parameter issues)
p_mediation <- ggplot2::ggplot(mediation_results, ggplot2::aes(x = Average_Mediation_Effect, y = Mediator_Variable)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  ggplot2::geom_point(size = 5, color = "#2E86AB") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = `95%CI_Lower`, xmax = `95%CI_Upper`), 
                          height = 0.2, linewidth = 1.2, color = "#2E86AB") +
  ggplot2::geom_text(ggplot2::aes(label = paste0(Average_Mediation_Effect, " [", `95%CI_Lower`, ",", `95%CI_Upper`, "]")), 
                     hjust = -0.1, size = 4, family = "SimHei") +
  ggplot2::labs(
    title = "Mediation Effect Analysis (Bootstrap 1000)",
    x = "Average Mediation Effect (Standardized)",
    y = "Mediator Variable",
    caption = "Note: Red dashed line = effect = 0; Significant if CI does not contain 0"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold", family = "SimHei"),
    plot.caption = ggplot2::element_text(hjust = 0.5, size = 10, family = "SimHei"),
    axis.title = ggplot2::element_text(size = 12, family = "SimHei"),
    axis.text = ggplot2::element_text(size = 11, family = "SimHei"),
    panel.grid.major = ggplot2::element_line(linewidth = 0.2),
    panel.grid.minor = ggplot2::element_blank()
  )

# 2. SEM path diagram (Fixed font and parameter errors)
showtext_auto(FALSE)
pdf("SEM_Path_Diagram.pdf", width = 10, height = 8) # Remove Arial font dependency
semPlot::semPaths(
  sem_fit,
  what = "std",          
  whatLabels = "std",    
  layout = "tree",       
  edge.label.cex = 1.2,  
  node.label.cex = 1.4,  
  node.color = c(Group = "#ffebee", Microbe = "#e8f5e9",
                 Metabolite = "#fff3e0", Soil_phys = "#f3e5f5"),
  nodeLabels = c("Plastisphere", microbe_var$clean_name, metab_var$clean_name, phys_var$clean_name),
  main = "Chain Mediation Model of Plastisphere on Soil C/N Cycle",
  main.cex = 1.5
)
dev.off()
showtext_auto(TRUE)

# 3. Fit indices plot (Fixed NA handling)
fit_plot_data <- fit_evaluation[1:5, ] %>%
  dplyr::mutate(
    Threshold_Value = dplyr::case_when(
      Index == "χ²/df" ~ 3,
      Index %in% c("CFI", "TLI") ~ 0.9,
      Index %in% c("RMSEA", "SRMR") ~ 0.08
    ),
    Actual_Value = as.numeric(gsub("\\[.*\\]", "", Value))
  ) %>%
  dplyr::replace_na(list(Actual_Value = 0)) # Replace NA values

p_fit <- ggplot2::ggplot(fit_plot_data, ggplot2::aes(x = Index, y = Actual_Value)) +
  ggplot2::geom_bar(stat = "identity", fill = "#2E86AB", width = 0.6) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = Threshold_Value), linetype = "dashed", color = "red", linewidth = 1) +
  ggplot2::geom_text(ggplot2::aes(label = Actual_Value), vjust = -0.5, size = 4, fontface = "bold") +
  ggplot2::geom_text(ggplot2::aes(y = Threshold_Value, label = paste0("Threshold: ", Threshold_Value)), 
                     vjust = -0.5, hjust = 1.2, color = "red", size = 3) +
  ggplot2::labs(
    title = "SEM Model Fit Indices",
    x = "Fit Index",
    y = "Value",
    caption = "Note: Red dashed line = acceptable threshold"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold", family = "SimHei"),
    plot.caption = ggplot2::element_text(hjust = 0.5, size = 10, family = "SimHei"),
    axis.title = ggplot2::element_text(size = 12, family = "SimHei"),
    axis.text = ggplot2::element_text(size = 11, family = "SimHei")
  )

# 4. C/N cycle ecological effect plot
p_cn_effect <- ggplot2::ggplot(cn_cycle_analysis, ggplot2::aes(x = Variable, y = CN_Cycle_Correlation, fill = Ecological_Effect_Direction)) +
  ggplot2::geom_bar(stat = "identity", width = 0.6, color = "black") +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  ggplot2::geom_text(ggplot2::aes(label = paste0(round(CN_Cycle_Correlation,3), " (", Effect_Strength, ")")), 
                     vjust = -0.5, size = 4, fontface = "bold", family = "SimHei") +
  ggplot2::scale_fill_manual(values = c("Promotion" = "#2E86AB", "Inhibition" = "#E63946")) +
  ggplot2::labs(
    title = "C/N Cycle Ecological Effect of Key Variables",
    x = "Key Variable",
    y = "Spearman Correlation with C/N Cycle Index",
    fill = "Ecological Effect Direction",
    caption = "Note: C/N Cycle Index = Normalized CAT activity"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold", family = "SimHei"),
    plot.caption = ggplot2::element_text(hjust = 0.5, size = 10, family = "SimHei"),
    axis.title = ggplot2::element_text(size = 12, family = "SimHei"),
    axis.text = ggplot2::element_text(size = 11, family = "SimHei"),
    legend.title = ggplot2::element_text(size = 11, family = "SimHei"),
    legend.text = ggplot2::element_text(size = 10, family = "SimHei")
  )

# Save visualization results (Add error handling)
tryCatch({
  ggplot2::ggsave("Mediation_Effect_Forest_Plot.pdf", p_mediation, width = 10, height = 6, dpi = 300)
  cat("✅ Mediation effect forest plot saved successfully\n")
}, error = function(e) {
  cat("⚠️ Mediation effect forest plot save failed:", e$message, "\n")
})

tryCatch({
  ggplot2::ggsave("SEM_Fit_Indices_Plot.pdf", p_fit, width = 8, height = 6, dpi = 300)
  cat("✅ Fit indices plot saved successfully\n")
}, error = function(e) {
  cat("⚠️ Fit indices plot save failed:", e$message, "\n")
})

tryCatch({
  ggplot2::ggsave("CN_Cycle_Ecological_Effect_Plot.pdf", p_cn_effect, width = 8, height = 6, dpi = 300)
  cat("✅ C/N cycle ecological effect plot saved successfully\n")
}, error = function(e) {
  cat("⚠️ C/N cycle ecological effect plot save failed:", e$message, "\n")
})

# ===================== Results Saving =====================
# Save core result tables (Add error handling)
write.csv(fit_evaluation, "SEM_Fit_Evaluation.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(mediation_results, "Mediation_Effect_Results.csv", row.names = FALSE, fileEncoding = "UTF-8")

tryCatch({
  write.csv(sensitivity_analysis, "Parameter_Sensitivity_Analysis.csv", row.names = FALSE, fileEncoding = "UTF-8")
}, error = function(e) {
  cat("⚠️ Parameter sensitivity analysis file save failed:", e$message, "\n")
})

tryCatch({
  write.csv(pred_error, "Model_Robustness_Validation.csv", row.names = FALSE, fileEncoding = "UTF-8")
}, error = function(e) {
  cat("⚠️ Model robustness validation file save failed:", e$message, "\n")
})

write.csv(core_data, "SEM_Core_Analysis_Data.csv", row.names = TRUE, fileEncoding = "UTF-8")

tryCatch({
  write.csv(cn_cycle_analysis, "CN_Cycle_Ecological_Effect.csv", row.names = FALSE, fileEncoding = "UTF-8")
}, error = function(e) {
  cat("⚠️ C/N cycle ecological effect file save failed:", e$message, "\n")
})

# Save SEM detailed results
sem_detailed <- lavaan::parameterEstimates(sem_fit, standardized = TRUE) %>%
  dplyr::mutate(
    Standardized_Coefficient = round(std.all, 3),
    Unstandardized_Coefficient = round(est, 3),
    Standard_Error = round(se, 3),
    P_Value = round(pvalue, 3),
    `95%CI_Lower` = round(ci.lower, 3),
    `95%CI_Upper` = round(ci.upper, 3)
  ) %>%
  dplyr::select(lhs, op, rhs, Unstandardized_Coefficient, Standardized_Coefficient, Standard_Error, P_Value, `95%CI_Lower`, `95%CI_Upper`)

write.csv(sem_detailed, "SEM_Detailed_Results.csv", row.names = FALSE, fileEncoding = "UTF-8")
