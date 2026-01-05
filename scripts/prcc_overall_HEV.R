library(epiR)
library(sensitivity)

# find directory where script is located
this_file_path <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) NULL
)

# fallback if NULL
if (is.null(this_file_path) && requireNamespace("rstudioapi", quietly = TRUE)) {
  this_file_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
}
if (is.null(this_file_path)) {
  this_file_path <- getwd()
}

# build data folder path relative to script
data_dir <- file.path(this_file_path, "..", "data")  # go up from 'scripts' to repo root, then into 'data'
csv_file <- file.path(data_dir, "sensitivity_HEV_data.csv")

# check if the file exists
if (!file.exists(csv_file)) {
  stop("Data file not found: ", csv_file, "\nPlease run the MATLAB scripts first to generate the 'data' folder and CSV files.")
}

df <- read.csv(csv_file, header = TRUE)

params <- c("q", "phi", "epsilon", "p_v", "nu_max", "wepsilon_max", "a", "c")

interventions <- 2:7
thresholds <- c(1, 10, 25)

params_for_intervention <- list(
  "2" = c("q","phi"),
  "3" = c("q","phi","epsilon","p_v","nu_max"),
  "4" = c("q","phi","wepsilon_max","a"),
  "5" = c("q","phi","c"),
  "6" = c("q","phi","wepsilon_max","a","c"),
  "7" = params  # all
)

# store PRCC results in a list
prcc_results <- list()

for (i in interventions) {
  these_params <- params_for_intervention[[as.character(i)]]
  for (t in thresholds) {
    rr_col <- paste0("RR_T", t, "_I", i)
    # subset only relevant parameters plus the response column
    current_df <- df[, c(these_params, rr_col)]
    prcc_obj   <- epi.prcc(current_df, sided.test = 2, conf.level = 0.95)
    prcc_results[[paste0("I", i, "_T", t)]] <- prcc_obj
  }
}

# print
for (i in interventions) {
  for (t in thresholds) {
    cat(sprintf("\nResults for Intervention %d, Threshold %d:\n", i, t))
    print(prcc_results[[paste0("I", i, "_T", t)]])
  }
}

intervention_labels <- c(
  "2" = "(A) Testing & Isolation",
  "3" = "(B) Vaccination",
  "4" = "(C) Sanitation",
  "5" = "(D) Water Treatment",
  "6" = "(E) WASH",
  "7" = "(F) Vaccination + WASH"
)

####################################### T=1 and T=25 side by side:
create_prcc_overall_plot <- function(prcc_results, outpath) {
  pdf(outpath, width = 12, height = 6)
  layout(matrix(1:2, nrow = 1))  # 1 row, 2 columns
  
  par(mar = c(10, 5, 5, 2),
      cex.axis = 1.5,
      cex.lab = 1.6,
      cex.main = 1.6)
  
  all_colors <- c(
    q = "#D3D3D3",
    phi = "#D3D3D3",
    epsilon = "#D3D3D3",
    p_v = "#D3D3D3",
    nu_max = "#D3D3D3",
    epsilon_omega_max = "#D3D3D3",
    a = "#D3D3D3",
    c = "#D3D3D3"
  )
  
  param_keys <- c("q", "phi", "epsilon", "p_v", "nu_max", "epsilon_omega_max", "a", "c")
  labels_expr <- expression(italic(q), italic(phi), italic(epsilon), italic(p[v]), italic(nu[max]),
                            italic(epsilon[omega]^max), italic(a), italic(c))
  
  labels_panel <- c("(A)", "(B)")
  thresholds <- c(1, 25)
  
  for (i in 1:2) {
    Tval <- thresholds[i]
    prcc_data <- prcc_results[[paste0("I7_T", Tval)]]
    rownames(prcc_data) <- param_keys
    
    param_order <- param_keys[order(prcc_data$est, decreasing = TRUE)]
    idx <- match(param_order, param_keys)
    ordered_labels <- labels_expr[idx]
    ordered_colors <- unname(all_colors[param_keys[idx]])
    
    bar_pos <- barplot(
      prcc_data[param_order, "est"],
      ylim = c(-1, 1.25),
      col = ordered_colors,
      ylab = expression(bold("PRCC")),
      names.arg = ordered_labels,
      las = 2,
      cex.names = 1.9
    )
    
    abline(h = 0, lty = 2)
    arrows(
      bar_pos, prcc_data[param_order, "lower"],
      bar_pos, prcc_data[param_order, "upper"],
      length = 0.05, angle = 90, code = 3
    )
    
    # Add (A) or (B)
    usr <- par("usr")
    mtext(text = labels_panel[i], side = 3, line = 3.5, at = usr[1], adj = 0, font = 2, cex = 2)
  }
  
  dev.off()
}

# detect current script directory (same logic as above)
this_file_path <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) NULL
)

if (is.null(this_file_path) && requireNamespace("rstudioapi", quietly = TRUE)) {
  this_file_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
}
if (is.null(this_file_path)) {
  this_file_path <- getwd()
}

# build figures folder path relative to script location
fig_dir <- file.path(this_file_path, "..", "figures")

# create 'figures' folder if it doesn't exist
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# build full path for Figure5
outname <- file.path(fig_dir, "Figure5.pdf")
create_prcc_overall_plot(prcc_results, outname)

