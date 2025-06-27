################################################################################
#
# AUTHOR: Dr. Tanuj Puri, DATE: 2021, updated: 2025
#
################################################################################
#
# FILE: dxa_data_analysis.R
#
# PURPOSE: 
#   This R script performs the core analysis used in the associated 
#   scientific study published in the British Journal of Radiology (BJR).
#   It validates multiple historical equations for estimating Lean Body Mass (LBM) 
#   against reference DXA-derived values, using both real and simulated datasets.
#
#   The script includes data loading, cleaning, transformation, statistical evaluation, 
#   correlation analysis, and Bland–Altman comparison plots.
#
#   The results generated from this script (applied to real data, not available in 
#   this repository) were published the BJR article (`Article_BJR_20210378.pdf` 
#   and its supplement).
#
# STRUCTURE:
#   The script is organized into three main parts:
#     1. Function Definitions: 
#        - Data cleaning, inspection, demographics summary
#        - Statistical summaries, correlation, and Bland–Altman plotting
#     2. Analysis Workflow:
#        - Load and inspect data (real or synthetic)
#        - Clean data by removing predefined error codes
#        - Subgroup analysis by sex (Male, Female, All)
#        - Compute LBM using 10 published equations in literature
#        - Compare computed values with reference DXA measurements
#     3. Output:
#        - Correlation plots (printed in the R console)
#        - Bland–Altman plots (saved as `.jpg` images in the working directory)
#
# DEPENDENCIES:
#   R Packages required (install if missing):
#     - BlandAltmanLeh: For Bland–Altman plot visualization
#     - ggplot2:        For correlation plots
#
#   To install, run:
#     install.packages("BlandAltmanLeh")
#     install.packages("ggplot2")
#
# INPUT:
#   - The primary dataset used is: `simulated_dxa_dataset.csv`. It is a modified 
#     version of the original data source, provided by Prof. Glen Blake, KCL.
#
#   - Required columns include: p_ID, sex, htn09 (height in cm), wtn09 (weight in kg), 
#     AGEN09 (age in yrs), bmi09, and DXA-derived body composition values (e.g., DXAWBBMC09).
#
# OUTPUT:
#   - Console: Summary statistics, sex counts, correlation plots.
#   - Files:   Bland–Altman plots as `.jpg` for All, Male, and Female groups.
#              These are saved in the current working directory.
#
# USAGE:
#   - Can be run in RStudio, base R, or Google Colab (auto-detects Colab path).
#   - Ensure `simulated_dxa_dataset.csv` is in your working directory.
#   - Modify the `read_dataset()` function if using real data or a different path.
#
# LICENSE:
#   - This script is released under the GNU General Public License v3.0
#     for transparency, reproducibility, and academic reuse.
#
#
################################################################################

#=============================START OF FUNCTIONS================================
#
# FUNCTION 1: clear all figures and variables-----------------------------------
clearAll <- function() {
  print("Before clear: ")
  print(ls(envir = .GlobalEnv))                    # list variables in global env

  rm(list=ls(envir=.GlobalEnv), envir=.GlobalEnv); # clear all variables in global env
  graphics.off();                                  # close all plots

  print("After clear: ")
  print(ls(envir = .GlobalEnv))                    # list variables in global env
  }


# TEST 1: clearAll FUNCTION. To run, change "FALSE" to "TRUE"
if (FALSE) {
  a <- 10;
  b <- 'hello';
  plot(1:10, 11:20);
  clearAll()
}
#-------------------------------------------------------------------------------

# FUNCTION 2: checkData: check the datastructure--------------------------------
checkData <- function(checkingdata) {
  col_names = names(checkingdata)

  return( list(
      col_names  = col_names,
      dimensions = dim(checkingdata),
      structure  = capture.output(str(checkingdata)),
      head       = head(checkingdata),
      tail       = tail(checkingdata)
  ))
}
# TEST 2: checkData FUNCTION. To run, change "FALSE" to "TRUE"
if (FALSE) {
  df <- mtcars #sample data
  summary_info <- checkData(df)
  print(summary_info)
}
#-------------------------------------------------------------------------------

# FUNCTION 3: countSex: count number of male & females in data------------------
countSex <- function(data) {
  if (all("sex" != names(data))) {
      stop("Error: 'sex' colum not found in data")
  }else{
      males   <- sum(data$sex == 1, na.rm = TRUE)
      females <- sum(data$sex == 2, na.rm = TRUE)
  }
  return(c(Males=males, Females=females))
}
# TEST 3: countSex FUNCTION. To run, change "FALSE" to "TRUE"
if (FALSE) {
  df <- data.frame(sex = c(1, 2, 1, 2, 1, 2, NA, 2)) #sample data
  sex_count <- countSex(df) # call the function
  print(sex_count)
}
#------------------------------------------------------------------------------
# FUNCTION 4: Function to summarize demographic data: height, weight, age
demoData <- function(demographicData, digits=2, print=TRUE){
    required_cols <- c("htn09","wtn09","AGEN09")
    missing       <- setdiff( required_cols, names(demographicData) )

    if (length(missing)>0){
        stop(paste("Missing columns: ", paste(missing, collapse=", ")))
    }
    summary_stats <- data.frame(
        Variable = c("Height", "Weight", "Age"),
        Mean     = round( c( mean(demographicData$htn09), mean(demographicData$wtn09), mean(demographicData$AGEN09) ) , digits),
        Min      = round( c( min(demographicData$htn09), min(demographicData$wtn09), min(demographicData$AGEN09) ) , digits),
        Max      = round( c( max(demographicData$htn09), max(demographicData$wtn09), max(demographicData$AGEN09) ) , digits)
        )
    if (print){
        print(summary_stats)
    }
    return(summary_stats)
}
# TEST 4: demoData FUNCTION. To run, change "FALSE" to "TRUE"
if(FALSE){
  test_df <- data.frame(
    htn09  = c(160.5, 170.2, 165.0),    # height in cm
    wtn09  = c(60.0, 72.5, 68.3),       # weight in kg
    AGEN09 = c(61, 63, 65),             # age in years
    temp09 = c(01, 39, 500)             # age in years
  )
  summary_output <- demoData(test_df, print=FALSE);
  print(summary_output)
}
#-------------------------------------------------------------------------------
# FUNCTION 5: Bland-Altman plot function
install.packages("BlandAltmanLeh")
library(BlandAltmanLeh)
baPlotFun <- function(A, B, eq, sex, sex_data, fig_name){
    # Input check
    if ( length(A) != length(B) ) {
        stop("A and B must be same length.")
    }
    # Correlations
    pear <- cor.test(x=A, y=B, method='pearson')
    sper <- cor.test(x=A, y=B, method='spearman')

    # Normality of differences
    diff           <- A - B;
    shapiro_result <- shapiro.test(diff);
    normality      <- if(shapiro_result$p.value <= 0.05) "N" else "Y";

    # Regression slope
    slope <- cor.test(diff, 0.5 * (A + B), method = "pearson")

    # Plotting
    print(paste("Saving plot to:", fig_name))
    graphics.off()
    jpeg(fig_name)  # Save to JPEG

    bland.altman.plot(
    A, B,
    conf.int = 0.95,
    pch = 19,
    col = sex,
    main = paste("Eq", eq, "vs LBM:", sex_data),
    sub = paste0(
         "p.r=",   round(pear$estimate, 4), ", p.p=", round(pear$p.value, 2),
      " | s.r=",   round(sper$estimate, 2), ", s.p=", round(sper$p.value, 2),
      " | bias=",  round(bland.altman.stats(A, B)$mean.diffs, 2),
      " | lim=",   round(2 * bland.altman.stats(A, B)$critical.diff, 2),
      " | norm=",  normality,
      " | slp.r=", round(slope$estimate, 2), ", slp.p=", round(slope$p.value, 2)
    )
  );
  legend("topleft", legend = c("Male", "Female"), fill = 1:2)

  # Calculate stats for labeling
  ba_stats  <- bland.altman.stats(A, B)
  mean_diff <- ba_stats$mean.diffs
  crit_diff <- ba_stats$critical.diff

  # Y-axis values for lines
  lower_limit <- mean_diff - crit_diff
  upper_limit <- mean_diff + crit_diff

  # X-position: pick a suitable spot, e.g., mean of X
  x_pos <- mean( (A + B) / 2 )

  # Add labels near lines (adjust y-offset if needed)
  text(x = x_pos, y = mean_diff,   labels = "Mean (Bias)", pos = 3, cex = 0.8, col = "blue")
  text(x = x_pos, y = upper_limit, labels = "+1.96 SD",   pos = 3,  cex = 0.8, col = "red")
  text(x = x_pos, y = lower_limit, labels = "−1.96 SD",   pos = 1,  cex = 0.8, col = "red")

  dev.off()

  # Clean up
  rm(diff, A, B, pear, sper, shapiro_result, slope)
}

# TEST 5: baPlotFun FUNCTION. To run, change "FALSE" to "TRUE"
if(FALSE){
  set.seed(42)
  A   <- rnorm(50, mean=69, sd=5)
  B   <- A + rnorm(50, mean=0, sd=2)
  sex <- sample(1:2, 50, replace=TRUE)
  baPlotFun(A, B, eq="1", sex=sex, sex_data="Both", fig_name="/content/bland_altman_test.jpg")
}
#-------------------------------------------------------------------------------
# FUNCTION 6: mean, standard error of the mean, 95% CI function
meanSemCiFun <- function(A, B, print_results = TRUE) {
  # Calculate differences
  diffs <- A - B

  # Mean difference
  mean_diff <- mean(diffs)

  # Standard error of the mean (SEM)
  sem <- sd(diffs) / sqrt(length(diffs))

  # 95% Confidence Interval (assuming normal distribution, approx 1.96 for 95%)
  ci_lower <- mean_diff - 1.96 * sem
  ci_upper <- mean_diff + 1.96 * sem
  ci       <- c(ci_lower, ci_upper)

  # Optional print output
  if (print_results) {
    cat("Mean difference:", mean_diff, "\n")
    cat("SEM:", sem, "\n")
    cat("95% CI:", ci_lower, "-", ci_upper, "\n")
  }

  # Return results as a named list
  return(list(mean_diff=mean_diff, sem=sem, ci=ci))
}
# TEST 6: meanSemCiFun FUNCTION. To run, change "FALSE" to "TRUE"
if(FALSE){
  A <- c(10, 20, 30, 40, 50)
  B <- c(12, 18, 29, 41, 49)
  results <- meanSemCiFun(A, B, print_results=FALSE); # Call function
  print(results)# Check returned values
}
#-------------------------------------------------------------------------------
# FUNCTION 7: read file "simulated_dxa_dataset.csv" from local environment
read_dataset <- function(filename = "simulated_dxa_dataset.csv") {
    # Helper to detect Google Colab
    is_colab <- Sys.getenv("COLAB_GPU") != "" || grepl("/content", getwd())

    # Decide base directory based on enviorment
    if (isTRUE(is_colab)) {
      file_path <- "/content/simulated_dxa_dataset.csv"
      } else {
      file_path <- file.path(getwd(), "simulated_dxa_dataset.csv")
    }

    # Check file exist and read
    if (file.exists(file_path)){
        message("Rading file from: ", file_path)
        df <- read.csv(file_path)
        return(df)
    }else{
        stop("File not found at: ", file_path)
    }
}
# TEST 7: read_dataset FUNCTION. To run, change "FALSE" to "TRUE"
if(FALSE){
  dataset <- read_dataset()
  head(dataset)
}
#-------------------------------------------------------------------------------
# FUNCTION 8: to remove error codes from specific columns
clean_dxa_dataset <- function(df, error_codes, columns) {
  for (col in columns) {
    before <- nrow(df)
    df     <- df[!df[[col]] %in% error_codes, ]
    after  <- nrow(df)
    cat("Removed", before - after, "rows with error codes from", col, "\n")
  }
  return(df)
}
# TEST 8: clean_dxa_data FUNCTION. To run, change "FALSE" to "TRUE"
if (FALSE){
    browser()  # Execution will pause here

    dataset <- read_dataset() # read data from file
    str(dataset)              # check the original data
    print(sapply(dataset, class))   # Inspect column types

    error_codes <- c(9, 7777, 8888, 9999, 777777, 999999);
    columns     <- setdiff(colnames(dataset), "p_ID") # all columns except p_ID
    print(columns)

    # Convert all selected columns to numeric (necessary??)
    dataset[columns] <- lapply( dataset[columns], as.numeric )

    dataset_clean <- clean_dxa_dataset(dataset, error_codes, columns)# Clean the dataset

    str(dataset_clean)# Optionally check the result
}
#-------------------------------------------------------------------------------
# FUNCTION 9: validate the data from equations
validate_lbm_columns <- function(df, lbm_cols, reference_vector = NULL) {
  cat("Checking LBM equation outputs...\n")

  # Column check
  missing <- setdiff(lbm_cols, colnames(df))
  if (length(missing) > 0) {
    warning("Missing columns: ", paste(missing, collapse=", "))
    return()
  }

  # NA check
  cat("NA counts:\n")
  print(sapply(df[lbm_cols], function(col) sum(is.na(col))))

  # Summary stats
  cat("\nSummary:\n")
  print(summary(df[lbm_cols]))

  # Correlation with provided reference (must be numeric vector)
  if ( missing(reference_vector) || is.null(reference_vector) || length(reference_vector) != nrow(df) )  {
    warning("Reference_vector must be a numeric vector of same length as the number of rows.")
    return()
  }

  cat("\n Correlation with reference LBM: \n")
  corrs <- sapply(df[lbm_cols], function(col) cor(col, reference_vector, use = "complete.obs"))
  print(round(corrs, 3))
}
#-------------------------------------------------------------------------------
# FUNCTION10: PLot correlations
install.packages("ggplot2")
library(ggplot2)
plot_lbm_correlations <- function(df, lbm_cols, reference_vector) {
  # Validate reference vector length
  if (length(reference_vector) != nrow(df)) {
    stop("Reference vector must have the same number of rows as the dataset.")
  }

  # Check that 'sex' column exists
  if (sum("sex" == colnames(df)) == 0) {
    stop("The dataset must contain a 'sex' column for coloring.")
  }

  # Convert sex to factor with labels
  df$sex <- factor(df$sex, levels = c(1, 2), labels = c("Male", "Female"))

  for (col_name in lbm_cols) {
    if (length(reference_vector) != nrow(df)) {
      warning(paste("Column", col_name, "not found in dataset. Skipping..."))
      next
    }

    plot_data <- data.frame(
      Estimate = df[[col_name]],
      Reference = reference_vector,
      Sex = df$sex
    )    

    p <- ggplot(plot_data, aes(x = Reference, y = Estimate, color = Sex)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = FALSE, color = "blue") +
      scale_color_manual(values = c("Male" = "red", "Female" = "black")) +
      labs(
        title = paste("Correlation Plot:", col_name),
        x = "Reference LBM",
        y = col_name,
        color = "Sex"
      ) +
      theme_minimal()

    print(p) # Print in the figure rather than saving it
  }
}
#===============================END OF FUNCTIONS================================

#==============================START OF ANALYSIS================================
print("--------------------LOAD DATA AND CHECK WHAT IT LOOKS LIKE---------------")
dataset <- read_dataset()
head(dataset)
print("--------------------------")
attr(dataset, "variable.labels")
print("--------------------------")
checkData(dataset)
print("--------------------------")
countSex(dataset)
print("--------------------------")

print("D--------------------ATA CLEAN UP STARTS---------------------------------")
str(dataset) # check the original data
print("--------------------------")
print(sapply(dataset, class)) # Inspect column types
print("--------------------------")
error_codes <- c(9, 7777, 8888, 9999, 777777, 999999);
columns     <- setdiff(colnames(dataset), "p_ID") # all columns except p_ID
print(columns)
print("ignore inf column")
print("--------------------------")
# Convert all selected columns to numeric (necessary??)
dataset[columns] <- lapply( dataset[columns], as.numeric )
dataset_clean    <- clean_dxa_dataset(dataset, error_codes, columns)# Clean dataset
print("--------------------------")
str(dataset_clean) # Check results
print("--------------------------")

print("--------------------RECHECK CLEANED DATA: FOR All, MALES AND FEMALES ----")
checkData(dataset_clean)
countSex(dataset_clean)
demoData(dataset_clean)
print("---------------------------")
dataset_males <- dataset_clean[dataset_clean$sex == 1,]
head(dataset_males)
checkData(dataset_males)
countSex(dataset_males)
demoData(dataset_males)
print("---------------------------")
dataset_females <- dataset_clean[dataset_clean$sex == 2,]
head(dataset_females)
checkData(dataset_females)
countSex(dataset_females)
demoData(dataset_females)
print("---------------------------")

print("--------------------CODE MODELS/EQUATIONS FROM THE PAPER-----------------")

# print("UNITS IN CSV FILE: NSHD_ID,AGEN09 (yrs),bmi09(kg/m^2),DXASBBMC09(g),DXASBFT09(g),DXASBLN09(g),DXASBMS09(g),DXAWBMC09(g),DXAWBFT09(g),DXAWBLN09(g),DXAWBMS09(g),htn09(cm),inf,sex,wtn09(kg)")

# EQ1: 1966 HUME,  W_kg  Height_cm  M=1 F=0 - CHECKED UNITS
dataset_clean$hume_1996 <- ifelse( dataset_clean$sex==1,
                     0.32810*dataset_clean$wtn09 + 0.33929*dataset_clean$htn09 - 29.5336,
                     0.29569*dataset_clean$wtn09 + 0.41813*dataset_clean$htn09 - 43.2933 )

# EQ2: 1971 HUME,  W_kg Height_cm M=1 F=0 - CHECKED UNITS
dataset_clean$hume_1971 <- ifelse( dataset_clean$sex==1,
                      (100/73)*(0.194786*dataset_clean$htn09 + 0.296785*dataset_clean$wtn09 - 14.012934),
                      (100/73)*(0.344547*dataset_clean$htn09 + 0.183809*dataset_clean$wtn09 - 35.270121) )

# EQ3: JAMES 1976 from book: W_kg Height_m - CHECKED UNITS
dataset_clean$james_1976 <- ifelse( dataset_clean$sex==1,
                      dataset_clean$wtn09 - 1.281*(dataset_clean$bmi09 - 10.13)*dataset_clean$wtn09/100,
                      dataset_clean$wtn09 - 1.480*(dataset_clean$bmi09 - 07.00)*dataset_clean$wtn09/100 )

# EQ4: James as used by 1981 HALLYCK, W_kg Height_cm - CHECKED UNITS
dataset_clean$hallynck_1981 <- ifelse( dataset_clean$sex==1,
                      1.10*dataset_clean$wtn09 - 128*(dataset_clean$wtn09/dataset_clean$htn09)^2,
                      1.07*dataset_clean$wtn09 - 148*(dataset_clean$wtn09/dataset_clean$htn09)^2 )

# EQ5: BOER 1984, W_kg Height_m - CHECKED UNITS
dataset_clean$boer_1984 <- ifelse( dataset_clean$sex==1,
                      0.407*dataset_clean$wtn09 + 26.7*(dataset_clean$htn09/100) - 19.2,
                      0.252*dataset_clean$wtn09 + 47.3*(dataset_clean$htn09/100) - 48.3 )

# EQ6: 1991 DEURENBERG, W_kg  Height_m age_yrs sex_num 1_M,0_F - CHECKED UNITS
dataset_clean$deurenberg_1991 <- ifelse( dataset_clean$sex==1,
                      dataset_clean$wtn09 - (1.2*dataset_clean$bmi09  + 0.23*dataset_clean$AGEN09 - 16.2)*dataset_clean$wtn09/100 ,
                      dataset_clean$wtn09 - (1.2*dataset_clean$bmi09  + 0.23*dataset_clean$AGEN09 - 5.4)*dataset_clean$wtn09/100 )

# EQ7: 1977 COMMITTEE OF AMERICAN DIABETES, W_kg Height_cm - CHECKED UNITS
dataset_clean$zasadny_1993 <- ifelse( dataset_clean$sex==1,
                      48.0 + 1.06*(dataset_clean$htn09 - 152 ),
                      45.5 + 0.91*(dataset_clean$htn09 - 152) )

# EQ8: James as used by MORGAN 1994, W_kg Height_cm - CHECKED UNITS (INCORRECT VERSION OF JAMES)
dataset_clean$morgan_1994 <- ifelse( dataset_clean$sex==1,
                      1.10*dataset_clean$wtn09 - 120*(dataset_clean$wtn09/dataset_clean$htn09)^2,
                      1.07*dataset_clean$wtn09 - 148*(dataset_clean$wtn09/dataset_clean$htn09)^2 )

# EQ9: 2000 Dympna Gallagher, W_kg Height_m age_yrs sex_num - CHECKED UNITS
dataset_clean$gallagher_2000 <- ifelse( dataset_clean$sex==1,
                       dataset_clean$wtn09 - (64.5 - 848/dataset_clean$bmi09 + 0.079*dataset_clean$AGEN09 - 16.4*(dataset_clean$sex-0) + 0.05*(dataset_clean$sex-0)*dataset_clean$AGEN09 + 39*(dataset_clean$sex-0)/dataset_clean$bmi09)*dataset_clean$wtn09/100,
                       dataset_clean$wtn09 - (64.5 - 848/dataset_clean$bmi09 + 0.079*dataset_clean$AGEN09 - 16.4*(dataset_clean$sex-2) + 0.05*(dataset_clean$sex-2)*dataset_clean$AGEN09 + 39*(dataset_clean$sex-2)/dataset_clean$bmi09)*dataset_clean$wtn09/100)

# EQ10: Janmahasatian 2005, W_kg Height_m - CHECKED UNITS
dataset_clean$janmahasatian_2005 <- ifelse( dataset_clean$sex==1,
                      9270*dataset_clean$wtn09/(6680 + 216*dataset_clean$bmi09),
                      9270*dataset_clean$wtn09/(8780 + 244*dataset_clean$bmi09) )

# EQ11: ORIGINAL/TRUE FROM DXA
#dataset$original <- (dataset$DXAWBBMC09 + dataset$DXAWBLN09)/1000
#dataset_clean$original <- dataset_clean$dxa_lbm_kg # To check if we get correlation =1


print("--------------------VALIDATE THE NEW DATA FROM EQUATIONS-----------------")

print("=> ALL: Compare 10 equations with the original DXA data")
head(dataset_clean)
expected_cols <- c("hume_1996","hume_1971","james_1976","hallynck_1981",
                   "boer_1984","deurenberg_1991","zasadny_1993","morgan_1994",
                   "gallagher_2000","janmahasatian_2005")

#reference_vector <- (dataset_clean$DXAWBBMC09 + dataset_clean$DXAWBLN09)/1000 ; #Original
reference_vector <- dataset_clean$dxa_lbm_kg ; #Original
validate_lbm_columns(dataset_clean, expected_cols, reference_vector)

print("--------------------Corelation Plots Printed (Not saved)-----------------")

plot_lbm_correlations(dataset_clean, expected_cols, reference_vector)

print("--------------------Bland-ALtman Plots (saved)---------------------------")
############ BOTH (m+f)
if (TRUE){
  # Detect Colab
  is_colab <- Sys.getenv("COLAB_GPU") != "" || grepl("/content", getwd())

  # Set base dir
  base_dir <- if (is_colab) "/content" else getwd()

  for (i in seq_along(expected_cols)) {
    col_name <- expected_cols[i]
    col_data <- dataset_clean[[col_name]]  # Get column by name
    eq_number <- as.character(i)           # For labeling

    fig_name <- file.path(base_dir, paste0(col_name, "_Both.jpg"))

    baPlotFun(A = col_data,
              B = reference_vector,
              eq = eq_number,
              sex = dataset_clean$sex,
              sex_data = "Both",
              fig_name = fig_name)
  }
}
############ MALES (m)
if (TRUE){
  # Detect Colab
  is_colab <- Sys.getenv("COLAB_GPU") != "" || grepl("/content", getwd())

  # Set base dir
  base_dir <- if (is_colab) "/content" else getwd()
  for (i in seq_along(expected_cols)) {
    col_name <- expected_cols[i]

    male_idx <- which(dataset_clean$sex == 1)

    col_data_male <- dataset_clean[[col_name]][male_idx]
    reference_male <- reference_vector[male_idx]
    sex_male <- dataset_clean$sex[male_idx]

    eq_number <- as.character(i)           # For labeling

    fig_name <- file.path(base_dir, paste0(col_name, "_Males.jpg"))

    baPlotFun(A = col_data_male,
              B = reference_male,
              eq = eq_number,
              sex = sex_male,
              sex_data = "Males",
              fig_name = fig_name)
  }
}
############ FEMALES (f)
if (TRUE){
  # Detect Colab
  is_colab <- Sys.getenv("COLAB_GPU") != "" || grepl("/content", getwd())

  # Set base dir
  base_dir <- if (is_colab) "/content" else getwd()

  for (i in seq_along(expected_cols)) {
    col_name <- expected_cols[i]

    female_idx <- which(dataset_clean$sex == 2)

    col_data_female <- dataset_clean[[col_name]][female_idx]
    reference_female <- reference_vector[female_idx]
    sex_female <- dataset_clean$sex[female_idx]

    eq_number <- as.character(i)           # For labeling

    fig_name <- file.path(base_dir, paste0(col_name, "_Females.jpg"))

    baPlotFun(A = col_data_female,
              B = reference_female,
              eq = eq_number,
              sex = sex_female,
              sex_data = "Females",
              fig_name = fig_name)
  }
}
#==============================END OF ANALYSIS==================================
