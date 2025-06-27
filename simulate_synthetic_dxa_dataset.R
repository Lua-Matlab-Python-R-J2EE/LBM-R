################################################################################
# R SCRIPT: To simulate synthetic DXA dataset with errors
#
#
# AUTHOR: Dr. Tanuj Puri, DATE: March 2025
#
#
# DESCRIPTION: This script generates a synthetic dataset simulating total body 
# DXA data for 2000 individuals. It inclues demographic and body composition
# variables, modeled using realistic distributions for males and females aged 
# 60-65 years. The script also introduces deliberate errors in atleast 500
# data cells to simulate real-world data entries or collection anomalies
# such as missing, unkown, zero, negatove, or non-numeric values.
# 
#
# VARIABLES USED:  
# p_ID:       Particupant ID (integer, unique)
# AGEN09:     Age of participant (years: range 60-65)
# sex:        Biological sex (1=Male, 2=Female)
# wtn09:      Body wight (kg)
# htn09:      Height (cm)
# dxa_lbm_kg: Lean body mass in kg (derived from DXAWBLN09)
# DXASBBMC09: Subtotal bone mineral content **excluding head** (g)
#             Codes: 7777 = Not attended Clinic/test not done, 9999 = Unknown
# DXAWBBMC09: Whole Body Bone Mineral Content (g)
#             Codes: 7777 = Not attended Clinic/test not done, 9999 = Unknown
# DXASBFT09:  Subtotal Fat mass **excluding head** (g)
#             Codes: 777777 = Not attended Clinic/test not done, 999999 = Unknown
# DXAWBFT09:  Whole Body Fat Mass (g)
#             Codes: 777777 = Not attended Clinic/test not done, 999999 = Unknown
# DXASBLN09:  Subtotal Lean mass **excluding head** (g)
#             Codes: 777777 = Not attended Clinic/test not done, 999999 = Unknown
# DXAWBLN09:  Whole Body Lean Mass (g)
#             Codes: 777777 = Not attended Clinic/test not done, 999999 = Unknown
# DXASBMS09:  Subtotal Body Mass (Fat + Lean, **excluding head)** (g)
# DXAWBMS09:  Whole Body Mass (Fat + Lean) (g)
# bmi09:      Body Mass Index (kg/m²), calculated from weight and height
# inf:        Random binary indicator (e.g., for infection status or flag)
#
#
# STATISTICAL ASSUMPTIONS TAKEN FROM:  
# Tanuj Puri, Glen M Blake, Comparison of ten predictive equations for estimating 
# lean body mass with dual-energy X-ray absorptiometry in older patients, British 
# Journal of Radiology, Volume 95, Issue 1133, 1 May 2022, 20210378, 
# https://doi.org/10.1259/bjr.20210378 
#
# Males:
#   Age      : Mean = 63 yrs (range 60–65)
#   Weight   : Mean = 85.3 kg (range 50.6–128.5)
#   Height   : Mean = 175.3 cm (range 156.7–195.7)
#
# Females:
#   Age      : Mean = 63 yrs (range 60–65)
#   Weight   : Mean = 72.3 kg (range 38.0–136.5)
#   Height   : Mean = 162.2 cm (range 144.3–179.0)
#
# Combined:
#   Age      : Mean = 63 yrs (range 60–65)
#   Weight   : Mean = 78.5 kg (range 38.0–136.5)
#   Height   : Mean = 168.5 cm (range 144.3–195.7)
#
#
# ERROR INJECTION: At least 500 data entries across random row/columns are 
# replaced with Zero values, negatove values, or Non-numeric strings. These 
# simulate real-world data quality issues, allowing for testing of data cleaning 
# and validation  procedures.
# 
#
# OUTPUT: 
#   A data frame with 2000 rows and 16 columns.
#   Optionally saved as a CSV file for use in further analysis.
#
#
################################################################################
#
#
# HELPER FUNCTION 1: To generate random values with bounds and mean
generate_values <- function(n, mean, min, max, sd=NULL) {
    if (is.null(sd))
    {
    # Estimate SD as (max-min)/6, assuming normal distribution covering 99.7% values
     sd <- (max-min)/6 
    }
    # Generate n random numbers from normal distribution with specified mean/SD
    vals <- rnorm(n, mean=mean, sd=sd)

    # Clamp vals to stay between [min, max] range
    vals <- pmin( pmax(vals, min) , max)

    # Rounds all values in the vals to 1 decimal place before returning
    invisible(round(vals,1)) 
}
#
#
#-------------------------------------------------------------------------------
#
#
# HELPER FUNCTION 2: To generate  dummy DXA values with approximate logic
generate_dxa_values <- function(n, min_val, max_val, mising_code, unknown_code){
    # Generate n values from uniform distribution between min/max
    #values <- runif( n, min=min_val, max=max_val )
    values <- rnorm( n, mean=0.5*(min_val+max_val), sd=(max_val-min_val)/6 )

    # Inject 2% of n to be missing enteries
    idx_missing <- sample( 1:n, size=round(0.02*n) ) 

    # Inject 2% of n to be unknown enteries, that are not missing entries
    idx_unknown <- sample( setdiff(1:n, idx_missing), size=round(0.02*n) )
    
    # Replace selected positions with respective missing/unknown codes
    values[idx_missing] <- mising_code
    values[idx_unknown] <- unknown_code
    
    invisible(round(values, 2)) # 2 deciaml places
}   
#
#
#-------------------------------------------------------------------------------
# 
#
# HELPER FUNCTION 3: To summarize rows with error codes in DXA columns
vew_error_rows <- function( df, dx_columns, error_codes, show_sample=TRUE, sample_n=3 ){
  # Find rows where any DXA column has an error code
  error_row_logical <- apply( df[dx_columns], 1, function(row) any(row %in% error_codes) )

  # count total rows with errors using two different ways
  total_rows1 <- nrow( df[error_row_logical, ] )
  total_rows2 <- sum( error_row_logical )

  # print counts
  cat("Total rows1: ",total_rows1, "\n")
  cat("Total rows2: ",total_rows2, "\n")

  # Print a sample of these rows
  if (show_sample && total_rows1>0)
  {
      cat("Sample rows with error codes: \n")
      print( head( df[error_row_logical,], sample_n) )
  } else if (total_rows1==0) 
  {
      cat("No rows with error codes found. \n")
  }
  invisible(df[error_row_logical, ])                                                                                       
}
#
#
#------------------------------------------------------------------------------
#
#
# HELPER FUNCTION 4: To save the dataframe as a CSV file, if filename is provided.
#                    It also returns the dataframe. Checks if the code is run on
#                    local enviorment or colab, and saves the file accordingly.                             
save_to_csv <- function(df, filename = NULL, row_names = FALSE) {  
  # Detect colab
  is_colab <- Sys.getenv("COLAB_GPU") != "" || grepl("/content", getwd()) 
  
  base_dir <- if (is_colab) "/content" else getwd() #set base dir based on env
  
  if (!is.null(filename)) {
    full_path <- file.path(base_dir, filename) 
    write.csv(df, file = full_path, row.names = row_names)# Save the file
    cat("Data saved to:", full_path, "\n")
  } else {
    cat("No filename provided; data not saved.\n")
  }

  invisible(df)
}
#
#
################################################################################
# Load package for data manupulation
library(dplyr)
#---------------------------------
# Set seed for reproducibility of random number generators
set.seed(123)

# Total number of rows required
n_total <- 2000

# Number or females/males based on approx proportions
n_females <- 811
n_males   <- 747
n_other   <- n_total - (n_males + n_females)

# Generate sex vector (1 = Male, 2 = Female)
sex <- c( rep(1,n_males), rep(2,n_females), sample( 1:2,n_other,replace=TRUE) )

# Generate age using generate_values() fn, given min/max/mean
age <- generate_values( n_total, mean=63, min=60, max=65 ) 

# Generate weight based on sex and using generate_values() fn, given min/max/mean
weight <- ifelse( sex==1, generate_values( n_total, mean=85.3, min=50.6, max=128.5 ), 
                          generate_values( n_total, mean=72.3, min=38.0, max=136.5 ))

# Generate height based on sex and using generate_values() fn, given min/max/mean
height <- ifelse( sex==1, generate_values( n_total, mean=175.3, min=156.7, max=195.7 ), 
                          generate_values( n_total, mean=162.2, min=144.3, max=179.0 ))

# Compute BMD = weight(kg)/ (height(m))^2, round to 2 decimal places
bmi09 <- round( weight / (height/100)^2 , 2) 

# Generate DXA field values for bone mineral content (g)
DXAWBBMC09 <- generate_dxa_values(n_total, 10, 30, 7777, 9999)#Whole
DXASBBMC09 <- generate_dxa_values(n_total,  8, 25, 7777, 9999)#No-head

# Generate DXA field values for fat mass (g)
DXAWBFT09 <- generate_dxa_values(n_total, 10000, 35000, 777777, 999999)#Whole
DXASBFT09 <- generate_dxa_values(n_total, 10000, 30000, 777777, 999999)#No-head

# Generate DXA field values for lean mass (g)
DXAWBLN09 <- generate_dxa_values(n_total, 25000, 70000, 777777, 999999)#Whole
DXASBLN09 <- generate_dxa_values(n_total, 20000, 60000, 777777, 999999)#No-head

# Generate DXA field values for body mass (g) - fat mass + lean mass
DXAWBMS09 <- DXAWBFT09 + DXAWBLN09 #whole
DXASBMS09 <- DXASBFT09 + DXASBLN09 #No-head

dxa_lbm_kg <- round(DXAWBLN09/1000, 2) #in kg , round to 2 decimal

# Generate wtn09 and htn09 from weight and height, respectively
wtn09 <- weight
htn09 <- height

# generate p_ID and AGEN09 from age and ID, respectively
p_ID   <- 1:n_total
AGEN09 <- age 

# Initialize the data frame so that the data can be saved as Excel/CSV file
df <- data.frame(
  p_ID       = p_ID,
  AGEN09     = AGEN09,
  sex        = sex,
  wtn09      = wtn09,
  htn09      = htn09,
  dxa_lbm_kg = dxa_lbm_kg,
  DXASBBMC09 = DXASBBMC09,
  DXASBFT09  = DXASBFT09,
  DXASBLN09  = DXASBLN09,
  DXASBMS09  = DXASBMS09,
  DXAWBBMC09 = DXAWBBMC09,
  DXAWBFT09  = DXAWBFT09,
  DXAWBLN09  = DXAWBLN09,
  DXAWBMS09  = DXAWBMS09,
  bmi09      = bmi09,
  inf        = sample(0:1, n_total, replace = TRUE)  
)

# View sample of dataset
# head(df)
# df[df$DXASBBMC09 == 7777, ]
# df[df$DXASBLN09 == 777777, ]
# head( df[, dx_columns], 15 )
# table(df$DXAWBBMC09 %in% error_codes)

# View other samples with error codes
dx_columns  <- c("DXAWBBMC09", "DXASBBMC09", "DXASBFT09", "DXASBLN09", 
                 "DXASBMS09", "DXAWBFT09", "DXAWBLN09")
error_codes <- c(7777, 9999, 777777, 999999)
error_rows  <- vew_error_rows( df, dx_columns, error_codes, show_sample=FALSE, sample_n=6)

# Save to CSV if filename is provided. It also returns the dataframe.
df_new <- save_to_csv(df, "simulated_dxa_dataset.csv")
cat("-------------------------------------------------------------------------")
head(df_new)
cat("-------------------------------------------------------------------------")
head(df)
cat("-------------------------------------------------------------------------")
# ------------------------------------ END ------------------------------------ 
