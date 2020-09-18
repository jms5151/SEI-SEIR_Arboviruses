# run linear regressions to compare correlations between model predictions
# and observations with socio-ecological factors 

# load data
SE_data <- read.csv("data/socio_ecological_factors.csv")

# create a vector of socio-ecological factors
SE_factors <- colnames(SE_data)[4:18]

# create a vector of response variables
Resp_vars <- c("vectors", "cases")

# create dataframe to store linear regression results
SE_analysis_df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(SE_analysis_df) <- c("SE_factor", "Response_variable", "Rsq", "Pvalue")

# run linear regressions
for(i in Resp_vars){
  df <- subset(SE_data, Response_variable == i)
  for(j in SE_factors){
    reg_out <- lm(df$Correlation_coefficient ~ df[,j])
    rsq <- format(round(summary(reg_out)$r.squared, 2), nsmall = 2)
    pval <- format(round(summary(reg_out)$coefficients[2,4], 2), nsmall = 2)
    SE_analysis_df[nrow(SE_analysis_df)+1, ] <- c(j, i, rsq, pval)
  }
}
