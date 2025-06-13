# Load necessary libraries
library(tseries)   # For time series analysis
library(zoo)       # For working with indexed totally ordered observations
library(forecast)  # For forecasting functions
require(fUnitRoots)  # For unit root tests
require(urca)        # For unit root and cointegration tests

# Define custom functions

# pierce.test function:
# This function implements the Pierce (or portmanteau) test to check for residual autocorrelation.
# It accepts four arguments:
# - x: The time series data or residuals (as a zoo object) to analyze.
# - m: The number of lags to include in the test.
# - s: The step size (typically 1).
# - para: The number of model parameters to exclude from the degrees of freedom.
#
# The function:
# 1. Computes the autocorrelation function (ACF) up to m*s lags.
# 2. Calculates the test statistic as per the Pierce test formula.
# 3. Returns a list containing:
#    - `statistic`: The calculated test statistic.
#    - `p.value`: The p-value indicating the significance of the test statistic.
pierce.test <- function(x = NULL, m = NULL, s = NULL, para = NULL) {
  if (is.null(x) || is.null(m) || is.null(s) || is.null(para)) {
    stop("All parameters (x, m, s, para) must be provided")
  }
  x <- coredata(x)  # Extract numeric data from zoo object
  kk <- acf(x, m * s, plot = FALSE)  # Calculate autocorrelation function
  n <- length(x)  # Number of observations
  statistics <- 0
  for (i in 1:m) {
    statistics <- statistics + (1 / (n - (i * s))) * kk$acf[(i * s) + 1]^2
  }
  statistics <- statistics * n * (n + 2)
  p.value <- 1 - pchisq(statistics, m - para)
  return(list(statistic = statistics, p.value = p.value))  # Return values as a list
}

# Ljung.Box2 function:
# This function performs the Ljung-Box test to check residual autocorrelation.
# It accepts three parameters:
# - x: The time series data or residuals (as a zoo object) to analyze.
# - lag: The number of lags to include in the test.
# - param: The number of model parameters to adjust for.
#
# The function:
# 1. Computes the ACF up to the specified `lag` (ignoring the first autocorrelation).
# 2. Calculates the test statistic using the Ljung-Box formula.
# 3. Returns a list containing:
#    - `statistic`: The test statistic.
#    - `gdl`: The degrees of freedom (lag minus the number of parameters).
#    - `pvalue`: The p-value indicating the significance of the test.
Ljung.Box2 <- function(x, lag = 15, param = 1) {
  if (NCOL(x) > 1) stop("x is not a vector or univariate time series")
  x <- coredata(x)  # Extract numeric data from zoo object
  cor <- acf(x, lag.max = lag, plot = FALSE)
  n <- length(x)
  PARAMETER <- lag - param
  obs <- cor$acf[2:(lag + 1)]
  STATISTIC <- n * (n + 2) * sum((1 / seq(n - 1, n - lag)) * obs^2)
  PVAL <- 1 - pchisq(STATISTIC, lag - param)
  return(list(statistic = STATISTIC, df = PARAMETER, pvalue = PVAL))
}

# Ljung.Box function:
# The Ljung.Box function provides an extended Ljung-Box test across multiple lag lengths.
# It accepts four parameters:
# - x: The time series data or residuals (as a zoo object) to analyze.
# - maxlag: The maximum number of lags to include in the analysis.
# - par: The number of model parameters to adjust for.
# - all: If TRUE, computes the test across all lags up to `maxlag` and plots the p-values.
#
# The function:
# 1. Calls `Ljung.Box2` for each lag up to `maxlag`, excluding the specified `par` parameters.
# 2. Constructs a matrix containing the test statistic, degrees of freedom, and p-value at each lag.
# 3. Optionally plots the p-values as a bar chart if `all` is TRUE.
# 4. Returns a matrix containing the statistics for each lag.
Ljung.Box <- function(x, maxlag = 24, par = 0, all = TRUE) {
  if (all == TRUE) {
    LB <- matrix(0, nrow = maxlag, ncol = 3)
    dimnames(LB) <- list(NULL, c("statistic", "df", "pvalue"))
    for (i in (par + 1):maxlag) {
      lbi <- Ljung.Box2(x, lag = i, param = par)
      if (length(lbi) == 3) {
        LB[i, ] <- c(round(lbi$statistic, 2), lbi$df, lbi$pvalue)
      } else {
        warning(paste("Ljung.Box2 did not return expected results for lag", i))
      }
    }
    plot(LB[, 3], type = "h", lwd = 4, ylim = c(0, 1), ylab = "P-Value", xlab = "Lag")
    abline(h = 0.05, lty = 2, col = "red")
    abline(h = 0)
  }
  if (all == FALSE) LB <- Ljung.Box2(x, lag = maxlag, param = par)
  
  return(LB)
}

# stat.mod function:
# The stat.mod function provides comprehensive model diagnostics for a fitted time series model.
# It takes three arguments:
# - fit.y: The fitted model object (e.g., ARIMA or other time series model).
# - maxlag: The maximum number of lags for the Ljung-Box test to check residual autocorrelation.
# - ordini: The number of model parameters.
# 
# The function calculates and prints:
# 1. Coefficients, standard errors, t-statistics, and p-values of the model parameters.
# 2. Model diagnostics, such as the residual variance (sigma2), Akaike Information Criterion (AIC), and log-likelihood.
# 3. Results from the Ljung-Box test to check for autocorrelation in the residuals up to `maxlag` lags.
# 
# The results are formatted and displayed for clear interpretation of the model's fit and adequacy.
stat.mod=function(fit.y,maxlag=15,ordini=0)
{
  Coef=fit.y$coef
  Std.Err=sqrt(diag(fit.y$var.coef))
  tstat=fit.y$coef/sqrt(diag(fit.y$var.coef))
  pval=2*pnorm(abs(tstat),lower.tail=F)
  tabella=data.frame(Coef,Std.Err,tstat,pval)
  sigma2=fit.y$sigma2
  names(sigma2)="sigma2"
  AIC=fit.y$aic
  names(AIC)="AIC"
  loglik=fit.y$loglik
  names(loglik)="loglik"
  cc1=Ljung.Box(fit.y$resid,maxlag,ordini)
  riga0="------------- Parameters estimation -----------"
  names(riga0)=" "
  riga="-----------------------------------------------"
  names(riga)=" "
  print(riga0, quote=F)
  print(tabella)
  print(riga, quote=F)
  riga2="------------- Ljung-Box Test ---------------"
  names(riga2)=" "
  print(c(sigma2,AIC,loglik))
  print(riga2, quote=F)
  print(t(cc1))
}

# Qtests function:
# This function performs a series of Ljung-Box tests up to lag `k` for the given series.
# It accepts three arguments:
# - series: The time series data to analyze.
# - k: The maximum lag to test.
# - fitdf: The number of model parameters to adjust for.
#
# The function:
# 1. Applies the Ljung-Box test at each lag up to `k`.
# 2. Returns a matrix containing the lag and the corresponding p-value.
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

# adfTest_valid function:
# This function performs the Augmented Dickey-Fuller (ADF) test for stationarity, increasing lags until no autocorrelation in residuals.
# It accepts three arguments:
# - series: The time series data to analyze.
# - kmax: The maximum number of lags to test.
# - type: The type of test ("c" for constant, "ct" for constant and trend, etc.).
#
# The function:
# 1. Performs the ADF test with increasing lags until residuals show no autocorrelation.
# 2. Prints whether the residuals are autocorrelated or not at each step.
# 3. Returns the final ADF test result.
adfTest_valid <-
  function(series,kmax,type){ #ADF tests until no more autocorrelated residuals
    k <- 0
    noautocorr <- 0
    while (noautocorr==0){
      cat(paste0("ADF with ",k, " lags: residuals OK? "))
      adf <- adfTest(series,lags=k,type=type)
      pvals <- Qtests(adf@test$lm$residuals,24,fitdf=length(adf@test$lm$coefficients))[,2]
      if (sum(pvals<0.05,na.rm=T) == 0) {
        noautocorr <- 1; cat("OK \n")}
      else cat("NO \n")
      k <- k + 1
    }
    return(adf)
  }

# modelchoice function:
# This function fits an ARIMA model and checks the significance of AR and MA terms and the residuals' autocorrelation.
# It accepts four arguments:
# - p: The AR order.
# - q: The MA order.
# - data: The time series data.
# - k: The maximum number of lags for the Ljung-Box test.
#
# The function:
# 1. Tries to fit an ARIMA model with the given orders.
# 2. Checks the significance of AR and MA terms and residuals' autocorrelation.
# 3. Returns a list with the model's coefficients, standard errors, and Ljung-Box test results.
modelchoice <- function(p, q, data, k = 24) {
  estim <- try(arima(data, order = c(p, 0, q), optim.control = list(maxit = 20000)), silent = TRUE)
  
  if (class(estim) == "try-error") {
    return(c("p" = p, "q" = q, "arsignif" = NA, "masignif" = NA, "resnocorr" = NA, "ok" = NA))
  }
  
  coefs <- coef(estim)
  se <- sqrt(diag(vcov(estim)))
  
  arsignif <- if (p == 0) NA else {
    ar_coefs <- coef(estim)[1:p]
    ar_se <- sqrt(diag(vcov(estim)))[1:p]
    sum(2 * (1 - pnorm(abs(ar_coefs / ar_se))) <= 0.05) == p
  }
  
  masignif <- if (q == 0) NA else {
    ma_coefs <- coef(estim)[(p + 1):(p + q)]
    ma_se <- sqrt(diag(vcov(estim)))[(p + 1):(p + q)]
    sum(2 * (1 - pnorm(abs(ma_coefs / ma_se))) <= 0.05) == q
  }
  
  resnocorr <- sum(Box.test(estim$residuals, lag = k, type = "Ljung-Box")$p.value <= 0.05, na.rm = TRUE) == 0
  checks <- c(arsignif, masignif, resnocorr)
  ok <- as.numeric(sum(checks, na.rm = TRUE) == (3 - sum(is.na(checks))))
  
  ljung_box <- Ljung.Box(estim$residuals, maxlag = k, par=p+q)
  
  # Print coefficients and their standard errors
  cat(paste0("ARMA(", p, ",", q, ") coefficients:\n"))
  print(coefs)
  cat("Standard errors:\n")
  print(se)
  
  # Print Ljung-Box test result
  cat("Ljung-Box test result:\n")
  print(ljung_box)
  
  return(list(p = p, q = q, arsignif = arsignif, masignif = masignif, resnocorr = resnocorr, ok = ok,
              coefficients = coefs, se = se, ljung_box = ljung_box))
}

# armamodelchoice function:
# This function finds the best ARMA model within given p and q ranges by checking each model's validity.
# It accepts three arguments:
# - pmax: The maximum AR order to test.
# - qmax: The maximum MA order to test.
# - y_df: The differenced time series data.
#
# The function:
# 1. Tests all combinations of AR and MA orders up to pmax and qmax.
# 2. Stores results for each model's coefficients, standard errors, and diagnostic checks.
# 3. Returns a data frame with the comparison of all tested models.
armamodelchoice <- function(pmax, qmax, y_df) {
  pqs <- expand.grid(p = 0:pmax, q = 0:qmax)
  results <- list()
  
  for (i in 1:nrow(pqs)) {
    p <- pqs$p[i]
    q <- pqs$q[i]
    cat(paste0("Computing ARMA(", p, ",", q, ") \n"))
    result <- modelchoice(p, q, coredata(y_df))
    results[[i]] <- result
  }
  
  comparison <- do.call(rbind, results)
  comparison <- as.data.frame(comparison, stringsAsFactors = FALSE)
  comparison[, c("p", "q", "arsignif", "masignif", "resnocorr", "ok")] <- 
    lapply(comparison[, c("p", "q", "arsignif", "masignif", "resnocorr", "ok")], as.numeric)
  
  return(comparison)
}

########
# PART 1
########

# Set working directory to the location of the data file
setwd("C:/Users/Marco Plazzogna/Desktop/MAGISTRALE/PRIMO ANNO/SECONDO SEMESTRE/TIME SERIES/Assignment")

# Read data from a CSV file, skipping the first 4 lines and assuming ";" as the separator
data <- read.csv("valeurs_mensuelles.csv", sep = ";", skip = 4, header = FALSE)

# Extract the second column (assuming it contains the time series data)
serie <- data[, 2]

# Convert the series into a zoo object starting from January 2000 with monthly frequency
y <- zoo(serie, order.by = as.Date(seq(from = as.Date("2000-01-01"), by = "month", length.out = length(serie))))

# Plot the original time series
plot(y, main = "Original Time Series", ylab = "Value", xlab = "Time")

# Perform Augmented Dickey-Fuller test to check for stationarity
adf <- adfTest_valid(y,24,"ct")
print(adf)

adf <- adfTest_valid(y,24,"c")
print(adf)

# Perron-Phillips test to check for stationarity
pp.test(y)
y.pp = ur.pp(y, type = "Z-tau", model = "trend")
summary(y.pp)

# Optionally, add ACF and PACF plots for the series to inspect autocorrelations
acf(coredata(y), main = "ACF of Series")
pacf(coredata(y), main = "PACF of Series")

# Difference the series to achieve stationarity
y_df = diff(y)

# Plot the differenced series
plot(y_df, main = "Differenced Time Series", ylab = "Value", xlab = "Time")

# Plot ACF and PACF for the differenced series
acf(coredata(y_df), main = "ACF of Differenced Series")
pacf(coredata(y_df), main = "PACF of Differenced Series")

# Set maximum AR and MA orders for model selection
pmax = 5
qmax = 1

########
# PART 2
########

# Find the best ARMA models within the specified p and q ranges
armamodels <- armamodelchoice(pmax,qmax, y_df)

# Select models that are well-adjusted and valid
selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),] #models well adjusted and valid
selec[, c("p", "q", "arsignif", "masignif", "resnocorr", "ok")]

# Fit several ARMA models
models <- list()
aic_values <- numeric()
bic_values <- numeric()

# Find the best model with AIC and BIC
for (i in 1:nrow(selec)) {
  p <- selec$p[i]
  q <- selec$q[i]
  fit <- tryCatch({
    Arima(coredata(y_df), order = c(p, 0, q))
  }, error = function(e) NULL)
  
  if (!is.null(fit)) {
    models[[paste0("ARMA(", p, ",", q, ")")]] <- fit
    aic_values <- c(aic_values, fit$aic)
    bic_values <- c(bic_values, BIC(fit))
  }
}

# Create a data frame to compare models
comparison <- data.frame(
  Model = names(models),
  AIC = aic_values,
  BIC = bic_values
)

# Print the model comparison table
print(comparison)

# Find the best models based on AIC and BIC
best_aic_model_index <- which.min(comparison$AIC)
best_bic_model_index <- which.min(comparison$BIC)

best_aic_model_name <- names(models)[best_aic_model_index]
best_bic_model_name <- names(models)[best_bic_model_index]

best_aic_model <- models[[best_aic_model_name]]
best_bic_model <- models[[best_bic_model_name]]

# Display summaries of both models
cat("Model with the best AIC:\n")
summary(best_aic_model)

cat("\nModel with the best BIC:\n")
summary(best_bic_model)

best_model = best_aic_model

# Check residuals of the best AIC model
cat("\nResidual diagnostics for the best AIC model:\n")
checkresiduals(best_model)

# Perform Ljung-Box test on residuals of the best AIC model
cat("\nLjung-Box test for residuals (best AIC model):\n")
Ljung.Box(best_model$residuals, maxlag = 24, par=1)

########
# PART 3
########

# Extract the last 30 observations from the original series
last_30_index <- tail(index(forecasts$x), 30)
last_30_data <- window(forecasts$x, start = last_30_index[1])

h <- 2

# Generate forecasts
forecasts <- forecast(best_model, h = h, level = 0.95)

# Extract forecast values and confidence intervals
forecast_values <- forecasts$mean
lower_bound <- forecasts$lower  # 95% lower bound
upper_bound <- forecasts$upper  # 95% upper bound

# Assuming you have the index of the last fitted value and the first forecast value
last_fitted_time <- tail(time(forecasts$fitted), 1)
first_forecast_time <- head(time(forecasts$mean), 1)

last_fitted_value <- tail(forecasts$fitted, 1)
first_forecast_value <- head(forecasts$mean, 1)

# Plot the last 30 observations with forecasts and confidence intervals
plot(last_30_data, main = "Forecast with 95% Confidence Intervals", xlab = "Time", ylab = "Value", type = "l", xlim=c(260,290))
lines(forecasts$fitted, col = "blue")
points(last_30_data, col="black")
points(forecasts$fitted, col="blue")
points(time(forecasts$mean), forecasts$mean, col = "red")
arrows(time(forecasts$mean), lower_bound, time(forecasts$mean), upper_bound, code = 3, angle = 90, length = 0.1, col = "red", lty=3)
lines(forecasts$mean, col="red")

# Draw a link between the last blue point and the first red point
segments(last_fitted_time, last_fitted_value, first_forecast_time, first_forecast_value, col="red")

# Add legend
legend("bottomleft", legend = c("Series", "Fitted", "Forecast", "95% CI"), col = c("black", "blue", "red", "red"), lty = c(1, 1, 1, 3), pch = c(NA, NA, 1, NA), bty = "o", cex = 0.7)

############
# EXTRA PART
############

# see the predictions on the original series and compare them with the actual observations

# Remove the last two observations from the series y
y_new <- window(y, end = index(y)[length(y) - 2])

# Estimate an ARIMA(0,1,1) for the new series y
model_arima <- Arima(y_new, order = c(0, 1, 1))

# Set the forecast horizon
h <- 2

# Generate forecasts for h=2
forecasts <- forecast(model_arima, h = h, level = 95)

# Extract forecast values and confidence intervals
forecast_values <- forecasts$mean
lower_bound <- forecasts$lower  # 95% lower bound
upper_bound <- forecasts$upper  # 95% upper bound
forecast_time <- seq(from = end(y_new), by = "month", length.out = h+1)[-1]

# Filter the original time series starting from 2021
start_date <- as.Date("2021-01-01")
y_filtered <- window(y, start = start_date)

# Plot the original series starting from 2021 with the forecasts
plot(y_filtered, main = "Original Series with Forecasts and 95% Confidence Intervals (Starting from 2021)", 
     ylab = "Value", xlab = "Time", col = "black", ylim=c(80,120))

# Add the forecast values to the plot
points(forecast_time, forecast_values, col = "red", type = "o")

# Add the confidence intervals
for (i in 1:length(forecast_values)) {
  arrows(forecast_time[i], lower_bound[i], forecast_time[i], upper_bound[i], code = 3, angle = 90, length = 0.05, col = "blue")
}

# Add legend
legend("topleft", legend = c("Original Series", "Forecast", "95% CI"), 
       col = c("black", "red", "blue"), lty = c(1, 1, NA), 
       pch = c(NA, NA, 16), bty = "o", cex = 0.7)

