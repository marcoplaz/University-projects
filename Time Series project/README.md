# Project Summary: Time Series Analysis of French Defense Industry Production

## Objective

This project analyzes the seasonally and working-day adjusted industrial production index for the manufacture of weapons and ammunition in France (NAF rev. 2, class 25.40), using monthly data from January 2000 to January 2024. The goal is to model and forecast this index using ARMA/ARIMA models, focusing on identifying the best-fitting model and its predictive performance.

## Methodology

### Data Preprocessing

- The original time series shows a general upward trend with spikes in volatility, notably around 2008–2010 and post-2020.
- Augmented Dickey-Fuller (ADF) and Phillips-Perron (PP) tests yielded mixed results.
- The series was differenced to induce stationarity.
- ACF and PACF plots of the differenced series showed no signs of strong seasonality.

### Model Selection

- A grid search over ARMA(p, q) models was conducted with \( p \leq 5 \) and \( q \leq 1 \).
- The selection process checked:
  - Significance of AR and MA parameters,
  - Residual independence using the Ljung-Box test,
  - Model fit via AIC and BIC.
- The best model identified was **ARIMA(0,1,1)**.

## Diagnostics

- Residuals from the ARIMA(0,1,1) model are approximately normal and centered at zero.
- Ljung-Box tests confirm no residual autocorrelation.
- The model adequately captures the underlying dynamics of the differenced series.

## Forecasting

- 1-step and 2-step forecasts were computed with 95% confidence intervals.
- Forecast behavior aligns with an MA(1) model:
  - The first forecast uses the last residual;
  - Later forecasts revert to the estimated process mean (≈ 0.1422).

## Further Considerations

- Structural breaks (2008 financial crisis, COVID-19) suggest the benefit of modeling sub-periods separately.
- An open research question is whether an exogenous series \( Y_t \) can improve forecasts for the production index:
  - Granger causality tests could assess whether \( Y_t \) helps predict \( X_t \).
  - Real-world example: Defense spending announcements might inform future production.
