# Synthesis of the Report: *Dynamic Models with Latent Variables*

## Objective
This report analyzes the paper "Stationarity and Ergodicity of Markov Switching Positive Conditional Mean Models" (Aknouche & Francq, 2022), focusing on the theoretical properties and numerical validation of MS-PCM models for non-negative time series.

## Theoretical Contributions

- The paper studies Markov Switching Positive Conditional Mean (MS-PCM) models, generalizing ACD, INGARCH, and Beta models through regime-switching mechanisms.
- It establishes conditions for:
  - Stationarity and ergodicity via spectral radius tests on transition-dependent matrices.
  - Existence of finite moments, especially under past-regime dependent switching.
- It distinguishes three formulations:
  - Past-regime dependent switching
  - Present-regime dependent switching
  - Present-regime mean-dependent switching
- All results extend to non-linear dynamics under Lipschitz assumptions.

## Critical Analysis

**Strengths:**
- Strong mathematical rigor and integration of spectral methods.
- Broad applicability to positive-valued time series.

**Limitations:**
- No empirical dataset used.
- Estimation procedures only briefly mentioned.
- Some assumptions not fully detailed (e.g., stochastic ordering).

## Numerical Validation

Simulations confirm the theoretical results:

- **Poisson MS-INGARCH**:
  - Empirical mean ≈ 6.13, variance ≈ 8.13.
  - Matches theoretical expectations.

- **Negative Binomial MS-INGARCH**:
  - Empirical mean ≈ 6.14, variance ≈ 12.25.
  - Captures overdispersion well.

These simulations validate stationarity, ergodicity, and the model’s flexibility.

## Conclusion

The report confirms the paper’s theoretical soundness and its relevance for modeling regime-switching in non-negative time series. It sets a foundation for future work on estimation techniques and extensions to higher-order or multivariate settings.
