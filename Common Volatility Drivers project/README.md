# Synthesis of the Project: *Global COVOL Risk Mitigation Portfolio*

## Objective

This project replicates and extends part of the analysis from the paper *"What are the events that shake our world? Measuring and hedging global COVOL"* by Engle & Campos-Martins (2023). The aim is to build a **risk mitigation portfolio** that hedges against global volatility shocks (global COVOL), particularly in an asset management context.

## Dataset and Preprocessing

- Monthly returns from a cross-section of **equity indices** and **factor portfolios**.
- Period: 2000–2023.
- Returns are cleaned, normalized, and converted to excess returns for GARCH modeling.

## Methodology

1. **Step 1: GARCH(1,1) Estimation**
   - Each series is fitted with an AR(1)-GARCH(1,1) model.
   - Extracts conditional variances as volatility proxies.

2. **Step 2: COVOL Factor Extraction**
   - Applies Principal Component Analysis (PCA) to the **log-volatilities**.
   - First principal component is interpreted as the **global COVOL factor**.
   - Loadings are normalized to define **COVOL betas**.

3. **Step 3: Portfolio Optimization**
   - Constructs **COVOL-mitigating portfolios** by minimizing exposure to the COVOL factor under return and volatility constraints.
   - Compared to:
     - **Equal-weighted portfolio**
     - **Minimum variance portfolio**
     - **COVOL-agnostic constrained portfolio**

## Results

- The COVOL-mitigating portfolio shows **lower drawdowns** and improved **Sharpe ratios**, especially during high-volatility periods (e.g., 2008 crisis, COVID-19).
- Outperforms benchmarks in terms of stability and downside protection.
- Confirms the **economic relevance** of global COVOL for asset allocation and stress-resilient investing.

## Conclusion

This project successfully replicates the core logic of the Engle & Campos-Martins methodology and validates the benefits of incorporating global COVOL information into portfolio construction. The resulting portfolio shows enhanced resilience during global volatility spikes, making it a valuable tool for institutional risk management.

*Note*: the code for this part is not available given the contract signed with Lombard Odier about the property of such coding materials.
