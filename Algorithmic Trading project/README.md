# Synthesis: Algorithmic Trading Project

## Title
**Algorithmic Trading: An Analysis of ‘A Mean-Field Game of Market-Making against Strategic Traders’**

## Overview

This project presents a structured review and critical analysis of a recent and innovative paper by Baldacci, Bergault, and Possamaï. The paper introduces a **mean-field game (MFG)** model designed to analyze interactions between a **major market-maker** and many **strategic minor market-takers**. The work is situated within the fields of algorithmic trading, market microstructure, and stochastic control theory.

## Objectives

The paper's central goal is to model optimal market-making strategies while accounting for the **strategic behavior of market-takers**. Traditional models become computationally infeasible when the number of traders increases, so the authors adopt a **major-minor MFG** framework to allow for tractable analysis and simulation.

## Key Contributions

- Introduction of a **Markovian MFG framework** to study market-making with strategic traders.
- Development of **existence and verification theorems** for equilibrium strategies.
- Numerical simulations illustrating the effect of inventory, private signals, and spread dynamics on trading behavior.
- Highlighting the role of **signal volatility** in disrupting inventory alignment and reducing market activity.

## Critical Analysis

The students identify the paper as a significant advancement over classical models (e.g., Avellaneda–Stoikov), particularly due to its **inclusion of trader strategies** that depend on private signals and inventory levels. They contextualize the work within a broader literature of MFGs in finance, especially French research institutions like École Polytechnique and Dauphine-PSL.

While the **existence of equilibrium** is rigorously proven, the **uniqueness** is only conjectured based on numerical stability. The students explore possible theoretical tools (e.g., monotonicity methods, master equation theory) for proving uniqueness but conclude that further research is needed given the jump-driven and discrete nature of the system.

## Conclusion

This project demonstrates a deep understanding of a cutting-edge application of mean-field game theory in finance. It highlights the theoretical and practical relevance of modeling strategic interactions between a market-maker and market-takers. Despite some limitations, notably the lack of uniqueness proof, the reviewed paper is positioned as a **milestone for future research** in MFG-based market modeling.
