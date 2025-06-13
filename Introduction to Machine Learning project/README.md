# Synthesis of the Project: *Classification Analysis of Bank Marketing Data*

## Objective

This project analyzes a Portuguese bank marketing dataset to predict whether a client subscribes to a term deposit (target variable: **`y`**). The dataset is used to evaluate classification models and address class imbalance, using both course-taught techniques and additional methods.

## Dataset Overview

- Contains client demographic and campaign-related attributes.
- Binary target variable: `y` = 'yes' or 'no' (term deposit subscription).
- No missing values.
- Notable **class imbalance**: 'no' responses heavily outnumber 'yes'.

## Methodological Approach

The analysis applies and compares several classification techniques:
- **Logistic Regression**
- **K-Nearest Neighbors (KNN)**
- **Support Vector Machines (SVM)**
- **Neural Networks**

Additional tools include:
- **Oversampling/Undersampling**: to correct class imbalance and improve model performance.
- **ROC curve & AUC**: to evaluate classifiers beyond accuracy, especially under imbalance.

## Goals

- Understand how client and campaign features influence subscription behavior.
- Assess the performance of each classifier in predicting the minority class ('yes').
- Demonstrate the impact of resampling and ROC-based evaluation on model selection.

## Conclusion

The project combines theoretical and practical methods to address a classic classification task with imbalanced data. By extending the standard workflow with resampling techniques and ROC/AUC analysis, it offers a more robust and nuanced evaluation of predictive performance.
