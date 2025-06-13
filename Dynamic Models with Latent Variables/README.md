# Synthesis: Deep Machine Learning Theory Project

## Title
**Learning Lenses for Inverse Graphics: A Critical Review**

## Overview

This project presents a structured review of the paper *“Learning Lenses for Inverse Graphics”* by Kulkarni et al., which proposes a novel deep generative model for vision tasks based on **analysis-by-synthesis**. The goal of the original paper is to build a model that can infer interpretable latent variables (such as pose, lighting, or shape) from images, and use them to reconstruct the observed data. This work lies at the intersection of **computer vision**, **probabilistic modeling**, and **deep generative learning**.

## Objectives

The key objective is to emulate the **inverse graphics pipeline**: given an image, the model should infer the underlying parameters that describe how the image was generated. The authors design a **compositional generative model** with interpretable latent codes and train it via **variational inference** combined with deep neural networks.

## Key Contributions

- Introduction of an **encoder-decoder architecture** with convolutional neural networks for the image encoder and a compositional decoder that mirrors graphics engines.
- Use of a **variational autoencoder (VAE)** framework to perform inference and learning, extended to support structured latent variables.
- Demonstration of the model’s performance on **face images** and **3D chair renderings**, showing that the learned representations are disentangled and interpretable.
- Visualization of latent space traversals to validate semantic alignment of latent dimensions (e.g., changing azimuth or elevation corresponds to smooth changes in the output).

## Critical Analysis

The students appreciate the **originality** of the approach in combining deep learning with probabilistic graphics models. They emphasize the modular structure of the generative model as a key advantage for interpretability and generalization.

However, they also highlight several **limitations**:
- The model assumes strong prior knowledge of rendering parameters and relies on synthetic datasets, which limits generalization to real-world images.
- The inference is only approximate and depends on the quality of the variational posterior, which might be too simple to capture the true posterior.
- The training process is computationally intensive and highly sensitive to network architecture and optimization.

## Observations

- The model successfully disentangles generative factors such as pose, lighting, and shape.
- It provides a meaningful test case for bridging **graphics and vision** through learned representations.
- While not directly scalable to high-resolution or natural images, the approach opens up research directions in **interpretable and compositional generative modeling**.

## Conclusion

This review reflects a strong understanding of the theoretical and architectural choices made in the paper. The students underline the model's **contribution to interpretable deep learning** and its role in advancing **inverse graphics** through deep generative models. Despite its reliance on synthetic datasets, the paper is considered a **pioneering effort** in combining structured latent spaces with deep inference frameworks.
