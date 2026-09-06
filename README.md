This notebook performs single-cell RNA sequencing (scRNA-seq) data integration and evaluates various clustering and embedding methods using `scib-metrics`.

## Table of Contents

- [Introduction](#introduction)
- [Data](#data)
- [Preprocessing](#preprocessing)
- [Integration Methods](#integration-methods)
- [Clustering Methods](#clustering-methods)
- [Benchmarking](#benchmarking)
- [Visualization](#visualization)
- [Results](#results)

## Introduction

This project aims to benchmark different data integration and clustering techniques on scRNA-seq data to assess their performance in correcting batch effects and preserving biological variation. We use the `scib-metrics` library for a comprehensive evaluation.

## Data

We use the `pbmc3k` dataset from `scanpy`. The data is artificially batched and cell types are simulated for demonstration purposes.

## Preprocessing

The raw count data undergoes standard scRNA-seq preprocessing steps:

1.  **Filtering:** Cells with less than 200 genes and genes expressed in less than 3 cells are removed.
2.  **Normalization:** Total-count normalization to 1e4 counts per cell.
3.  **Log-transformation:** Log1p transformation of the normalized data.
4.  **Highly Variable Genes (HVG) Selection:** Selection of 1000 most highly variable genes.
5.  **PCA:** Principal Component Analysis for dimensionality reduction.

## Integration Methods

Two popular integration methods are applied to correct for batch effects:

-   **Harmony:** A method that iteratively adjusts cell embeddings to remove batch effects while preserving biological variation.
-   **scVI:** A deep generative model for scRNA-seq data that learns a low-dimensional latent representation.

## Clustering Methods

After integration, several clustering algorithms are applied to the data's low-dimensional embeddings:

-   **KMeans:** A centroid-based clustering algorithm.
-   **t-SNE:** A dimensionality reduction technique often used for visualization and sometimes for clustering in conjunction with other methods.
-   **Diffusion Maps:** A non-linear dimensionality reduction technique particularly suited for uncovering underlying manifold structures.
-   **DBSCAN:** A density-based spatial clustering of applications with noise algorithm.
-   **Hierarchical Clustering:** A method that builds a hierarchy of clusters.
-   **Louvain/Leiden:** Graph-based clustering algorithms widely used in scRNA-seq for community detection.

## Benchmarking

The `scib-metrics` `Benchmarker` is used to evaluate the performance of different integration and clustering methods. Metrics include:

-   **BioConservation Metrics:** Silhouette score (label), NMI, ARI for clustering quality, and isolated label scores.
-   **BatchCorrection Metrics:** Graph connectivity to assess the removal of batch effects.

## Visualization

UMAP (Uniform Manifold Approximation and Projection) is used for visualizing the high-dimensional data in 2D space. Embeddings from different integration methods and clustering results are plotted, colored by cell type, technical batch, and predicted clusters to visually assess their performance.

## Results

The benchmarking results are presented in a DataFrame, providing a quantitative comparison of each method across various biological conservation and batch correction metrics. UMAP visualizations allow for qualitative assessment of batch mixing and cell type separation.
