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

## Draft Summary of Results for README:

The benchmarking results, as evaluated by `scib-metrics`, provide insights into the performance of different integration and clustering techniques. The `df_extended` DataFrame summarizes various biological conservation and batch correction metrics for unintegrated PCA, Harmony, scVI, t-SNE, and Diffusion Maps embeddings, coupled with different clustering methods (KMeans, DBSCAN, Hierarchical, Louvain/Leiden).

Key observations include:

*   **Integration Effectiveness:** Harmony and scVI generally show improved batch correction scores (e.g., graph connectivity) compared to the unintegrated PCA, indicating their success in mitigating batch effects.
*   **Biological Conservation:** Metrics like Silhouette score (label) and NMI/ARI for various clustering methods help assess how well biological cell type information is preserved or recovered after integration and clustering. Initial KMeans NMI and ARI scores are low across most embeddings, suggesting that these simple clustering methods on their own might not fully capture the complex biological structure or that the default parameters need tuning. The `X_diffmap` embedding shows a notably higher KMeans NMI and ARI, suggesting it might provide a more distinct separation for these clusters with KMeans.
*   **New Embeddings:** The `X_tsne` embedding shows similar performance to unintegrated PCA across the displayed metrics, while `X_diffmap` appears to offer different characteristics, particularly in clustering quality as seen with KMeans NMI/ARI.
*   **Isolated Labels:** The `Isolated labels` metric indicates the ability to detect distinct cell populations, with scores generally consistent across most embeddings but a slight decrease for `X_diffmap`.

Overall, the benchmark provides a quantitative framework to compare methods. Further detailed analysis of specific metrics and their trade-offs is crucial for selecting the most appropriate method for a given biological question.

## How to Interpret Benchmarking Scores

When reviewing the benchmarking results, it's important to understand what each metric signifies:

*   **Integration Effectiveness (Batch Correction Metrics):** These metrics quantify how well the integration method removes batch effects.
    *   **Graph Connectivity:** Measures the connectivity of the graph between cells of different batches. Higher values (closer to 1) indicate better mixing of batches, suggesting successful batch effect removal. A perfectly integrated dataset would have high graph connectivity.
    *   **LISI (Batch):** Local Inverse Simpson Index for batch. Higher values indicate that cells from different batches are well-mixed locally, implying effective batch correction. Values closer to the number of batches suggest ideal mixing.
    *   **kBET:** k-Nearest Neighbor Batch Effect Test. Low kBET scores (closer to 0) suggest that cells from different batches are well-mixed in local neighborhoods, indicating effective batch correction.
    *   **PCR comparison:** Principal Component Regression comparison. A higher value (closer to 1) indicates that the explained variance by batch decreases significantly after integration, suggesting successful batch effect removal.

*   **Biological Conservation (BioConservation Metrics):** These metrics evaluate how well the biological signal (e.g., cell type distinctions) is preserved or enhanced after integration.
    *   **Silhouette Score (Label):** Measures how similar a cell is to its own cell type cluster compared to other cell type clusters. Higher values (closer to 1) indicate better separation of true biological cell types. A score near 0 means overlapping clusters, and negative scores suggest misclassified cells.
    *   **NMI (Normalized Mutual Information) & ARI (Adjusted Rand Index) for Clustering:** These scores compare the predicted clusters (from methods like KMeans, Louvain, etc.) to the ground truth cell type labels. Higher values (closer to 1) indicate a better agreement between the predicted clusters and the true biological cell types, suggesting that the clustering method on that embedding effectively recovers biological structure.
    *   **cLISI (Cell Type):** Local Inverse Simpson Index for cell type. Higher values indicate that cells of the same cell type are close together locally, implying good preservation of biological distinction.

*   **Isolated Labels:** This metric quantifies the presence of isolated cell types or rare populations that might be negatively impacted by integration. A higher score is generally better, indicating that distinct cell populations remain distinguishable after integration.

**General Interpretation:**

*   Ideally, an effective integration method will show **high scores for batch correction metrics** (e.g., Graph Connectivity, LISI Batch) and **high scores for biological conservation metrics** (e.g., Silhouette, NMI/ARI for true labels). There is often a trade-off between these two categories.
*   When comparing different embeddings (`X_pca_unintegrated`, `X_harmony`, `X_scvi`, `X_tsne`, `X_diffmap`), look for methods that strike a good balance. For example, `X_harmony` and `X_scvi` are expected to outperform `X_pca_unintegrated` in batch correction, while ideally maintaining or improving biological conservation.
*   Clustering metrics like KMeans NMI and ARI on specific embeddings (e.g., `X_diffmap`) can reveal which embeddings are most amenable to particular clustering algorithms for discerning biological groups.
