# RISE-X

Robust Image-on-scalar Regression with Simultaneous Activation Region and Extreme Signal Detection


## Installation

```r
remotes::install_github("YutingDuan22/RISEX")
```

## Reproducible Workflow

```r
library(RISEX)
adj <- coordinates_to_adjacency_matrix.nf(voxels, connectivity = 6)
voxel_wise_fit <- voxel.wise.fit(img, cov_mat, voxels, cut_alpha = 0.95)
RISEX <- RISEX(img=img, voxels=voxels, cov_mat=cov_mat,
              voxel_wise_fit$psi_mat, voxel_wise_fit$lambda, adj = adj,
              beta0 = voxel_wise_fit$beta,
              sigma = voxel_wise_fit$sigma1,
              sigma0 = voxel_wise_fit$sigma0)
```

## Evaluation
```r
metrics(true_beta, RISEX$beta)
```

![RISE-X Example](Simu_RISEX.png)
