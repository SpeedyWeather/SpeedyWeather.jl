# Sensitivity Examples

This folder contains examples for sensitivity analysis using Enzyme and Checkpointing, similar to those featured in [Moses et al. (2025): "DJ4Earth: Differentiable, and Performance-portable Earth System Modeling via Program Transformations"](https://essopenarchive.org/users/1000301/articles/1360053-dj4earth-differentiable-and-performance-portable-earth-system-modeling-via-program-transformations). Please consult the paper for more details on the methodology and the results. 

Be aware that first time compilation of the gradient / `autodiff` can take a while, subsequent executions will be very fast. The code is only tested on Julia 1.10 so far. 