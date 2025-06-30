# glmnet-cox-1se

Compare the `lambda.1se` rule behavior for Cox regressions in glmnet 4.1-8 vs. 4.1-9.

## Reproducibility

First, restore the renv environment:

```r
renv::activate()
renv::restore()
```

Then, run `benchmark.R`.

Note that the script installs different versions of glmnet at runtime.
