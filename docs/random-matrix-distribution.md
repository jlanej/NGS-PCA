# Random Matrix Distribution

## Background

The [Halko, Martinsson, and Tropp (2011)](https://doi.org/10.1137/090771806) randomized SVD algorithm begins by projecting the input matrix onto a low-dimensional random subspace. The quality of this initial projection influences the approximation, though subsequent power (subspace) iterations substantially reduce any dependence on the specific distribution used.

## Supported distributions

NGS-PCA supports two distributions for the initial random test matrix (selected via `-distribution`):

| Option | Description |
|---|---|
| `UNIFORM` | Entries drawn uniformly from [0, 1) using a Mersenne Twister PRNG. This is the default, chosen for backward compatibility. |
| `GAUSSIAN` | Entries drawn from a standard normal distribution N(0, 1) using the same PRNG. This matches the distribution specified in the original Halko et al. paper. |

## Effect on results

In practice, the choice of distribution has minimal effect on the resulting principal components when power iterations are used (the default). This was validated directly using the bundled 1000 Genomes chr1 example (n = 18 samples, 2 PCs):

- **PC1:** absolute Pearson *r* = 1.000, absolute Spearman *r* = 1.000
- **PC2:** absolute Pearson *r* = 0.9999, absolute Spearman *r* = 1.000

The full per-PC comparison (mean, median, SD, Pearson and Spearman correlations) is available in [`example/exampleOutput_1000G_chr1/distribution_comparison.txt`](../example/exampleOutput_1000G_chr1/distribution_comparison.txt) and is regenerated and checksum-verified on every CI run.

This behavior is consistent with the theoretical result that the power iteration scheme rapidly converges the estimated column space regardless of the initial random matrix, as long as the matrix has sufficient rank diversity in the test directions.

## Implementation details

Random values are generated using [Apache Commons Math `MersenneTwister`](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/random/MersenneTwister.html). The seed is controlled by the `-randomSeed` parameter (default: 42), ensuring deterministic and reproducible output across runs.

```java
MersenneTwister twister = new MersenneTwister(randomSeed);
// UNIFORM: twister.nextDouble()
// GAUSSIAN: twister.nextGaussian()
```

For strict reproducibility, fix both `-randomSeed` and `-distribution` in your analysis pipeline.
