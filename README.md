# Tail-PPML

Stata package `ppmlmn` for PPML gravity estimation with two-way (exporter and importer) fixed effects and heavy-tail-robust inference. `ppmlmn` reports a country-level $m$-out-of-$n$ bootstrap confidence interval that remains valid when bilateral shocks are Pareto-tailed with index $\alpha \in (1, 2)$, a regime in which the conventional HC0 sandwich interval under-covers.

## Install

```stata
net install ppmlmn, from("https://raw.githubusercontent.com/yulongwang06/Tail-PPML/main/") replace
```

Then

```stata
help ppmlmn
```

No other Stata packages are required. Full installation notes, including a manual-install path, are in [`INSTALL.md`](INSTALL.md).

## Minimal example

```stata
ppmlmn trade ln_dist contig comlang colony, exporter(iso_o) importer(iso_d)
```

This estimates

$$
\mathbb{E}[y_{ij} \mid x] \;=\; \exp\!\bigl(\phi^{\text{ex}}_i + \phi^{\text{im}}_j + x_{ij}'\beta\bigr),
$$

profiles the $2N$ country fixed effects out by Sinkhorn (matrix balancing), updates $\beta$ by Newton, and reports two 95% CIs for each slope: the HC0 (or HC1) sandwich, and the country-level $m$-out-of-$n$ bootstrap. Frequently used options:

- `reps(#)` number of bootstrap draws $M$ (default `999`)
- `g(#)` countries per subsample $G$ (default $\lfloor N^{0.7} \rfloor$)
- `seed(#)` random-number seed
- `studentizer(hc0 | hc1)` studentizer for the self-normalized statistic (default `hc1`)
- `level(#)` confidence level (default `95`)

Countries that are perfectly separated (a zero exporter or importer marginal) are dropped, as in standard PPML.

## Why the bootstrap

When $y_{ij}$ is heavy-tailed with index $\alpha \in (1,2)$, the PPML score has infinite variance. The point estimate $\hat\beta$ converges to a multivariate $\alpha$-stable distribution at rate $n^{1 - 1/\alpha}$ rather than to a Gaussian at rate $n^{-1/2}$, and the conventional sandwich variance estimator diverges in probability at rate $n^{2/\alpha - 1}$. Cluster-robust intervals do not restore coverage; they remain anchored to the wrong limiting law.

The country-level $m$-out-of-$n$ bootstrap draws $G$ countries without replacement, keeps the induced subnetwork of dyads, refits the two-way FE PPML, and self-normalizes each subsample slope by its own sandwich SE. Confidence limits invert the empirical quantiles of these self-normalized statistics. The procedure requires no estimation of $\alpha$ or of the stable scale and is asymptotically valid in both the Gaussian ($\alpha \ge 2$) and stable ($\alpha \in (1,2)$) regimes.

## Repository contents

```
ppmlmn.ado, ppmlmn.sthlp          the command and its help file
ppmlmn.pkg, stata.toc              net-install manifest
INSTALL.md                          installation notes
```

Python and MATLAB/Octave implementations of the same procedure, together with the replication scripts for the simulations and empirical applications in the paper, will be added under a `replication/` subfolder.

## Citation

Egger, P., T. Ji, and Y. Wang (2026), "PPML and Heavy-Tailed Trade Flows: Why Standard Inference Fails and How to Fix It."

## Author

Yulong Wang, Lehigh University, https://sites.google.com/site/yulongwanghome/. Please file questions and bug reports on the [issue tracker](https://github.com/yulongwang06/Tail-PPML/issues).
