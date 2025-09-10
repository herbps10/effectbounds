# `effectbounds`: estimating non-overlap bounds for causal effects

The identification of causal effects typically relies on the _overlap assumption_ (also known as _positivity_), which requires that all units have a positive probability of being in either the treatment or control group. When overlap is structurally violated, with some units having zero probability of receiving the treatment (or control), then causal effects are un-identified. When overlap is practically violated, with some units having very small probability of receiving the treatment (or control), then traditional causal inference estimators may have poor finite-sample performance.

Non-overlap bounds are an approach for estimating causal effects even when non-overlap is violated, by focusing on estimating _bounds_ on the effect rather than its precise value. 

The `R` package `effectbounds` provides tools for estimating non-overlap bounds.
