# affine-wall-straightening

Exact finite straightening of rational Hurwitz–Lerch multiple polylogarithms into cyclotomic multiple polylogarithms.

This repository contains the supplementary script for the paper

> Ramanuja, *Finite affine-wall straightening for rational Hurwitz–Lerch multiple polylogarithms*.

## What the script does

`example_7_2_reduction.py` implements the straightening algorithm of Definition 4.2 of the paper in exact arithmetic. Coefficients live in the group ring Q[μ_M], with M = qN. It then carries out the complete reduction of Example 7.2:

```
X = Σ_{n1 > n2 > n3 > 0} 1 / ( n1² (n2 − 1/2) (n3 + 1) )
```

that is, s = (2, 1, 1), t = (0, 1/2, −1), q = 2, with unit colours. The output is X written as an explicit rational linear combination of alternating multiple zeta values Li_h(±1, …).

The strict convention used throughout is

```
Li_{h1,…,hd}(y1,…,yd) = Σ_{n1 > … > nd > 0} y1^n1 ⋯ yd^nd / (n1^h1 ⋯ nd^hd).
```

## Requirements and usage

Python 3.8 or later. The script uses only the standard library.

```
python3 example_7_2_reduction.py
```

The output lists one line per symbol, giving its coefficient and then the symbol. For example:

```
     2  Li_{2,1,1}(1,1,1)
     8  Li_{1,1}(-1,1)
     8  Li_{1}(-1)
```

No symbol of the form Li_{1,…}(1, …) appears. This matches the paper's regularization-free boundary theorem (Theorem 6.6).

## Using the algorithm on other sums

`straighten_H(R, s, t)` straightens

```
H_s(ξ; t) = Σ_{n1 > … > nr > 0} Π_j ξ_j^{n_j} / (n_j − t_j)^{s_j}
```

for rational shifts t_j that are not positive integers. Here `R = CycRing(a, N, q)` encodes the colours ξ_j = exp(2πi a_j / N), and q is a common denominator of the t_j.

The result is a dictionary that maps symbols to coefficients:

- **Symbols** are keys `(h, args)`. Each entry `((lo, hi), l)` of `args` stands for the argument (α_lo ⋯ α_hi) · ω^l, where α_j = exp(2πi a_j / (qN)) and ω = exp(2πi/q). Convert it to an exponent with `R.arg_exp((lo, hi), l)`.
- **Coefficients** are elements of Q[μ_{qN}], stored as `{exponent: Fraction}`.
