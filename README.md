# affine-wall-straightening

Exact finite straightening of rational Hurwitz–Lerch multiple polylogarithms into cyclotomic multiple polylogarithms.

This repository accompanies the paper

> Ramanuja, *Finite affine-wall straightening for rational Hurwitz–Lerch multiple polylogarithms*.

## Contents

| Path | What it is |
|---|---|
| `paper/affine_wall_straightening.pdf` | The paper. |
| `paper/affine_wall_straightening.tex` | LaTeX source of the paper (`amsart`). |
| `example_7_2_reduction.py` | The supplementary script cited in the paper. It carries out the complete reduction of Example 7.2 in exact arithmetic. |
| `verification/` | Scripts that check the paper's statements by computation (see below). They are not cited in the paper, whose assertions rest on their written proofs. |

The strict convention used throughout is

```
Li_{h1,…,hd}(y1,…,yd) = Σ_{n1 > … > nd > 0} y1^n1 ⋯ yd^nd / (n1^h1 ⋯ nd^hd).
```

## The supplementary script

`example_7_2_reduction.py` implements the straightening algorithm of Definition 4.2. Coefficients are exact elements of the group ring Q[μ_M], with M = qN. The script reduces

```
X = Σ_{n1 > n2 > n3 > 0} 1 / ( n1² (n2 − 1/2) (n3 + 1) ),
```

that is, s = (2, 1, 1), t = (0, 1/2, −1), q = 2, with unit colours. The output is X as an explicit rational linear combination of alternating multiple zeta values Li_h(±1, …), one line per symbol, coefficient first:

```
python3 example_7_2_reduction.py
```

No symbol of the form Li_{1,…}(1, …) appears, as the paper's regularization-free boundary theorem (Theorem 6.6) predicts. The script needs only the Python standard library (Python 3.8 or later).

### Using the algorithm on other sums

`straighten_H(R, s, t)` straightens

```
H_s(ξ; t) = Σ_{n1 > … > nr > 0} Π_j ξ_j^{n_j} / (n_j − t_j)^{s_j}
```

for rational shifts t_j that are not positive integers. Here `R = CycRing(a, N, q)` encodes the colours ξ_j = exp(2πi a_j / N), and q is a common denominator of the t_j.

The result is a dictionary that maps symbols to coefficients:

- **Symbols** are keys `(h, args)`. Each entry `((lo, hi), l)` of `args` stands for the argument (α_lo ⋯ α_hi) · ω^l, where α_j = exp(2πi a_j / (qN)) and ω = exp(2πi/q). Convert it to an exponent with `R.arg_exp((lo, hi), l)`.
- **Coefficients** are elements of Q[μ_{qN}], stored as `{exponent: Fraction}`.

## Verification scripts

These need NumPy. Run them from anywhere, for example `python3 verification/check_boundary_cancellation.py`.

| Script | Checks | Expected outcome |
|---|---|---|
| `straighten.py` | The algorithm, with a complex floating-point ring (`NumRing`, for 0 < \|z_j\| < 1) as well as the exact ring. | Library module, used by the scripts below. |
| `formulas.py` | Independent implementations of the closed formulas of Theorems 5.5 (depth r−1) and 5.6 (top weight). | Library module. |
| `check_algorithm_vs_direct_sum.py` | Theorem 5.1 on 60 random cases with \|z_j\| < 1, against direct summation; also (7.1). | Worst relative error about 3·10⁻¹⁰; depth ≤ r and weight ≤ w in every case. |
| `check_depth_rm1_and_top_weight.py` | Theorems 5.5 and 5.6 on 80 random cases, in exact arithmetic. | 0 mismatches as formal elements. |
| `check_boundary_cancellation.py` | Theorem 6.6 on 300 random cases at roots of unity. | Dangerous coefficients vanish in all 192 cases with (s₁, ξ₁) ≠ (1,1), and survive in all 108 cases with (s₁, ξ₁) = (1,1). |
| `check_closed_forms.py` | The constants (7.3)–(7.6), using Richardson extrapolation. | Agreement to about 10⁻¹⁵. |
| `check_depth3_and_staircase.py` | Identity (7.2) at finite truncation, and the staircase values S₃, S₄ of Theorem A.1. | (7.2) holds exactly; the S₃, S₄ errors decay like log²M / M. |
| `example72_reduction.py` | Same reduction as the supplementary script, via the full library. | The 17 symbols of Example 7.2. |
| `check_example72_numerics.py` | The complete reduction of Example 7.2, numerically, with tail extrapolation. | X ≈ 0.7696049264; the two sides agree to about 10⁻¹⁰. Takes a minute or so. |
| `output_sizes.py` | Output sizes against the bounds of Theorem 3.5 and Proposition 5.8, and the top-weight part of the staircase. | Sizes far below the bounds; the staircase top-weight part is ζ*(2,1,1,1) term by term. |
