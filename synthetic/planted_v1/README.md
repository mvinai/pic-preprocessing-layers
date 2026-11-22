# Planted-spectrum recovery checks (SYNTHETIC ONLY, seed 6310)
Script: planted_spectra_check.py (run with .venv/bin/python). Results: planted_results.json.
Why: synthetic_v1 has latent factors but no known singular values, so it cannot test whether the metrics recover a planted rank or decay (revised_plan_v1 Section 4).
Cases at shapes 9505x82 and 9505x59: exact rank 10 + noise 0.01; power law alpha=1.5; exponential beta=0.15.
Findings:
- k95 (sigma^2 energy), permutation null and MP edge all recover planted rank 10.
- Decay fits recover alpha=1.50 and beta=0.15 exactly and pick the right family; on the rank-10 step spectrum they wrongly pick "power", so the decay model needs a step/rank-plus-noise alternative.
- Entropy effective rank and k95 on sigma energy are inflated by the noise floor (34 and 70 for true rank 10 at 82 cols); report alongside, not as the headline.
- MP edge over-counts on noise-free decaying spectra (5 vs ~2 by permutation); keep it as a reference line only, as planned.
- CUR error ratio vs SVD at k=10 is always >=1; CPQR (about 1.5-1.6) and DEIM (1.5-1.8) beat leverage (2.0-3.6) and uniform (2.2-3.6).
Not yet tested: bootstrap stability, randomized SVD, Belsley, U=A(I,J)^+ middle matrix, classifiers.
