"""
Four-patch SIR model with gravity-model spatial coupling.
LASER-framework style — runs with numpy only (no extra installs needed).

Patch layout (linear, 75 km spacing):
    Patch 0 ──75km── Patch 1 ──75km── Patch 2 ──75km── Patch 3

Gravity parameters: k=0.01, a=1, b=1, c=1.5
Row-normalisation cap: 15 % export fraction per patch per step.
"""

import numpy as np

# ── Patch geometry ────────────────────────────────────────────────────────────
NUM_PATCHES = 4
SPACING_KM  = 75.0
positions   = np.arange(NUM_PATCHES) * SPACING_KM   # [0, 75, 150, 225] km

# Pairwise Euclidean distances (km); diagonal is 0
dist = np.abs(positions[:, None] - positions[None, :])   # shape (4, 4)

# ── Population (initial) ──────────────────────────────────────────────────────
N = np.array([10_000, 8_000, 12_000, 6_000], dtype=float)

# ── Gravity-model spatial coupling ────────────────────────────────────────────
# Gravity flow:  G[i,j] = k * N[i]^a * N[j]^b / d[i,j]^c   (i ≠ j)
k_grav = 0.01
a_exp  = 1.0
b_exp  = 1.0
c_exp  = 1.5
MAX_EXPORT_FRAC = 0.15   # row-normalisation cap

def build_gravity_coupling(populations, distances, k, a, b, c, max_export_frac):
    """
    Return a (n_patches × n_patches) matrix `coupling` where
    coupling[i, j] is the *fraction* of patch i's population
    that moves to patch j each time step.

    Steps
    -----
    1. Compute raw gravity flows G[i,j] = k * N[i]^a * N[j]^b / d[i,j]^c
    2. Convert to movement fractions relative to source population.
    3. Row-normalise: if the total export fraction for patch i exceeds
       `max_export_frac`, scale every off-diagonal entry in that row down
       uniformly so the row sum equals exactly `max_export_frac`.
    """
    n = len(populations)

    # Raw gravity matrix — avoid division by zero on diagonal (d=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        G = np.where(
            distances > 0,
            k * (populations[:, None] ** a) * (populations[None, :] ** b)
              / (distances ** c),
            0.0,
        )

    # Fraction of source population N[i] that moves to patch j
    frac = G / populations[:, None]      # frac[i,j] = G[i,j] / N[i]

    # Row-normalise so total export ≤ max_export_frac
    row_export = frac.sum(axis=1, keepdims=True)   # shape (n, 1)
    scale = np.where(row_export > max_export_frac,
                     max_export_frac / row_export,
                     1.0)
    frac = frac * scale                            # coupling matrix

    return frac


coupling = build_gravity_coupling(N, dist, k_grav, a_exp, b_exp, c_exp, MAX_EXPORT_FRAC)

# ── Report coupling matrix ────────────────────────────────────────────────────
print("=" * 60)
print("Gravity coupling matrix  (fraction of source pop → dest)")
print("=" * 60)
header = "        " + "".join(f"  Patch {j}" for j in range(NUM_PATCHES))
print(header)
for i in range(NUM_PATCHES):
    row_str = "".join(f"  {coupling[i, j]:.5f}" for j in range(NUM_PATCHES))
    print(f"Patch {i}{row_str}")
print()
print("Row sums (total export fraction per patch):")
for i in range(NUM_PATCHES):
    print(f"  Patch {i}: {coupling[i].sum():.5f}  "
          f"({'capped at 15%' if coupling[i].sum() >= MAX_EXPORT_FRAC - 1e-9 else 'below cap'})")
print()

# ── SIR disease parameters ────────────────────────────────────────────────────
beta  = 0.30   # transmission rate  (day⁻¹)
gamma = 0.10   # recovery rate      (day⁻¹)
dt    = 1.0    # time step          (days)
T     = 180    # simulation length  (days)
steps = int(T / dt)

# ── Initial conditions — seed patch 0 ────────────────────────────────────────
I = np.array([10.0,  0.0,  0.0,  0.0])
R = np.zeros(NUM_PATCHES)
S = N - I - R

# ── Storage ───────────────────────────────────────────────────────────────────
S_hist = np.empty((steps + 1, NUM_PATCHES))
I_hist = np.empty((steps + 1, NUM_PATCHES))
R_hist = np.empty((steps + 1, NUM_PATCHES))
S_hist[0] = S
I_hist[0] = I
R_hist[0] = R

# ── Main simulation loop ──────────────────────────────────────────────────────
for step in range(steps):

    # 1. Within-patch SIR transitions (Euler step)
    new_inf = beta * S * I / N * dt     # force of infection × S
    new_rec = gamma * I * dt

    S_new = S - new_inf
    I_new = I + new_inf - new_rec
    R_new = R + new_rec

    # 2. Gravity-model spatial movement
    #    For each compartment, outflow[i,j] = coupling[i,j] * X[i]
    #    Net change for patch i: sum_j inflow[j→i] − sum_j outflow[i→j]
    def move(X):
        outflow = coupling * X[:, None]          # (n, n): outflow[i,j] from i
        net     = outflow.sum(axis=0) - outflow.sum(axis=1)   # inflow − outflow
        return X + net

    S_new = move(S_new)
    I_new = move(I_new)
    R_new = move(R_new)

    # Guard against floating-point negatives
    S = np.maximum(S_new, 0.0)
    I = np.maximum(I_new, 0.0)
    R = np.maximum(R_new, 0.0)

    S_hist[step + 1] = S
    I_hist[step + 1] = I
    R_hist[step + 1] = R

# ── Summary statistics ────────────────────────────────────────────────────────
print("=" * 60)
print("Simulation summary (gravity-coupled 4-patch SIR)")
print("=" * 60)
print(f"{'Patch':>7}  {'Pop (N)':>9}  {'Peak I':>9}  "
      f"{'Peak day':>9}  {'Attack rate':>12}")
for p in range(NUM_PATCHES):
    peak_val = I_hist[:, p].max()
    peak_day = int(I_hist[:, p].argmax() * dt)
    attack   = R_hist[-1, p] / N[p]
    print(f"{p:>7}  {N[p]:>9.0f}  {peak_val:>9.1f}  "
          f"{peak_day:>9}  {attack:>11.1%}")

print()
print("Final compartment totals across all patches:")
print(f"  S = {S_hist[-1].sum():.1f}")
print(f"  I = {I_hist[-1].sum():.1f}")
print(f"  R = {R_hist[-1].sum():.1f}")
print(f"  N = {(S_hist[-1] + I_hist[-1] + R_hist[-1]).sum():.1f}  "
      f"(initial N = {N.sum():.1f})")