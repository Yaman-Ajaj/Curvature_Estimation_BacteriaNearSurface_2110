# Trajectory Curvature Analysis Pipeline

A MATLAB toolkit for estimating the curvature of 2D bacterial (or other particle) swimming trajectories, validating those estimates against known ground truth using synthetic data, and comparing different curvature-estimation algorithms under realistic measurement noise.

## What this repo does

Given tracked (x, y, z, frame) trajectories, the pipeline:

1. **Isolates the near-surface portions** of each track (using a z-height window), since curvature is only computed on contiguous runs where the particle stays near a boundary.
2. **Estimates curvature** along each near-surface run using one of several interchangeable algorithms (circle fitting, tangent-angle differencing, spline smoothing + differentiation, local polynomial fitting, high-order finite differences).
3. **Validates** those estimators against **synthetic trajectories with known analytical curvature** (circles, ellipses, arcs, and tracks with spatially/temporally varying curvature), with and without added positional noise.
4. **Sweeps method parameters** (smoothing tolerance, window size, polynomial order, step size) to find the settings that minimize error against ground truth.
5. **Visualizes and compares** results across methods, tracks, and parameter choices.

## How it works (pipeline overview)

```
create_synthetic_trajectories.m
        │  generates clean + noisy tracks with known true curvature
        ▼
analyzeAllTracks.m
        │  loads a trajectory .mat file, splits into near-surface runs,
        │  runs the selected curvature_estimation_*.m method on each run
        ▼
curvatureResults_*.mat  (per-track, per-section curvature results)
        │
        ├── groundTruthComparison.m   → error vs. known true curvature, by shape type
        ├── sweep_*.m                 → parameter sweeps (tolerance, window, order, step)
        ├── compareCurvatureVariationAcrossResultFiles.m → compare methods/runs
        └── plotTrackByIndexComparative.m / PlotWithInset.m → per-track visualizations
```

### Core building blocks

- **`analyzeAllTracks.m`** — main driver. Prompts you to pick a curvature method, loads a trajectory file, finds near-surface runs per track (splitting on z-boundary crossings so curvature is never computed across a gap), and saves one results `.mat` file per run.
- **`create_synthetic_trajectories.m`** — generates five families of test trajectories (`Circular`, `GentleArc`, `Elliptical`, `SpatiallyVarying`, `FastVarying`), each with an analytically known curvature, and saves both a clean and a jittered (noisy) version for testing estimator robustness.
- **`curvature_estimation_*.m`** — the interchangeable curvature estimators:
  - `curvature_estimation_OsculatingCircle.m` — exact 3-point circumscribed-circle fit
  - `curvature_estimation_leastSquaresCircle.m` — over-determined (Kasa) least-squares circle fit over a window
  - `curvatureEstimationTangentAngle.m` — signed curvature from tangent-angle differencing
  - `curvature_estimation_centralDiff.m` — high-order central-difference derivatives plugged into the parametric curvature formula
  - `curvature_estimation_localPoly.m` — local polynomial fit, differentiated analytically
  - `curvature_estimation_splineDerivative.m` / `curvature_estimation_splineSmoothed*.m` — smoothing-spline based methods (smooth the path first, then apply one of the above formulas)
- **`smoothTrajectorySpline.m`** — shared helper: fits a tolerance-based smoothing spline (`spaps`) to noisy (x,y) data before differentiating.
- **`smoothCurvature.m`** — post-hoc smoothing of a computed curvature signal (`movmean`, `median`, or `sgolay`).
- **`sweep_*.m`** — grid-search helpers that run a chosen method across a parameter grid (spline tolerance, step size, window/order) and report/plot RMSE against ground truth.
- **`groundTruthComparison.m`** — loads one or more results files, computes absolute/relative error against the analytically known curvature per synthetic shape, trimming a percentage of points at each run's edges (`edgeTrimPercent.m`) to avoid edge-effect bias.
- **Plotting scripts** (`Plot_Trajectories.m`, `plotTrackByIndexComparative.m`, `PlotWithInset.m`) — visualize raw 3D tracks, and compare curvature-vs-frame curves across methods or against ground truth, with an inset of the trajectory itself.

## Typical workflow

1. Run `create_synthetic_trajectories.m` once to generate `Fake_synthetic_trajectories_Clean.mat` and `Fake_synthetic_trajectories_Noise.mat`.
2. Run `run_sweep_tests.m` to tune a method's parameters (tolerance multiplier, window size, polynomial order) against the noisy synthetic set, using the clean set's known curvature as ground truth.
3. Run `analyzeAllTracks.m` on your real (or synthetic) data with the tuned method/parameters to produce a `curvatureResults_*.mat` file.
4. Run `groundTruthComparison.m` (on synthetic data) or the plotting scripts (on real data) to inspect and validate the results.

## Requirements

- MATLAB with the **Curve Fitting Toolbox** (`spaps`, `fnval`, `fnder` are used by the spline-based methods).
- No other external dependencies.

## Repository contents

See the accompanying Word document (`Curvature_Analysis_Pipeline_Documentation.docx`) for a detailed, function-by-function explanation of every file, its inputs/outputs, and the underlying math.
