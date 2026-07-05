# Scenario-DeePC

**Scenario-based Data-Enabled Predictive Control: Robustification via the Scenario Approach**

MATLAB implementation accompanying the paper by Sebastian Zieglmeier, Nikolas Recke, and Mathias Hudoba de Badyn (Department of Technology Systems, University of Oslo).

Scenario-DeePC integrates the scenario approach into Data-Enabled Predictive Control (DeePC) to provide a distribution-free probabilistic guarantee on output-constraint satisfaction. The uncertainty description is built directly from the controller's own closed-loop prediction errors, keeping the method fully data-driven and free of distributional assumptions. This repository contains the code to reproduce the numerical examples in the paper.

## Requirements

- MATLAB (tested with R2024b)
- A convex optimization solver / modeling toolbox (preferably MOSEK with a license)

## Repository structure

| Path | Description |
|------|-------------|
| `simulate_Scenario_DeePC_Linear_Boeing_Online.m` | Boeing 747 example: nominal, robust, and adaptive Scenario-DeePC (online buffer). |
| `simulate_Scenario_DeePC_Linear_Boeing_real.m` | Boeing 747 example: closed-loop run used for the reported results. |
| `simulate_Scenario_DeePC_LPV.m` | Nonlinear two-tank (LPV) example. |
| `Controller/DeePC_fast.m` | Baseline DeePC controller. |
| `Controller/Scenario_DeePC_multi_Sc_cost.m` | Scenario-DeePC controller (scenario-averaged cost). |
| `System/Linear_Boeing_747.m` | Boeing 747 longitudinal model. |
| `System/LPV_2_Tank.m` | Two-tank LPV model. |
| `Functions/make_hankel_MIMO.m` | Hankel-matrix construction for MIMO data. |
| `Functions/discretize_LPV.m` | LPV discretization. |
| `Functions/get_ref.m`, `Functions/get_ref2.m` | Reference-trajectory generators. |
| `Functions/Smooth_Step.m` | Smoothed step reference. |
| `Functions/system_boundaries.m`, `Functions/system_boundaries_multi.m` | Input/output constraint definitions. |

## Running the examples

Open MATLAB in the repository root and run one of the top-level scripts:

```matlab
% Linear Boeing 747 example
simulate_Scenario_DeePC_Linear_Boeing_Online

% Nonlinear two-tank (LPV) example
simulate_Scenario_DeePC_LPV
```

Each script sets up the system, collects data, builds the scenario buffer from prediction errors, and runs the closed-loop comparison between standard DeePC and Scenario-DeePC.

## Citation

If you use this code, please cite the paper:

```bibtex
@article{Zieglmeier_ScenarioDeePC,
  title   = {Scenario-based Data-Enabled Predictive Control: Robustification via the Scenario Approach},
  author  = {Zieglmeier, Sebastian and Recke, Nikolas and Hudoba de Badyn, Mathias},
  year    = {2026},
  note    = {Preprint. Submitted for possible publication.}
}
```

## License

See the repository license file. The accompanying manuscript is a preprint submitted for possible publication; copyright may be transferred without notice, after which the published version may differ from this one.
