# Figure scripts

These are the MATLAB scripts used to generate the main figures of the ketamine
manuscript. They are provided for **reproducibility and archival** purposes: they
read aggregated, cross-animal source data through the lab's internal data paths
(`getDataRoot()`, `getCADatasets`, etc.) and are not intended to be run by
external users out of the box. For a runnable, self-contained walkthrough of the
key analyses on a single-animal sample dataset, see `../AnalysisTutorial/`.

| Figure | Script | Topic |
|--------|--------|-------|
| Fig 1 | `Setup_BehLFPImpact_v5.m`      | Behavioral & LFP impact of ketamine; experimental setup |
| Fig 2 | `DissociationBehHPC_v9.m`      | Dissociation between behavior and hippocampal coding (theta sequences by velocity) |
| Fig 3 | `WakingRestDissociation_v6.m`  | Waking vs. rest dissociation (population-vector cosine-similarity decoding) |
| Fig 4 | `WeakerFunctionalConnect_v7.m` | Weaker functional connectivity (spike mutual information, dimensionality) |
| Fig 5 | `SleepStateIntrusion_v12.m`    | Sleep-state intrusion (bispectrum, theta wave asymmetry) |
| Fig 6 | `LDSREMandRun_v5.m`            | Linear dynamical system model — REM and Run |

MATLAB version: developed and tested on R2017b.
