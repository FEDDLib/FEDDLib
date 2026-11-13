# Artery drug transport (amlodipine study)

A structure-chemical interaction run on the `Artery_dan_SCI` mesh: the artery
is pressurised, and from t = 860 s a drug enters through both walls — inner
flags 5, 7, 8 and outer flags 4, 6, 9 — and diffuses through the tissue.

All seven material regions carry `ActiveBool`, `growthBool` and
`reorientationBool` at 0, and `M`, the parameter the study's reaction sweep
varied, is 0 in every region.

Build target `paper_amlodipine`. It needs Interface2 (the AceGen elements);
without it `main` returns `EXIT_FAILURE` after the MPI session, like every
other AceGen example here. It is a research run — 1500 s of simulated time,
11 time-step segments, the finest of them 12000 steps of 0.025 s — and is
therefore built but not registered as a ctest.

```
mpirun -n <ranks> ./paper_amlodipine
```

Every input file is a command-line option (`--simulationsParameters`,
`--materialParameters`, `--solverParameters`,
`--preconditionerParametersStructure`, `--preconditionerParametersChemistry`),
so a different parameter set needs no new executable and no new directory.

## This case is the unperturbed one

The study behind it swept four parameters around this configuration. Only the
centre point is kept here; each sweep is reproduced by copying the two XML
files, editing the one parameter, and passing them on the command line.

| Swept | Where | Here | The study's values |
|---|---|---|---|
| `D0` | `materialParameters_dan.xml` | `7.e-3 … 3.5e-3` | entry 1 ∈ {0.0007, 0.007, 0.07} × entry 3 ∈ {0.00035, 0.0035, 0.035}, nine cases |
| `M` | `materialParameters_dan.xml` | all `0.` (off) | entries 2 and 3 from `-3.e-1 / -1.5e-1` down to `-3.e-6 / -1.5e-6`, six cases |
| `Kappa` | `materialParameters_dan.xml` | `103.783 / 69.5349` | three sets, in directories named 15, 30, and 45; this is the 30 set |
| `Pressure Reduction Amount mmHg` | `simulationParameters.xml` | `0.0` | 0, 10, 20, 30 |

The drug itself is a fifth switch: `Inflow Start Time` is 860.0 here, and the
study's drug-free cases set it to 1.e7, i.e. past the end of the run.

The study also ran four further arteries — kim_guzman, narula, phinikaridou
and plasschaert — each with its own mesh and `materialParameters_*.xml`. Their
`main.cpp` differed from this one only in the default material file name, so
they need no code of their own; point `--materialParameters` at the other file
and copy the matching mesh.
