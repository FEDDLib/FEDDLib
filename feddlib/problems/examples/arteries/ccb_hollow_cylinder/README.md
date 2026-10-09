# ccb_hollow_cylinder

A hollow cylinder (Ri=1.0, Ro=1.3, L=2.0) with the structure-chemistry
interaction (SCI) problem and the non-CMM SMC element
(`SCI_SMC_Active_Growth_Reorientation`), solved with FEDDLib's
NOX/Belos/FROSch solver stack. `../ccb_hollow_cylinder_cmm` is the same case
with the constrained-mixture (CMM) element.

- Boundary conditions: bottom and top faces fixed axially, three
  circumferential pins for static determinacy, pressure on the inner wall
  ramped 0 -> 20 (plain pressure units) over t = [0, 10].
- All `*Bool` material flags are 0 (purely passive, load-driven response).
- `solverParameters.xml` prints NOX's outer iterations (`Outer Iteration`,
  `Outer Iteration StatusTest`) and Belos's details of every linear solve
  (`Linear Solver Details`, `Inner Iteration`).

## Build

From a configured Trilinos+FEDDLib(+Interface2/AceGenInterface) build
directory:

```
make ccb_hollow_cylinder
# or: make -j <N>   (builds everything, this target included)
```

## Run

```
cd <build_dir>/feddlib/problems/examples/arteries/ccb_hollow_cylinder
mpirun -np 4 ./ccb_hollow_cylinder.exe 2>&1 | tee run.log
```

(The binary may be called `problems_ccb_hollow_cylinder.exe` depending on
TriBITS target aliasing in this build -- check
`ctest -N | grep ccb_hollow_cylinder` if `ccb_hollow_cylinder.exe` isn't
found next to the copied XML/mesh files.)

## Mesh / BC flags

See `meshes/ccb_hollow_cylinder/hollow_cylinder_p1.mesh` (P1; FEDDLib builds
the P2 mesh) and the comment block at the top of `main.cpp` for the
face/vertex flags (bottom/top axial fix, 3 circumferential pins, inner-wall
pressure).

## Adaptive time stepping

With adaptive time stepping a time step that fails is repeated from the state
it started from with a smaller size, and the size grows again after a run of
converged time steps. A time step fails when an element cannot compute its
state (its Gauss-point iteration returns an error code), when the Newton
iteration does not converge in `MaxNonLinIts` or diverges, or when the solver
breaks down. The parameters, in `Timestepping Parameter`:

| Parameter | Default | Meaning |
|---|---|---|
| `Adaptive Time Stepping` | false | default of the segments' `Adaptive` |
| `Time Step Reduction Factor` | 0.5 | a failed time step is repeated with its size times this |
| `Time Step Increase Factor` | 2 | the size grows by this ... |
| `Converged Time Steps Before Increase` | 5 | ... after this many converged time steps in a row |
| `Newton Divergence Factor` | 1e3 | a Newton iteration whose residual exceeds this times its first one fails the time step (0: off) |
| `Accept Element Failures At Minimum dt` | false | default of the segments' parameter of that name |

and in each segment of `Timestepping Intervalls`:

| Parameter | Default | Meaning |
|---|---|---|
| `dt` | | the size the segment starts with (and, without `Adaptive`, takes throughout) |
| `Adaptive` | `Adaptive Time Stepping` | whether the segment adapts its size; a segment that does not, e.g. a pressure ramp, takes `dt` |
| `Maximum dt` | `dt` | the largest size the segment grows to |
| `Minimum dt` | `dt`/1000 | below it the run stops with an error naming the failure |
| `Accept Element Failures At Minimum dt` | false | a time step that elements fail at the smallest size is tried once more with their failures only counted, and accepted if the Newton iteration converges |

The log shows every repeated time step (`[adaptive] t = ...: the time step of
... failed (...); repeating it with ...`) and, at the end, how many there were.

`simulationParameters_adaptive*.xml` are the inputs of the ctest
`ccb_hollow_cylinder_adaptive`. Without failures the run is that of the fixed
time step size; a repeated time step gives the result of a run that took the
accepted sizes directly, to the last bit, unless the preconditioner keeps its
coarse basis (`Reuse: Coarse Basis`), which it then took from the failed
attempt's first matrix (differences at the Newton tolerance).

Load functions that compute the load at `t + Load Step Size` (as the artery
cases do) assume a fixed time step size: use them with adaptive time stepping
only in segments that do not adapt. This case's load function uses the time
of the end of the time step and works with any size.
