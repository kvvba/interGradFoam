# interGradFoam: changes from mulesGradFoam (theory)

interGradFoam is mulesGradFoam with a different treatment of the phase change
mass source. The thermal model (interfacial heat flux from the
alpha1 = 0.5 iso-surface with Tsat as a Dirichlet condition, separate TL/TV
equations), surface tension and momentum are unchanged.

Notation: alpha1 liquid volume fraction, rho1 liquid and rho2 vapour density,
U one-fluid velocity, mDot volumetric mass source [kg/m^3/s] (> 0 for
evaporation), qDot interfacial heat flux per unit volume, hlv latent heat.


## 1. Sharp mass source instead of Hardt & Wondra smoothing

**mulesGradFoam.** The source qDot/hlv, which is non-zero only in the cells
cut by the interface, is
1. smoothed with a Helmholtz equation, `diff*lap(mDot) - mDot = -qDot/hlv`,
2. removed from the interface region (only alpha1 < 0.001 and
   alpha1 > 0.999 are kept),
3. rescaled so the vapour source and the liquid sink each integrate to the
   evaporation rate.

The bubble growth then depends on the numerical parameters `diff` and the
cut-off values, and the source does not sit where phase change takes place.

**interGradFoam.** As in PSI-BOIL (Sato & Niceno 2013; Lal, Sato & Niceno
2015):

    mDot = qDot/hlv     only in the cells cut by the interface, zero elsewhere

There is no smoothing and no shifting, and therefore no `diff` parameter.


## 2. Consistent volume and liquid-fraction equations

With the source located in the interface cells, the governing equations are

    div(U) = mDot (1/rho2 - 1/rho1)                        (B)  pressure eq.
    d(alpha1)/dt + div(alpha1 U) = -mDot/rho1              (A)  VOF eq.

(A) + (B) conserve mixture mass exactly. Two changes follow from this.

- **Pressure equation.** The volume source is `mDot (1/rho2 - 1/rho1)`.
  mulesGradFoam used `mDot/rho`, which is only correct when the source sits
  in pure vapour and the sink in pure liquid. Dividing by the mixture density
  in an interface cell gives the wrong expansion.
- **VOF equation.** The liquid sink `-mDot/rho1` now appears explicitly.
  mulesGradFoam had no source in the alpha equation. Phase change entered
  only through the velocity divergence, which Hardt & Wondra had moved out of
  the interface.

The interFoam-type alpha equation contains a term `alpha1 * div(U)`. In the
incompressible case it compensates the continuity error. Here, the
phase-change part of div(U) must not act on alpha1, because (A) has no such
term. It is therefore removed from the divergence:

    divU_alpha = div(phi) - mDot (1/rho2 - 1/rho1)

What remains is the continuity error only, which is the meaning the term has
in interFoam.

*Why the interPhaseChangeFoam pattern was rejected.* Proposed an
implicit coefficient `Sp = -mDot (1/rho2 - 1/rho1)` to cancel div(U). In
explicit MULES, Sp acts on the new alpha1 and div(U) on the old one, so the
two only cancel to O(dt). Because rho1/rho2 is about 1600, the error relative
to the evaporated liquid volume is about (rho1/rho2) x (change of alpha1 per
step). On the spherical test this created spurious liquid of the order of 60x
the evaporated volume per step, and the mass continuity error rose from
1e-10 to 4e-5.


## 3. Cap and redistribution of the source

mDot follows from heat flux times interface area. It does not know how much
liquid the cell holds, so the source alone can empty a cell (or overfill it,
for condensation) within one step. As in PSI-BOIL (`mdot_cut`):

    evaporation:   mDot <= rho1 * alpha1 / dt
    condensation:  mDot >= -rho1 * (1 - alpha1) / dt

The excess goes to the face neighbours: towards the liquid for evaporation,
towards the vapour for condensation. Each face is weighted with the square
of the cosine between the interface normal and the face normal. The
receivers are capped again, and any remainder is dropped and logged.

mDot is used by the alpha step of the *next* time step. With an adaptive
time step that step can be up to 1.2x longer, so the cap uses 1.2 x dt.


## 4. Face flux of the volume fraction out of source cells (displacedPhaseFlux)

**Problem.** All of the expansion `mDot (1/rho2 - 1/rho1)` happens in the
cut cell. An algebraic VOF scheme takes the face value of outgoing flow
mainly from the donor cell, so the source cell pushes its own vapour/liquid
mixture into the liquid cell next to it. That neighbour then loses liquid
while the source cell still holds some. The alpha1 = 0.5 iso-surface, and
with it the source, moves on early. The source cell is left with
alpha1 ~ 0.15. Nothing removes this residual: there is no source there any
more, the vapour inside the bubble does not move, and MULES interface
compression scales with |U| ~ 0.

**Model.** Fluid leaving a source cell through a face is the phase on the
other side of that face: liquid towards the liquid, vapour towards the
vapour. The face value is therefore the alpha1 of the *receiving* cell. The
vapour produced stays in the source cell and displaces the liquid. The cell
empties completely before the interface moves to the next cell. These are
the kinematics of a sharp interface, which geometric VOF gets automatically.

Boundedness is still enforced by the MULES limiter. The rule is related to
the donor-acceptor scheme of Hirt & Nichols (1981). It is exact for pure
expansion, and an approximation for flow along the interface through source
cells.


## 5. Required conditions on the volume fraction field

Two numerical effects can still leave liquid inside the bubble.

- **MULES local-extremum bound.** The limiter does not let a cell fall below
  the minimum alpha1 of its neighbours, except by the small explicit sink
  `mDot dt/rho1`. If vapour-side neighbours hold residual liquid, the source
  cell stalls at that level. `extremaCoeff 0.1` relaxes the bound towards the
  global [0, 1] range.
- **Initial interface.** Cells on the vapour side of alpha1 = 0.5 are never
  source cells, so any liquid they hold at t = 0 stays there. With Hardt &
  Wondra the vapour produced inside the bubble flushes it out. The initial
  alpha1 must therefore be sharp, e.g. exact volume fractions, not a profile
  smeared over several cells.

With a sharp initial interface the residual vanishes completely, and
`extremaCoeff` makes no difference. It is kept as a safeguard for situations
that create partially filled vapour-side cells. The residual liquid also
drives spurious vapour currents inside the bubble, through CSF acting on its
alpha gradients; removing it reduced them by about 7x.


## 6. Alpha sub-cycling

Allowed. The sink and divU act in each sub-cycle, so they act over the full
time step in total, as do the cap and the volume balance. mDot, and with it
the set of source cells, is computed once per time step and stays fixed over
the sub-cycles. On the spherical case, 3 sub-cycles gave the same growth,
balance and interface as 1. They cost 11 % more and allowed small
undershoots (min alpha1 ~ -2e-5), so 1 is recommended.


## 7. Results so far (7.5 um spherical bubble growth, Ja ~ 15, to 0.1-0.15 ms)

| | HW (diff 2.25e-10) | sharp (final variant) |
|---|---|---|
| R/R_analytical at 0.1 ms | -18 % | -17 % |
| mass continuity error | ~3e-6 | ~1e-10 |
| liquid volume imbalance / source | n/a | ~0.002 % |
| liquid left inside bubble | 0 | 0 |

Both under-predict mostly in the first 0.05 ms. At t0 the thermal boundary
layer (~3 um) is thinner than one cell, so this is likely a resolution
effect.

The interfacial heat flux, and therefore mDot, varies ~13 % from cell to
cell along the interface, and the total rate jumps up to ~4 % between steps.
The smoothing in mulesGradFoam hid this; the sharp source passes it on to
the pressure and velocity fields.


## 8. Open points

- Flow boiling (Maity): bubble foot / contact line, where the cap is expected
  to be active; tangential flow through source cells (approximation in 4.).
- Implicit MULES predictor (MULESCorr) not yet tested.
- Grid convergence (3.75 um) for the radius deficit.
- Density-scaled CSF (the vapour side currently gets
  ~rho1/rho2 larger capillary acceleration), smoother mDot (one-sided
  directional gradients), liquid-velocity extension.
