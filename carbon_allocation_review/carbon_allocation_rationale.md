# Step-by-step explanation of `carbon_allocation_offline.f90`

## 1. Purpose of this module

This module implements an offline carbon allocation solver for one average woody individual over one allocation period.

The allocation period is intentionally generic. It can represent one year, one month, one season, or any other period chosen by the host model. The solver itself does not know whether the carbon increment came from daily, monthly, or annual accumulation.

The central question solved by the module is:

> Given the current carbon pools of a woody individual and the carbon available for growth during an allocation period, how should this carbon be distributed among leaves, fine roots, and sapwood while satisfying allometric constraints?

The module does **not** compute photosynthesis, maintenance respiration, phenology, mortality, water balance, or daily carbon fluxes. It only receives an already-computed amount of available carbon and allocates it.

---

## 2. Main biological pools and variables

The woody individual is represented by four carbon pools:

\[
L = \text{leaf carbon mass}
\]

\[
R = \text{fine-root carbon mass}
\]

\[
S = \text{sapwood carbon mass}
\]

\[
H = \text{heartwood carbon mass}
\]

During normal allocation, carbon is allocated to the three living pools:

\[
\Delta B = \Delta L + \Delta R + \Delta S
\]

where:

\[
\Delta B = \text{carbon available for allocation}
\]

\[
\Delta L = \text{leaf carbon increment}
\]

\[
\Delta R = \text{fine-root carbon increment}
\]

\[
\Delta S = \text{sapwood carbon increment}
\]

Heartwood is not directly allocated new carbon during normal allocation. It can increase when sapwood is converted to heartwood, especially in abnormal allocation or in a separate turnover process.

---

## 3. Main allometric relationships

The solver is based on a set of allometric constraints. These constraints determine whether a proposed allocation is structurally consistent.

### 3.1 Leaf-to-root functional balance

The model assumes a functional balance between leaf mass and fine-root mass:

\[
L = lm2rm \times R
\]

where:

\[
lm2rm = \text{leaf-to-fine-root mass ratio}
\]

After allocation, the final pools should satisfy:

\[
L_{new} = lm2rm \times R_{new}
\]

If the candidate leaf increment is:

\[
x = \Delta L
\]

then:

\[
L_{new} = L + x
\]

Therefore:

\[
R_{new} = \frac{L+x}{lm2rm}
\]

and the fine-root increment is:

\[
\Delta R = R_{new} - R
\]

\[
\Delta R = \frac{L+x}{lm2rm} - R
\]

This is why the code computes:

```fortran
root_new   = leaf_new / params%leaf_to_root_ratio
delta_root = root_new - state%root_mass
```

This relationship is important because once `delta_leaf` is chosen, `delta_root` is no longer independent.

---

### 3.2 Carbon conservation

The available carbon over the allocation period must be distributed among living tissues:

\[
\Delta B = \Delta L + \Delta R + \Delta S
\]

Solving for sapwood increment:

\[
\Delta S = \Delta B - \Delta L - \Delta R
\]

Substituting:

\[
\Delta R = \frac{L+x}{lm2rm} - R
\]

we get:

\[
\Delta S = \Delta B - x - \left(\frac{L+x}{lm2rm} - R\right)
\]

\[
\Delta S = \Delta B - x - \frac{L+x}{lm2rm} + R
\]

The final sapwood mass is:

\[
S_{new} = S + \Delta S
\]

Therefore:

\[
S_{new} = S + \Delta B - x - \frac{L+x}{lm2rm} + R
\]

This is the expression that often looks confusing in the original LPJ-style code. It is simply the final sapwood mass written as a function of the candidate leaf increment.

---

### 3.3 Pipe model: leaf area and sapwood area

The pipe model links leaf area to sapwood cross-sectional area:

\[
LA = latosa \times SA
\]

where:

\[
LA = \text{leaf area}
\]

\[
SA = \text{sapwood cross-sectional area}
\]

Leaf area is computed from leaf mass and specific leaf area:

\[
LA = L \times SLA
\]

Therefore:

\[
L \times SLA = latosa \times SA
\]

and:

\[
SA = \frac{L \times SLA}{latosa}
\]

After allocation:

\[
SA_{new} = \frac{L_{new} \times SLA}{latosa}
\]

This relationship prevents the model from producing too much leaf area without enough conducting sapwood.

---

### 3.4 Sapwood mass, sapwood area, and height

Sapwood mass is related to sapwood volume and wood density:

\[
S = wooddens \times V_{sapwood}
\]

Sapwood volume is approximated as:

\[
V_{sapwood} = h \times SA
\]

Therefore:

\[
S = wooddens \times h \times SA
\]

Solving for height:

\[
h = \frac{S}{wooddens \times SA}
\]

After allocation:

\[
h_{pipe} = \frac{S_{new}}{wooddens \times SA_{new}}
\]

Using the pipe model:

\[
SA_{new} = \frac{L_{new} \times SLA}{latosa}
\]

so:

\[
h_{pipe} = \frac{S_{new}}{wooddens \times \frac{L_{new} \times SLA}{latosa}}
\]

This is the height implied by the combination of final leaf mass, final sapwood mass, SLA, wood density, and the pipe model.

---

### 3.5 Height-diameter allometry and total stem geometry

The model also uses a height-diameter relationship:

\[
h = allom2 \times D^{allom3}
\]

where:

\[
D = \text{stem diameter}
\]

The total stem volume is approximated as a cylinder:

\[
V_{stem} = h \times \frac{\pi D^2}{4}
\]

Total stem carbon includes sapwood and heartwood:

\[
C_{stem} = S + H
\]

and:

\[
C_{stem} = wooddens \times V_{stem}
\]

or:

\[
V_{stem} = \frac{S+H}{wooddens}
\]

The goal is to eliminate diameter and express stem geometry in terms of height.

From:

\[
h = allom2 \times D^{allom3}
\]

we solve for diameter:

\[
D = \left(\frac{h}{allom2}\right)^{1/allom3}
\]

Substitute this into stem volume:

\[
V_{stem} = \frac{\pi}{4} \times h \times \left(\frac{h}{allom2}\right)^{2/allom3}
\]

Rearranging gives:

\[
h^{1 + 2/allom3} = \frac{V_{stem} \times allom2^{2/allom3}}{\pi/4}
\]

Since:

\[
V_{stem} = \frac{S_{new}+H}{wooddens}
\]

we get:

\[
h_{stem}^{1 + 2/allom3}
=
allom2^{2/allom3}
\times
\frac{(S_{new}+H)/wooddens}{\pi/4}
\]

This is what the code calls:

```fortran
height_power_from_stem_geometry
```

It is not exactly height. It is height raised to the power:

\[
1 + \frac{2}{allom3}
\]

---

## 4. Why the normal allocation problem is solved with bisection

Normal allocation has three unknown increments:

\[
\Delta L, \Delta R, \Delta S
\]

But the code reduces the problem to one unknown:

\[
x = \Delta L
\]

Once \(x\) is chosen:

\[
\Delta R = \frac{L+x}{lm2rm} - R
\]

and:

\[
\Delta S = \Delta B - x - \Delta R
\]

Therefore the only remaining question is:

> Which value of \(x = \Delta L\) makes the final plant structure consistent with both the pipe model and total stem geometry?

The code defines a residual function:

\[
f(x) = height\_power\_from\_stem\_geometry - height\_power\_from\_pipe\_model
\]

The correct allocation satisfies:

\[
f(x) = 0
\]

When:

\[
f(x) = 0
\]

then:

\[
height\_power\_from\_stem\_geometry = height\_power\_from\_pipe\_model
\]

Because plant height is positive, matching height raised to this power is equivalent to matching height itself.

The equation is nonlinear and difficult to solve analytically, so the code uses the bisection method.

---

## 5. Structure of the module

The module is organized into three main parts:

1. Type definitions and constants.
2. Functions that compute derived quantities and residuals.
3. The main allocation routine and helper routines.

---

## 6. Module header and visibility

The module starts with:

```fortran
module carbon_allocation_offline_kernel
```

This defines a Fortran module. A module groups related types, constants, functions, and subroutines.

The line:

```fortran
use, intrinsic :: iso_fortran_env, only: real64
```

imports the `real64` kind from the intrinsic Fortran environment. This allows all real variables to be declared consistently as 64-bit real numbers:

```fortran
real(real64)
```

The line:

```fortran
implicit none
```

prevents Fortran from creating undeclared variables automatically. This is essential for avoiding hidden bugs caused by misspelled variable names.

The line:

```fortran
private
```

means that everything in the module is private by default. Only the items explicitly declared as `public` can be accessed from outside the module.

The public interface includes:

```fortran
public :: AllocationParameters
public :: PlantCarbonState
public :: AllocationResult
public :: allocate_woody_individual
public :: allocation_residual
public :: leaf_requirement_for_existing_sapwood
public :: height_from_total_stem_carbon
public :: height_from_pipe_model_state
```

This means that external code can create parameter objects, plant-state objects, allocation-result objects, and can call the main allocation and diagnostic functions.

---

## 7. Numerical constants and solver settings

The module defines:

```fortran
real(real64), parameter :: pi = 3.1415926535897932384626433832795_real64
```

This is used in the stem-volume equation:

\[
V_{stem} = h \times \frac{\pi D^2}{4}
\]

The bisection solver uses:

```fortran
real(real64), parameter :: default_x_tolerance = 1.0e-8_real64
real(real64), parameter :: default_f_tolerance = 1.0e-10_real64
integer, parameter :: default_max_iterations = 200
integer, parameter :: default_scan_segments  = 100
```

These are numerical settings, not biological parameters.

`default_x_tolerance` controls how small the interval in \(\Delta L\) must become before bisection accepts a solution.

`default_f_tolerance` controls how close the residual must be to zero.

`default_max_iterations` prevents the bisection loop from running forever.

`default_scan_segments` controls how many subintervals are used when scanning for a sign-changing bracket.

---

## 8. `AllocationParameters`

The type:

```fortran
type :: AllocationParameters
```

defines all parameters needed by the allocation solver.

It includes:

- `sla`: specific leaf area.
- `latosa`: leaf-area-to-sapwood-area ratio.
- `wood_density`: wood density.
- `leaf_to_root_ratio`: final leaf-to-root mass ratio for the allocation period.
- `allom2`: coefficient of the height-diameter allometry.
- `allom3`: exponent of the height-diameter allometry.

This type does not assign values. It only defines the structure. Values are assigned by the calling program, for example:

```fortran
params%sla = 12.0_real64
params%latosa = 8000.0_real64
```

---

## 9. `PlantCarbonState`

The type:

```fortran
type :: PlantCarbonState
```

stores the plant state before the current allocation event.

It includes:

- `leaf_mass`
- `root_mass`
- `sapwood_mass`
- `heartwood_mass`
- `height`

This is the old/current structural state used as input.

If allocation is annual, this is the state before the annual allocation event.

If allocation is monthly, this is the state before the monthly allocation event.

If the host model does not store height, the module provides functions to compute height from existing carbon pools before allocation.

---

## 10. `AllocationResult`

The type:

```fortran
type :: AllocationResult
```

stores the output of the allocation routine.

It includes flags indicating which allocation pathway was used:

```fortran
logical :: normal_allocation
logical :: converged
integer :: iterations
```

It also stores the carbon increments:

```fortran
real(real64) :: delta_leaf
real(real64) :: delta_root
real(real64) :: delta_sapwood
```

and final carbon pools:

```fortran
real(real64) :: leaf_mass_new
real(real64) :: root_mass_new
real(real64) :: sapwood_mass_new
real(real64) :: heartwood_mass_new
```

It also stores structural outputs:

```fortran
real(real64) :: sapwood_area_new
real(real64) :: height_new
real(real64) :: stem_diameter_new
```

Finally, it stores diagnostic residuals:

```fortran
real(real64) :: carbon_balance_error
real(real64) :: leaf_root_residual
real(real64) :: pipe_model_residual
real(real64) :: allocation_residual_final
```

These residuals are not biological processes. They are checks to verify whether the solution satisfies the intended equations.

---

## 11. `contains`

The keyword:

```fortran
contains
```

marks the beginning of the functions and subroutines defined inside the module.

Everything before `contains` defines types, constants, and public/private visibility.

Everything after `contains` defines executable procedures.

---

## 12. Function `height_from_total_stem_carbon`

This function computes plant height from sapwood mass, heartwood mass, wood density, and height-diameter allometry.

It uses:

\[
C_{stem} = S + H
\]

\[
V_{stem} = \frac{C_{stem}}{wooddens}
\]

and:

\[
V_{stem} = h \times \frac{\pi D^2}{4}
\]

with:

\[
h = allom2 \times D^{allom3}
\]

The function eliminates diameter and solves:

\[
h = \left(
\frac{V_{stem} \times allom2^{2/allom3}}{\pi/4}
\right)^{1/(1+2/allom3)}
\]

This function is useful when the host model does not store height explicitly.

It returns zero height if total stem carbon is non-positive.

What it does not do:

- It does not update any carbon pool.
- It does not allocate carbon.
- It does not check the pipe model.
- It only computes a structural height from existing woody carbon.

---

## 13. Function `height_from_pipe_model_state`

This function computes height from the current leaf and sapwood pools using the pipe model.

The pipe model gives:

\[
SA = \frac{L \times SLA}{latosa}
\]

Sapwood mass gives:

\[
S = wooddens \times h \times SA
\]

Therefore:

\[
h = \frac{S}{wooddens \times SA}
\]

This function is valid only if the current leaf and sapwood pools are already approximately consistent with the pipe model.

What it does not do:

- It does not use heartwood.
- It does not enforce height-diameter allometry.
- It does not allocate carbon.

For the offline allocation module, `height_from_total_stem_carbon` is usually the safer first choice because it uses total woody structure directly.

---

## 14. Function `leaf_requirement_for_existing_sapwood`

This function computes the minimum leaf mass required to maintain the current sapwood mass at the current height.

Starting from:

\[
LA = latosa \times SA
\]

and:

\[
LA = L \times SLA
\]

also:

\[
S = wooddens \times h \times SA
\]

so:

\[
SA = \frac{S}{wooddens \times h}
\]

Substitute into the pipe model:

\[
L \times SLA = latosa \times \frac{S}{wooddens \times h}
\]

Solving for leaf mass:

\[
L_{required} = \frac{latosa \times S}{wooddens \times h \times SLA}
\]

This is used to compute the lower bound for normal allocation:

\[
\Delta L_{min} = L_{required} - L
\]

What this function does not do:

- It does not decide whether allocation is normal or abnormal.
- It does not compute root or sapwood increments.
- It only computes the leaf mass required by the current sapwood-height state.

---

## 15. Function `allocation_residual`

This function computes the nonlinear residual used by the bisection solver.

Inputs:

- current plant state,
- allocation parameters,
- available carbon,
- candidate leaf increment `delta_leaf`.

Output:

- residual value \(f(\Delta L)\).

The function asks:

> If this candidate `delta_leaf` were used, would the final plant structure be consistent with both total stem geometry and the pipe model?

### 15.1 Candidate leaf increment

The candidate is:

\[
x = \Delta L
\]

The new leaf mass is:

\[
L_{new} = L + x
\]

### 15.2 Root increment

Functional balance requires:

\[
R_{new} = \frac{L_{new}}{lm2rm}
\]

so:

\[
\Delta R = R_{new} - R
\]

### 15.3 Sapwood increment

Carbon conservation gives:

\[
\Delta S = \Delta B - \Delta L - \Delta R
\]

and:

\[
S_{new} = S + \Delta S
\]

### 15.4 Guard against invalid candidate values

If the candidate produces:

\[
L_{new} \leq 0
\]

or:

\[
S_{new} \leq 0
\]

then the function returns a huge value:

```fortran
f = huge(1.0_real64)
```

This prevents invalid candidate points from being accepted during root search.

### 15.5 Compute transformed height from stem geometry

The function computes:

\[
h_{stem}^{1 + 2/allom3}
\]

from total stem carbon:

\[
S_{new} + H
\]

This is stored as:

```fortran
height_power_from_stem_geometry
```

### 15.6 Compute transformed height from pipe model

The function also computes:

\[
h_{pipe}^{1 + 2/allom3}
\]

from final leaf and sapwood masses:

\[
L_{new}, S_{new}
\]

This is stored as:

```fortran
height_power_from_pipe_model
```

### 15.7 Residual

The function returns:

\[
f(\Delta L) = h_{stem}^{1+2/allom3} - h_{pipe}^{1+2/allom3}
\]

In code:

```fortran
f = height_power_from_stem_geometry - height_power_from_pipe_model
```

The correct normal allocation satisfies:

\[
f(\Delta L) = 0
\]

What this function does not do:

- It does not update plant pools.
- It does not decide allocation pathway.
- It only evaluates whether a proposed leaf increment is structurally consistent.

---

## 16. Subroutine `allocate_woody_individual`

This is the main routine. It receives the old plant state, the parameters, and the available carbon. It returns the allocation result.

### 16.1 Initialize output

The line:

```fortran
result = AllocationResult()
```

resets the output object to default values.

### 16.2 Basic parameter checks

The routine checks whether the allometric parameters are positive:

- SLA must be positive.
- `latosa` must be positive.
- wood density must be positive.
- leaf-to-root ratio must be positive.
- `allom2` and `allom3` must be positive.

If any of these are invalid, the routine returns with an error message.

It also checks whether input height is positive. If not, it returns with an error message.

### 16.3 Compute minimum leaf increment

The routine computes:

\[
L_{required} = \frac{latosa \times S}{wooddens \times h \times SLA}
\]

Then:

\[
\Delta L_{min} = L_{required} - L
\]

This is the minimum leaf increment required to maintain current sapwood under the pipe model.

### 16.4 Compute minimum root increment

The required root mass associated with `leaf_required` is:

\[
R_{required} = \frac{L_{required}}{lm2rm}
\]

Therefore:

\[
\Delta R_{min} = R_{required} - R
\]

### 16.5 Decide whether normal allocation is feasible

Normal allocation is attempted only if:

\[
\Delta L_{min} > 0
\]

\[
\Delta R_{min} > 0
\]

and:

\[
\Delta L_{min} + \Delta R_{min} \leq \Delta B
\]

In words:

> There must be enough carbon available to pay for the minimum leaf and root increments needed to maintain the current sapwood.

If these conditions fail, the routine switches to abnormal allocation.

---

## 17. Normal allocation pathway

If normal allocation is feasible, the solver searches for \(\Delta L\) using bisection.

### 17.1 Lower bound

The lower bound is:

\[
lower = \Delta L_{min}
\]

This is the minimum acceptable leaf increment.

### 17.2 Upper bound

The upper bound is the maximum possible leaf increment if no carbon goes to new sapwood:

\[
\Delta S = 0
\]

Then:

\[
\Delta B = \Delta L + \Delta R
\]

with:

\[
\Delta R = \frac{L+\Delta L}{lm2rm} - R
\]

Substitute:

\[
\Delta B = \Delta L + \frac{L+\Delta L}{lm2rm} - R
\]

Rearrange:

\[
\Delta B + R - \frac{L}{lm2rm}
= \Delta L \left(1 + \frac{1}{lm2rm}\right)
\]

Therefore:

\[
\Delta L_{max} =
\frac{\Delta B - L/lm2rm + R}{1 + 1/lm2rm}
\]

This is the upper bound.

### 17.3 Invalid interval

If:

\[
upper \leq lower
\]

then the bisection interval is invalid. There is no meaningful range in which to search for a normal allocation solution. The routine switches to abnormal allocation.

### 17.4 Find a sign-changing bracket

Bisection requires two points, `left` and `right`, such that:

\[
f(left) \times f(right) \leq 0
\]

This means the residual changes sign between the two points, so a root lies between them.

The routine first evaluates:

\[
f(lower)
\]

and:

\[
f(upper)
\]

If the full interval does not show a sign change, the routine scans the interval in smaller subintervals. The number of subintervals is controlled by:

```fortran
default_scan_segments
```

The scan step is:

\[
scan\_step = \frac{upper - lower}{default\_scan\_segments}
\]

During the scan, the code stores the previous tested point and residual:

```fortran
previous_x
previous_f
```

and compares them with the current point. When a sign change is found, that subinterval becomes the initial bracket for bisection.

### 17.5 Bisection loop

Once a valid bracket is found, the routine repeatedly computes:

\[
mid = \frac{left + right}{2}
\]

and evaluates:

\[
f(mid)
\]

If:

\[
|f(mid)| \leq \text{default\_f\_tolerance}
\]

or:

\[
|right-left| \leq \text{default\_x\_tolerance}
\]

the solution is accepted.

If not, the routine keeps the half of the interval where the sign change remains.

### 17.6 Finalize normal allocation

Once `delta_leaf` is found, the routine calls:

```fortran
finalize_allocation_from_delta_leaf
```

This computes:

\[
\Delta R = \frac{L+\Delta L}{lm2rm} - R
\]

\[
\Delta S = \Delta B - \Delta L - \Delta R
\]

and updates the final pools.

---

## 18. Subroutine `finalize_allocation_from_delta_leaf`

This routine takes the selected `delta_leaf` and computes the final allocation result.

It computes:

\[
\Delta R
\]

from functional balance.

It computes:

\[
\Delta S
\]

from carbon conservation.

Then it updates:

\[
L_{new} = L + \Delta L
\]

\[
R_{new} = R + \Delta R
\]

\[
S_{new} = S + \Delta S
\]

For normal allocation:

\[
H_{new} = H
\]

because heartwood is not directly increased by new carbon allocation in this routine.

The routine also computes:

\[
SA_{new} = \frac{L_{new} \times SLA}{latosa}
\]

Then:

\[
h_{new} = \frac{S_{new}}{wooddens \times SA_{new}}
\]

and:

\[
D_{new} = \left(\frac{h_{new}}{allom2}\right)^{1/allom3}
\]

It also computes diagnostic residuals.

---

## 19. Diagnostic residuals

### 19.1 Carbon balance error

\[
carbon\_balance\_error
=
\Delta B - (\Delta L + \Delta R + \Delta S)
\]

For normal allocation, this should be close to zero.

### 19.2 Leaf-root residual

\[
leaf\_root\_residual
=
L_{new} - lm2rm \times R_{new}
\]

This should be close to zero if functional balance is satisfied.

### 19.3 Pipe-model residual

\[
pipe\_model\_residual
=
LA_{new} - latosa \times SA_{new}
\]

This should be close to zero if the pipe model is satisfied.

### 19.4 Allocation residual final

\[
allocation\_residual\_final = f(\Delta L)
\]

This is the residual of the nonlinear equation solved by bisection.

---

## 20. Abnormal allocation pathway

Abnormal allocation is used when normal allocation is not feasible.

This can happen when the available carbon is not enough to increase leaf, root, and sapwood while satisfying the allometric constraints.

In abnormal allocation, the routine does not solve the full nonlinear normal-allocation problem. Instead, it applies a corrective allocation to restore functional balance and adjust sapwood.

### 20.1 First step: restore leaf-root balance with available carbon

The abnormal routine first assumes:

\[
\Delta B = \Delta L + \Delta R
\]

and:

\[
L + \Delta L = lm2rm \times (R + \Delta R)
\]

From the second equation:

\[
\Delta R = \frac{L+\Delta L}{lm2rm} - R
\]

Substitute into the first:

\[
\Delta B = \Delta L + \frac{L+\Delta L}{lm2rm} - R
\]

Solving for \(\Delta L\):

\[
\Delta L =
\frac{\Delta B - L/lm2rm + R}{1 + 1/lm2rm}
\]

This is the same expression used for the upper bound in normal allocation.

### 20.2 If `delta_leaf > 0`

If the calculated leaf increment is positive, the remaining available carbon is assigned to roots:

\[
\Delta R = \Delta B - \Delta L
\]

If this produces negative root allocation, the routine instead allocates all available carbon to leaves and adjusts roots to satisfy final leaf-root balance.

### 20.3 If `delta_leaf <= 0`

If the calculated leaf increment is negative, the routine assigns all available carbon to roots:

\[
\Delta R = \Delta B
\]

Then it reduces leaves to satisfy:

\[
L_{new} = lm2rm \times R_{new}
\]

so:

\[
\Delta L = lm2rm \times (R + \Delta R) - L
\]

### 20.4 Sapwood adjustment

After leaf and root pools are adjusted, the routine computes the sapwood mass required by the pipe model at the current height:

\[
SA_{required} = \frac{L_{new} \times SLA}{latosa}
\]

\[
S_{required} = wooddens \times h_{old} \times SA_{required}
\]

Therefore:

\[
S_{required} = \frac{L_{new} \times SLA}{latosa} \times wooddens \times h_{old}
\]

The sapwood increment is:

\[
\Delta S = S_{required} - S
\]

In abnormal allocation, this is expected to be negative in many cases, meaning that excess sapwood is converted to heartwood.

### 20.5 Heartwood update

If:

\[
\Delta S < 0
\]

then:

\[
H_{new} = H + |\Delta S|
\]

The code implements this as:

```fortran
result%heartwood_mass_new = state%heartwood_mass + max(-result%delta_sapwood, 0.0_real64)
```

This means reduced sapwood is added to heartwood.

What abnormal allocation does not yet do:

- It does not explicitly return litter fluxes from reduced leaves or roots.
- It does not explicitly return the sapwood-to-heartwood flux as a named output.
- It does not use a storage carbon pool.
- It does not solve the full nonlinear bisection problem.

For full integration, litter and sapwood-to-heartwood fluxes should be tracked explicitly.

---

## 21. What this module currently does

This module:

1. Receives current woody carbon pools.
2. Receives an amount of available carbon for allocation.
3. Receives allometric parameters.
4. Determines whether normal allocation is feasible.
5. If feasible, solves for the leaf increment by bisection.
6. Computes root and sapwood increments from functional balance and carbon conservation.
7. Computes final pools and structural variables.
8. Provides diagnostic residuals.
9. Falls back to abnormal allocation when normal allocation is not feasible.

---

## 22. What this module currently does not do

This module does **not**:

1. Compute GPP.
2. Compute maintenance respiration.
3. Compute growth respiration.
4. Compute daily water stress.
5. Compute phenology.
6. Compute mortality.
7. Compute establishment.
8. Compute tissue turnover.
9. Represent a storage or reserve carbon pool.
10. Accumulate daily NPP.
11. Decide whether allocation should be annual, monthly, seasonal, or daily.
12. Track litter fluxes explicitly in abnormal allocation.
13. Track sapwood-to-heartwood conversion as a named flux output.

These processes must be handled outside this module or added in later versions.

---

## 23. Important caution about timestep

The module can be called for any allocation period, but this does not mean it should be called every day.

If it is called daily without a carbon storage pool, `c_available` may be too small to support structural allocation. This can cause frequent abnormal allocation, artificial sapwood reduction, or unrealistic day-to-day structural changes.

A safer implementation strategy is:

1. Compute photosynthesis, respiration, water stress, and carbon balance daily.
2. Accumulate available carbon over a longer period.
3. Call this allocation module monthly, seasonally, or annually.
4. Test annual allocation first because it is closest to the original LPJ-style logic.
5. Test monthly or seasonal allocation only after the annual offline solver is validated.

---

## 24. Recommended offline testing strategy

Before integrating with the full model, the module should be tested offline.

### Test 1: compilation

Compile with strong checks:

```bash
gfortran -Wall -Wextra -fcheck=all -ffpe-trap=invalid,zero,overflow -c carbon_allocation_offline.f90
```

### Test 2: simple normal allocation

Use a plant state and carbon increment large enough to allow normal allocation.

Expected behavior:

- `normal_allocation = .true.`
- `converged = .true.`
- `carbon_balance_error` close to zero.
- `leaf_root_residual` close to zero.
- `pipe_model_residual` close to zero.

### Test 3: forced abnormal allocation

Use a very small `c_available` to force abnormal allocation.

Expected behavior:

- `normal_allocation = .false.`
- `converged = .true.`
- sapwood may decrease.
- heartwood may increase.
- diagnostic message should report abnormal allocation.

### Test 4: sensitivity to timestep

Run the same total carbon input as:

- one annual call,
- twelve monthly calls,
- many daily calls.

Compare:

- final leaf mass,
- final root mass,
- final sapwood mass,
- final heartwood mass,
- frequency of abnormal allocation.

If daily allocation produces many abnormal events, this is likely a timestep artifact rather than a biological signal.

---

## 25. Recommended future improvements

The current module is a strong first offline solver, but several improvements are recommended before final model integration.

### 25.1 Rename `allocation_residual`

A clearer name would be:

```fortran
allometric_residual_for_delta_leaf
```

because the function specifically computes the allometric residual for a candidate leaf increment.

### 25.2 Use relative tolerances

The current tolerances are absolute. For robust model integration, tolerances should scale with carbon-pool magnitude.

A better approach is:

\[
x_{tol} = \max(x_{abs\_tol}, x_{rel\_tol} \times carbon\_scale)
\]

where `carbon_scale` could be the maximum of available carbon and current living pools.

### 25.3 Add explicit flux outputs

The abnormal allocation routine should eventually return:

- leaf-to-litter flux,
- root-to-litter flux,
- sapwood-to-heartwood flux.

This is important for carbon conservation at the ecosystem level.

### 25.4 Add storage if using daily allocation

If structural allocation is attempted daily, a carbon storage pool should be added:

\[
C_{storage,t+1} = C_{storage,t} + C_{available,t} - C_{allocated,t}
\]

Without this, daily structural allocation may be numerically unstable and biologically unrealistic.

### 25.5 Improve diagnostic outputs

Diagnostics should eventually include:

- allocation pathway,
- bisection bounds,
- bisection iterations,
- final residuals,
- whether abnormal allocation occurred,
- carbon lost to litter,
- sapwood converted to heartwood.

---

## 26. Conceptual summary

The module solves carbon allocation by turning a three-variable allocation problem into a one-variable root-finding problem.

Instead of independently solving for:

\[
\Delta L, \Delta R, \Delta S
\]

it chooses:

\[
x = \Delta L
\]

Then computes:

\[
\Delta R = \frac{L+x}{lm2rm} - R
\]

and:

\[
\Delta S = \Delta B - x - \Delta R
\]

The correct \(x\) is the one that makes stem geometry and pipe-model geometry consistent.

Normal allocation means the model found a positive-growth solution for leaf, root, and sapwood.

Abnormal allocation means the normal solution was not feasible, so the model applies a corrective adjustment that may reduce pools to restore allometry.

The key point is:

> The bisection method does not allocate all carbon by itself. It only finds the leaf increment that makes the allometric constraints close. Once leaf increment is known, root and sapwood increments follow from functional balance and carbon conservation.

