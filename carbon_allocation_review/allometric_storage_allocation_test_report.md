# Technical report: allometric carbon allocation with labile storage and validation tests

## 1. Purpose of this report

This report documents the current allometric carbon allocation scheme implemented in `carbon_allocation_offline_storage_with_checks.f90` and summarizes the validation tests performed before adding growth respiration. The objective was to verify whether the new gradual allocation routine with a labile storage pool is numerically stable, conserves carbon, and behaves consistently under a broad range of physiological and allometric conditions.

The report focuses on the routine currently used in the storage-based offline tests:

```fortran
allocate_gradual_with_storage(...)
```

The legacy routines `allocate(...)` and `abnormal_allocation(...)` are still present in the module, but they are not called by the current storage-based test drivers. Therefore, the tests described here validate the gradual storage allocation pathway, not the legacy rigid allocation pathway.

---

## 2. Conceptual rationale of the allocation scheme

The new allocation scheme was designed to represent carbon allocation at a daily timestep without forcing the plant to instantaneously satisfy all allometric constraints. In the previous rigid allocation logic, all available carbon was allocated immediately to leaf, fine root, and sapwood compartments in a way that attempted to close the allometric equations at the end of the timestep. This type of formulation is mathematically convenient, but it is biologically strong for daily allocation because real plants do not instantaneously reorganize their architecture every day.

The gradual storage-based scheme follows a different logic:

1. Daily NPP enters a labile carbon storage pool.
2. Structural growth demand is computed from two components:
   - allometric correction demand, which moves the plant toward leaf-root and pipe-model balance;
   - background balanced-growth demand, which allows continued structural growth even when the plant is already close to allometric equilibrium.
3. Carbon is withdrawn from storage and allocated to structural pools only up to the minimum of:
   - available storage after daily NPP input;
   - total structural demand;
   - a maximum daily allocation limit.
4. The allocated carbon is distributed proportionally among leaf, fine root, and sapwood according to their relative demands.
5. Height is updated from total stem carbon rather than by forcing the pipe-model residual to zero at each daily timestep.

This formulation makes allocation sink-limited rather than purely source-driven. It allows temporary carbon accumulation in storage and gradual convergence toward allometric consistency.

---

## 3. Carbon pools and state variables

The current test implementation tracks the following plant carbon compartments:

| Compartment | Role in the allocation scheme |
|---|---|
| `leaf_mass` | Structural carbon in leaves. Determines leaf area through SLA. |
| `root_mass` | Structural carbon in fine roots. Constrained by the target leaf-root ratio. |
| `sapwood_mass` | Conductive stem carbon. Used in the pipe-model relationship. |
| `heartwood_mass` | Non-living or inactive stem carbon. It is kept constant in the current gradual allocation routine. |
| `carbon_storage` | Labile carbon pool that receives daily NPP and pays for structural allocation. |
| `height` | Recomputed from total stem carbon using the height-diameter/allometric relationship. |

At this stage, storage is a labile carbon pool but does not yet include explicit storage respiration, mobilization costs, or reserve turnover. Growth respiration is also not included yet.

---

## 4. Main allometric constraints

### 4.1 Leaf-root balance

The target leaf-root relationship is:

```text
leaf_mass = leaf_to_root_ratio × root_mass
```

In the base test case, `leaf_to_root_ratio = 1`, so the target is approximately:

```text
leaf_mass = root_mass
```

The diagnostic residual is:

```text
leaf_root_residual = leaf_mass - leaf_to_root_ratio × root_mass
```

A positive residual indicates that leaf mass is high relative to root mass. A negative residual indicates that root mass is high relative to leaf mass.

### 4.2 Pipe-model relationship

The pipe model links leaf area and sapwood area:

```text
leaf_area = latosa × sapwood_area
```

Leaf area is computed from leaf carbon and SLA:

```text
leaf_area = leaf_mass × SLA
```

Sapwood area is derived from sapwood mass, wood density, and height:

```text
sapwood_area = sapwood_mass / (wood_density × height)
```

The pipe-model residual is therefore:

```text
pipe_residual = leaf_area - latosa × sapwood_area
```

A negative pipe residual indicates that leaf area is too low relative to the existing sapwood area. A positive pipe residual indicates that leaf area is high relative to sapwood area.

### 4.3 Height update

In the gradual allocation routine, height is not adjusted to force the pipe model to close instantly. Instead, height is recomputed from total stem carbon:

```text
total_stem_carbon = sapwood_mass + heartwood_mass
```

This prevents unrealistic daily height jumps and allows pipe-model residuals to decline gradually through structural allocation.

---

## 5. Daily allocation algorithm

### 5.1 Daily NPP input to storage

Annual or rate-based NPP input is converted to daily carbon input using:

```text
npp_daily = npp_rate × dt_years
```

With a daily timestep:

```text
dt_years = 1 / 365
```

Daily NPP is first added to labile storage:

```text
storage_after_npp = carbon_storage_before + npp_daily
```

If daily NPP is negative and exceeds available storage, storage is clamped at zero and the remaining deficit is reported as:

```text
unmet_storage_deficit
```

This makes the carbon deficit explicit rather than hiding it inside a storage clamp.

### 5.2 Allometric correction demand

The routine computes target structural states based on the current allometric imbalance. The plant is not forced to reach these targets immediately. Instead, the difference between current and target values is converted into a daily correction demand using:

```text
daily_correction_demand = allometric_deficit / allometric_adjustment_days
```

The parameter `allometric_adjustment_days` therefore controls how rapidly the plant moves toward allometric balance.

### 5.3 Background balanced-growth demand

If the plant is already close to allometric equilibrium, pure correction-based allocation can become very small. In that case, positive NPP may accumulate indefinitely in storage even though a biologically realistic plant should continue growing in a balanced way.

To avoid this behavior, the routine includes background balanced-growth demand:

```text
leaf_background_demand = leaf_target × dt_years / leaf_background_timescale_years
root_background_demand = root_target × dt_years / root_background_timescale_years
sapwood_background_demand = sapwood_target × dt_years / sapwood_background_timescale_years
```

In the tested background-growth mode, the timescales were:

| Pool | Background growth timescale |
|---|---:|
| Leaf | 3 years |
| Fine root | 3 years |
| Sapwood | 15 years |

When background growth is disabled, these timescales are set to zero and allocation depends only on explicit allometric correction deficits.

### 5.4 Daily allocation limiter

The total amount of carbon that can be allocated to structural growth in one timestep is limited by:

```text
carbon_to_allocate = min(storage_after_npp, total_demand_daily, max_daily_allocation)
```

where:

```text
max_daily_allocation = max_allocation_fraction × living_carbon
```

and:

```text
living_carbon = leaf_mass + root_mass + sapwood_mass
```

This prevents unrealistically large structural jumps in a single daily timestep.

### 5.5 Proportional distribution among structural pools

Once the total daily structural allocation is determined, carbon is distributed proportionally to the relative demand of each pool:

```text
frac_leaf    = leaf_demand_daily / total_demand_daily
frac_root    = root_demand_daily / total_demand_daily
frac_sapwood = sapwood_demand_daily / total_demand_daily
```

Then:

```text
delta_leaf    = frac_leaf × carbon_to_allocate
delta_root    = frac_root × carbon_to_allocate
delta_sapwood = frac_sapwood × carbon_to_allocate
```

The structural pools are updated by adding these increments. Heartwood is unchanged in the current gradual allocation routine.

---

## 6. Carbon-accounting diagnostics added to the routine

Before adding growth respiration, explicit carbon-accounting diagnostics were added to the allocation output. These diagnostics are essential because growth respiration will introduce an additional carbon sink, and the allocation scheme must be demonstrably conservative before adding that process.

### 6.1 Structural carbon balance

This diagnostic checks whether all carbon assigned to structural growth was actually distributed among leaf, root, and sapwood increments:

```text
structural_balance_error = carbon_to_allocate - (delta_leaf + delta_root + delta_sapwood)
```

Expected value:

```text
structural_balance_error ≈ 0
```

### 6.2 Storage carbon balance

This diagnostic checks whether final storage follows from initial storage, daily NPP, any explicit unmet deficit, and structural allocation:

```text
storage_balance_error = carbon_storage_after -
  (carbon_storage_before + npp_daily + unmet_storage_deficit - carbon_to_allocate)
```

Expected value:

```text
storage_balance_error ≈ 0
```

### 6.3 Whole-plant carbon balance

This diagnostic checks the combined change in structural carbon plus storage:

```text
whole_plant_balance_error =
  (whole_carbon_after - whole_carbon_before) -
  (npp_daily + unmet_storage_deficit)
```

Expected value:

```text
whole_plant_balance_error ≈ 0
```

The `unmet_storage_deficit` term is included because when negative NPP exceeds available storage, the allocation routine reports the unpaid carbon demand rather than allowing storage to become negative.

---

## 7. One-year reference test

A one-year reference simulation was run with:

```text
npp_rate = 3.5
initial leaf = 1.0
initial root = 0.8
initial sapwood = 10.0
initial heartwood = 20.0
initial storage = 0.0
```

The main result after 365 days was:

| Variable | Initial | Final |
|---|---:|---:|
| Leaf mass | 1.000 | 2.121 |
| Root mass | 0.800 | 2.045 |
| Sapwood mass | 10.000 | 10.670 |
| Height | 13.136 | 13.190 |
| Storage | 0.000 | 0.466 |
| Leaf-root residual | 0.200 | 0.076 |
| Pipe-model residual | -12.360 | -0.425 |

This behavior is consistent with the intended gradual allocation logic:

- leaf, root, and sapwood increased;
- root increased proportionally more than leaf, reducing the initial leaf-root imbalance;
- the pipe-model residual became much less negative, indicating movement toward pipe-model consistency;
- storage began to accumulate once daily structural demand became lower than daily NPP input;
- carbon-accounting errors remained near numerical precision, approximately 10^-15.

The annual carbon balance also closed: annual NPP was approximately equal to structural growth plus the increase in storage.

---

## 8. Unit tests performed

A separate Fortran unit-test driver was created:

```text
test_storage_allocation_unit_cases.f90
```

The goal was to test expected behavior in simple, interpretable cases. All tests passed.

| Test | Scenario | Expected behavior | Result |
|---|---|---|---|
| A | Zero NPP and zero storage | No structural growth should occur; storage remains zero. | PASS |
| B | Positive NPP with imbalanced plant | Structural pools grow and allometric residuals improve. | PASS |
| C1 | Equilibrated plant with background growth enabled | The plant can continue balanced structural growth. | PASS |
| C2 | Equilibrated plant with background growth disabled | Structural growth is near zero and carbon accumulates in storage. | PASS |
| D | Positive initial storage and zero NPP | Storage can support structural growth when demand exists. | PASS |
| E1 | Negative NPP with enough storage | Storage is consumed without creating an unmet deficit. | PASS |
| E2 | Negative NPP exceeding storage | Storage goes to zero and the remaining deficit is reported. | PASS |
| F | Very high storage and demand | Daily allocation is constrained by `max_allocation_fraction`. | PASS |

Each unit test also checked:

- no negative structural pools;
- no negative storage;
- no negative structural increments;
- structural balance error near zero;
- storage balance error near zero;
- whole-plant balance error near zero.

---

## 9. Sensitivity-test design

A broader sensitivity-test driver was then created:

```text
test_storage_allocation_sensitivity.f90
```

The test matrix included:

| Dimension | Values tested |
|---|---:|
| NPP rate | -0.5, 0.0, 0.5, 3.5, 8.0 |
| Initial storage | 0.0, 0.5, 5.0 |
| Allometric adjustment time | 30, 365, 730 days |
| Maximum daily allocation fraction | 0.001, 0.005, 0.02 |
| Background-growth mode | enabled, disabled |
| Trait combinations | base, low SLA, high SLA, dense wood |
| Initial plant states | original imbalanced, equilibrated, leaf-rich, root-rich, sapwood-rich |

This resulted in:

```text
5 × 3 × 3 × 3 × 2 × 4 × 5 = 5400 scenarios
```

Each scenario was run for:

```text
5 years = 1825 daily timesteps
```

The output was saved in:

```text
storage_allocation_sensitivity_summary.csv
```

---

## 10. Sensitivity-test results

All 5400 sensitivity scenarios completed successfully.

| Diagnostic | Result |
|---|---:|
| Total scenarios | 5400 |
| Passed scenarios | 5400 |
| Failed scenarios | 0 |
| Maximum absolute structural balance error | 1.11 × 10^-16 |
| Maximum absolute storage balance error | 0.00 |
| Maximum absolute whole-plant balance error | 2.59 × 10^-14 |
| Minimum living pool observed | 0.8 |

These results indicate that the gradual storage allocation routine is numerically stable across the tested parameter space. No NaN, Inf, negative pool, negative storage, or carbon-accounting failure was detected.

### 10.1 Storage accumulation

Although all scenarios passed numerically, the sensitivity tests revealed an important biological/modeling behavior: storage can accumulate substantially when structural demand is low relative to NPP input.

Across the 5400 scenarios:

| Storage criterion | Number of scenarios |
|---|---:|
| Maximum storage fraction > 0.30 | 680 |
| Maximum storage fraction > 0.50 | 216 |

This pattern was strongly associated with the absence of background balanced-growth demand:

| Background-growth mode | Scenarios with storage fraction > 0.30 | Scenarios with storage fraction > 0.50 |
|---|---:|---:|
| Enabled | 45 | 0 |
| Disabled | 635 | 216 |

This confirms that if allocation depends only on allometric correction deficits, an already equilibrated plant can accumulate carbon in storage without converting it into structural biomass. This is numerically valid but may be biologically undesirable unless storage accumulation is explicitly intended and constrained by additional processes such as reserve turnover, storage respiration, or sink limitation.

### 10.2 Maximum-storage case

The maximum storage fraction observed was approximately:

```text
max_storage_fraction = 0.649
```

This occurred under the following conditions:

| Parameter | Value |
|---|---:|
| Background growth | Disabled |
| NPP rate | 8.0 |
| Initial storage | 5.0 |
| Allometric adjustment time | 365 days |
| Maximum daily allocation fraction | 0.001 |
| Initial state | Approximately equilibrated |
| Trait case | Low SLA strategy |
| Cumulative NPP over 5 years | 40.0 |
| Cumulative structural allocation | 0.0 |
| Final storage | 45.0 |

This is an expected result for an equilibrated plant with background growth disabled: because there is no allometric deficit, there is no structural demand, so incoming NPP remains in storage.

### 10.3 Interpretation of structural allocation fractions

Some scenarios produced structural allocation fractions greater than 1 when expressed relative to positive NPP alone. This is not necessarily an error because positive initial storage can also finance structural growth. Therefore, `structural_fraction_of_positive_npp` should be interpreted as a diagnostic index, not as a strict carbon-use efficiency.

---

## 11. Interpretation of the tests

The tests support three main conclusions.

First, the gradual storage allocation routine is numerically robust under the tested conditions. It conserves carbon to numerical precision, avoids negative pools, respects the daily allocation limiter, and handles negative NPP by explicitly reporting unmet storage deficits.

Second, the background balanced-growth term is functionally important. Without it, plants that are close to allometric equilibrium may stop structural growth even under positive NPP, causing storage to accumulate. This behavior is consistent with the equations, but it may not be biologically desirable for normal growth conditions. Therefore, this term should be interpreted as a balanced structural growth demand rather than an arbitrary basal growth term.

Third, the tests validate the numerical basis of the allocation scheme, but they do not by themselves prove that the chosen parameter values are biologically optimal. Parameters such as `allometric_adjustment_days`, `max_allocation_fraction`, and the background growth timescales should still be treated as sensitive model parameters and evaluated against expected plant growth behavior or empirical benchmarks when possible.

---

## 12. Implications for adding growth respiration

The current scheme treats `carbon_to_allocate` as carbon converted into structural growth. Growth respiration has not yet been included. Because the allocation routine is now carbon-conservative without growth respiration, the next step can be implemented and tested cleanly.

The recommended formulation is to define structural growth as the net biomass increment and compute the total carbon cost of that growth using a growth efficiency. For example, if growth respiration is 25% of the total carbon cost, then growth efficiency is 0.75:

```text
growth_efficiency = 0.75
```

For a requested structural increment:

```text
structural_growth = delta_leaf + delta_root + delta_sapwood
```

The total carbon withdrawn from storage should be:

```text
total_growth_cost = structural_growth / growth_efficiency
```

Growth respiration would then be:

```text
growth_respiration = total_growth_cost - structural_growth
```

After this modification, the storage balance should become:

```text
carbon_storage_after = carbon_storage_before + npp_daily - total_growth_cost
```

rather than:

```text
carbon_storage_after = carbon_storage_before + npp_daily - structural_growth
```

This distinction is essential. With growth respiration included, the carbon withdrawn from storage will be larger than the carbon incorporated into biomass.

---

## 13. Recommended next validation step

After implementing growth respiration, the same unit tests and the same 5400-scenario sensitivity matrix should be rerun. The expected changes are:

1. Structural growth should be lower for the same amount of available storage.
2. Storage should decline faster because it pays both structural growth and growth respiration.
3. A new diagnostic should verify that:

```text
total_carbon_withdrawn_from_storage = structural_growth + growth_respiration
```

4. Whole-plant carbon balance should include growth respiration as an explicit carbon loss:

```text
change_in_structural_carbon + change_in_storage + growth_respiration = npp_daily
```

or, depending on the model-level definition of NPP:

```text
change_in_structural_carbon + change_in_storage = NPP_after_growth_respiration
```

The key conceptual decision before coding is whether the input to this allocation routine represents carbon before or after growth respiration. The recommended approach is to pass carbon available before growth respiration, let the allocation routine compute structural growth and growth respiration, and report both explicitly.

---

## 14. Summary

The gradual storage allocation scheme provides a more biologically plausible daily allocation framework than rigid instantaneous allocation because it allows carbon storage, sink limitation, gradual allometric correction, and bounded daily growth. The validation tests show that the current implementation conserves carbon and behaves robustly across a broad range of parameter combinations.

The most important modeling insight from the tests is that background balanced-growth demand is necessary to avoid excessive storage accumulation when plants are already close to allometric equilibrium. Therefore, this term should be retained or replaced by a more mechanistic sink-demand formulation before moving to full growth respiration.

The code is now in a suitable state to proceed to the next development step: adding growth respiration with explicit carbon-accounting diagnostics.

---

## Appendix A. Diagnostic figures

The following plots were generated from the sensitivity-test output and can be used for visual inspection of the storage and allometric responses.

![Final storage fraction by NPP, with background balanced-growth demand enabled.](/mnt/data/storage_fraction_by_npp_background_1.png)

![Final storage fraction by NPP, with background balanced-growth demand disabled.](/mnt/data/storage_fraction_by_npp_background_2.png)

![Fraction of positive NPP converted into structural biomass, with background balanced-growth demand enabled.](/mnt/data/structural_fraction_by_npp_background_1.png)

![Fraction of positive NPP converted into structural biomass, with background balanced-growth demand disabled.](/mnt/data/structural_fraction_by_npp_background_2.png)

![Final pipe-model residual by initial state case.](/mnt/data/pipe_residual_by_state_case.png)
