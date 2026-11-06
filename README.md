# GIS busbar: thermal, short circuit force and enclosure pressure study

Personal project by Saad Ur Rehman. **Every input is ASSUMED.** No value comes from a datasheet, a measurement, a standard
table or an employer. Results describe this model only and are not design values.

## What it does
`gis_busbar.py` checks one metre of a three phase gas insulated switchgear (GIS) busbar segment with a common enclosure:
1. **Electrical and thermal.** AC resistance of an aluminium tube conductor with temperature dependent resistivity and a simple
   skin effect model (current in the outer skin depth of the wall). Joule loss per conductor plus an assumed enclosure eddy loss.
   Steady state temperatures from a coupled network: conductor to enclosure by gas convection and radiation, enclosure to ambient
   by natural convection and radiation, solved with SciPy (brentq) and iterated on resistance.
2. **Short circuit.** Peak current, electromagnetic force per metre for a triangular phase arrangement, bending stress of the
   conductor tube as a simply supported span between support insulators (static, no dynamic factor) and the support force.
3. **Enclosure pressure.** Thin wall hoop stress of the enclosure at the pressure test level.
4. **Sweep.** Conductor wall thickness from 4 to 12 mm, reporting the smallest wall that meets all assumed limits.

Run `python gis_busbar.py` (writes `results/wall_sweep.csv`, `results/summary.json`, `results/study.png`) and `python -m pytest`.

## Results from the run
| wall mm | conductor C | enclosure C | loss W/m per conductor | short circuit stress MPa | passes |
|---|---|---|---|---|---|
| 8 | 132.4 | 71.9 | 180.4 | 219.6 | no |
| 10 | 113.9 | 65.2 | 140.0 | 186.8 | no |
| 12 | 102.0 | 61.0 | 115.2 | 165.5 | yes |

- Smallest wall meeting all limits: **12 mm** (conductor 102.0 C, enclosure 61.0 C, 8.96 kg/m of conductor).
- Short circuit: 100 kA peak, 9.62 kN/m force, 28.9 kN on each support insulator.
- Enclosure: 1.05 MPa test pressure, 39.4 MPa hoop stress, margin 2.03 against the assumed proof stress.
- Skin depth in the conductor at 80 C and 50 Hz: 13.9 mm, so the whole wall of the swept range carries current.

## ASSUMED inputs
- 3150 A rms per phase, 50 Hz, ambient 40 C.
- Conductor: aluminium tube, 100 mm outer diameter, resistivity 2.9e-8 ohm m at 20 C, temperature coefficient 0.004 1/K, emissivity 0.30.
- Enclosure: 600 mm inner diameter, 8 mm wall, emissivity 0.80 inside and 0.90 outside, eddy loss 40% of the conductor losses.
- Convection 6 W/m2K inside the gas and 5 W/m2K outside in air.
- Acceptance limits of this study: conductor 105 C, enclosure 70 C, conductor proof stress 170 MPa, enclosure proof stress 80 MPa.
- Short circuit 40 kA rms with peak factor 2.5, phase centre distance 0.18 m, 3 m span between support insulators.
- Enclosure design pressure 0.70 MPa gauge with a test factor of 1.5.

## Limits
- One metre slice, no axial heat flow, no contacts, joints or plug in connections, no solar load.
- Skin effect is a step model, not a Bessel solution. Proximity effect between phases is ignored.
- Enclosure losses are a fixed fraction, not computed from induced currents.
- Short circuit stress is static (no IEC 60865 dynamic factors, no conductor oscillation), so it is a first estimate.
- Insulating gas properties are folded into the assumed convection coefficient. No dielectric (electric field) design is done.
- Not validated against a measurement. The tests check formulas and the energy balance, not model accuracy.
