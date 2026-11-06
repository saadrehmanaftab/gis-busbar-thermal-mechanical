"""Thermal, short circuit force and enclosure pressure study of a three phase GIS busbar segment.

Personal project by Saad Ur Rehman. Every input below is ASSUMED (typical order of magnitude values),
none comes from a datasheet, a measurement or an employer. See README.md for the full list and the limits.
"""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SIGMA = 5.670374419e-8      # Stefan Boltzmann, W/m2K4
MU0 = 4e-7 * math.pi        # H/m
T0 = 273.15

# ---------------------------------------------------------------- ASSUMED inputs
A = dict(
    I_rated=3150.0,          # A rms per phase
    f=50.0,                  # Hz
    D_c=0.100,               # conductor tube outer diameter, m
    rho20=2.9e-8,            # aluminium conductor resistivity at 20 C, ohm m
    alpha=0.004,             # 1/K
    eps_c=0.30,              # conductor emissivity
    D_e=0.600,               # enclosure inner diameter, m
    t_e=0.008,               # enclosure wall, m
    eps_e_in=0.80,           # enclosure inner emissivity (painted)
    eps_e_out=0.90,          # enclosure outer emissivity (painted)
    h_gas=6.0,               # convection conductor to gas to enclosure, W/m2K
    h_air=5.0,               # natural convection enclosure to ambient, W/m2K
    T_amb=40.0,              # C
    enc_loss_frac=0.40,      # enclosure eddy loss as fraction of total conductor loss
    T_c_max=105.0,           # conductor limit, C (acceptance criterion of this study)
    T_e_max=70.0,            # accessible enclosure limit, C (acceptance criterion of this study)
    I_k=40e3,                # short circuit current, A rms
    kappa=2.5,               # peak factor, ip = kappa * I_k
    a=0.18,                  # phase centre distance, triangular arrangement, m
    L_span=3.0,              # span between support insulators, m
    Rp02_c=170e6,            # conductor proof stress, Pa
    p_design=0.70e6,         # enclosure design gauge pressure, Pa
    test_factor=1.5,         # pressure test factor applied to the design pressure
    Rp02_e=80e6,             # enclosure alloy proof stress, Pa
)


def skin_depth(rho, f):
    return math.sqrt(2 * rho / (2 * math.pi * f * MU0))


def conductor_resistance(t_c, T, a=A):
    """AC resistance per metre. Current is taken to flow in the outer min(t, skin depth) of the wall."""
    rho = a["rho20"] * (1 + a["alpha"] * (T - 20.0))
    t_eff = min(t_c, skin_depth(rho, a["f"]))
    area = math.pi / 4 * (a["D_c"] ** 2 - (a["D_c"] - 2 * t_eff) ** 2)
    return rho / area


def heat_out_enclosure(T_e, a=A):
    """Heat leaving 1 m of enclosure to ambient by convection and radiation, W/m."""
    Ao = math.pi * (a["D_e"] + 2 * a["t_e"])
    Te, Ta = T_e + T0, a["T_amb"] + T0
    return Ao * (a["h_air"] * (Te - Ta) + a["eps_e_out"] * SIGMA * (Te ** 4 - Ta ** 4))


def heat_conductor_to_enclosure(T_c, T_e, a=A):
    """Heat from one conductor to the enclosure by gas convection and radiation, W/m."""
    Ac = math.pi * a["D_c"]
    Ae_share = math.pi * a["D_e"] / 3.0
    Tc, Te = T_c + T0, T_e + T0
    rad = SIGMA * Ac * (Tc ** 4 - Te ** 4) / (1 / a["eps_c"] + (Ac / Ae_share) * (1 / a["eps_e_in"] - 1))
    return a["h_gas"] * Ac * (T_c - T_e) + rad


def thermal(t_c, a=A):
    """Steady state conductor and enclosure temperature for wall thickness t_c (iterates on resistance)."""
    T_c = 80.0
    for _ in range(100):
        P_c = a["I_rated"] ** 2 * conductor_resistance(t_c, T_c, a)          # per conductor, W/m
        P_tot = 3 * P_c * (1 + a["enc_loss_frac"])
        T_e = brentq(lambda T: heat_out_enclosure(T, a) - P_tot, a["T_amb"], 400.0)
        T_c_new = brentq(lambda T: heat_conductor_to_enclosure(T, T_e, a) - P_c, T_e, 600.0)
        if abs(T_c_new - T_c) < 1e-6:
            break
        T_c = T_c_new
    return dict(T_c=T_c_new, T_e=T_e, P_c=P_c, P_tot=P_tot)


def short_circuit(t_c, a=A):
    ip = a["kappa"] * a["I_k"]
    F = MU0 / (2 * math.pi) * (math.sqrt(3) / 2) * ip ** 2 / a["a"]            # N/m, triangular arrangement
    D, d = a["D_c"], a["D_c"] - 2 * t_c
    Z = math.pi * (D ** 4 - d ** 4) / (32 * D)
    sigma = F * a["L_span"] ** 2 / (8 * Z)                                     # simply supported span, static
    return dict(ip=ip, F_per_m=F, sigma=sigma, support_force=F * a["L_span"])


def enclosure_pressure(a=A):
    p_test = a["p_design"] * a["test_factor"]
    sigma = p_test * a["D_e"] / (2 * a["t_e"])                                 # thin wall hoop stress
    return dict(p_test=p_test, sigma_hoop=sigma, margin=a["Rp02_e"] / sigma)


def sweep(walls_mm=np.arange(4, 13, 1)):
    rows = []
    for w in walls_mm:
        t = w / 1000.0
        th, sc = thermal(t), short_circuit(t)
        ok = th["T_c"] <= A["T_c_max"] and th["T_e"] <= A["T_e_max"] and sc["sigma"] <= A["Rp02_c"]
        rows.append(dict(wall_mm=w, T_conductor_C=round(th["T_c"], 1), T_enclosure_C=round(th["T_e"], 1),
                         loss_per_conductor_W_m=round(th["P_c"], 1), sc_stress_MPa=round(sc["sigma"] / 1e6, 1),
                         mass_kg_m=round(2700 * math.pi / 4 * (A["D_c"] ** 2 - (A["D_c"] - 2 * t) ** 2), 2), passes=ok))
    return pd.DataFrame(rows)


def main(out=Path(__file__).parent / "results"):
    out.mkdir(exist_ok=True)
    df = sweep()
    df.to_csv(out / "wall_sweep.csv", index=False)
    passing = df[df.passes]
    best = passing.iloc[0] if len(passing) else None
    sc = short_circuit(best.wall_mm / 1000.0) if best is not None else None
    pr = enclosure_pressure()
    summary = dict(
        skin_depth_mm=round(skin_depth(A["rho20"] * (1 + A["alpha"] * 80), A["f"]) * 1000, 1),
        min_wall_mm=int(best.wall_mm) if best is not None else None,
        at_min_wall=best.to_dict() if best is not None else None,
        peak_current_kA=sc["ip"] / 1e3 if sc else None,
        force_kN_per_m=round(sc["F_per_m"] / 1e3, 2) if sc else None,
        support_force_kN=round(sc["support_force"] / 1e3, 1) if sc else None,
        enclosure_test_pressure_MPa=round(pr["p_test"] / 1e6, 2),
        enclosure_hoop_stress_MPa=round(pr["sigma_hoop"] / 1e6, 1),
        enclosure_margin=round(pr["margin"], 2),
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(df.wall_mm, df.T_conductor_C, "o-", label="conductor")
    ax[0].plot(df.wall_mm, df.T_enclosure_C, "s-", label="enclosure")
    ax[0].axhline(A["T_c_max"], ls="--", c="C0"); ax[0].axhline(A["T_e_max"], ls="--", c="C1")
    ax[0].set_xlabel("conductor wall, mm"); ax[0].set_ylabel("temperature, C"); ax[0].legend()
    ax[1].plot(df.wall_mm, df.sc_stress_MPa, "o-"); ax[1].axhline(A["Rp02_c"] / 1e6, ls="--")
    ax[1].set_xlabel("conductor wall, mm"); ax[1].set_ylabel("short circuit bending stress, MPa")
    fig.tight_layout(); fig.savefig(out / "study.png", dpi=150)
    print(df.to_string(index=False)); print(json.dumps(summary, indent=2, default=str))
    return df, summary


if __name__ == "__main__":
    main()
