import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gis_busbar as g


def test_hoop_stress_matches_hand_calc():
    r = g.enclosure_pressure()
    hand = 0.70e6 * 1.5 * 0.600 / (2 * 0.008)
    assert math.isclose(r["sigma_hoop"], hand, rel_tol=1e-12)


def test_force_per_length_matches_hand_calc():
    ip = 2.5 * 40e3
    hand = 2e-7 * (math.sqrt(3) / 2) * ip ** 2 / 0.18
    assert math.isclose(g.short_circuit(0.008)["F_per_m"], hand, rel_tol=1e-9)


def test_energy_balance_closes():
    th = g.thermal(0.008)
    assert math.isclose(g.heat_out_enclosure(th["T_e"]), th["P_tot"], rel_tol=1e-6)
    assert math.isclose(g.heat_conductor_to_enclosure(th["T_c"], th["T_e"]), th["P_c"], rel_tol=1e-6)


def test_thicker_wall_runs_cooler_and_stiffer():
    a, b = g.thermal(0.005), g.thermal(0.010)
    assert b["T_c"] <= a["T_c"]
    assert g.short_circuit(0.010)["sigma"] < g.short_circuit(0.005)["sigma"]
