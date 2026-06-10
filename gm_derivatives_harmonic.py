#!/usr/bin/env python3
"""
gm_derivatives_harmonic.py — gm2 / gm3 vs. Vgs by single-tone harmonic
======================================================================
distortion, using the MXG (CW source) + UXA (harmonic receiver) + a
gate-bias SMU, on-wafer.

Principle
---------
Drive the gate with a small tone of amplitude A (volts, at the gate
plane). The drain-current power series gives harmonic amplitudes

    i(f0)  = gm1 * A
    i(2f0) = (gm2 / 4)  * A^2
    i(3f0) = (gm3 / 24) * A^3

so measuring the three harmonic powers delivered to the load yields the
derivatives directly — no numerical differentiation, no noise blow-up.

Amplitude calibration (two options, both implemented)
-----------------------------------------------------
  NOMINAL : A from the MXG set level minus the input-path loss,
            assuming a high-impedance gate at f0 (voltage doubling at
            an unterminated port):  A = sqrt(8 * Z0 * P_avail_tip).
  SELFCAL : if you pass gm1(Vgs) from an independent source (Y21 of
            pulsed S-parameters, or the SMU Id-Vgs derivative), the
            measured fundamental current pins A at every bias point:
            A = i(f0) / gm1.  This removes the input-path and
            gate-impedance assumptions entirely and is the
            recommended mode.

Then:
    gm2 = 4  * i(2f0) / A^2          [A/V^2]
    gm3 = 24 * i(3f0) / A^3          [A/V^3]

Current from measured power: the harmonic power P_n measured by the UXA
(de-embedded to the probe tip by loss_out_db) is delivered into Z0:

    i_n = sqrt(2 * P_n / Z0) * (1 + Z0/rds)        # rds current-divider
                                                   # correction, optional

Validity guards built in
------------------------
  * SLOPE CHECK : at the first bias point the drive is stepped over a
    small range and the HD2/HD3 slopes are fit; they must be ~2 and
    ~3 dB/dB (polynomial regime). Deviation => drive too high (higher-
    order terms) or too low (noise floor).
  * FLOOR CHECK : the noise floor is read next to each harmonic; any
    harmonic within `floor_margin_db` of the floor is flagged NaN
    rather than reported as data.
  * Choose f0 low enough (default 10 MHz) that Cgs/Cgd nonlinearity and
    feedback are negligible and the bias tees still pass it; keep the
    drive small (default -20 dBm at the MXG).

Hardware
--------
  MXG  : CW at f0, ALC ON (CW => ALC is correct here)
  UXA  : SA mode, narrow-span peak-marker reads at f0, 2f0, 3f0
  SMU  : gate bias + Id readback. A generic SCPI SMU class is provided
         (override `GateSMU` for your instrument's dialect); a manual
         prompt fallback is included. Drain supply assumed fixed and
         externally controlled (it is part of the bias condition, not
         swept here).

Output: gm_derivatives.csv with
  vgs, id_ma, p1/p2/p3 (dBm at tip), i1/i2/i3 (A), A_used (V),
  gm1_est (S), gm2 (A/V^2), gm3 (A/V^3)
plus the located gm3 zero crossing (sweet-spot bias), if any.
"""

import logging
import math
import time
from dataclasses import dataclass, field

import numpy as np
import pyvisa

log = logging.getLogger("gm_harm")
logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s")

Z0 = 50.0


# ----------------------------------------------------------------------
@dataclass
class GmConfig:
    # Instruments
    mxg_visa: str = "TCPIP0::192.168.0.20::hislip0::INSTR"
    uxa_visa: str = "TCPIP0::192.168.0.10::hislip0::INSTR"
    smu_visa: str | None = "TCPIP0::192.168.0.30::inst0::INSTR"  # None => manual

    # Tone plan
    f0_hz: float = 10e6            # low: below reactive-NL territory,
                                   # above bias-tee cutoff — verify yours
    p_mxg_dbm: float = -20.0       # small-signal; slope check validates

    # On-wafer path losses at f0 (NOT your 7 GHz values — re-run
    # loss_cal.py at f0; cables/probes are much kinder at 10 MHz)
    loss_in_db: float = 0.8        # MXG -> gate probe tip @ f0
    loss_out_db: float = 21.2      # drain probe tip -> UXA @ f0 (incl. pad)

    # Gate bias sweep
    vgs_start: float = -3.5
    vgs_stop: float = -1.0
    vgs_step: float = 0.025
    vgs_abs_min: float = -6.0      # hard clamp — never exceed
    vgs_abs_max: float = 0.0
    gate_compliance_a: float = 1e-3
    settle_s: float = 0.3          # bias settle (traps want more; raise it)

    # Optional rds correction and self-calibration
    rds_ohm: float | None = None   # e.g. 400.0; None => skip divider corr.
    gm1_table: list[tuple[float, float]] = field(default_factory=list)
    # [(vgs, gm1_S), ...] from Y21 / pulsed S-par. Empty => NOMINAL mode.

    # UXA receiver settings
    span_hz: float = 50e3
    rbw_hz: float = 1e3
    avg_count: int = 10
    floor_margin_db: float = 10.0

    # Slope check
    slope_check: bool = True
    slope_span_db: float = 6.0
    slope_points: int = 4
    slope_tol: float = 0.25        # accept 2.0+/-tol and 3.0+/-tol dB/dB

    csv_out: str = "gm_derivatives.csv"
    timeout_ms: int = 60_000


# ----------------------------------------------------------------------
# Instrument wrappers
# ----------------------------------------------------------------------
def open_visa(addr: str, timeout_ms: int):
    rm = pyvisa.ResourceManager()
    inst = rm.open_resource(addr)
    inst.timeout = timeout_ms
    inst.read_termination = "\n"
    inst.write_termination = "\n"
    log.info("%s -> %s", addr, inst.query("*IDN?").strip())
    return inst


class GateSMU:
    """Generic SCPI SMU (SOUR:VOLT / MEAS:CURR? dialect).

    Subclass and override set_voltage()/read_current()/output() for
    Keithley 24xx (:SOUR:VOLT:LEV, :READ?) or 4200A-SCS setups.
    """

    def __init__(self, addr: str | None, compliance_a: float, timeout_ms: int):
        self.manual = addr is None
        if not self.manual:
            self.inst = open_visa(addr, timeout_ms)
            self.inst.write(f":SENSe:CURRent:PROTection {compliance_a}")

    def output(self, on: bool):
        if self.manual:
            return
        self.inst.write(f":OUTPut:STATe {'ON' if on else 'OFF'}")

    def set_voltage(self, v: float):
        if self.manual:
            input(f"  >> Set Vgs = {v:+.3f} V on the gate supply, then Enter...")
            return
        self.inst.write(f":SOURce:VOLTage {v}")

    def read_current(self) -> float:
        if self.manual:
            try:
                return float(input("  >> Enter measured Ig or Id (A), blank=NaN: ") or "nan")
            except ValueError:
                return float("nan")
        return float(self.inst.query(":MEASure:CURRent?"))


class UXAReceiver:
    """Narrow-span peak-marker power reads in SA mode."""

    def __init__(self, addr: str, cfg: GmConfig):
        self.u = open_visa(addr, cfg.timeout_ms)
        self.cfg = cfg
        u = self.u
        u.write("*CLS;*RST")
        u.write(":INSTrument:SELect SA")
        u.query("*OPC?")
        u.write(":INITiate:CONTinuous OFF")
        u.write(f":SENSe:FREQuency:SPAN {cfg.span_hz}")
        u.write(f":SENSe:BANDwidth:RESolution {cfg.rbw_hz}")
        u.write(":SENSe:AVERage:STATe ON")
        u.write(f":SENSe:AVERage:COUNt {cfg.avg_count}")
        u.write(":DISPlay:WINDow:TRACe:Y:RLEVel -10")
        u.write(":SENSe:POWer:RF:ATTenuation 6")

    def read_power_dbm(self, freq_hz: float) -> float:
        u = self.u
        u.write(f":SENSe:FREQuency:CENTer {freq_hz}")
        u.write(":INITiate:IMMediate")
        u.query("*OPC?")
        u.write(":CALCulate:MARKer1:MAXimum")
        return float(u.query(":CALCulate:MARKer1:Y?"))

    def read_floor_dbm(self, freq_hz: float) -> float:
        """Noise floor estimate just off the harmonic."""
        return self.read_power_dbm(freq_hz + 4 * self.cfg.span_hz)


# ----------------------------------------------------------------------
# Math
# ----------------------------------------------------------------------
def dbm_to_w(p_dbm: float) -> float:
    return 10 ** (p_dbm / 10) / 1e3


def harmonic_current(p_tip_dbm: float, rds_ohm: float | None) -> float:
    """Drain-current harmonic amplitude from tip-referred power into Z0."""
    i_load = math.sqrt(2 * dbm_to_w(p_tip_dbm) / Z0)
    if rds_ohm:
        i_load *= (1 + Z0 / rds_ohm)        # current-divider correction
    return i_load


def nominal_gate_amplitude(p_mxg_dbm: float, loss_in_db: float) -> float:
    """High-Z gate: open-port voltage from available power at the tip."""
    p_avail = dbm_to_w(p_mxg_dbm - loss_in_db)
    return math.sqrt(8 * Z0 * p_avail)


def interp_gm1(table: list[tuple[float, float]], vgs: float) -> float | None:
    if not table:
        return None
    t = sorted(table)
    v = [x[0] for x in t]
    g = [x[1] for x in t]
    if vgs < v[0] or vgs > v[-1]:
        return None
    return float(np.interp(vgs, v, g))


# ----------------------------------------------------------------------
# Measurement
# ----------------------------------------------------------------------
def slope_check(mxg, rx: UXAReceiver, cfg: GmConfig):
    """Verify HD2/HD3 grow at 2 and 3 dB/dB — confirms polynomial regime."""
    levels = np.linspace(cfg.p_mxg_dbm - cfg.slope_span_db,
                         cfg.p_mxg_dbm, cfg.slope_points)
    h2, h3 = [], []
    for p in levels:
        mxg.write(f":SOURce:POWer {p}")
        time.sleep(0.1)
        h2.append(rx.read_power_dbm(2 * cfg.f0_hz))
        h3.append(rx.read_power_dbm(3 * cfg.f0_hz))
    s2 = np.polyfit(levels, h2, 1)[0]
    s3 = np.polyfit(levels, h3, 1)[0]
    log.info("Slope check: HD2 %.2f dB/dB (want 2), HD3 %.2f dB/dB (want 3)",
             s2, s3)
    ok = (abs(s2 - 2) <= cfg.slope_tol) and (abs(s3 - 3) <= cfg.slope_tol)
    if not ok:
        log.warning("Slopes off nominal: drive too high (compression/higher "
                    "orders) if slopes sag, too low (floor) if they flatten. "
                    "Adjust p_mxg_dbm and re-run.")
    mxg.write(f":SOURce:POWer {cfg.p_mxg_dbm}")
    return ok


def measure_point(rx: UXAReceiver, cfg: GmConfig) -> dict:
    """Read the three harmonics + floors at the current bias."""
    out = {}
    for n in (1, 2, 3):
        f = n * cfg.f0_hz
        p = rx.read_power_dbm(f)
        floor = rx.read_floor_dbm(f)
        if p - floor < cfg.floor_margin_db:
            log.warning("  H%d only %.1f dB above floor — flagged invalid",
                        n, p - floor)
            p = float("nan")
        out[f"p{n}_dbm"] = p
    return out


def run(cfg: GmConfig):
    mxg = open_visa(cfg.mxg_visa, cfg.timeout_ms)
    rx = UXAReceiver(cfg.uxa_visa, cfg)
    smu = GateSMU(cfg.smu_visa, cfg.gate_compliance_a, cfg.timeout_ms)

    # --- MXG: clean CW tone, ALC ON (correct for CW)
    mxg.write("*RST")
    mxg.write(":SOURce:RADio:ARB:STATe OFF")
    mxg.write(":OUTPut:MODulation:STATe OFF")
    mxg.write(f":SOURce:FREQuency {cfg.f0_hz}")
    mxg.write(f":SOURce:POWer {cfg.p_mxg_dbm}")
    mxg.write(":SOURce:POWer:ALC:STATe ON")

    # --- Bias plan sanity
    lo, hi = sorted((cfg.vgs_start, cfg.vgs_stop))
    if lo < cfg.vgs_abs_min or hi > cfg.vgs_abs_max:
        raise ValueError("Vgs sweep violates absolute gate limits")
    n_pts = int(round((hi - lo) / cfg.vgs_step)) + 1
    vgs_list = np.linspace(cfg.vgs_start, cfg.vgs_stop, n_pts)

    rows = []
    try:
        smu.output(True)
        smu.set_voltage(vgs_list[0])
        input("Probes down, drain bias applied, gate at start value?  "
              "Enter to enable RF...")
        mxg.write(":OUTPut:STATe ON")
        time.sleep(0.3)

        if cfg.slope_check:
            slope_check(mxg, rx, cfg)

        a_nom = nominal_gate_amplitude(cfg.p_mxg_dbm, cfg.loss_in_db)
        log.info("Nominal gate amplitude A = %.4f V (%.1f mVpk)",
                 a_nom, 1e3 * a_nom)

        for vgs in vgs_list:
            smu.set_voltage(float(vgs))
            time.sleep(cfg.settle_s)
            i_d = smu.read_current()

            m = measure_point(rx, cfg)
            # de-embed to the drain probe tip
            p1 = m["p1_dbm"] + cfg.loss_out_db
            p2 = m["p2_dbm"] + cfg.loss_out_db
            p3 = m["p3_dbm"] + cfg.loss_out_db

            i1 = harmonic_current(p1, cfg.rds_ohm) if not math.isnan(p1) else float("nan")
            i2 = harmonic_current(p2, cfg.rds_ohm) if not math.isnan(p2) else float("nan")
            i3 = harmonic_current(p3, cfg.rds_ohm) if not math.isnan(p3) else float("nan")

            gm1_ext = interp_gm1(cfg.gm1_table, float(vgs))
            if gm1_ext and not math.isnan(i1):
                A = i1 / gm1_ext                  # SELFCAL: pins A per bias
                gm1 = gm1_ext
                mode = "selfcal"
            else:
                A = a_nom                          # NOMINAL
                gm1 = i1 / A if not math.isnan(i1) else float("nan")
                mode = "nominal"

            gm2 = 4.0 * i2 / A**2 if not math.isnan(i2) else float("nan")
            gm3 = 24.0 * i3 / A**3 if not math.isnan(i3) else float("nan")

            log.info("Vgs %+6.3f V  Id %7.3f mA | gm1 %7.2f mS  "
                     "gm2 %8.2f mA/V^2  gm3 %9.2f mA/V^3  [%s]",
                     vgs, 1e3 * i_d, 1e3 * gm1, 1e3 * gm2, 1e3 * gm3, mode)
            rows.append(dict(vgs=float(vgs), id_a=i_d,
                             p1_dbm=p1, p2_dbm=p2, p3_dbm=p3,
                             i1=i1, i2=i2, i3=i3, A=A,
                             gm1=gm1, gm2=gm2, gm3=gm3, cal=mode))
    finally:
        mxg.write(":OUTPut:STATe OFF")
        # leave the gate biased — lifting bias under probes is its own SOP

    # --- gm3 zero crossing (sweet spot). gm3 from HD3 power is a
    # magnitude; the sign flips where |gm3| dips to a sharp minimum.
    g3 = np.array([r["gm3"] for r in rows])
    v = np.array([r["vgs"] for r in rows])
    valid = ~np.isnan(g3)
    if valid.sum() > 4:
        k = int(np.nanargmin(np.abs(g3)))
        log.info("Sharpest |gm3| minimum (candidate sweet spot): "
                 "Vgs = %+.3f V", v[k])
        log.info("NOTE: HD3 measures |gm3|; confirm it is a true null "
                 "(sign change) by the depth of the dip and/or an IM3-vs-"
                 "bias null at the same Vgs.")

    with open(cfg.csv_out, "w") as f:
        cols = list(rows[0].keys())
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(str(r[c]) for c in cols) + "\n")
    log.info("Written %s (%d bias points)", cfg.csv_out, len(rows))
    return rows


if __name__ == "__main__":
    run(GmConfig())
