#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MA Millionaire Surtax: Panel ETL + DiD/Event-Study Analysis (auto-controls version)
"""

import os
import sys
import warnings
from typing import List, Dict, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

# ---------------- Configuration ----------------

STATES = ["MA", "NH", "VT", "CT", "RI", "ME", "NY"]
YEARS = list(range(2013, 2025))  # 2013–2024 inclusive
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
OUT_DIR  = os.path.join(os.path.dirname(__file__), "data_out")
os.makedirs(OUT_DIR, exist_ok=True)

# Expected filenames
FN_POP = "population.csv"
FN_REV = "revenue.csv"
FN_CPI = "cpi_u.csv"
FN_TAX = "tax_rates.csv"
FN_ACS = "acs_migration.csv"     # optional
FN_IRS = "irs_migration.csv"     # optional
FN_CTL = "controls.csv"          # optional

# Column name fallbacks
POP_COLS = dict(state="state", year="year", population="population")
REV_COLS = dict(state="state", year="year", revenue="total_tax_revenue_nominal")
CPI_COLS = dict(year="year", cpi="cpi")
TAX_COLS = dict(state="state", year="year", top_rate="top_marginal_rate")
ACS_COLS = dict(state="state", year="year", net_mig="net_domestic_migration")
IRS_COLS = dict(state="state", year="year", net_returns="net_migration_returns", net_agi="net_migration_agi")

BASE_CPI_YEAR = 2024  # real dollars base year

def _template_path(name: str) -> str:
    return os.path.join(DATA_DIR, name.replace(".csv", "_template.csv"))

def _ensure_csv_or_template(path: str, headers: List[str]) -> Optional[pd.DataFrame]:
    if os.path.exists(path):
        try:
            return pd.read_csv(path)
        except Exception as e:
            print(f"[ERROR] Failed to read {path}: {e}")
            sys.exit(1)
    else:
        tpath = _template_path(os.path.basename(path))
        if not os.path.exists(tpath):
            pd.DataFrame(columns=headers).to_csv(tpath, index=False)
        print(f"[MISSING] {os.path.basename(path)} not found. Created template at {tpath}")
        return None

def _harmonize(df: pd.DataFrame, mapping: Dict[str, str]) -> pd.DataFrame:
    inv = {}
    for logical, filecol in mapping.items():
        if filecol in df.columns:
            inv[filecol] = logical
    if inv:
        df = df.rename(columns=inv)
    return df

def _state_year_grid(states: List[str], years: List[int]) -> pd.DataFrame:
    return pd.MultiIndex.from_product([states, years], names=["state", "year"]).to_frame(index=False)

def _require_columns(df: pd.DataFrame, required: List[str], tag: str):
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{tag} is missing required columns: {missing}. Found columns: {list(df.columns)}")

def load_population() -> Optional[pd.DataFrame]:
    path = os.path.join(DATA_DIR, FN_POP)
    df = _ensure_csv_or_template(path, [POP_COLS["state"], POP_COLS["year"], POP_COLS["population"]])
    if df is None: return None
    df = _harmonize(df, POP_COLS)
    _require_columns(df, ["state","year","population"], "population.csv")
    return df[df["state"].isin(STATES) & df["year"].between(min(YEARS), max(YEARS))]

def load_revenue() -> Optional[pd.DataFrame]:
    path = os.path.join(DATA_DIR, FN_REV)
    df = _ensure_csv_or_template(path, [REV_COLS["state"], REV_COLS["year"], REV_COLS["revenue"]])
    if df is None: return None
    df = _harmonize(df, REV_COLS)
    _require_columns(df, ["state","year","revenue"], "revenue.csv")
    return df[df["state"].isin(STATES) & df["year"].between(min(YEARS), max(YEARS))]

def load_cpi() -> Optional[pd.DataFrame]:
    path = os.path.join(DATA_DIR, FN_CPI)
    df = _ensure_csv_or_template(path, [CPI_COLS["year"], CPI_COLS["cpi"]])
    if df is None: return None
    df = _harmonize(df, CPI_COLS)
    _require_columns(df, ["year","cpi"], "cpi_u.csv")
    return df[df["year"].between(min(YEARS), max(YEARS))]

def load_tax_rates() -> Optional[pd.DataFrame]:
    path = os.path.join(DATA_DIR, FN_TAX)
    df = _ensure_csv_or_template(path, [TAX_COLS["state"], TAX_COLS["year"], TAX_COLS["top_rate"]])
    if df is None: return None
    df = _harmonize(df, TAX_COLS)
    _require_columns(df, ["state","year","top_rate"], "tax_rates.csv")
    return df[df["state"].isin(STATES) & df["year"].between(min(YEARS), max(YEARS))]

def load_acs_migration() -> Optional[pd.DataFrame]:
    path = os.path.join(DATA_DIR, FN_ACS)
    df = _ensure_csv_or_template(path, [ACS_COLS["state"], ACS_COLS["year"], ACS_COLS["net_mig"]])
    if df is None: return None
    df = _harmonize(df, ACS_COLS)
    _require_columns(df, ["state","year","net_mig"], "acs_migration.csv")
    return df[df["state"].isin(STATES) & df["year"].between(min(YEARS), max(YEARS))]

def load_irs_migration() -> Optional[pd.DataFrame]:
    path = os.path.join(DATA_DIR, FN_IRS)
    df = _ensure_csv_or_template(path, [IRS_COLS["state"], IRS_COLS["year"], IRS_COLS["net_returns"], IRS_COLS["net_agi"]])
    if df is None: return None
    df = _harmonize(df, IRS_COLS)
    _require_columns(df, ["state","year","net_returns","net_agi"], "irs_migration.csv")
    return df[df["state"].isin(STATES) & df["year"].between(min(YEARS), max(YEARS))]

def load_controls() -> Optional[pd.DataFrame]:
    path = os.path.join(DATA_DIR, FN_CTL)
    if not os.path.exists(path):
        tpath = _template_path(FN_CTL)
        if not os.path.exists(tpath):
            pd.DataFrame(columns=["state","year"]).to_csv(tpath, index=False)
        print(f"[OPTIONAL] {FN_CTL} not found. Created template at {tpath}")
        return None
    df = pd.read_csv(path)
    _require_columns(df, ["state","year"], "controls.csv")
    return df[df["state"].isin(STATES) & df["year"].between(min(YEARS), max(YEARS))]

def build_panel(pop, rev, cpi, tax, acs=None, irs=None, controls=None) -> pd.DataFrame:
    panel = _state_year_grid(STATES, YEARS)
    panel = panel.merge(pop, on=["state","year"], how="left")
    panel = panel.merge(rev, on=["state","year"], how="left")
    panel = panel.merge(cpi, on=["year"], how="left", suffixes=("","_cpi"))
    panel = panel.merge(tax, on=["state","year"], how="left")
    if acs is not None:
        panel = panel.merge(acs, on=["state","year"], how="left")
    if irs is not None:
        panel = panel.merge(irs, on=["state","year"], how="left")
    if controls is not None:
        dup = [c for c in controls.columns if c in panel.columns and c not in ["state","year"]]
        if dup:
            controls = controls.rename(columns={c: f"{c}_ctl" for c in dup})
        panel = panel.merge(controls, on=["state","year"], how="left")

    panel["treatment"] = (panel["state"] == "MA").astype(int)
    panel["post"] = (panel["year"] >= 2023).astype(int)
    panel["did"] = panel["treatment"] * panel["post"]

    # Deflate to base
    base_vals = panel.loc[panel["year"]==BASE_CPI_YEAR, "cpi"].dropna().unique()
    if base_vals.size == 0:
        warnings.warn(f"No CPI for {BASE_CPI_YEAR}. Using max available year.")
        latest = panel["year"].dropna().max()
        base_vals = panel.loc[panel["year"]==latest, "cpi"].dropna().unique()
    base_cpi = float(base_vals[0])

    panel["revenue_real"] = panel["revenue"] * (base_cpi / panel["cpi"])
    panel["rev_per_capita_real"] = panel["revenue_real"] / panel["population"]
    panel["rev_per_capita_nominal"] = panel["revenue"] / panel["population"]
    if "net_mig" in panel.columns:
        panel["net_mig_per_1k"] = (panel["net_mig"] / panel["population"]) * 1000.0
    return panel

def detect_controls(panel: pd.DataFrame) -> list:
    exclude = {
        "state","year","population","revenue","cpi","top_rate","net_mig","net_returns","net_agi",
        "treatment","post","did","revenue_real","rev_per_capita_real","rev_per_capita_nominal",
        "net_mig_per_1k"
    }
    ctrl = []
    for c in panel.columns:
        if c in exclude: 
            continue
        if pd.api.types.is_numeric_dtype(panel[c]):
            if panel.groupby("state")[c].nunique().max() > 1:
                ctrl.append(c)
    return ctrl

def run_did(panel: pd.DataFrame, outcome: str, controls: list, out_txt: str, with_controls: bool):
    df = panel.dropna(subset=[outcome, "did"])
    if with_controls and controls:
        rhs = "did + " + " + ".join(controls) + " + C(state) + C(year)"
    else:
        rhs = "did + C(state) + C(year)"
    formula = f"{outcome} ~ {rhs}"
    model = smf.ols(formula, data=df).fit(cov_type="cluster", cov_kwds={"groups": df["state"]})
    with open(out_txt, "w") as f:
        f.write(model.summary().as_text())
        f.write("\n\nControls used: " + (", ".join(controls) if with_controls and controls else "(none)"))
    print(("[DiD + controls]" if with_controls else "[DiD]"), "Results written to", out_txt)

def run_event_study(panel: pd.DataFrame, outcome: str, controls: list,
                    out_txt: str, out_png: str, window_pre: int=6, window_post: int=2,
                    with_controls: bool=False):
    df = panel.copy()
    df["event_time"] = df["year"] - 2023
    df = df[(df["event_time"] >= -window_pre) & (df["event_time"] <= window_post)].copy()
    for k in range(-window_pre, window_post+1):
        if k == -1: continue
        df[f"et_{k}"] = ((df["event_time"] == k).astype(int) * df["treatment"]).astype(int)
    et_terms = " + ".join([f"et_{k}" for k in range(-window_pre, window_post+1) if k != -1])
    if with_controls and controls:
        rhs = f"{et_terms} + " + " + ".join(controls) + " + C(state) + C(year)"
    else:
        rhs = f"{et_terms} + C(state) + C(year)"
    formula = f"{outcome} ~ {rhs}"
    model = smf.ols(formula, data=df).fit(cov_type="cluster", cov_kwds={"groups": df["state"]})
    with open(out_txt, "w") as f:
        f.write(model.summary().as_text())
        f.write("\n\nControls used: " + (", ".join(controls) if with_controls and controls else "(none)"))
    print(("[EventStudy + controls]" if with_controls else "[EventStudy]"), "Results written to", out_txt)
    import numpy as np
    coefs, lower, upper, es = [], [], [], []
    for k in range(-window_pre, window_post+1):
        if k == -1: continue
        term = f"et_{k}"
        if term in model.params.index:
            b = model.params[term]; se = model.bse[term]
            coefs.append(b); lower.append(b - 1.96*se); upper.append(b + 1.96*se); es.append(k)
    plt.figure(figsize=(7,4.5))
    plt.axhline(0, linestyle="--")
    plt.axvline(-1, linestyle=":")
    if es:
        plt.errorbar(es, coefs, yerr=[np.array(coefs)-np.array(lower), np.array(upper)-np.array(coefs)],
                     fmt='o', capsize=3)
    plt.title(f"Event-Study: Effect on {outcome}{' (with controls)' if with_controls else ''}")
    plt.xlabel("Event time (years, 0 = 2023)")
    plt.ylabel("Treatment effect (relative to t=-1)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    print("[Event Study] Plot saved to", out_png)

def main():
    pop = load_population()
    rev = load_revenue()
    cpi = load_cpi()
    tax = load_tax_rates()
    acs = load_acs_migration()
    irs = load_irs_migration()
    controls = load_controls()

    missing = [name for name, df in [(FN_POP,pop),(FN_REV,rev),(FN_CPI,cpi),(FN_TAX,tax)] if df is None]
    if missing:
        print("\n[SETUP NEEDED] Missing required files in ./data:")
        for m in missing: print(" -", m)
        print("Templates created. Fill and re-run.\n")
        _state_year_grid(STATES, YEARS).to_csv(os.path.join(OUT_DIR, "panel_shell.csv"), index=False)
        return

    panel = build_panel(pop, rev, cpi, tax, acs, irs, controls)

    out_csv = os.path.join(OUT_DIR, "panel_final.csv")
    panel.to_csv(out_csv, index=False)
    print("[OK] Panel saved to", out_csv)

    try:
        pivot = panel.pivot_table(index="year", columns="state", values="rev_per_capita_real")
        plt.figure(figsize=(8,5))
        for s in ["MA","NH","VT","CT","RI","ME","NY"]:
            if s in pivot.columns:
                plt.plot(pivot.index, pivot[s], label=s, linewidth=1.8)
        plt.title("Real Revenue per Capita (Pre & Post)")
        plt.xlabel("Year")
        plt.ylabel(f"Real $ per person (base {BASE_CPI_YEAR})")
        plt.legend(ncol=3, fontsize=8)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "pretrends_revpc_real.png"), dpi=180)
        print("[Plot] Pretrend figure saved.")
    except Exception as e:
        print("[WARN] Could not plot pretrends:", e)

    ctrl_vars = detect_controls(panel)
    with open(os.path.join(OUT_DIR, "detected_controls.txt"), "w") as f:
        f.write("\n".join(ctrl_vars) if ctrl_vars else "(none)")

    # DiD
    run_did(panel, "rev_per_capita_real", ctrl_vars,
            os.path.join(OUT_DIR, "did_results_nocontrols.txt"), with_controls=False)
    run_did(panel, "rev_per_capita_real", ctrl_vars,
            os.path.join(OUT_DIR, "did_results_withcontrols.txt"), with_controls=True)

    # Event-study
    run_event_study(panel, "rev_per_capita_real", ctrl_vars,
                    os.path.join(OUT_DIR, "event_study_nocontrols.txt"),
                    os.path.join(OUT_DIR, "event_study_nocontrols.png"),
                    window_pre=6, window_post=2, with_controls=False)
    run_event_study(panel, "rev_per_capita_real", ctrl_vars,
                    os.path.join(OUT_DIR, "event_study_withcontrols.txt"),
                    os.path.join(OUT_DIR, "event_study_withcontrols.png"),
                    window_pre=6, window_post=2, with_controls=True)

if __name__ == "__main__":
    main()
