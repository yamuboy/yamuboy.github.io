"""OpenCellID + FCC ETL & Simple Propagation Modeling
=====================================================
This module provides a minimal end‑to‑end pipeline:

1. Download / ingest OpenCellID & FCC data (assumes CSV extracts already on disk).
2. Load into PostGIS (or build a local GeoParquet cache if PostGIS not available).
3. Cluster cell observations into sites (DBSCAN).
4. Assign operator + candidate 5G band(s) using a license table.
5. Create empirical coverage polygons (concave hull of measurements).
6. (Optional) Create modeled coverage via a simple Hata/3GPP Urban Macro pathloss.
7. Expose a simple FastAPI endpoint to serve GeoJSON.

Dependencies (install via pip):
    geopandas, pandas, shapely, sqlalchemy, psycopg2-binary, scikit-learn, numpy,
    fastapi, uvicorn, pyproj, rasterio, shapely[vectorized]

PostGIS connection string example:
    postgresql+psycopg2://user:password@localhost:5432/cellmap

NOTE: Real FCC ULS license footprints involve county/PEA polygons. For brevity we
use a simplified placeholder license CSV with columns:
    operator, band, low_mhz, high_mhz, geometry (WKT polygon)

OpenCellID CSV (subset) columns used:
    radio,mcc,mnc,nci,lat,lon,range,samples,updated_at

This script focuses on NR cells. Extend for LTE similarly.
"""
from __future__ import annotations
import os
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Point
from shapely.ops import unary_union
from sklearn.cluster import DBSCAN
from sqlalchemy import create_engine, text

# --------------------------------------------------
# Configuration
# --------------------------------------------------
POSTGIS_URL = os.getenv("POSTGIS_URL", "postgresql+psycopg2://postgres:postgres@localhost:5432/cellmap")
OCID_CSV = Path("./opencellid_nr_sample.csv")            # Provide your path
FCC_LICENSE_CSV = Path("./fcc_licenses_simplified.csv")  # Provide license polygons
OUTPUT_CACHE = Path("./cache")
OUTPUT_CACHE.mkdir(exist_ok=True)

# DBSCAN parameters for site clustering (approx 150 m)
CLUSTER_EPS_METERS = 150
MIN_SAMPLES = 1

# Coverage modeling
RSRP_THRESHOLD_DBM = -110  # coverage edge
TX_POWER_DBM = 46          # assumed per site EIRP
CABLE_LOSS_DB = 2
ANTENNA_GAIN_DB = 15
NOISE_FIGURE_DB = 5

# --------------------------------------------------
# Helper utilities
# --------------------------------------------------
def load_opencellid(csv_path: Path) -> gpd.GeoDataFrame:
    df = pd.read_csv(csv_path)
    df = df[df["radio"].isin(["NR", "LTE"])]  # allow extension
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat), crs="EPSG:4326")
    return gdf


def load_fcc_licenses(csv_path: Path) -> gpd.GeoDataFrame:
    lic = pd.read_csv(csv_path)
    lic_gdf = gpd.GeoDataFrame(lic, geometry=gpd.GeoSeries.from_wkt(lic["geometry"]), crs="EPSG:4326")
    return lic_gdf


def db_connect():
    engine = create_engine(POSTGIS_URL)
    return engine


# --------------------------------------------------
# Band inference: Map MCC/MNC → operator and intersect license polygons
# --------------------------------------------------
MCC_MNC_OPERATOR = {
    (310, 260): "T-Mobile",
    (311, 480): "Verizon",
    (310, 410): "AT&T",
    # Extend as needed...
}


def assign_operator(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    gdf["operator"] = gdf.apply(lambda r: MCC_MNC_OPERATOR.get((r.get("mcc"), r.get("mnc")), "Unknown"), axis=1)
    return gdf


def infer_bands(cells: gpd.GeoDataFrame, licenses: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    # Spatial join each cell with operator's license polygons
    # Filter license rows by operator first (performance)
    out_bands: List[str] = []
    for idx, row in cells.iterrows():
        op = row.operator
        geom = row.geometry
        lic_subset = licenses[licenses.operator == op]
        match = lic_subset[lic_subset.contains(geom)]
        if match.empty:
            out_bands.append(None)
        else:
            # choose first band arbitrarily; could keep list
            out_bands.append(";".join(sorted(match.band.unique())))
    cells["bands"] = out_bands
    return cells


# --------------------------------------------------
# Site clustering using DBSCAN in a projected CRS (meters)
# --------------------------------------------------
def cluster_sites(cells: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    # Project to local equal area (approx) for distance clustering
    proj = cells.to_crs(3857)
    coords = np.vstack([proj.geometry.x, proj.geometry.y]).T
    eps = CLUSTER_EPS_METERS  # meters ~ WebMercator
    db = DBSCAN(eps=eps, min_samples=MIN_SAMPLES).fit(coords)
    cells["site_id"] = db.labels_
    # Site centroids
    site_geoms = proj.groupby("site_id").geometry.apply(lambda g: unary_union(list(g)).centroid)
    sites = gpd.GeoDataFrame({"site_id": site_geoms.index}, geometry=site_geoms.values, crs=proj.crs)
    sites = sites.to_crs(4326)
    # Aggregate operator + bands
    agg = cells.groupby("site_id").agg({"operator": lambda s: s.mode().iat[0] if not s.mode().empty else "Unknown", "bands": lambda s: ";".join(sorted({b for x in s.dropna() for b in x.split(";")}))})
    sites = sites.merge(agg, on="site_id", how="left")
    return cells, sites


# --------------------------------------------------
# Empirical coverage: concave hull (alpha shape) per site+band
# --------------------------------------------------
from shapely.ops import polygonize
from shapely.geometry import MultiPoint


def alpha_shape(points: List[Point], alpha: float) -> Optional[object]:
    """Compute an alpha shape (concave hull). alpha smaller → more detail."""
    if len(points) < 4:
        return MultiPoint(points).convex_hull
    try:
        import shapely
        from shapely import Delaunay
    except Exception:
        # fallback using scipy
        from scipy.spatial import Delaunay as SciDelaunay
        coords = np.array([[p.x, p.y] for p in points])
        tri = SciDelaunay(coords)
        edges = set()
        for ia, ib, ic in tri.simplices:
            pa, pb, pc = coords[[ia, ib, ic]]
            a = np.linalg.norm(pa - pb)
            b = np.linalg.norm(pb - pc)
            c = np.linalg.norm(pc - pa)
            s = (a + b + c) / 2.0
            area = max(s * (s - a) * (s - b) * (s - c), 0) ** 0.5
            if area == 0: continue
            circum_r = a * b * c / (4.0 * area)
            if circum_r < 1.0 / alpha:
                edges.update({tuple(sorted((ia, ib))), tuple(sorted((ib, ic))), tuple(sorted((ic, ia)))})
        from collections import defaultdict
        edge_segments = []
        for (i, j) in edges:
            edge_segments.append(((coords[i][0], coords[i][1]), (coords[j][0], coords[j][1])))
        m = shapely.geometry.MultiLineString(edge_segments)
        polygons = list(polygonize(m))
        if not polygons:
            return MultiPoint(points).convex_hull
        return unary_union(polygons)
    # If modern shapely has Delaunay convenience (Shapely 2.0), implement separately


def build_empirical_coverage(cells: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    coverage_polys = []
    for (site_id, bands), group in cells.groupby(["site_id", "bands"]):
        if pd.isna(bands):
            continue
        pts = list(group.to_crs(3857).geometry)
        hull = alpha_shape(pts, alpha=0.0015)  # tune alpha
        if hull is None:
            continue
        poly = gpd.GeoDataFrame({"site_id": [site_id], "bands": [bands]}, geometry=[hull], crs=3857).to_crs(4326)
        coverage_polys.append(poly)
    if coverage_polys:
        return pd.concat(coverage_polys, ignore_index=True)
    return gpd.GeoDataFrame(columns=["site_id", "bands", "geometry"], crs=4326)


# --------------------------------------------------
# Simple pathloss / RSRP model (COST231-Hata Urban) for FR1 demonstration
# --------------------------------------------------
def cost231_hata(fc_mhz: float, d_km: np.ndarray, hb_m: float = 30, hm_m: float = 1.5, city: str = "large") -> np.ndarray:
    # fc: 150-2000 MHz; for higher NR bands this is approximate
    a_hm = 3.2 * (np.log10(11.75 * hm_m) ** 2) - 4.97 if city == "large" else (1.1 * np.log10(fc_mhz) - 0.7) * hm_m - (1.56 * np.log10(fc_mhz) - 0.8)
    L = 46.3 + 33.9 * np.log10(fc_mhz) - 13.82 * np.log10(hb_m) - a_hm + (44.9 - 6.55 * np.log10(hb_m)) * np.log10(np.maximum(d_km, 0.001)) + 3  # +3 dB metropolitan
    return L


def rsrp_model(site_point: Point, band_freq_mhz: float, radius_km: float = 5.0, grid_res_m: int = 250) -> gpd.GeoDataFrame:
    # Build grid around site (WebMercator)
    site_merc = gpd.GeoSeries([site_point], crs=4326).to_crs(3857).iloc[0]
    x0, y0 = site_merc.x, site_merc.y
    half = radius_km * 1000
    xs = np.arange(x0 - half, x0 + half, grid_res_m)
    ys = np.arange(y0 - half, y0 + half, grid_res_m)
    xx, yy = np.meshgrid(xs, ys)
    d = np.sqrt((xx - x0) ** 2 + (yy - y0) ** 2) / 1000.0  # km
    pl = cost231_hata(band_freq_mhz, d_km=d, hb_m=30, hm_m=1.5)  # pathloss dB
    tx_eirp_dbm = TX_POWER_DBM + ANTENNA_GAIN_DB - CABLE_LOSS_DB
    rsrp = tx_eirp_dbm - pl  # crude (ignores sectorization)
    mask = rsrp >= RSRP_THRESHOLD_DBM
    if not mask.any():
        return gpd.GeoDataFrame()
    import rasterio.features
    # Build polygon from mask
    transform = rasterio.transform.from_origin(xs[0], ys[0], grid_res_m, grid_res_m)
    shapes = rasterio.features.shapes(mask.astype("uint8"), transform=transform)
    polys = [gpd.GeoSeries.from_wkt(str(g)) for g, val in shapes if val == 1]
    # Simpler: collect cell centers
    points = [Point(x, y) for x, y, m in zip(xx.ravel(), yy.ravel(), mask.ravel()) if m]
    if not points:
        return gpd.GeoDataFrame()
    cov = alpha_shape(points, alpha=0.0008)
    gdf = gpd.GeoDataFrame({"frequency_mhz": [band_freq_mhz]}, geometry=[cov], crs=3857).to_crs(4326)
    return gdf


# --------------------------------------------------
# FastAPI service to serve coverage
# --------------------------------------------------
from fastapi import FastAPI
app = FastAPI(title="Cell Coverage API")

@app.get("/sites")
def list_sites():
    engine = db_connect()
    with engine.connect() as conn:
        df = gpd.read_postgis("SELECT site_id, operator, bands, geom FROM sites", conn, geom_col="geom")
    return json.loads(df.to_crs(4326).to_json())

@app.get("/coverage/empirical")
def empirical():
    engine = db_connect()
    with engine.connect() as conn:
        df = gpd.read_postgis("SELECT site_id, bands, geom FROM coverage_empirical", conn, geom_col="geom")
    return json.loads(df.to_crs(4326).to_json())

# --------------------------------------------------
# 3GPP TR 38.901 Propagation Model (UMa FR1 + FR2 with LOS prob & shadow fading)
# --------------------------------------------------
"""
Adds:
  * LOS probability (p_LOS) per TR 38.901 UMa (Table 7.4.2-1)
  * Median pathloss for LOS/NLOS (FR1 UMa) and LOS/NLOS (FR2 simplified)
  * Shadow fading (log-normal) with standard deviations from Table 7.4.1-1
  * Stochastic combination: PL = p_LOS * (PL_LOS + X_LOS) + (1-p_LOS) * (PL_NLOS + X_NLOS)
    For coverage mapping we draw one realization (set RANDOM_SEED for determinism).

FR1 UMa (sub-6 GHz):
  LOS: see earlier equations.
  NLOS: PL_NLOS = 13.54 + 39.08 log10(d_3D) + 20 log10(fc) - 0.6(h_UT-1.5)
  σ_LOS=4 dB, σ_NLOS=6 dB.

p_LOS_UMa(d_2D):
  For d_2D <= 18 m: p_LOS=1
  Else: p_LOS = exp(-(d_2D-18)/63) * (1 - exp(-d_2D/63)) + exp(-d_2D/63)

FR2 (simplified using 38.901 UMa LOS, NLOS equations):
  LOS: PL_LOS = 32.4 + 21 log10(d_3D) + 20 log10(fc)
  NLOS: PL_NLOS = 38.3 log10(d_3D) + 17.3 + 24.9 log10(fc)
  σ_LOS=4 dB, σ_NLOS=7.82 dB (approx from spec).

NOTE: Full spec distinguishes street canyon, O2I, etc. This is a pragmatic subset.
"""
C = 3e8
RANDOM_SEED = 42
rng = np.random.default_rng(RANDOM_SEED)


def p_los_uma(d2d_m: np.ndarray) -> np.ndarray:
    p = np.ones_like(d2d_m, dtype=float)
    gt = d2d_m > 18
    p[gt] = np.exp(-(d2d_m[gt] - 18) / 63.0) * (1 - np.exp(-d2d_m[gt] / 63.0)) + np.exp(-d2d_m[gt] / 63.0)
    return p


def pl_38901_uma(fr_ghz: float, d2d_m: np.ndarray, h_bs: float = 30.0, h_ut: float = 1.5) -> tuple[np.ndarray, np.ndarray]:
    d3d = np.sqrt(d2d_m**2 + (h_bs - h_ut) ** 2)
    d_bp = 4 * h_bs * h_ut * fr_ghz * 1e9 / C
    pl_los = np.where(
        d2d_m <= d_bp,
        28.0 + 22 * np.log10(np.maximum(d3d, 1.0)) + 20 * np.log10(fr_ghz),
        28.0 + 40 * np.log10(np.maximum(d3d, 1.0)) + 20 * np.log10(fr_ghz) - 9 * np.log10(d_bp**2 + (h_bs - h_ut) ** 2),
    )
    pl_nlos = 13.54 + 39.08 * np.log10(np.maximum(d3d, 1.0)) + 20 * np.log10(fr_ghz) - 0.6 * (h_ut - 1.5)
    return pl_los, pl_nlos


def pl_38901_fr2(fr_ghz: float, d2d_m: np.ndarray, h_bs: float = 30.0, h_ut: float = 1.5) -> tuple[np.ndarray, np.ndarray]:
    d3d = np.sqrt(d2d_m**2 + (h_bs - h_ut) ** 2)
    pl_los = 32.4 + 21 * np.log10(np.maximum(d3d, 1.0)) + 20 * np.log10(fr_ghz)
    pl_nlos = 38.3 * np.log10(np.maximum(d3d, 1.0)) + 17.3 + 24.9 * np.log10(fr_ghz)
    return pl_los, pl_nlos


def modeled_coverage_38901(site_point: Point, band_freq_mhz: float, radius_km: float = 5.0, grid_res_m: int = 200) -> gpd.GeoDataFrame:
    """Generate coverage polygon using TR 38.901 with LOS probability & shadow fading.
    Chooses FR1 UMa if freq < 7 GHz else FR2 simplified. Draws one log-normal sample.
    """
    fr_ghz = band_freq_mhz / 1000.0
    site_merc = gpd.GeoSeries([site_point], crs=4326).to_crs(3857).iloc[0]
    x0, y0 = site_merc.x, site_merc.y
    half = radius_km * 1000
    xs = np.arange(x0 - half, x0 + half, grid_res_m)
    ys = np.arange(y0 - half, y0 + half, grid_res_m)
    xx, yy = np.meshgrid(xs, ys)
    d2d = np.sqrt((xx - x0) ** 2 + (yy - y0) ** 2)

    if fr_ghz < 7:  # FR1 UMa
        pl_los, pl_nlos = pl_38901_uma(fr_ghz, d2d)
        sigma_los, sigma_nlos = 4.0, 6.0
        p_los = p_los_uma(d2d)
    else:  # FR2 simplified
        pl_los, pl_nlos = pl_38901_fr2(fr_ghz, d2d)
        # LOS probability for FR2 UMa approximated using same as FR1 (spec has variant); simplify here
        p_los = p_los_uma(d2d)
        sigma_los, sigma_nlos = 4.0, 7.82

    # Shadow fading draws
    sf_los = rng.normal(0, sigma_los, size=d2d.shape)
    sf_nlos = rng.normal(0, sigma_nlos, size=d2d.shape)

    # Expected pathloss via mixing: For stochastic map choose scenario by Bernoulli(p_LOS)
    bern = rng.random(size=d2d.shape) < p_los
    pl = np.where(bern, pl_los + sf_los, pl_nlos + sf_nlos)

    tx_eirp_dbm = TX_POWER_DBM + ANTENNA_GAIN_DB - CABLE_LOSS_DB
    rsrp = tx_eirp_dbm - pl
    mask = rsrp >= RSRP_THRESHOLD_DBM

    points = [Point(x, y) for x, y, m in zip(xx.ravel(), yy.ravel(), mask.ravel()) if m]
    if not points:
        return gpd.GeoDataFrame()
    cov = alpha_shape(points, alpha=0.0008)
    gdf = gpd.GeoDataFrame(
        {"frequency_mhz": [band_freq_mhz], "model": ["38901_FR1" if fr_ghz < 7 else "38901_FR2"]},
        geometry=[cov], crs=3857
    ).to_crs(4326)
    return gdf

# --------------------------------------------------
# Dockerfile (save as Dockerfile in project root)
# --------------------------------------------------
"""Dockerfile
---------------------------------
FROM python:3.11-slim
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
RUN apt-get update && apt-get install -y build-essential gdal-bin libgdal-dev postgresql-client && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
EXPOSE 8000
CMD ["uvicorn", "this_module:app", "--host", "0.0.0.0", "--port", "8000"]

# requirements.txt should include:
# geopandas pandas shapely sqlalchemy psycopg2-binary scikit-learn fastapi uvicorn rasterio numpy
"""

# --------------------------------------------------
# Main ETL orchestration
# --------------------------------------------------
def run_pipeline():
    print("Loading OpenCellID…")
    cells = load_opencellid(OCID_CSV)
    print(f"Loaded {len(cells)} cells")
    cells = assign_operator(cells)

    print("Loading FCC licenses…")
    licenses = load_fcc_licenses(FCC_LICENSE_CSV)

    print("Inferring bands…")
    cells = infer_bands(cells, licenses)

    print("Clustering sites…")
    cells, sites = cluster_sites(cells)

    print("Building empirical coverage hulls…")
    empirical_cov = build_empirical_coverage(cells[cells.radio == "NR"])

    # Persist to PostGIS
    engine = db_connect()
    with engine.begin() as conn:
        cells.to_postgis("ocid_cells", conn, if_exists="replace")
        sites.to_postgis("sites", conn, if_exists="replace")
        empirical_cov.rename(columns={"geometry": "geom"}).to_postgis("coverage_empirical", conn, if_exists="replace")
    print("Pipeline complete.")


if __name__ == "__main__":
    run_pipeline()
    # To launch API separately:
    # uvicorn this_module:app --reload --port 8000

# --------------------------------------------------
# Example usage snippet (outside of ETL) for propagation models
# --------------------------------------------------
"""Example
from shapely.geometry import Point
import geopandas as gpd

# Assume you have a site centroid (lon, lat) and want a modeled FR1 and FR2 coverage polygon
site = Point(-73.9857, 40.7484)  # Example: NYC

# FR1 example at 3700 MHz (n77)
fr1_cov = modeled_coverage_38901(site, band_freq_mhz=3700.0, radius_km=3.0, grid_res_m=100)
fr1_cov.to_file("coverage_fr1.geojson", driver="GeoJSON")

# FR2 example at 28000 MHz (n257/n261 area)
fr2_cov = modeled_coverage_38901(site, band_freq_mhz=28000.0, radius_km=1.5, grid_res_m=50)
fr2_cov.to_file("coverage_fr2.geojson", driver="GeoJSON")

print("Saved FR1 + FR2 modeled coverage.")
"""

# --------------------------------------------------
# FCC ID Inspection Route for Vendor Attribution
# --------------------------------------------------
"""Workflow
1. Field survey collects (site_id, fcc_id_string, photo_path) -> CSV "fcc_id_observations.csv".
   fcc_id format usually: GRANTEECODE-PRODUCTCODE or GRANTEECODEPRODUCT (dash optional).
2. Script normalizes IDs and queries FCC Equipment Authorization search to extract vendor (grantee) name.
3. Store (site_id, vendor, fcc_id, confidence) into PostGIS 'site_vendors' and update 'sites.vendor'.

NOTE: FCC site uses anti‑automation protections. Heavy scraping may fail or require manual captcha.
Provide a manual override CSV (site_id,vendor) for failed lookups.
"""
import re
import time
import requests
from bs4 import BeautifulSoup

FCC_SEARCH_URL = "https://apps.fcc.gov/oetcf/eas/reports/FccIdSearch.cfm"
HEADERS = {"User-Agent": "Mozilla/5.0 (compatible; CellMapperBot/0.1)"}

FCC_ID_RE = re.compile(r"([A-Za-z0-9]{3})[- ]?([A-Za-z0-9]{1,14})")

def normalize_fcc_id(raw: str) -> str | None:
    if not isinstance(raw, str):
        return None
    m = FCC_ID_RE.search(raw.upper().strip())
    if not m: return None
    return f"{m.group(1)}-{m.group(2)}"  # canonical with dash

def fetch_vendor_for_fcc_id(fcc_id: str, pause: float = 1.0) -> str | None:
    """Return vendor (grantee) name from FCC site or None.
    Simple HTML scrape; may break if FCC layout changes.
    """
    time.sleep(pause)  # be polite
    try:
        resp = requests.get(FCC_SEARCH_URL, params={"fcc_id": fcc_id}, headers=HEADERS, timeout=15)
        if resp.status_code != 200:
            return None
        soup = BeautifulSoup(resp.text, "html.parser")
        # Find table cell containing 'Grantee Code' then read following sibling cells until 'Applicant' or similar
        tables = soup.find_all("table")
        for table in tables:
            tds = table.find_all("td")
            for i, td in enumerate(tds):
                if "Grantee Code" in td.get_text():
                    # vendor often appears in same row under 'Applicant' or 'Grantee' label; search nearby
                    # Heuristic: next non-empty td after this row containing 'Name'/'Applicant'
                    for j in range(i, min(i+10, len(tds))):
                        text = tds[j].get_text(strip=True)
                        if text.lower().startswith("applicant") or text.lower().startswith("grantee"):  # label cell
                            if j+1 < len(tds):
                                vendor = tds[j+1].get_text(strip=True)
                                if vendor:
                                    return vendor
        return None
    except Exception:
        return None


def ingest_vendor_observations(observation_csv: Path, manual_override_csv: Path | None = None):
    obs = pd.read_csv(observation_csv)
    obs["fcc_id_norm"] = obs["fcc_id"].apply(normalize_fcc_id)
    obs = obs.dropna(subset=["fcc_id_norm"])  # remove invalid

    vendors = []
    for fcc_id in obs["fcc_id_norm"].unique():
        vendor = fetch_vendor_for_fcc_id(fcc_id)
        vendors.append({"fcc_id_norm": fcc_id, "vendor": vendor})
    vendor_df = pd.DataFrame(vendors)

    # Merge back
    obs = obs.merge(vendor_df, on="fcc_id_norm", how="left")

    # Apply manual overrides (site_id,vendor)
    if manual_override_csv and manual_override_csv.exists():
        overrides = pd.read_csv(manual_override_csv)
        obs = obs.merge(overrides, on="site_id", how="left", suffixes=("", "_override"))
        obs["vendor_final"] = obs["vendor_override"].fillna(obs["vendor"])
    else:
        obs["vendor_final"] = obs["vendor"]

    # Confidence: 1 if scraped or overridden, else 0
    obs["confidence"] = obs["vendor_final"].apply(lambda v: 1.0 if isinstance(v, str) and v.strip() else 0.0)

    # Load into PostGIS
    engine = db_connect()
    with engine.begin() as conn:
        # Upsert table
        obs[["site_id", "fcc_id_norm", "vendor_final", "confidence"]].rename(columns={"vendor_final": "vendor"}).to_sql("site_vendors", conn, if_exists="replace", index=False)
        # Update sites table vendor column (add if missing)
        conn.execute(text("ALTER TABLE sites ADD COLUMN IF NOT EXISTS vendor text"))
        conn.execute(text("ALTER TABLE sites ADD COLUMN IF NOT EXISTS vendor_confidence numeric"))
        conn.execute(text("UPDATE sites s SET vendor = v.vendor, vendor_confidence = v.confidence FROM site_vendors v WHERE s.site_id = v.site_id"))
    print("Vendor ingestion complete.")

# Example invocation:
# ingest_vendor_observations(Path("fcc_id_observations.csv"), manual_override_csv=Path("vendor_manual.csv"))

# --------------------------------------------------
# Command-line wrapper
# --------------------------------------------------
"""Usage
python this_module.py pipeline \
    --ocid ./opencellid_nr_sample.csv \
    --licenses ./fcc_licenses_simplified.csv

python this_module.py vendor --observations fcc_id_observations.csv --overrides vendor_manual.csv

python this_module.py model --lon -73.9857 --lat 40.7484 --freq 3700 --out coverage_fr1.geojson
"""
import argparse

def cli():
    parser = argparse.ArgumentParser(description="CellMapper ETL / Vendor / Modeling CLI")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_pipe = sub.add_parser("pipeline", help="Run full ETL pipeline")
    p_pipe.add_argument("--ocid", type=Path, default=OCID_CSV)
    p_pipe.add_argument("--licenses", type=Path, default=FCC_LICENSE_CSV)

    p_vendor = sub.add_parser("vendor", help="Ingest FCC ID observations for vendor attribution")
    p_vendor.add_argument("--observations", type=Path, required=True)
    p_vendor.add_argument("--overrides", type=Path, default=None)

    p_model = sub.add_parser("model", help="Generate a modeled coverage polygon GeoJSON")
    p_model.add_argument("--lon", type=float, required=True)
    p_model.add_argument("--lat", type=float, required=True)
    p_model.add_argument("--freq", type=float, required=True, help="Frequency MHz")
    p_model.add_argument("--radius-km", type=float, default=3.0)
    p_model.add_argument("--grid", type=int, default=100, help="Grid resolution meters")
    p_model.add_argument("--out", type=Path, required=True)

    args = parser.parse_args()
    if args.cmd == "pipeline":
        global OCID_CSV, FCC_LICENSE_CSV
        OCID_CSV = args.ocid
        FCC_LICENSE_CSV = args.licenses
        run_pipeline()
    elif args.cmd == "vendor":
        ingest_vendor_observations(args.observations, manual_override_csv=args.overrides)
    elif args.cmd == "model":
        from shapely.geometry import Point
        site = Point(args.lon, args.lat)
        cov = modeled_coverage_38901(site, band_freq_mhz=args.freq, radius_km=args.radius_km, grid_res_m=args.grid)
        if cov.empty:
            print("No coverage polygon generated (likely threshold too high).")
        else:
            cov.to_file(args.out, driver="GeoJSON")
            print(f"Saved modeled coverage to {args.out}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] in {"pipeline", "vendor", "model"}:
        cli()
    else:
        run_pipeline()

