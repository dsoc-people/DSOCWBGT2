# app.py
# -*- coding: utf-8 -*-
# Streamlit App — Kentucky WBGT + WeatherSTEM Map (stable county focus, no flicker)

import requests
import pandas as pd
import numpy as np
import folium
from streamlit_folium import st_folium
import streamlit as st
import inspect
from streamlit_autorefresh import st_autorefresh
import geopandas as gpd
from shapely.geometry import Point
import branca.colormap as cm

# --- Custom Station Coordinates Provided by User ---
station_coords_text = """FARM,36.93,-86.47
RSVL,36.85,-86.92
MRHD,38.22,-83.48
MRRY,36.61,-88.34
PCWN,37.28,-84.96
HTFD,37.45,-86.89
CMBA,37.12,-85.31
CRMT,37.94,-85.67
LXGN,37.93,-84.53
BLRK,37.47,-86.33
SCTV,36.74,-86.21
PRNC,37.09,-87.86
BMBL,36.86,-83.84
PGHL,36.94,-87.48
LSML,38.08,-84.90
ERLN,37.32,-87.49
OLIN,37.36,-83.96
QKSD,37.54,-83.32
SWON,38.53,-84.77
LGNT,37.54,-84.63
MROK,36.95,-85.99
PVRT,37.54,-87.28
BNGL,37.36,-85.49
CRRL,38.67,-85.15
HRDB,37.77,-84.82
FRNY,37.72,-87.90
GRDR,36.79,-85.45
RPTN,37.36,-88.07
ELST,37.71,-84.18
DRFN,36.88,-88.32
BTCK,37.01,-88.96
WLBT,37.83,-85.96
WSHT,37.97,-82.50
WNCH,38.01,-84.13
CCLA,36.67,-88.67
BNVL,37.28,-84.67
RNDH,37.45,-82.99
HCKM,36.85,-88.34
RBSN,37.42,-83.02
HHTS,36.96,-85.64
PRYB,36.83,-83.17
CADZ,36.83,-87.86
ALBN,36.71,-85.14
HUEY,38.97,-84.72
VEST,37.41,-82.99
GRHM,37.82,-87.51
MQDY,37.71,-86.50
CLSL,38.28,-84.10
CHTR,38.58,-83.42
FLRK,36.77,-84.48
DORT,37.28,-82.52
FCHV,38.16,-85.38
LGRN,38.46,-85.47
HDYV,37.26,-85.78
LUSA,38.10,-82.60
PRST,38.09,-83.76
BRND,37.95,-86.22
LRTO,37.63,-85.37
HDGV,37.57,-85.70
WTBG,37.13,-82.84
SWZR,36.67,-86.61
CCTY,37.29,-87.16
ZION,36.76,-87.21
PSPG,37.01,-86.37
BMTN,36.92,-82.91
WDBY,37.18,-86.65
DANV,37.62,-84.82
CROP,38.33,-85.17
HARD,37.76,-86.46
GAMA,36.66,-85.80
DABN,37.18,-84.56
DIXO,37.52,-87.69
WADD,38.09,-85.14
EWPK,37.04,-86.35
RFVC,37.46,-83.16
RFSM,37.43,-83.18
CARL,38.32,-84.04
MONT,36.87,-84.90
BAND,37.13,-88.95
WOOD,36.99,-84.97
DCRD,37.87,-83.65
SPIN,38.13,-84.50
GRBG,37.21,-85.47
PBDY,37.14,-83.58
BLOM,37.96,-85.31
LEWP,37.92,-86.85
STAN,37.85,-83.88
BEDD,38.63,-85.32
"""


# ---------------- Streamlit Setup ----------------
st.set_page_config(page_title="Kentucky WBGT Monitor", layout="wide")
st.title("🌡️ Kentucky WBGT / Weather Map Dashboard")

# Version-safe autorefresh
sig = inspect.signature(st_autorefresh)
if "rerun" in sig.parameters:
    refresh_counter = st_autorefresh(
        interval=5 * 60 * 1000,
        limit=None,
        key="wbgt_refresh",
        rerun=False
    )
else:
    refresh_counter = st_autorefresh(
        interval=5 * 60 * 1000,
        limit=None,
        key="wbgt_refresh"
    )

if "last_map" not in st.session_state:
    st.session_state["last_map"] = None

year = "2025"

# ---------------- Sidebar ----------------
st.sidebar.header("Map Controls")
selected_var = st.sidebar.selectbox(
    "Variable to Display:",
    ["WBGT (°F)", "Temperature (°F)", "Dewpoint (°F)", "Wind Speed (mph)"]
)

# ---------------- County Loader ----------------
@st.cache_data
def load_ky_counties():
    """Return reliable Kentucky county polygons."""
    try:
        url = "https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json"
        j = requests.get(url, timeout=10).json()
        feats = [f for f in j["features"] if f["properties"]["STATE"] == "21"]
        for f in feats:
            f["properties"]["NAME"] = f["properties"]["NAME"].title()
        gdf = gpd.GeoDataFrame.from_features(feats, crs="EPSG:4326")
        return gdf
    except Exception:
        # minimal fallback dataset
        fallback = {
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "properties": {"NAME": "Warren"},
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [[
                            [-86.6, 36.8],
                            [-86.2, 36.8],
                            [-86.2, 37.1],
                            [-86.6, 37.1],
                            [-86.6, 36.8]
                        ]]
                    }
                },
                {
                    "type": "Feature",
                    "properties": {"NAME": "Hardin"},
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [[
                            [-86.2, 37.5],
                            [-85.7, 37.5],
                            [-85.7, 37.9],
                            [-86.2, 37.9],
                            [-86.2, 37.5]
                        ]]
                    }
                },
                {
                    "type": "Feature",
                    "properties": {"NAME": "Daviess"},
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [[
                            [-87.4, 37.6],
                            [-86.9, 37.6],
                            [-86.9, 37.9],
                            [-87.4, 37.9],
                            [-87.4, 37.6]
                        ]]
                    }
                }
            ]
        }
        return gpd.GeoDataFrame.from_features(fallback["features"], crs="EPSG:4326")

counties_gdf = load_ky_counties()
county_list = sorted(counties_gdf["NAME"].unique())
selected_county = st.sidebar.selectbox("Select a Kentucky County:", county_list)

# Mapping for compact IDs to WeatherSTEM URL keys
# Moved up here for wider scope and logical grouping
name_variants = {
    "WKUChaos": "WKU Chaos",
    "WKUCHAOS": "WKU Chaos",
    "Etown": "E'town",
    "WKUIMFields": "WKU IM Fields",
    "Owensboro": "Owensboro",
    "Glasgow": "Glasgow",
    "WKU": "WKU",
}

# ---------------- WeatherSTEM URLs ----------------
urls = {
    "WKU": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wku/latest.json",
    "WKU Chaos": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wkuchaos/latest.json",
    "WKU IM Fields": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wkuimfields/latest.json",
    "E'town": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/hardin/wswelizabethtown/latest.json", # Fixed typo 'wswuelizabethtown'
    "Owensboro": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/daviess/wswowensboro/latest.json",
    "Glasgow": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/barren/wswglasgow/latest.json",
    "Maker's Mark Warehouse": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/marion-ky/makersmarkwarehouse/latest.json",
    "Maker's Mark St Mary": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/marion-ky/makersmarkstmary/latest.json",
    "Maker's Mark Lebanon": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/marion-ky/makersmarklebanon/latest.json",
    "Maker's Mark Innovation Garden": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/marion-ky/makersmark/latest.json",
    "Jim Beam Booker Noe": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/nelson/jbbookernoe/latest.json",
    "Jim Beam Bardstown": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/nelson/jbbardstown/latest.json",
    "Jim Beam Clermont": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/bullitt/jbclermont/latest.json",
    "Jim Beam Old Crow": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/franklin-ky/jboldcrow/latest.json",
    "Jim Beam Grand Dad": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/franklin-ky/jbgranddad/latest.json",
    "Woodford County Courthouse": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/woodford/courthouse/latest.json",
    "Adair County High School": "https://cdn/weatherstem.com/dashboard/data/dynamic/model/adair/achs/latest.json",
    "Clinton County High School": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/clinton/clintonhs/latest.json",
    "Novelis Guthrie": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/todd/novelis/latest.json"
}


# ---------------- Helper Functions ----------------
def extract_value(records, target):
    for r in records:
        if target.lower() in r.get("sensor_name", "").lower():
            return r.get("value")
    return None

def farenheit_to_celsius(f):
    return (f - 32) * 5 / 9

def celsius_to_farenheit(c):
    return c * 9 / 5 + 32

def dbdp2wb(tc, dc, p):
    return (tc + dc) / 2

def wbgt(tempF, mph, rad, bar, dpF):
    tempC = farenheit_to_celsius(tempF)
    mps = mph * 0.44704
    tempK = tempC + 273.15
    tempG = np.nan if rad is None or np.isnan(rad) else (
        tempK + (rad - 30) / (0.0252 * rad + 10.5 * mps + 22.5 + 1e-9) - 273.15
    )
    p = bar * 3.38639
    dpC = farenheit_to_celsius(dpF)
    wbc = dbdp2wb(tempC, dpC, p) if not (
        np.isnan(tempC) or np.isnan(dpC) or np.isnan(p)
    ) else np.nan
    wbgt_c = (
        0.7 * wbc + 0.2 * tempG + 0.1 * tempC
        if not (np.isnan(wbc) or np.isnan(tempG) or np.isnan(tempC))
        else np.nan
    )
    return celsius_to_farenheit(wbgt_c)

# ---------------- Station coords loader (from your text) ----------------
@st.cache_data
def load_station_coords():
    station_data = []
    for line in station_coords_text.strip().split("\n"):
        parts = line.split(",")
        if len(parts) == 3:
            station_data.append(
                {"abbrev": parts[0], "lat": float(parts[1]), "lon": float(parts[2])}
            )
    df = pd.DataFrame(station_data)
    coords = {r["abbrev"]: (r["lat"], r["lon"]) for _, r in df.iterrows()}
    return df, coords

stations_df, station_coords = load_station_coords()

# ---------------- Data Fetch ----------------
@st.cache_data(ttl=300)
def fetch_weatherstem_data(): # Renamed from fetch_weatherstem
    data = []
    for site, url in urls.items():
        try:
            j = requests.get(url, timeout=10).json()
            records = j.get("records", [])
            wbgt_val = extract_value(records, "Wet Bulb Globe Temperature")
            temp = extract_value(records, "Thermometer")
            dew = extract_value(records, "Dewpoint")
            wind = extract_value(records, "Anemometer")
            t = j.get("time", "N/A")
            data.append({
                "name": site,
                "observation_time": t,
                "WBGT (°F)": wbgt_val,
                "Temperature (°F)": temp,
                "Dewpoint (°F)": dew,
                "Wind Speed (mph)": wind,
                "source": "White Squirrel Weather"
            })
        except Exception:
            data.append({
                "name": site,
                "observation_time": "Error",
                "WBGT (°F)": None,
                "Temperature (°F)": None,
                "Dewpoint (°F)": None,
                "Wind Speed (mph)": None,
                "source": "White Squirrel Weather"
            })
    return pd.DataFrame(data)

@st.cache_data(ttl=300)
def fetch_mesonet_data(station_ids, station_coords_map):
    """Fetches and processes data for a list of Mesonet station IDs."""
    mesonet_data_rows = []
    for station_id in station_ids:
        mesonet_data_rows.append(process_station_data(station_id, station_coords_map))
    return pd.DataFrame(mesonet_data_rows)

@st.cache_data(ttl=300)
def process_station_data(station_id, coords):
    """Always return a row for this station, even if data download/parsing fails."""
    lat, lon = coords.get(station_id, (None, None))

    try:
        murl = f"https://d266k7wxhw6o23.cloudfront.net/data/{station_id}/{year}/manifest.json"
        manifest = requests.get(murl, timeout=15).json()
        if not manifest:
            raise ValueError("Empty manifest")

        latest_day = max(manifest.keys())
        key = manifest[latest_day]["key"]
        data = requests.get(
            f"https://d266k7wxhw6o23.cloudfront.net/{key}",
            timeout=15
        ).json()
        df = pd.DataFrame(data["rows"], columns=data["columns"])
        cols = ["TAIR", "DWPT", "WSPD", "SRAD", "PRES", "UTCTimestampCollected"]
        if not all(c in df.columns for c in cols):
            raise ValueError("Missing required columns")

        tair_c, dwpt_c, wspd_mps, srad, pres_hpa = [
            df[c].dropna().iloc[-1] for c in cols[:-1]
        ]
        pres_inhg = pres_hpa * 0.02953
        obs_time = df["UTCTimestampCollected"].dropna().iloc[-1]

        wbgt_f = wbgt(
            celsius_to_farenheit(tair_c),
            wspd_mps * 2.23694,
            srad,
            pres_inhg,
            celsius_to_farenheit(dwpt_c),
        )

        return {
            "name": station_id,
            "latitude": lat,
            "longitude": lon,
            "WBGT (°F)": wbgt_f,
            "Temperature (°F)": celsius_to_farenheit(tair_c),
            "Dewpoint (°F)": celsius_to_farenheit(dwpt_c),
            "Wind Speed (mph)": wspd_mps * 2.23694,
            "observation_time": obs_time,
            "source": "Mesonet",
        }

    except Exception:
        # Return row with None values so marker still shows
        return {
            "name": station_id,
            "latitude": lat,
            "longitude": lon,
            "WBGT (°F)": None,
            "Temperature (°F)": None,
            "Dewpoint (°F)": None,
            "Wind Speed (mph)": None,
            "observation_time": "Error",
            "source": "Mesonet",
        }

if refresh_counter:
    fetch_weatherstem_data.clear()
    fetch_mesonet_data.clear()
    process_station_data.clear() # Clear cache for process_station_data directly if it's cached

with st.spinner("Fetching latest WBGT data..."):
    # Identify Mesonet station abbreviations from the combined list
    mesonet_abbrevs_for_processing = []
    for abbrev in stations_df["abbrev"].tolist():
        normalized_name = name_variants.get(abbrev, abbrev)
        if normalized_name not in urls: # If it's NOT a WeatherSTEM station based on URL list
            mesonet_abbrevs_for_processing.append(abbrev)

    mesonet_df = fetch_mesonet_data(mesonet_abbrevs_for_processing, station_coords)

    ws_df = fetch_weatherstem_data()

# --- WeatherSTEM coordinates from station list ---
# This block assigns latitude/longitude to the ws_df
known_coords_for_ws = {} # Renamed this variable to avoid any potential scope issues
for line in station_coords_text.strip().split("\n"):
    parts = line.split(",")
    if len(parts) != 3:
        continue
    raw_name = parts[0]
    normalized_name = name_variants.get(raw_name, raw_name)
    if normalized_name in urls:
        known_coords_for_ws[normalized_name] = (float(parts[1]), float(parts[2]))

for i, row in ws_df.iterrows():
    lat, lon = known_coords_for_ws.get(row["name"], (None, None))
    ws_df.loc[i, ["latitude", "longitude"]] = lat, lon

# Combine AFTER coordinates are attached
combined = pd.concat([mesonet_df, ws_df], ignore_index=True)

# ---------------- Map Color ----------------
def variable_color(val, var):
    if pd.isna(val):
        return "#808080"
    if var in ["Temperature (°F)", "Dewpoint (°F)"]:
        cmap = cm.LinearColormap(["#0000FF", "#00FF00", "#FF0000"], vmin=30, vmax=100)
        return cmap(val)
    elif var == "WBGT (°F)":
        if val < 66:
            return "#008000"
        elif val < 74:
            return "#FEF200"
        elif val < 83:
            return "#FF0000"
        else:
            return "#000000"
    elif var == "Wind Speed (mph)":
        cmap = cm.LinearColormap(["#FFFFFF", "#00FFFF", "#0000FF"], vmin=0, vmax=20)
        return cmap(val)
    return "#808080"

# ---------------- Main Map ----------------
center_lat = combined["latitude"].dropna().mean()
center_lon = combined["longitude"].dropna().mean()
m = folium.Map(location=[center_lat, center_lon], zoom_start=7, control_scale=True)

for _, row in combined.iterrows():
    lat, lon = row.get("latitude"), row.get("longitude")
    if pd.isna(lat) or pd.isna(lon):
        continue
    val = row.get(selected_var)
    color = variable_color(val, selected_var)
    popup = (
        f"<b>{row['name']} ({row['source']})</b><br>"
        f"{selected_var}: {val if pd.notna(val) else 'N/A'}"
    )
    folium.CircleMarker(
        location=[lat, lon],
        radius=7,
        color="black" if row["source"] == "White Squirrel Weather" else color,
        weight=2 if row["source"] == "White Squirrel Weather" else 1,
        fill=True,
        fill_color=color,
        fill_opacity=0.85,
        popup=popup,
    ).add_to(m)

st.session_state["last_map"] = m
st_folium(st.session_state["last_map"], width=1000, height=650)

# ---------------- County Focus ----------------
st.markdown("### 🧭 County Focus View")
county_geom = counties_gdf[counties_gdf["NAME"] == selected_county].geometry.iloc[0]
county_bounds = county_geom.bounds
county_map = folium.Map(
    location=[
        (county_bounds[1] + county_bounds[3]) / 2,
        (county_bounds[0] + county_bounds[2]) / 2,
    ],
    zoom_start=9,
    control_scale=True,
)

folium.GeoJson(
    county_geom.__geo_interface__,
    style_function=lambda x: {
        "fillColor": "#ff7800",
        "color": "black",
        "weight": 2,
        "fillOpacity": 0.25,
    },
).add_to(county_map)

# Find stations within polygon
points = [
    Point(lon, lat)
    for lon, lat in zip(combined["longitude"], combined["latitude"])
]
combined["in_county"] = [county_geom.contains(p) for p in points]
subset = combined[combined["in_county"]]

for _, row in subset.iterrows():
    val = row.get(selected_var)
    color = variable_color(val, selected_var)
    popup = (
        f"<b>{row['name']} ({row['source']})</b><br>"
        f"{selected_var}: {val if pd.notna(val) else 'N/A'}"
    )
    folium.CircleMarker(
        location=[row.latitude, row.longitude],
        radius=8,
        color="black" if row["source"] == "White Squirrel Weather" else color,
        weight=2 if row["source"] == "White Squirrel Weather" else 1,
        fill=True,
        fill_color=color,
        fill_opacity=0.9,
        popup=popup,
    ).add_to(county_map)

st_folium(county_map, width=850, height=450)
