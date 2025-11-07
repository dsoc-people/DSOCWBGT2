# app.py
# -*- coding: utf-8 -*-
# Streamlit App — Smooth-refresh Mesonet + WeatherSTEM WBGT Map (no flashing, black-bordered WS markers)

import requests
import pandas as pd
import numpy as np
import folium
from streamlit_folium import st_folium
import streamlit as st
from streamlit_autorefresh import st_autorefresh
import geopandas as gpd
import branca.colormap as cm
from io import BytesIO

# ---------------- Streamlit Setup ----------------
st.set_page_config(page_title="Kentucky WBGT Monitor", layout="wide")
st.title("🌡️ Kentucky WBGT / Weather Map Dashboard")

# 🔁 Auto-refresh every 5 minutes (300,000 ms)
# Compatible with streamlit-autorefresh==1.0.1 (no rerun argument)
refresh_counter = st_autorefresh(interval=5 * 60 * 1000, limit=None, key="wbgt_refresh")

year = "2025"

# Maintain previous map between refreshes
if "last_map" not in st.session_state:
    st.session_state["last_map"] = None

# ---------------- Sidebar Controls ----------------
st.sidebar.header("Map Controls")

selected_var = st.sidebar.selectbox(
    "Variable to Display:",
    ["WBGT (°F)", "Temperature (°F)", "Dewpoint (°F)", "Wind Speed (mph)"]
)

# ---------------- Safe County Loader ----------------
@st.cache_data
def load_ky_counties():
    """Load Kentucky counties from Plotly datasets (fallback if offline)."""
    try:
        url = "https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json"
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        gdf = gpd.read_file(BytesIO(resp.content))
        gdf = gdf[gdf["STATE"] == "21"]  # Kentucky FIPS = 21
        gdf["NAME"] = gdf["NAME"].str.title()
        return gdf
    except Exception:
        fallback = {
            "type": "FeatureCollection",
            "features": [
                {"type": "Feature", "properties": {"NAME": "Warren"},
                 "geometry": {"type": "Polygon", "coordinates": [[[-86.6, 36.8], [-86.2, 36.8], [-86.2, 37.1], [-86.6, 37.1], [-86.6, 36.8]]]}},
                {"type": "Feature", "properties": {"NAME": "Hardin"},
                 "geometry": {"type": "Polygon", "coordinates": [[[-86.2, 37.5], [-85.7, 37.5], [-85.7, 37.9], [-86.2, 37.9], [-86.2, 37.5]]]}},
                {"type": "Feature", "properties": {"NAME": "Daviess"},
                 "geometry": {"type": "Polygon", "coordinates": [[[-87.4, 37.6], [-86.9, 37.6], [-86.9, 37.9], [-87.4, 37.9], [-87.4, 37.6]]]}}
            ]
        }
        return gpd.GeoDataFrame.from_features(fallback["features"])

counties_gdf = load_ky_counties()
county_list = sorted(counties_gdf["NAME"].unique())
selected_county = st.sidebar.selectbox("Select a Kentucky County:", county_list)

# ---------------- WeatherSTEM URLs ----------------
urls = {
    "WKU": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wku/latest.json",
    "WKU Chaos": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wkuchaos/latest.json",
    "WKU IM Fields": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wkuimfields/latest.json",
    "E'town": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/hardin/wswelizabethtown/latest.json",
    "Owensboro": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/daviess/wswowensboro/latest.json",
    "Glasgow": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/barren/wswglasgow/latest.json"
}

# ---------------- Helper Functions ----------------
def extract_value(records, target):
    for r in records:
        if target.lower() in r.get("sensor_name", "").lower():
            return r.get("value")
    return None

def farenheit_to_celsius(f): return (f - 32) * 5 / 9
def celsius_to_farenheit(c): return c * 9 / 5 + 32
def dbdp2wb(tc, dc, p): return (tc + dc) / 2

def wbgt(tempF, mph, rad, bar, dpF):
    tempC = farenheit_to_celsius(tempF)
    mps = mph * 0.44704
    tempK = tempC + 273.15
    tempG = np.nan if rad is None or np.isnan(rad) else tempK + (rad - 30) / (0.0252 * rad + 10.5 * mps + 22.5 + 1e-9) - 273.15
    p = bar * 3.38639
    dpC = farenheit_to_celsius(dpF)
    wbc = dbdp2wb(tempC, dpC, p) if not (np.isnan(tempC) or np.isnan(dpC) or np.isnan(p)) else np.nan
    wbgt_c = 0.7 * wbc + 0.2 * tempG + 0.1 * tempC if not (np.isnan(wbc) or np.isnan(tempG) or np.isnan(tempC)) else np.nan
    return celsius_to_farenheit(wbgt_c)

# ---------------- Data Fetch ----------------
@st.cache_data(ttl=300)
def fetch_weatherstem():
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
                "name": site, "observation_time": "Error", "WBGT (°F)": None,
                "Temperature (°F)": None, "Dewpoint (°F)": None, "Wind Speed (mph)": None,
                "source": "White Squirrel Weather"
            })
    return pd.DataFrame(data)

@st.cache_data(ttl=300)
def process_station_data(station_id, coords):
    try:
        murl = f"https://d266k7wxhw6o23.cloudfront.net/data/{station_id}/{year}/manifest.json"
        manifest = requests.get(murl, timeout=15).json()
        latest_day = max(manifest.keys())
        key = manifest[latest_day]["key"]
        data = requests.get(f"https://d266k7wxhw6o23.cloudfront.net/{key}", timeout=15).json()
        df = pd.DataFrame(data["rows"], columns=data["columns"])
        cols = ["TAIR", "DWPT", "WSPD", "SRAD", "PRES", "UTCTimestampCollected"]
        if not all(c in df.columns for c in cols): return None
        tair_c, dwpt_c, wspd_mps, srad, pres_hpa = [df[c].dropna().iloc[-1] for c in cols[:-1]]
        pres_inhg = pres_hpa * 0.02953
        obs_time = df["UTCTimestampCollected"].dropna().iloc[-1]
        wbgt_f = wbgt(celsius_to_farenheit(tair_c), wspd_mps * 2.23694, srad, pres_inhg, celsius_to_farenheit(dwpt_c))
        lat, lon = coords.get(station_id, (None, None))
        return {
            "name": station_id, "latitude": lat, "longitude": lon,
            "WBGT (°F)": wbgt_f,
            "Temperature (°F)": celsius_to_farenheit(tair_c),
            "Dewpoint (°F)": celsius_to_farenheit(dwpt_c),
            "Wind Speed (mph)": wspd_mps * 2.23694,
            "observation_time": obs_time, "source": "Mesonet"
        }
    except Exception:
        return None

@st.cache_data
def load_station_coords():
    url = "https://d266k7wxhw6o23.cloudfront.net/metadata/stations_468eb55962c18d1fc333160925381b9d6fb5eb86cd6fbbfbfc285b1d6fcfe7a0.json"
    df = pd.DataFrame(requests.get(url).json())
    coords = {r["abbrev"]: (r["lat"], r["lon"]) for _, r in df.iterrows()}
    return df, coords

stations_df, station_coords = load_station_coords()

# Refresh only data cache every cycle
if refresh_counter:
    fetch_weatherstem.clear()
    process_station_data.clear()

# ---------------- Build Dataset ----------------
with st.spinner("Fetching latest WBGT data..."):
    ws_df = fetch_weatherstem()
    mesonet_df = pd.DataFrame([r for r in (process_station_data(s, station_coords)
                                           for s in stations_df["abbrev"].tolist()) if r])
    combined = pd.concat([mesonet_df, ws_df], ignore_index=True)

# Add known WS coordinates
known_coords = {
    "WKU": (36.9855, -86.4551),
    "WKU Chaos": (36.9855, -86.4551),
    "WKU IM Fields": (36.9809, -86.4614),
    "E'town": (37.6939, -85.8594),
    "Owensboro": (37.7719, -87.1112),
    "Glasgow": (36.9959, -85.9119),
}
for i, row in ws_df.iterrows():
    lat, lon = known_coords.get(row["name"], (None, None))
    ws_df.loc[i, "latitude"], ws_df.loc[i, "longitude"] = lat, lon
combined.update(ws_df)

# ---------------- Color Function ----------------
def variable_color(val, var):
    if pd.isna(val): return "#808080"
    if var in ["Temperature (°F)", "Dewpoint (°F)"]:
        cmap = cm.LinearColormap(["#0000FF", "#00FF00", "#FF0000"], vmin=30, vmax=100)
        return cmap(val)
    elif var == "WBGT (°F)":
        if val < 66: return "#008000"
        elif val < 74: return "#FEF200"
        elif val < 83: return "#FF0000"
        else: return "#000000"
    elif var == "Wind Speed (mph)":
        cmap = cm.LinearColormap(["#FFFFFF", "#00FFFF", "#0000FF"], vmin=0, vmax=20)
        return cmap(val)
    else:
        return "#808080"

# ---------------- Build Map ----------------
center_lat = combined["latitude"].dropna().mean()
center_lon = combined["longitude"].dropna().mean()
m = folium.Map(location=[center_lat, center_lon], zoom_start=7, control_scale=True)

mesonet_layer = folium.FeatureGroup(name="Mesonet")
ws_layer = folium.FeatureGroup(name="White Squirrel Weather")

for _, row in combined.iterrows():
    lat, lon = row.get("latitude"), row.get("longitude")
    if pd.isna(lat) or pd.isna(lon):
        continue
    val = row.get(selected_var)
    popup = f"<b>{row['name']} ({row['source']})</b><br>{selected_var}: {val if pd.notna(val) else 'N/A'}<br>Obs: {row.get('observation_time','N/A')}"
    color = variable_color(val, selected_var)

    if row["source"] == "White Squirrel Weather":
        # Black-bordered WeatherSTEM marker
        folium.CircleMarker(
            location=[lat, lon],
            radius=7,
            color="black",
            weight=2,
            fill=True,
            fill_color=color,
            fill_opacity=0.85,
            popup=folium.Popup(popup, max_width=250),
