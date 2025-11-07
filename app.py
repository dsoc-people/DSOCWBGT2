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
import shapely.geometry as geom
from shapely.geometry import Point, shape
import branca.colormap as cm

# ---------------- Streamlit Setup ----------------
st.set_page_config(page_title="Kentucky WBGT Monitor", layout="wide")
st.title("🌡️ Kentucky WBGT / Weather Map Dashboard")

# Version-safe autorefresh
sig = inspect.signature(st_autorefresh)
if "rerun" in sig.parameters:
    refresh_counter = st_autorefresh(interval=5 * 60 * 1000, limit=None, key="wbgt_refresh", rerun=False)
else:
    refresh_counter = st_autorefresh(interval=5 * 60 * 1000, limit=None, key="wbgt_refresh")

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
                {"type": "Feature", "properties": {"NAME": "Warren"},
                 "geometry": {"type": "Polygon",
                              "coordinates": [[[-86.6,36.8],[-86.2,36.8],[-86.2,37.1],[-86.6,37.1],[-86.6,36.8]]]}},
                {"type": "Feature", "properties": {"NAME": "Hardin"},
                 "geometry": {"type": "Polygon",
                              "coordinates": [[[-86.2,37.5],[-85.7,37.5],[-85.7,37.9],[-86.2,37.9],[-86.2,37.5]]]}},
                {"type": "Feature", "properties": {"NAME": "Daviess"},
                 "geometry": {"type": "Polygon",
                              "coordinates": [[[-87.4,37.6],[-86.9,37.6],[-86.9,37.9],[-87.4,37.9],[-87.4,37.6]]]}}
            ]
        }
        return gpd.GeoDataFrame.from_features(fallback["features"], crs="EPSG:4326")

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
    tempG = np.nan if rad is None or np.isnan(rad) else tempK + (rad - 30) / \
        (0.0252 * rad + 10.5 * mps + 22.5 + 1e-9) - 273.15
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
                "name": site, "observation_time": "Error",
                "WBGT (°F)": None, "Temperature (°F)": None,
                "Dewpoint (°F)": None, "Wind Speed (mph)": None,
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
        return {"name": station_id, "latitude": lat, "longitude": lon,
                "WBGT (°F)": wbgt_f, "Temperature (°F)": celsius_to_farenheit(tair_c),
                "Dewpoint (°F)": celsius_to_farenheit(dwpt_c),
                "Wind Speed (mph)": wspd_mps * 2.23694,
                "observation_time": obs_time, "source": "Mesonet"}
    except Exception:
        return None

@st.cache_data
def load_station_coords():
    url = "https://d266k7wxhw6o23.cloudfront.net/metadata/stations_468eb55962c18d1fc333160925381b9d6fb5eb86cd6fbbfbfc285b1d6fcfe7a0.json"
    df = pd.DataFrame(requests.get(url).json())
    coords = {r["abbrev"]: (r["lat"], r["lon"]) for _, r in df.iterrows()}
    return df, coords

stations_df, station_coords = load_station_coords()

if refresh_counter:
    fetch_weatherstem.clear()
    process_station_data.clear()

with st.spinner("Fetching latest WBGT data..."):
    ws_df = fetch_weatherstem()
    mesonet_df = pd.DataFrame([r for r in (process_station_data(s, station_coords)
                                           for s in stations_df["abbrev"].tolist()) if r])
    combined = pd.concat([mesonet_df, ws_df], ignore_index=True)

# Add known WS coordinates
known_coords = {"WKU": (36.9855, -86.4551), "WKU Chaos": (36.9855, -86.4551),
                "WKU IM Fields": (36.9809, -86.4614), "E'town": (37.6939, -85.8594),
                "Owensboro": (37.7719, -87.1112), "Glasgow": (36.9959, -85.9119)}
for i, row in ws_df.iterrows():
    lat, lon = known_coords.get(row["name"], (None, None))
    ws_df.loc[i, ["latitude", "longitude"]] = lat, lon
combined.update(ws_df)

# ---------------- Map Color ----------------
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
    return "#808080"

# ---------------- Main Map ----------------
center_lat = combined["latitude"].dropna().mean()
center_lon = combined["longitude"].dropna().mean()
m = folium.Map(location=[center_lat, center_lon], zoom_start=7, control_scale=True)

for _, row in combined.iterrows():
    lat, lon = row.get("latitude"), row.get("longitude")
    if pd.isna(lat) or pd.isna(lon): continue
    val = row.get(selected_var)
    color = variable_color(val, selected_var)
    popup = f"<b>{row['name']} ({row['source']})</b><br>{selected_var}: {val if pd.notna(val) else 'N/A'}"
    folium.CircleMarker(
        location=[lat, lon],
        radius=7,
        color="black" if row["source"] == "White Squirrel Weather" else color,
        weight=2 if row["source"] == "White Squirrel Weather" else 1,
        fill=True,
        fill_color=color,
        fill_opacity=0.85,
        popup=popup
    ).add_to(m)

st.session_state["last_map"] = m
st_folium(st.session_state["last_map"], width=1000, height=650)

# ---------------- County Focus ----------------
st.markdown("### 🧭 County Focus View")
county_geom = counties_gdf[counties_gdf["NAME"] == selected_county].geometry.iloc[0]
county_bounds = county_geom.bounds
county_map = folium.Map(location=[(county_bounds[1]+county_bounds[3])/2,
                                  (county_bounds[0]+county_bounds[2])/2],
                        zoom_start=9, control_scale=True)

folium.GeoJson(county_geom.__geo_interface__,
               style_function=lambda x: {"fillColor": "#ff7800", "color": "black",
                                         "weight": 2, "fillOpacity": 0.25}).add_to(county_map)

# Find stations within polygon
points = [Point(lon, lat) for lon, lat in zip(combined["longitude"], combined["latitude"])]
combined["in_county"] = [county_geom.contains(p) for p in points]
subset = combined[combined["in_county"]]

for _, row in subset.iterrows():
    val = row.get(selected_var)
    color = variable_color(val, selected_var)
    folium.CircleMarker(
        location=[row.latitude, row.longitude],
        radius=8,
        color="black" if row["source"] == "White Squirrel Weather" else color,
        weight=2 if row["source"] == "White Squirrel Weather" else 1,
        fill=True,
        fill_color=color,
        fill_opacity=0.9,
        popup=f"<b>{row['name']} ({row['source']})</b><br>{selected_var}: {val if pd.notna(val) else 'N/A'}"
    ).add_to(county_map)

st_folium(county_map, width=850, height=450)
