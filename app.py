# app.py
# -*- coding: utf-8 -*-
import requests
import pandas as pd
import numpy as np
import folium
from streamlit_folium import st_folium
import streamlit as st
from branca.element import MacroElement
from jinja2 import Template

# --- Streamlit Page Setup ---
st.set_page_config(page_title="Kentucky WBGT Map", layout="wide")

st.title("🌡️ Kentucky Mesonet — Wet Bulb Globe Temperature (WBGT)")
st.markdown("""
Displays the **current estimated WBGT (°F)** across all available Kentucky Mesonet stations.  
Data are pulled live from the Mesonet API (`d266k7wxhw6o23.cloudfront.net`).
""")

# --- Config ---
BASE = "https://d266k7wxhw6o23.cloudfront.net/"
YEAR = "2025"
GOOGLE_API_KEY = "YOUR_GOOGLE_MAPS_API_KEY"  # <-- replace with your real key

# --- Helper Functions ---
def farenheit_to_celsius(temp_f):
    return (temp_f - 32) * 5 / 9

def celsius_to_farenheit(temp_c):
    return temp_c * 9 / 5 + 32

def dbdp2wb(tempC, dpC, p):
    """Simplified wet-bulb approximation."""
    return (tempC + dpC) / 2

def wbgt(tempF, mph, rad, bar, dpF):
    tempC = farenheit_to_celsius(tempF)
    mps = mph * 0.44704
    tempK = tempC + 273.15

    # Globe temperature approximation
    if rad is None or np.isnan(rad):
        tempG = np.nan
    else:
        tempG = tempK + (rad - 30) / (0.0252 * rad + 10.5 * mps + 22.5 + 1e-9)
        tempG -= 273.15

    # Convert pressure to kPa
    p = bar * 3.38639
    dpC = farenheit_to_celsius(dpF)

    if np.isnan(tempC) or np.isnan(dpC) or np.isnan(p):
        wbc = np.nan
    else:
        wbc = dbdp2wb(tempC, dpC, p)

    if np.isnan(wbc) or np.isnan(tempG) or np.isnan(tempC):
        wbgt_c = np.nan
    else:
        wbgt_c = 0.7 * wbc + 0.2 * tempG + 0.1 * tempC

    return celsius_to_farenheit(wbgt_c)

# --- Load Station Metadata ---
@st.cache_data
def load_station_metadata():
    url = f"{BASE}metadata/manifest.json"
    manifest = requests.get(url).json()
    stations_key = manifest["stations"]["key"]
    stations_url = f"{BASE}{stations_key}"
    data = requests.get(stations_url).json()
    df = pd.DataFrame(data)
    df.rename(columns={"abbrev": "Station", "lat": "Latitude", "lon": "Longitude"}, inplace=True)
    return df

stations_df = load_station_metadata()
station_coords = {row["Station"]: (row["Latitude"], row["Longitude"]) for _, row in stations_df.iterrows()}

# --- Process Single Station ---
@st.cache_data(ttl=900)
def process_station(station_id):
    try:
        manifest_url = f"{BASE}data/{station_id}/{YEAR}/manifest.json"
        manifest = requests.get(manifest_url).json()
        latest_day = max(manifest.keys())
        key = manifest[latest_day]["key"]
        data_url = f"{BASE}{key}"
        data = requests.get(data_url).json()
        df = pd.DataFrame(data["rows"], columns=data["columns"])

        req_cols = ["TAIR", "DWPT", "WSPD", "SRAD", "PRES", "UTCTimestampCollected"]
        if not all(col in df.columns for col in req_cols):
            return None

        tair_c = df["TAIR"].dropna().iloc[-1]
        dwpt_c = df["DWPT"].dropna().iloc[-1]
        wspd_mps = df["WSPD"].dropna().iloc[-1]
        srad = df["SRAD"].dropna().iloc[-1]
        pres_hpa = df["PRES"].dropna().iloc[-1]
        pres_inhg = pres_hpa * 0.02953
        obs_time = df["UTCTimestampCollected"].dropna().iloc[-1]

        tair_f = celsius_to_farenheit(tair_c)
        dwpt_f = celsius_to_farenheit(dwpt_c)
        wspd_mph = wspd_mps * 2.23694

        wbgt_f = wbgt(tair_f, wspd_mph, srad, pres_inhg, dwpt_f)
        lat, lon = station_coords.get(station_id, (None, None))

        return {
            "Station": station_id,
            "WBGT (°F)": wbgt_f,
            "Latitude": lat,
            "Longitude": lon,
            "Time": obs_time
        }
    except Exception:
        return None

# --- Compute All WBGTs ---
with st.spinner("Loading latest Mesonet station data…"):
    results = [process_station(station) for station in stations_df["Station"].tolist()]
    results = [r for r in results if r is not None]
    wbgt_df = pd.DataFrame(results)

if wbgt_df.empty:
    st.error("No valid WBGT data could be retrieved.")
    st.stop()

# --- Map Rendering ---
avg_lat = wbgt_df["Latitude"].mean()
avg_lon = wbgt_df["Longitude"].mean()
m = folium.Map(location=[avg_lat, avg_lon], zoom_start=7, tiles=None)

# Google Maps Satellite + Labels
folium.TileLayer(
    tiles=f"https://mt1.google.com/vt/lyrs=y&x={{x}}&y={{y}}&z={{z}}&key={GOOGLE_API_KEY}",
    attr="Google Maps",
    name="Google Satellite",
    overlay=False,
    control=True
).add_to(m)

# County + State Boundaries (US Census)
folium.GeoJson(
    "https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json",
    name="U.S. Counties",
    style_function=lambda x: {
        "color": "#ffffff",
        "weight": 0.5,
        "fillOpacity": 0
    }
).add_to(m)

folium.GeoJson(
    "https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/us-states.json",
    name="U.S. States",
    style_function=lambda x: {
        "color": "#ffff00",
        "weight": 1.0,
        "fillOpacity": 0
    }
).add_to(m)

def wbgt_color(val):
    if pd.isna(val):
        return "gray"
    elif val < 85:
        return "#8BC34A"  # Low
    elif val < 88:
        return "#FFEB3B"  # Moderate
    elif val < 90:
        return "#F44336"  # High
    else:
        return "#212121"  # Extreme

for _, row in wbgt_df.iterrows():
    if pd.notna(row["Latitude"]) and pd.notna(row["Longitude"]):
        color = wbgt_color(row["WBGT (°F)"])
        folium.CircleMarker(
            location=[row["Latitude"], row["Longitude"]],
            radius=7,
            color=color,
            fill=True,
            fill_color=color,
            fill_opacity=0.9,
            popup=folium.Popup(
                f"<b>{row['Station']}</b><br>"
                f"WBGT: {row['WBGT (°F)']:.1f} °F<br>"
                f"Observed: {row['Time']}",
                max_width=250
            ),
        ).add_to(m)

# --- Add Horizontal Color Bar ---
legend_html = """
<div style="
    position: fixed; 
    bottom: 25px; left: 50%; transform: translateX(-50%);
    width: 400px; height: 45px; 
    background: linear-gradient(to right, #8BC34A 25%, #FFEB3B 25%, #FFEB3B 50%, #F44336 50%, #F44336 75%, #212121 75%);
    border: 2px solid grey;
    border-radius: 6px;
    z-index:9999; font-size:14px; color:white; text-align:center;">
    <div style="padding-top:6px;">
        <b>WBGT (°F)</b> — 
        <span style="color:#8BC34A;">80–85 Low</span> |
        <span style="color:#FFEB3B;">85–88 Moderate</span> |
        <span style="color:#F44336;">88–90 High</span> |
        <span style="color:#FFFFFF;">>90 Extreme</span>
    </div>
</div>
"""
macro = MacroElement()
macro._template = Template(legend_html)
m.get_root().add_child(macro)

folium.LayerControl().add_to(m)
st_folium(m, width=1000, height=650)

# --- Table Below Map ---
st.markdown("### Station WBGT Table")
st.dataframe(
    wbgt_df.sort_values("WBGT (°F)", ascending=False)
            .style.format({"WBGT (°F)": "{:.1f}"}),
    use_container_width=True,
    hide_index=True,
)

st.caption("WBGT is a simplified estimate based on air temperature, dew point, wind speed, solar radiation, and pressure.")
