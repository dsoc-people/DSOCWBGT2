# app.py
# -*- coding: utf-8 -*-
# Streamlit App — Combined Mesonet + WeatherSTEM WBGT Map (auto-refresh every 5 min)

import requests
import pandas as pd
import numpy as np
import folium
from streamlit_folium import st_folium
import streamlit as st
from streamlit_autorefresh import st_autorefresh

# ---------------- Streamlit Setup ----------------
st.set_page_config(page_title="Kentucky WBGT Monitor", layout="wide")
st.title("TESTY")
st.caption("Live Wet Bulb Globe Temperature (°F) across Kentucky Mesonet and White Squirrel Weather (WeatherSTEM) stations. Auto-updates every 5 minutes.")

# 🔁 Auto-refresh every 5 minutes (300,000 ms)
st_autorefresh(interval=5 * 60 * 1000, key="wbgt_autorefresh")

year = "2025"

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

# ---------------- WeatherSTEM Fetch ----------------
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
                "wbgt_f": wbgt_val,
                "Temperature (°F)": temp,
                "Dewpoint (°F)": dew,
                "Wind Speed (mph)": wind,
                "source": "White Squirrel Weather"
            })
        except Exception:
            data.append({
                "name": site, "observation_time": "Error", "wbgt_f": None,
                "Temperature (°F)": None, "Dewpoint (°F)": None, "Wind Speed (mph)": None,
                "source": "White Squirrel Weather"
            })
    return pd.DataFrame(data)

# ---------------- Mesonet Fetch ----------------
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
        return {"name": station_id, "latitude": lat, "longitude": lon, "wbgt_f": wbgt_f,
                "observation_time": obs_time, "source": "Mesonet"}
    except Exception:
        return None

# ---------------- Station Coordinates ----------------
@st.cache_data
def load_station_coords():
    url = "https://d266k7wxhw6o23.cloudfront.net/metadata/stations_468eb55962c18d1fc333160925381b9d6fb5eb86cd6fbbfbfc285b1d6fcfe7a0.json"
    df = pd.DataFrame(requests.get(url).json())
    coords = {r["abbrev"]: (r["lat"], r["lon"]) for _, r in df.iterrows()}
    return df, coords

stations_df, station_coords = load_station_coords()

# ---------------- Build Combined Dataset ----------------
with st.spinner("Fetching latest WBGT data..."):
    ws_df = fetch_weatherstem()
    mesonet_df = pd.DataFrame([r for r in (process_station_data(s, station_coords)
                                           for s in stations_df["abbrev"].tolist()) if r])
    combined = pd.concat([mesonet_df, ws_df], ignore_index=True)

# ---------------- Add Known Coordinates for WeatherSTEM ----------------
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

# ---------------- Map Rendering ----------------
def wbgt_color(w):
    """Apply WBGT color scale."""
    if w is None or pd.isna(w):
        return "#808080"  # gray for missing data
    if 40 <= w <= 64:
        return "#008000"  # green
    elif 65 <= w <= 72:
        return "#FEF200"  # Yellow
    elif 73 <= w <= 81:
        return "#FF0000"  # red
    else:
        return "#000000"  # black for ≥82

center_lat = combined["latitude"].dropna().mean()
center_lon = combined["longitude"].dropna().mean()
m = folium.Map(location=[center_lat, center_lon], zoom_start=7, control_scale=True)

mesonet_layer = folium.FeatureGroup(name="Mesonet")
ws_layer = folium.FeatureGroup(name="White Squirrel Weather")

for _, row in combined.iterrows():
    lat, lon = row.get("latitude"), row.get("longitude")
    if pd.isna(lat) or pd.isna(lon): continue
    wbgt_val = row.get("wbgt_f")
    popup = f"<b>{row['name']} ({row['source']})</b><br>WBGT: {wbgt_val if pd.notna(wbgt_val) else 'N/A'} °F<br>Obs: {row.get('observation_time','N/A')}"
    folium.CircleMarker(
        location=[lat, lon], radius=7, color=wbgt_color(wbgt_val),
        fill=True, fill_opacity=0.8,
        popup=folium.Popup(popup, max_width=250),
        tooltip=f"{row['name']}: WBGT {wbgt_val:.1f}°F" if pd.notna(wbgt_val) else f"{row['name']}: N/A"
    ).add_to(mesonet_layer if row["source"] == "Mesonet" else ws_layer)

mesonet_layer.add_to(m)
ws_layer.add_to(m)
folium.LayerControl(collapsed=False).add_to(m)

# ---------------- Legend ----------------
legend_html = """
<div style='position: fixed; bottom: 30px; left: 30px; z-index:9999;
 background: white; padding: 10px; border-radius:8px; font-size:12px'>
 <b>WBGT (°F)</b><br>
 <div><span style='background:#008000;width:12px;height:12px;display:inline-block;'></span> 40–64 (Safe)</div>
 <div><span style='background:#808080;width:12px;height:12px;display:inline-block;'></span> 65–72 (Caution)</div>
 <div><span style='background:#FF0000;width:12px;height:12px;display:inline-block;'></span> 73–81 (Danger)</div>
 <div><span style='background:#000000;width:12px;height:12px;display:inline-block;border:1px solid #999;'></span> ≥82 (Extreme)</div>
 <div><span style='background:#808080;width:12px;height:12px;display:inline-block;'></span> N/A (No Data)</div>
</div>
"""
m.get_root().html.add_child(folium.Element(legend_html))
st_folium(m, width=1000, height=650)
