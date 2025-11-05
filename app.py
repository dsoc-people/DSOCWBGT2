# app.py
# -*- coding: utf-8 -*-
# Streamlit App — Combined Mesonet + WeatherSTEM WBGT Map

import requests
import pandas as pd
import numpy as np
import folium
from streamlit_folium import st_folium
import streamlit as st
from geopy.geocoders import Nominatim
from geopy.extra.rate_limiter import RateLimiter

# ---------------- Streamlit Setup ----------------
st.set_page_config(page_title="Kentucky WBGT Monitor", layout="wide")
st.title("🌡️ Kentucky Mesonet & WeatherSTEM — WBGT Monitor")
st.caption("Live Wet Bulb Globe Temperature (°F) across Kentucky Mesonet and White Squirrel Weather (WeatherSTEM) stations.")

# ---------------- Configuration ----------------
year = "2025"

# URLs for WeatherSTEM stations
urls = {
    "WKU": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wku/latest.json",
    "WKU Chaos": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wkuchaos/latest.json",
    "WKU IM Fields": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wkuimfields/latest.json",
    "E'town": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/hardin/wswelizabethtown/latest.json",
    "Owensboro": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/daviess/wswowensboro/latest.json",
    "Glasgow": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/barren/wswglasgow/latest.json"
}

target_sensors = [
    "Wet Bulb Globe Temperature",
    "Thermometer",
    "Dewpoint",
    "Anemometer"
]

# ---------------- Helper Functions ----------------
def extract_value(records, target):
    for r in records:
        if target.lower() in r.get("sensor_name", "").lower():
            return r.get("value")
    return None

def farenheit_to_celsius(temp_f):
    return (temp_f - 32) * 5 / 9

def celsius_to_farenheit(temp_c):
    return temp_c * 9 / 5 + 32

def dbdp2wb(tempC, dpC, p):
    return (tempC + dpC) / 2

def wbgt(tempF, mph, rad, bar, dpF):
    tempC = farenheit_to_celsius(tempF)
    mps = mph * 0.44704
    tempK = tempC + 273.15
    if rad is None or np.isnan(rad):
        tempG = np.nan
    else:
        tempG = tempK + (rad - 30) / (0.0252 * rad + 10.5 * mps + 22.5 + 1e-9)
        tempG -= 273.15
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

# ---------------- WeatherSTEM Fetch ----------------
def fetch_weatherstem():
    data = []
    for site, url in urls.items():
        try:
            r = requests.get(url, timeout=10)
            r.raise_for_status()
            j = r.json()
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
            data.append({"name": site, "observation_time": "Error", "wbgt_f": None,
                         "Temperature (°F)": None, "Dewpoint (°F)": None, "Wind Speed (mph)": None,
                         "source": "White Squirrel Weather"})
    return pd.DataFrame(data)

# ---------------- Mesonet Data Fetch ----------------
@st.cache_data(ttl=900)
def process_station_data(station_id, coords):
    try:
        manifest_url = f"https://d266k7wxhw6o23.cloudfront.net/data/{station_id}/{year}/manifest.json"
        manifest = requests.get(manifest_url).json()
        latest_day = max(manifest.keys())
        key = manifest[latest_day]["key"]
        data_url = f"https://d266k7wxhw6o23.cloudfront.net/{key}"
        data = requests.get(data_url).json()
        df = pd.DataFrame(data["rows"], columns=data["columns"])
        cols = ["TAIR", "DWPT", "WSPD", "SRAD", "PRES", "UTCTimestampCollected"]
        if not all(c in df.columns for c in cols):
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
    coords = {row["abbrev"]: (row["lat"], row["lon"]) for _, row in df.iterrows()}
    return df, coords

stations_df, station_coords = load_station_coords()
station_abbrevs = stations_df["abbrev"].tolist()

# ---------------- Build Combined Dataset ----------------
with st.spinner("Fetching latest WBGT data..."):
    ws_df = fetch_weatherstem()
    mesonet_results = [process_station_data(s, station_coords) for s in station_abbrevs]
    mesonet_results = [r for r in mesonet_results if r]
    mesonet_df = pd.DataFrame(mesonet_results)
    combined = pd.concat([mesonet_df, ws_df], ignore_index=True)

# ---------------- Add Coordinates for WeatherSTEM ----------------
known_coords = {
    "WKU": (36.9855, -86.4551),
    "WKU Chaos": (36.9855, -86.4551),
    "WKU IM Fields": (36.9809, -86.4614),
    "E'town": (37.6939, -85.8594),
    "Owensboro": (37.7719, -87.1112),
    "Glasgow": (36.9959, -85.9119),
}
geolocator = Nominatim(user_agent="wbgt-map")
geocode = RateLimiter(geolocator.geocode, min_delay_seconds=1)
for i, row in ws_df.iterrows():
    if row["name"] in known_coords:
        lat, lon = known_coords[row["name"]]
    else:
        loc = geocode(f"{row['name']}, Kentucky")
        lat, lon = (loc.latitude, loc.longitude) if loc else (None, None)
    ws_df.loc[i, "latitude"] = lat
    ws_df.loc[i, "longitude"] = lon
combined.update(ws_df)

# ---------------- Map Rendering ----------------
def wbgt_color(w):
    if w is None or pd.isna(w): return "#808080"
    elif w < 60: return "#2ca02c"
    elif w < 70: return "#ff7f0e"
    elif w < 80: return "#d62728"
    else: return "#800026"

center_lat = combined["latitude"].dropna().mean()
center_lon = combined["longitude"].dropna().mean()
m = folium.Map(location=[center_lat, center_lon], zoom_start=7, control_scale=True)

mesonet_layer = folium.FeatureGroup(name="Mesonet")
ws_layer = folium.FeatureGroup(name="White Squirrel Weather")
for _, row in combined.iterrows():
    lat, lon = row.get("latitude"), row.get("longitude")
    if pd.isna(lat) or pd.isna(lon):
        continue
    wbgt_val = row.get("wbgt_f")
    popup_html = f"""
    <div style='font-family:system-ui'>
    <b>{row['name']} ({row['source']})</b><br>
    WBGT: {wbgt_val if pd.notna(wbgt_val) else 'N/A'} °F<br>
    Obs: {row.get('observation_time','N/A')}
    </div>"""
    marker = folium.CircleMarker(
        location=[lat, lon],
        radius=7,
        color=wbgt_color(wbgt_val),
        fill=True,
        fill_opacity=0.8,
        popup=popup_html,
        tooltip=f"{row['name']}: WBGT {wbgt_val:.1f}°F" if pd.notna(wbgt_val) else f"{row['name']}: N/A"
    )
    if row["source"] == "Mesonet":
        marker.add_to(mesonet_layer)
    else:
        marker.add_to(ws_layer)

mesonet_layer.add_to(m)
ws_layer.add_to(m)
folium.LayerControl(collapsed=False).add_to(m)

legend_html = """
<div style='position: fixed; bottom: 30px; left: 30px; z-index:9999;
 background: white; padding: 10px; border-radius:8px; font-size:12px'>
 <b>WBGT (°F)</b><br>
 <div><span style='background:#2ca02c;width:12px;height:12px;display:inline-block;'></span> <60</div>
 <div><span style='background:#ff7f0e;width:12px;height:12px;display:inline-block;'></span> 60–69.9</div>
 <div><span style='background:#d62728;width:12px;height:12px;display:inline-block;'></span> 70–79.9</div>
 <div><span style='background:#800026;width:12px;height:12px;display:inline-block;'></span> ≥80</div>
 <div><span style='background:#808080;width:12px;height:12px;display:inline-block;'></span> N/A</div>
</div>
"""
m.get_root().html.add_child(folium.Element(legend_html))

st_folium(m, width=1000, height=650)

# ---------------- Data Table ----------------
st.markdown("### WBGT Table (All Stations)")
st.dataframe(
    combined[["name", "source", "wbgt_f", "observation_time"]]
    .sort_values("wbgt_f", ascending=False)
    .reset_index(drop=True)
    .style.format({"wbgt_f": "{:.1f}"})
)
