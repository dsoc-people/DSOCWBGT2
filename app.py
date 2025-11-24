# app.py
# -*- coding: utf-8 -*-
# Unified Streamlit App — WeatherSTEM + Mesonet + USGS + County Focus

import streamlit as st
import requests
import pandas as pd
import numpy as np
import json
import folium
from streamlit_folium import st_folium
from shapely.geometry import Point
import geopandas as gpd
import branca.colormap as cm
import inspect
from streamlit_autorefresh import st_autorefresh

# ------------------------------------------------------
# Streamlit Setup
# ------------------------------------------------------
st.set_page_config(page_title="KY Weather + WBGT + River Gauges", layout="wide")
st.title("🌡️ Kentucky Weather / WBGT / River Gauge Dashboard")

# Autorefresh (5 min)
sig = inspect.signature(st_autorefresh)
if "rerun" in sig.parameters:
    refresh_counter = st_autorefresh(interval=5 * 60 * 1000, limit=None, key="auto_refresh", rerun=False)
else:
    refresh_counter = st_autorefresh(interval=5 * 60 * 1000, limit=None, key="auto_refresh")

# Map persistence
if "last_map" not in st.session_state:
    st.session_state["last_map"] = None

year = "2025"

# ------------------------------------------------------
# Sidebar Controls
# ------------------------------------------------------
st.sidebar.header("Map Controls")

selected_var = st.sidebar.selectbox(
    "Variable to Display:",
    ["WBGT (°F)", "Temperature (°F)", "Dewpoint (°F)", "Wind Speed (mph)"]
)

# ------------------------------------------------------
# WeatherSTEM URLs
# ------------------------------------------------------
weatherstem_urls = {
    "WKU": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wku/latest.json",
    "WKU Chaos": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wkuchaos/latest.json",
    "WKU IM Fields": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wkuimfields/latest.json",
    "E'town": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/hardin/wswelizabethtown/latest.json",
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
    "Woodford Courthouse": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/woodford/courthouse/latest.json",
    "Adair County High School": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/adair/achs/latest.json",
    "Clinton County High School": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/clinton/clintonhs/latest.json",
    "Novelis Guthrie": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/todd/novelis/latest.json"
}

# ------------------------------------------------------
# Sensor extraction helper
# ------------------------------------------------------
def extract_value(records, target):
    for r in records:
        if target.lower() in r.get("sensor_name", "").lower():
            return r.get("value")
    return None

# ------------------------------------------------------
# Load Kentucky counties
# ------------------------------------------------------
@st.cache_data
def load_ky_counties():
    try:
        url = "https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json"
        j = requests.get(url, timeout=10).json()
        feats = [f for f in j["features"] if f["properties"]["STATE"] == "21"]
        for f in feats:
            f["properties"]["NAME"] = f["properties"]["NAME"].title()
        return gpd.GeoDataFrame.from_features(feats, crs="EPSG:4326")
    except:
        return gpd.GeoDataFrame()

counties_gdf = load_ky_counties()
county_list = sorted(counties_gdf["NAME"].unique())
selected_county = st.sidebar.selectbox("Select County:", county_list)

# ------------------------------------------------------
# Station coordinate overrides
# ------------------------------------------------------
KNOWN_COORDS = {
    "WKU": (36.9855, -86.4551),
    "WKU Chaos": (36.9855, -86.4551),
    "WKU IM Fields": (36.9809, -86.4614),
    "E'town": (37.6939, -85.8594),
    "Owensboro": (37.7719, -87.1112),
    "Glasgow": (36.9959, -85.9119),
    "Maker's Mark Warehouse": (37.6333457845, -85.4075842212),
    "Maker's Mark St Mary": (37.5707524233, -85.3743790708),
    "Maker's Mark Lebanon": (37.5758692691, -85.2736659636),
    "Maker's Mark Innovation Garden": (37.64686, -85.34895),
    "Jim Beam Booker Noe": (37.8127589004, -85.6849316392),
    "Jim Beam Bardstown": (37.8344634433, -85.4711423977),
    "Jim Beam Clermont": (37.9317945798, -85.6520369416),
    "Jim Beam Old Crow": (38.1463823354, -84.8415031586),
    "Jim Beam Grand Dad": (38.215725282, -84.8093261477),
    "Woodford Courthouse": (38.052717, -84.73067),
    "Adair County High School": (37.107667, -85.32824),
    "Clinton County High School": (36.708211, -85.131276),
    "Novelis Guthrie": (36.6025431022, -87.7186136559)
}

# ------------------------------------------------------
# WeatherSTEM fetch
# ------------------------------------------------------
@st.cache_data(ttl=300)
def fetch_weatherstem():
    out = []
    for site, url in weatherstem_urls.items():
        try:
            j = requests.get(url, timeout=10).json()
            records = j.get("records", [])
            out.append({
                "name": site,
                "WBGT (°F)": extract_value(records, "Wet Bulb Globe Temperature"),
                "Temperature (°F)": extract_value(records, "Thermometer"),
                "Dewpoint (°F)": extract_value(records, "Dewpoint"),
                "Wind Speed (mph)": extract_value(records, "Anemometer"),
                "observation_time": j.get("time", "N/A"),
                "latitude": KNOWN_COORDS.get(site, (None, None))[0],
                "longitude": KNOWN_COORDS.get(site, (None, None))[1],
                "source": "WeatherSTEM"
            })
        except:
            out.append({
                "name": site, "WBGT (°F)": None, "Temperature (°F)": None,
                "Dewpoint (°F)": None, "Wind Speed (mph)": None,
                "observation_time": "Error",
                "latitude": KNOWN_COORDS.get(site, (None, None))[0],
                "longitude": KNOWN_COORDS.get(site, (None, None))[1],
                "source": "WeatherSTEM"
            })
    return pd.DataFrame(out)

# ------------------------------------------------------
# Mesonet WBGT physics
# ------------------------------------------------------
def f2c(f): return (f - 32) * 5/9
def c2f(c): return c * 9/5 + 32
def dbdp2wb(cT, cD, p): return (cT + cD) / 2

def wbgt_calc(tempF, mph, rad, bar_inHg, dpF):
    tempC = f2c(tempF)
    dpC = f2c(dpF)
    mps = mph * 0.44704
    tempK = tempC + 273.15

    if rad is None or np.isnan(rad):
        globC = np.nan
    else:
        globC = tempK + (rad - 30) / (0.0252*rad + 10.5*mps + 22.5 + 1e-9) - 273.15

    p = bar_inHg * 3.38639
    if any(np.isnan(v) for v in [tempC, dpC, p]):
        return np.nan

    wbc = dbdp2wb(tempC, dpC, p)
    if np.isnan(wbc) or np.isnan(globC):
        return np.nan

    return c2f(0.7*wbc + 0.2*globC + 0.1*tempC)

# ------------------------------------------------------
# Load Mesonet WBGT via manifest
# ------------------------------------------------------
STATION_LIST = list(KNOWN_COORDS.keys())

@st.cache_data(ttl=300)
def fetch_mesonet():
    out = []
    for sid in STATION_LIST:
        try:
            manifest = requests.get(
                f"https://d266k7wxhw6o23.cloudfront.net/data/{sid}/{year}/manifest.json",
                timeout=10
            ).json()
            latest_day = max(manifest.keys())
            key = manifest[latest_day]["key"]
            d = requests.get(f"https://d266k7wxhw6o23.cloudfront.net/{key}", timeout=10).json()
            df = pd.DataFrame(d["rows"], columns=d["columns"])

            req = ["TAIR", "DWPT", "WSPD", "SRAD", "PRES"]
            if not all(x in df.columns for x in req):
                continue

            tair_c = df["TAIR"].dropna().iloc[-1]
            dwpt_c = df["DWPT"].dropna().iloc[-1]
            wspd_mps = df["WSPD"].dropna().iloc[-1]
            srad = df["SRAD"].dropna().iloc[-1]
            pres_hpa = df["PRES"].dropna().iloc[-1]
            obs_t = df["UTCTimestampCollected"].dropna().iloc[-1]

            tair_f = c2f(tair_c)
            dwpt_f = c2f(dwpt_c)
            wspd_mph = wspd_mps * 2.23694
            pres_inHg = pres_hpa * 0.02953

            wb = wbgt_calc(tair_f, wspd_mph, srad, pres_inHg, dwpt_f)

            out.append({
                "name": sid,
                "WBGT (°F)": wb,
                "Temperature (°F)": tair_f,
                "Dewpoint (°F)": dwpt_f,
                "Wind Speed (mph)": wspd_mph,
                "observation_time": obs_t,
                "latitude": KNOWN_COORDS.get(sid, (None, None))[0],
                "longitude": KNOWN_COORDS.get(sid, (None, None))[1],
                "source": "Mesonet"
            })

        except:
            out.append({
                "name": sid, "WBGT (°F)": None,
                "Temperature (°F)": None, "Dewpoint (°F)": None,
                "Wind Speed (mph)": None,
                "observation_time": "Error",
                "latitude": KNOWN_COORDS.get(sid, (None, None))[0],
                "longitude": KNOWN_COORDS.get(sid, (None, None))[1],
                "source": "Mesonet"
            })
    return pd.DataFrame(out)

# ------------------------------------------------------
# USGS NWIS IV fetch
# ------------------------------------------------------
@st.cache_data(ttl=300)
def fetch_usgs():
    url = "https://waterservices.usgs.gov/nwis/iv/?format=json&stateCd=ky&siteStatus=active"
    try:
        j = requests.get(url, timeout=10).json()
        ts = j.get("value", {}).get("timeSeries", [])
    except:
        return pd.DataFrame()

    out = []
    for s in ts:
        try:
            site_info = s["sourceInfo"]
            var_info = s["variable"]
            site_id = site_info["siteCode"][0]["value"]
            latest = "N/A"
            vals = s["values"][0]["value"]
            if vals:
                latest = vals[-1]["value"]

            out.append({
                "Site ID": site_id,
                "Site Name": site_info["siteName"],
                "Latitude": site_info["geoLocation"]["geogLocation"]["latitude"],
                "Longitude": site_info["geoLocation"]["geogLocation"]["longitude"],
                "Parameter Name": var_info["variableName"],
                "Latest Value": latest,
                "Value Unit": var_info.get("unit", {}).get("unitAbbreviation", ""),
            })
        except:
            continue
    return pd.DataFrame(out)

# ------------------------------------------------------
# Fetch All Data
# ------------------------------------------------------
with st.spinner("Fetching live WeatherSTEM, Mesonet, and USGS data..."):
    df_ws = fetch_weatherstem()
    df_mesonet = fetch_mesonet()
    df_usgs = fetch_usgs()

combined_df = pd.concat([df_ws, df_mesonet], ignore_index=True)

# ------------------------------------------------------
# Color Logic (from your working app)
# ------------------------------------------------------
def variable_color(value, var):
    if pd.isna(value):
        return "#808080"
    if var in ["Temperature (°F)", "Dewpoint (°F)"]:
        cmap = cm.LinearColormap(["#0000FF", "#00FF00", "#FF0000"], vmin=30, vmax=100)
        return cmap(value)
    if var == "WBGT (°F)":
        if value < 66: return "#008000"
        if value < 74: return "#FEF200"
        if value < 83: return "#FF0000"
        return "#000000"
    if var == "Wind Speed (mph)":
        cmap = cm.LinearColormap(["#FFFFFF", "#00FFFF", "#0000FF"], vmin=0, vmax=20)
        return cmap(value)
    return "#808080"

# ------------------------------------------------------
# Main Map
# ------------------------------------------------------
center_lat = combined_df["latitude"].dropna().mean()
center_lon = combined_df["longitude"].dropna().mean()
m = folium.Map(location=[center_lat, center_lon], zoom_start=7, control_scale=True)

# Tile layers
tiles = {
    "CartoDB Voyager": "https://{s}.basemaps.cartocdn.com/rastertiles/voyager/{z}/{x}/{y}{r}.png",
    "CartoDB Dark Matter": "https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png",
    "Esri World Imagery": "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
    "OpenTopoMap": "https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png",
    "CartoDB Positron": "https://{s}.basemaps.cartocdn.com/rastertiles/light_all/{z}/{x}/{y}{r}.png",
}
for name, url in tiles.items():
    folium.TileLayer(tiles=url, name=name, control=True, overlay=False).add_to(m)

mesonet_layer = folium.FeatureGroup("Mesonet")
ws_layer = folium.FeatureGroup("WeatherSTEM")
usgs_layer = folium.FeatureGroup("USGS")
mesonet_layer.add_to(m)
ws_layer.add_to(m)
usgs_layer.add_to(m)

# ------------------------------------------------------
# Plot Mesonet + WeatherSTEM
# ------------------------------------------------------
for _, row in combined_df.iterrows():
    lat, lon = row["latitude"], row["longitude"]
    if pd.isna(lat) or pd.isna(lon): continue
    val = row.get(selected_var)
    color = variable_color(val, selected_var)

    popup = (
        f"<b>{row['name']} ({row['source']})</b><br>"
        f"{selected_var}: {val if pd.notna(val) else 'N/A'}<br>"
        f"Observed: {row['observation_time']}"
    )

    marker = folium.CircleMarker(
        location=[lat, lon],
        radius=7,
        color="black" if row["source"] == "WeatherSTEM" else color,
        weight=2 if row["source"] == "WeatherSTEM" else 1,
        fill=True, fill_color=color, fill_opacity=0.85,
        popup=popup
    )
    if row["source"] == "WeatherSTEM":
        marker.add_to(ws_layer)
    else:
        marker.add_to(mesonet_layer)

# ------------------------------------------------------
# Plot USGS (Always blue)
# ------------------------------------------------------
for site in df_usgs["Site ID"].unique():
    block = df_usgs[df_usgs["Site ID"] == site]
    r = block.iloc[0]
    lat, lon = r["Latitude"], r["Longitude"]
    popup = f"<b>{r['Site Name']}</b><br>Site ID: {site}<br>"
    for _, pr in block.iterrows():
        popup += f"{pr['Parameter Name']}: {pr['Latest Value']} {pr['Value Unit']}<br>"

    folium.Marker(
        location=[lat, lon],
        tooltip=f"{r['Site Name']} ({site})",
        popup=popup,
        icon=folium.Icon(color="blue", icon="tint", prefix="fa")
    ).add_to(usgs_layer)

folium.LayerControl(collapsed=False).add_to(m)

# ------------------------------------------------------
# Main Map Display
# ------------------------------------------------------
st.session_state["last_map"] = m
st_folium(m, width=1000, height=650)

# ------------------------------------------------------
# COUNTY FOCUS VIEW
# ------------------------------------------------------
st.markdown("### 🧭 County Focus")

geom = counties_gdf[counties_gdf["NAME"] == selected_county].geometry.iloc[0]
b = geom.bounds
cmx = folium.Map(location=[(b[1]+b[3])/2, (b[0]+b[2])/2], zoom_start=9, control_scale=True)

folium.GeoJson(
    geom.__geo_interface__,
    style_function=lambda x: {"fillColor": "#ff7800", "color": "black", "fillOpacity": 0.25, "weight": 2}
).add_to(cmx)

pts = [Point(lon, lat) for lon, lat in zip(combined_df["longitude"], combined_df["latitude"])]
combined_df["in_county"] = [geom.contains(p) for p in pts]

usgs_pts = [Point(lon, lat) for lon, lat in zip(df_usgs["Longitude"], df_usgs["Latitude"])]
df_usgs["in_county"] = [geom.contains(p) for p in usgs_pts]

# Weather/Mesonet
for _, row in combined_df[combined_df["in_county"]].iterrows():
    val = row.get(selected_var)
    color = variable_color(val, selected_var)
    popup = f"<b>{row['name']} ({row['source']})</b><br>{selected_var}: {val}"
    folium.CircleMarker(
        location=[row.latitude, row.longitude],
        radius=8,
        color="black" if row["source"]=="WeatherSTEM" else color,
        weight=2 if row["source"]=="WeatherSTEM" else 1,
        fill=True, fill_color=color, fill_opacity=0.9,
        popup=popup,
    ).add_to(cmx)

# USGS inside county
for _, row in df_usgs[df_usgs["in_county"]].iterrows():
    popup = f"<b>{row['Site Name']}</b><br>{row['Parameter Name']}: {row['Latest Value']} {row['Value Unit']}"
    folium.Marker(
        location=[row["Latitude"], row["Longitude"]],
        icon=folium.Icon(color="blue", icon="tint", prefix="fa"),
        popup=popup,
    ).add_to(cmx)

st_folium(cmx, width=850, height=450)
