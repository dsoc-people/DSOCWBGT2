# app.py
# -*- coding: utf-8 -*-
# Streamlit App — Consolidated WBGT + WeatherSTEM + Mesonet + River Gauge Map

import streamlit as st
from streamlit_folium import st_folium
from streamlit_autorefresh import st_autorefresh
import requests
import pandas as pd
import numpy as np
import folium
import json
import subprocess
import sys
import math
from geopy.geocoders import Nominatim
from geopy.extra.rate_limiter import RateLimiter


# ======================================================================================
# STREAMLIT SETUP
# ======================================================================================

st.set_page_config(page_title="Consolidated WBGT Map", layout="wide")
st.title("🌡️ Consolidated WBGT + WeatherSTEM + Mesonet + River Gauge Map")

# Auto-refresh every 5 minutes (300,000 ms)
st_autorefresh(interval=300000, key="wbgt_refresh")


# Ensure all dependencies exist
subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "folium", "geopy"])


# ======================================================================================
# WEATHERSTEM LATEST.JSON URLS (HARDCODED, KEEPING EXACTLY AS DISCUSSED)
# ======================================================================================

urls = {
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


# ======================================================================================
# HELPER FUNCTIONS
# ======================================================================================

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


def wbgt_calc(tempF, mph, rad, bar, dpF):
    tempC = farenheit_to_celsius(tempF)
    mps = mph * 0.44704
    tempK = tempC + 273.15

    # Compute globe temperature
    if rad is None or np.isnan(rad):
        tempG = np.nan
    else:
        tempG = tempK + (rad - 30) / (0.0252 * rad + 10.5 * mps + 22.5)
        tempG -= 273.15

    p = bar * 3.38639
    dpC = farenheit_to_celsius(dpF)

    if any([np.isnan(tempC), np.isnan(dpC), np.isnan(p)]):
        wb = np.nan
    else:
        wb = dbdp2wb(tempC, dpC, p)

    if any([np.isnan(wb), np.isnan(tempG), np.isnan(tempC)]):
        wbgt_c = np.nan
    else:
        wbgt_c = 0.7 * wb + 0.2 * tempG + 0.1 * tempC

    return celsius_to_farenheit(wbgt_c)


# ======================================================================================
# WEATHERSTEM LIVE → WHITE SQUIRREL WEATHER LAYER
# ======================================================================================

ws_data = []

for site, url in urls.items():
    try:
        js = requests.get(url, timeout=10).json()
        records = js.get("records", [])
        time_obs = js.get("time", "N/A")

        ws_data.append({
            "name": site,
            "observation_time": time_obs,
            "wbgt_f": extract_value(records, "Wet Bulb Globe Temperature"),
            "Temperature (°F)": extract_value(records, "Thermometer"),
            "Dewpoint (°F)": extract_value(records, "Dewpoint"),
            "Wind Speed (mph)": extract_value(records, "Anemometer"),
            "source": "White Squirrel Weather",
            "latitude": None,
            "longitude": None
        })

    except:
        ws_data.append({
            "name": site,
            "observation_time": "ERROR",
            "wbgt_f": np.nan,
            "Temperature (°F)": np.nan,
            "Dewpoint (°F)": np.nan,
            "Wind Speed (mph)": np.nan,
            "source": "White Squirrel Weather",
            "latitude": None,
            "longitude": None
        })

df_whitesquirrel = pd.DataFrame(ws_data)


# ======================================================================================
# MESONET STATION LIST — PASTE YOUR FULL BLOCK HERE
# ======================================================================================

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
WKUCHAOS,36.98582726072027,-86.44967208166477
WKUChaos,36.98582632147347,-86.44968800031974
Etown,37.69563805082102,-85.88387790284976
Glasgow,36.9774781561,-85.916651431
WKUIMFields,36.9774781561,-85.9166514315
Owensboro,36.9774781561,-85.9166514315
WKU,36.9774781561,-85.9166514315

"""

station_coords = {}
station_ids = []

for line in station_coords_text.strip().splitlines():
    parts = line.split(",")
    if len(parts) == 3:
        sid, lat, lon = parts
        try:
            station_coords[sid] = (float(lat), float(lon))
            station_ids.append(sid)
        except:
            pass


# ======================================================================================
# MESONET ARCHIVE WBGT PROCESSING
# ======================================================================================

def process_station_data(station_id, year, station_coords):
    try:
        manifest_url = f"https://d266k7wxhw6o23.cloudfront.net/data/{station_id}/{year}/manifest.json"
        manifest = requests.get(manifest_url, timeout=10).json()

        latest_day = max(manifest.keys())
        latest_key = manifest[latest_day]["key"]

        data_url = f"https://d266k7wxhw6o23.cloudfront.net/{latest_key}"
        data = requests.get(data_url, timeout=10).json()

        df = pd.DataFrame(data["rows"], columns=data["columns"])

        tair = df["TAIR"].dropna().iloc[-1]
        dwpt = df["DWPT"].dropna().iloc[-1]
        wspd = df["WSPD"].dropna().iloc[-1]
        srad = df["SRAD"].dropna().iloc[-1]
        pres = df["PRES"].dropna().iloc[-1]
        ob_time = df["UTCTimestampCollected"].dropna().iloc[-1]

        tair_f = celsius_to_farenheit(tair)
        dwpt_f = celsius_to_farenheit(dwpt)
        wspd_mph = wspd * 2.23694
        pres_inhg = pres * 0.02953

        wbgt_f = wbgt_calc(tair_f, wspd_mph, srad, pres_inhg, dwpt_f)

        lat, lon = station_coords.get(station_id, (None, None))

        return {
            "name": station_id,
            "observation_time": ob_time,
            "wbgt_f": wbgt_f,
            "Temperature (°F)": tair_f,
            "Dewpoint (°F)": dwpt_f,
            "Wind Speed (mph)": wspd_mph,
            "latitude": lat,
            "longitude": lon,
            "source": "Mesonet"
        }

    except:
        return None


mesonet_data = []
year = "2025"

for sid in station_ids:
    r = process_station_data(sid, year, station_coords)
    if r:
        mesonet_data.append(r)

df_mesonet = pd.DataFrame(mesonet_data)


# ======================================================================================
# USGS RIVER DATA (GAUGE HEIGHT, ETC.)
# ======================================================================================

river_url = "https://waterservices.usgs.gov/nwis/iv/?format=json&stateCd=ky&siteStatus=active"

river_rows = []

try:
    js = requests.get(river_url, timeout=10).json()
    ts_list = js.get("value", {}).get("timeSeries", [])

    for ts in ts_list:
        info = ts.get("sourceInfo", {})
        geo = info.get("geoLocation", {}).get("geogLocation", {})

        lat = geo.get("latitude")
        lon = geo.get("longitude")
        name = info.get("siteName", "Unknown")
        code = info.get("siteCode", [{}])[0].get("value", "N/A")

        vals = ts.get("values", [{}])[0].get("value", [])
        latest = vals[-1] if vals else {"value": None, "dateTime": None}

        river_rows.append({
            "name": name,
            "observation_time": latest.get("dateTime"),
            "reading": latest.get("value"),
            "wbgt_f": np.nan,
            "Temperature (°F)": np.nan,
            "Dewpoint (°F)": np.nan,
            "Wind Speed (mph)": np.nan,
            "latitude": lat,
            "longitude": lon,
            "source": "River Gauge"
        })

except:
    pass

df_river = pd.DataFrame(river_rows)


# ======================================================================================
# COMBINE ALL DATA SOURCES
# ======================================================================================

combined_df = pd.concat([df_whitesquirrel, df_mesonet, df_river], ignore_index=True)


# ======================================================================================
# COLOR + RADIUS FUNCTIONS
# ======================================================================================

def wbgt_color(w):
    if pd.isna(w): return "#808080"
    if w < 60: return "#2ca02c"
    if w < 70: return "#ff7f0e"
    if w < 80: return "#d62728"
    return "#800026"

def marker_radius(wind):
    try:
        return max(6, min(18, 6 + 1.2 * float(wind)))
    except:
        return 8


# ======================================================================================
# FOLIUM MAP
# ======================================================================================

valid = combined_df.dropna(subset=["latitude", "longitude"])
if valid.empty:
    center = [37.1, -85.9]
else:
    center = [valid["latitude"].mean(), valid["longitude"].mean()]

m = folium.Map(location=center, zoom_start=8, control_scale=True)

# Base layers
folium.TileLayer("cartodbpositron").add_to(m)
folium.TileLayer("cartodbdark_matter").add_to(m)
folium.TileLayer("openstreetmap").add_to(m)

# Data layers
layer_ws = folium.FeatureGroup(name="White Squirrel Weather").add_to(m)
layer_mesonet = folium.FeatureGroup(name="Mesonet").add_to(m)
layer_river = folium.FeatureGroup(name="USGS River Gauges").add_to(m)

# Plot markers
for _, r in combined_df.iterrows():
    lat, lon = r["latitude"], r["longitude"]
    if pd.isna(lat) or pd.isna(lon):
        continue

    color = wbgt_color(r["wbgt_f"])
    radius = marker_radius(r["Wind Speed (mph)"])

    popup_html = f"""
    <b>{r['name']}</b><br>
    Source: {r['source']}<br>
    Observed: {r['observation_time']}<br>
    WBGT: {r['wbgt_f']}<br>
    Temp: {r['Temperature (°F)']}<br>
    Dewpoint: {r['Dewpoint (°F)']}<br>
    Wind: {r['Wind Speed (mph)']} mph<br>
    Reading: {r.get('reading', 'N/A')}
    """

    marker = folium.CircleMarker(
        location=[lat, lon],
        radius=radius,
        color=color,
        fill=True,
        fill_color=color,
        fill_opacity=0.85,
        weight=1,
        popup=folium.Popup(popup_html, max_width=280)
    )

    if r["source"] == "Mesonet":
        marker.add_to(layer_mesonet)
    elif r["source"] == "White Squirrel Weather":
        marker.add_to(layer_ws)
    else:
        marker.add_to(layer_river)

folium.LayerControl().add_to(m)


# ======================================================================================
# DISPLAY MAP
# ======================================================================================

st_folium(m, width=1400, height=800)
