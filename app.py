# app.py
# -*- coding: utf-8 -*-
# Consolidated Weather Station Data Collection and Visualization — Streamlit Version

import requests
import pandas as pd
import numpy as np
import json
import folium
from folium.plugins import FloatImage
import streamlit as st
from streamlit_folium import st_folium
from geopy.geocoders import Nominatim
from geopy.extra.rate_limiter import RateLimiter
import math
import subprocess
import sys

# ---------------------------------------------------------------
# Streamlit Setup
# ---------------------------------------------------------------
st.set_page_config(page_title="Consolidated WBGT Dashboard", layout="wide")
st.title("Consolidated Weather Station Data & WBGT Visualization")

# ---------------------------------------------------------------
# (KEEPING this exactly as in your script)
# Installing packages dynamically
# ---------------------------------------------------------------
sys.path.append('path/to/your/module')
subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "folium", "geopy"])


# --- URLs for each WeatherSTEM station ---
urls = {
    "WKU": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wku/latest.json",
    "WKU Chaos": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wkuchaos/latest.json",
    "WKU IM Fields": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wkuimfields/latest.json",
    "E'town": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/hardin/wswelizabethtown/latest.json",
    "Owensboro": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/daviess/wswowensboro/latest.json",
    "Glasgow": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/barren/wswglasgow/latest.json",
    "Maker's Mark Warehouse": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/marion-ky/makersmarkwarehouse/latest.json",
    "Maker's Mark St Mary": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/marion-ky/makersmarkstmary/latest.json",
    "Maker's Mark Lebanon": "https://cdn/weatherstem.com/dashboard/data/dynamic/model/marion-ky/makersmarklebanon/latest.json",
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


# ---------------------------------------------------------------
# WeatherSTEM helper
# ---------------------------------------------------------------
def extract_value(records, target):
    for r in records:
        if target.lower() in r.get("sensor_name", "").lower():
            return r.get("value")
    return None


# ---------------------------------------------------------------
# Collect WeatherSTEM data
# ---------------------------------------------------------------
data = []
for site, url in urls.items():
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        site_json = response.json()

        records = site_json.get("records", [])
        time_obs = site_json.get("time", "N/A")

        wbgt = extract_value(records, "Wet Bulb Globe Temperature")
        temp = extract_value(records, "Thermometer")
        dewpt = extract_value(records, "Dewpoint")
        wind = extract_value(records, "Anemometer")

        data.append({
            "Site": site,
            "Observation Time": time_obs,
            "WBGT (°F)": wbgt,
            "Temperature (°F)": temp,
            "Dewpoint (°F)": dewpt,
            "Wind Speed (mph)": wind,
            "source": "WeatherSTEM Latest"
        })

    except Exception as e:
        data.append({
            "Site": site,
            "Observation Time": "Error",
            "WBGT (°F)": None,
            "Temperature (°F)": None,
            "Dewpoint (°F)": None,
            "Wind Speed (mph)": None,
            "source": "WeatherSTEM Latest"
        })

df = pd.DataFrame(data)


# ---------------------------------------------------------------
# Display WeatherSTEM table in Streamlit
# ---------------------------------------------------------------
st.subheader("Latest WeatherSTEM Observations")
st.dataframe(df, use_container_width=True)

df_whitesquirrel = df.copy()


# ---------------------------------------------------------------
# Hard-coded station coordinates
# ---------------------------------------------------------------
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
...
WKUIMFields,36.97747,-85.91665
Owensboro,36.97747,-85.91665
WKU,36.97747,-85.91665
"""


# ---------------------------------------------------------------
# USGS IV Kentucky
# ---------------------------------------------------------------
usgs_iv_ky_url = 'https://waterservices.usgs.gov/nwis/iv/?format=json&stateCd=ky&siteStatus=active'
usgs_json = None

try:
    r = requests.get(usgs_iv_ky_url)
    r.raise_for_status()
    usgs_json = r.json()
except:
    pass


# Parse USGS data
site_data = []
if usgs_json and "value" in usgs_json and "timeSeries" in usgs_json["value"]:
    for series in usgs_json["value"]["timeSeries"]:
        site_info = series["sourceInfo"]

        site_id = site_info["siteCode"][0]["value"]
        site_name = site_info["siteName"]
        lat = site_info["geoLocation"]["geogLocation"]["latitude"]
        lon = site_info["geoLocation"]["geogLocation"]["longitude"]

        variable_info = series["variable"]
        param_name = variable_info["variableName"]

        latest = "N/A"
        if series["values"] and series["values"][0]["value"]:
            latest = series["values"][0]["value"][-1]["value"]

        unit = variable_info["unit"].get("unitAbbreviation", "N/A")

        site_data.append({
            "Site ID": site_id,
            "Site Name": site_name,
            "Latitude": lat,
            "Longitude": lon,
            "Parameter Name": param_name,
            "Latest Value": latest,
            "Unit": unit
        })

ky_sites_df = pd.DataFrame(site_data)

st.subheader("USGS River Gauges — Active Kentucky Sites")
st.dataframe(ky_sites_df, use_container_width=True)


# ---------------------------------------------------------------
# Parse station coordinates
# ---------------------------------------------------------------
station_coords = {}
for line in station_coords_text.strip().splitlines():
    parts = line.split(',')
    if len(parts) == 3:
        abbrev = parts[0].strip()
        try:
            lat = float(parts[1].strip())
            lon = float(parts[2].strip())
            station_coords[abbrev] = (lat, lon)
        except:
            pass


# ---------------------------------------------------------------
# WBGT + Mesonet processing EXACTLY AS IN YOUR SCRIPT
# ---------------------------------------------------------------

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
        tempG = tempK + (rad - 30) / (0.0252 * rad + 10.5*mps + 22.5)
        tempG -= 273.15

    dpC = farenheit_to_celsius(dpF)
    p = bar * 3.38639

    wbc = dbdp2wb(tempC, dpC, p)

    if np.isnan(wbc) or np.isnan(tempG):
        return np.nan

    return celsius_to_farenheit(0.7*wbc + 0.2*tempG + 0.1*tempC)


year = "2025"


def process_station_data(station_id, year, station_coords):
    try:
        # EXACT ORIGINAL LOGIC PRESERVED
        manifest_url = f"https://d266k7wxhw6o23.cloudfront.net/data/{station_id}/{year}/manifest.json"
        manifest = requests.get(manifest_url).json()

        latest_day = max(manifest.keys())
        key = manifest[latest_day]["key"]

        data = requests.get(f"https://d266k7wxhw6o23.cloudfront.net/{key}").json()

        df = pd.DataFrame(data["rows"], columns=data["columns"])
        cols = ["TAIR","DWPT","WSPD","SRAD","PRES","UTCTimestampCollected"]

        if not all(c in df.columns for c in cols):
            return None

        tair = df["TAIR"].dropna().iloc[-1]
        dwpt = df["DWPT"].dropna().iloc[-1]
        wspd = df["WSPD"].dropna().iloc[-1]
        srad = df["SRAD"].dropna().iloc[-1]
        pres = df["PRES"].dropna().iloc[-1]
        timestamp = df["UTCTimestampCollected"].dropna().iloc[-1]

        wbgt_f = wbgt(
            celsius_to_farenheit(tair),
            wspd*2.23694,
            srad,
            pres*0.02953,
            celsius_to_farenheit(dwpt)
        )

        lat,lon = station_coords.get(station_id,(None,None))

        return {
            "name": station_id,
            "latitude": lat,
            "longitude": lon,
            "wbgt_f": wbgt_f,
            "observation_time": timestamp,
            "source":"Mesonet"
        }
    except:
        return None


# Run Mesonet WBGT
mesonet_results = []
for sid in station_coords.keys():
    result = process_station_data(sid, year, station_coords)
    if result:
        mesonet_results.append(result)

df_mesonet = pd.DataFrame(mesonet_results)


# ---------------------------------------------------------------
# Combine Mesonet + WeatherSTEM
# ---------------------------------------------------------------
df_whitesquirrel = df_whitesquirrel.rename(columns={"WBGT (°F)":"wbgt_f", "Site":"name"})

combined_df = pd.concat([df_mesonet, df_whitesquirrel], ignore_index=True)


# ---------------------------------------------------------------
# Folium Map
# ---------------------------------------------------------------
st.subheader("Combined WBGT Map")

m = folium.Map(location=[37.5,-85.5], zoom_start=7, control_scale=True)

mesonet_layer = folium.FeatureGroup("Mesonet").add_to(m)
ws_layer = folium.FeatureGroup("WeatherSTEM").add_to(m)
usgs_layer = folium.FeatureGroup("USGS").add_to(m)

def wbgt_color(w):
    if w is None or pd.isna(w):
        return "gray"
    if w < 60:
        return "green"
    elif w < 70:
        return "orange"
    elif w < 80:
        return "red"
    else:
        return "black"

# Add Mesonet + WeatherSTEM markers
for _,row in combined_df.iterrows():
    lat=row["latitude"]
    lon=row["longitude"]
    if pd.isna(lat) or pd.isna(lon):
        continue

    w=row["wbgt_f"]
    src=row["source"]
    popup=f"{row['name']}<br>WBGT: {w}"

    circle=folium.CircleMarker(
        location=[lat,lon],
        radius=7,
        color=wbgt_color(w),
        fill=True,
        fill_color=wbgt_color(w),
        popup=popup
    )

    if src=="Mesonet":
        circle.add_to(mesonet_layer)
    else:
        circle.add_to(ws_layer)


# Add USGS markers
for _,row in ky_sites_df.iterrows():
    folium.Marker(
        location=[row["Latitude"],row["Longitude"]],
        popup=f"{row['Site Name']}<br>{row['Parameter Name']}: {row['Latest Value']} {row['Unit']}",
        icon=folium.Icon(color="blue")
    ).add_to(usgs_layer)

folium.LayerControl().add_to(m)

st_data = st_folium(m, width=1200, height=700)
