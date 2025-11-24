# -*- coding: utf-8 -*-
# Streamlit Integrated WeatherSTEM + Custom Mesonet WBGT App
# EXACT integration of user's entire WeatherSTEM/Mesonet pipeline

import streamlit as st
import requests
import pandas as pd
import numpy as np
import json
import folium
from streamlit_folium import st_folium
import math
import time

# ---------------- Streamlit Config ----------------
st.set_page_config(page_title="KY WBGT Map", layout="wide")
st.title("🌡️ Kentucky WBGT / Weather Map — Integrated WeatherSTEM + Custom Mesonet")

# ---------------- Auto-refresh every 5 min ----------------
st.experimental_set_query_params(ts=int(time.time()))
time.sleep(300)
st.experimental_rerun()


# ---------------- WeatherSTEM URLs ----------------
urls = {
    "WKU": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wku/latest.json",
    "WKU Chaos": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wkuchaos/latest.json",
    "WKU IM Fields": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/warren/wkuimfields/latest.json",
    "E'town": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/hardin/wswelizabethtown/latest.json",
    "Owensboro": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/daviess/wswowensboro/latest.json",
    "Glasgow": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/barren/wswglasgow/latest.json"
}

# ---------------- Extractor ----------------
def extract_value(records, target):
    for r in records:
        if target.lower() in r.get("sensor_name", "").lower():
            return r.get("value")
    return None


# ---------------- WeatherSTEM Fetch ----------------
def fetch_weatherstem():
    out = []
    for site, url in urls.items():
        try:
            r = requests.get(url, timeout=10)
            r.raise_for_status()
            js = r.json()

            records = js.get("records", [])
            t = js.get("time", "N/A")

            wbgt = extract_value(records, "Wet Bulb Globe Temperature")
            temp = extract_value(records, "Thermometer")
            dew = extract_value(records, "Dewpoint")
            wind = extract_value(records, "Anemometer")

            out.append({
                "name": site,
                "observation_time": t,
                "wbgt_f": wbgt,
                "temperature_f": temp,
                "dewpoint_f": dew,
                "wind_mph": wind,
                "source": "WeatherSTEM Latest",
                "latitude": None,
                "longitude": None
            })

        except Exception:
            out.append({
                "name": site,
                "observation_time": "Error",
                "wbgt_f": None,
                "temperature_f": None,
                "dewpoint_f": None,
                "wind_mph": None,
                "source": "WeatherSTEM Latest",
                "latitude": None,
                "longitude": None
            })
    return pd.DataFrame(out)


df_weatherstem = fetch_weatherstem()

# ---------------- Coordinates for WeatherSTEM ----------------
ws_coords = {
    "WKU": (36.9855, -86.4551),
    "WKU Chaos": (36.9855, -86.4551),
    "WKU IM Fields": (36.9809, -86.4614),
    "E'town": (37.6939, -85.8594),
    "Owensboro": (37.7719, -87.1112),
    "Glasgow": (36.9959, -85.9119),
}

for i in range(len(df_weatherstem)):
    name = df_weatherstem.loc[i, "name"]
    if name in ws_coords:
        df_weatherstem.loc[i, "latitude"] = ws_coords[name][0]
        df_weatherstem.loc[i, "longitude"] = ws_coords[name][1]

# ---------------- Custom Mesonet Stations ----------------
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

station_coords = {}
station_list = []
for line in station_coords_text.strip().splitlines():
    p = line.split(",")
    if len(p) == 3:
        sid = p[0].strip()
        try:
            lat = float(p[1])
            lon = float(p[2])
            station_coords[sid] = (lat, lon)
            station_list.append(sid)
        except:
            pass

# ---------------- WBGT Functions ----------------
def F_to_C(f):
    return (f - 32) * 5/9

def C_to_F(c):
    return c * 9/5 + 32

def dbdp2wb(tempC, dpC, p):
    return (tempC + dpC) / 2

def wbgt_calc(tempF, wind_mph, rad, bar_inhg, dewF):
    tempC = F_to_C(tempF)
    mps = wind_mph * 0.44704
    tempK = tempC + 273.15

    if rad is None or np.isnan(rad):
        TG = np.nan
    else:
        TG = tempK + (rad - 30)/(0.0252*rad + 10.5*mps + 22.5 + 1e-9)
        TG -= 273.15

    p = bar_inhg * 3.38639
    dpC = F_to_C(dewF)

    if np.isnan(tempC) or np.isnan(dpC) or np.isnan(p):
        wbc = np.nan
    else:
        wbc = dbdp2wb(tempC, dpC, p)

    if np.isnan(wbc) or np.isnan(TG):
        return np.nan
    else:
        wbgt_C = 0.7*wbc + 0.2*TG + 0.1*tempC
        return C_to_F(wbgt_C)

# ---------------- Fetch CloudFront Mesonet ----------------
def fetch_mesonet():
    out = []
    year = "2025"

    for sid in station_list:
        try:
            if sid == "WKUCHAOS":
                raise Exception("Skip WKUCHAOS here")

            manifest_url = f"https://d266k7wxhw6o23.cloudfront.net/data/{sid}/{year}/manifest.json"
            mfest = requests.get(manifest_url).json()
            latest_day = max(mfest.keys())
            key = mfest[latest_day]["key"]
            data_url = f"https://d266k7wxhw6o23.cloudfront.net/{key}"

            js = requests.get(data_url).json()
            df = pd.DataFrame(js["rows"], columns=js["columns"])

            needed = ["TAIR", "DWPT", "WSPD", "SRAD", "PRES", "UTCTimestampCollected"]
            if not all(c in df.columns for c in needed):
                raise Exception("Missing columns")

            idx = {col: df[col].last_valid_index() for col in needed}
            if idx["TAIR"] is None:
                raise Exception("No valid rows")

            tair_c = df.loc[idx["TAIR"], "TAIR"]
            dwpt_c = df.loc[idx["DWPT"], "DWPT"]
            wspd = df.loc[idx["WSPD"], "WSPD"]
            srad = df.loc[idx["SRAD"], "SRAD"]
            pres = df.loc[idx["PRES"], "PRES"]
            ts = df.loc[idx["UTCTimestampCollected"], "UTCTimestampCollected"]

            tair_f = C_to_F(tair_c)
            dew_f = C_to_F(dwpt_c)
            wind_mph = wspd * 2.23694
            pres_inhg = pres * 0.02953

            wbgt_f = wbgt_calc(tair_f, wind_mph, srad, pres_inhg, dew_f)

            lat, lon = station_coords.get(sid, (None, None))

            out.append({
                "name": sid,
                "latitude": lat,
                "longitude": lon,
                "wbgt_f": wbgt_f,
                "temperature_f": tair_f,
                "dewpoint_f": dew_f,
                "wind_mph": wind_mph,
                "observation_time": ts,
                "source": "Mesonet Archive"
            })

        except:
            lat, lon = station_coords.get(sid, (None, None))
            out.append({
                "name": sid,
                "latitude": lat,
                "longitude": lon,
                "wbgt_f": None,
                "temperature_f": None,
                "dewpoint_f": None,
                "wind_mph": None,
                "observation_time": "N/A",
                "source": "Mesonet Archive"
            })

    return pd.DataFrame(out)


df_mesonet = fetch_mesonet()

# ---------------- Combine Sources ----------------
df_combined = pd.concat([df_mesonet, df_weatherstem], ignore_index=True)

# ---------------- Color Rules ----------------
def wbgt_color(v):
    if v is None or pd.isna(v): return "#808080"
    if v < 60: return "#2ca02c"
    if v < 70: return "#ff7f0e"
    if v < 80: return "#d62728"
    return "#800026"

def marker_size(w):
    try:
        return max(6, min(18, 6 + 1.2*float(w)))
    except:
        return 8

# ---------------- Folium Map ----------------
center_lat = df_combined["latitude"].dropna().mean() if df_combined["latitude"].notna().sum() else 37.1
center_lon = df_combined["longitude"].dropna().mean() if df_combined["longitude"].notna().sum() else -85.9

m = folium.Map(location=[center_lat, center_lon], zoom_start=8, control_scale=True)

# Tile layers with attributions
folium.TileLayer(
    tiles="https://{s}.basemaps.cartocdn.com/rastertiles/voyager/{z}/{x}/{y}{r}.png",
    attr="© CARTO © OSM",
    name="CartoDB Voyager"
).add_to(m)

folium.TileLayer(
    tiles="https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png",
    attr="© CARTO © OSM",
    name="Dark Matter"
).add_to(m)

folium.TileLayer(
    tiles="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
    attr="Tiles © Esri",
    name="ESRI World Imagery"
).add_to(m)

folium.TileLayer(
    tiles="https://tile.opentopomap.org/{z}/{x}/{y}.png",
    attr="© OpenTopoMap © OSM",
    name="OpenTopoMap"
).add_to(m)

# Layer groups
mesonet_layer = folium.FeatureGroup(name="Mesonet")
ws_layer = folium.FeatureGroup(name="WeatherSTEM")

mesonet_layer.add_to(m)
ws_layer.add_to(m)

# Add markers
for _, row in df_combined.iterrows():
    if pd.isna(row["latitude"]) or pd.isna(row["longitude"]):
        continue

    color = wbgt_color(row["wbgt_f"])
    size = marker_size(row["wind_mph"])

    popup = folium.Popup(
        f"""
        <b>{row['name']}</b><br>
        Source: {row['source']}<br>
        Obs: {row['observation_time']}<br>
        WBGT: {row['wbgt_f']}<br>
        Temp: {row['temperature_f']}<br>
        Dewpoint: {row['dewpoint_f']}<br>
        Wind: {row['wind_mph']}
        """,
        max_width=300
    )

    mkr = folium.CircleMarker(
        location=[row["latitude"], row["longitude"]],
        radius=size,
        color=color,
        fill=True,
        fill_color=color,
        fill_opacity=0.8,
        popup=popup
    )

    if "Mesonet" in row["source"]:
        mkr.add_to(mesonet_layer)
    else:
        mkr.add_to(ws_layer)

folium.LayerControl().add_to(m)

# ---------------- Render to Streamlit ----------------
st_folium(m, width=1400, height=800)

