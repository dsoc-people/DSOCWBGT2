# app.py
# -*- coding: utf-8 -*-
# Consolidated Weather Station Data Collection and Visualization — Streamlit MAP ONLY version

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
st.title("Consolidated WBGT Weather Map")

# (KEEPING your original dynamic installer — but Streamlit Cloud ignores it.)
sys.path.append('path/to/your/module')
subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "folium", "geopy"])

# ---------------------------------------------------------------
# WeatherSTEM URLs
# ---------------------------------------------------------------
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
    "Jim Beam Bardstown": "https://cdn/weatherstem.com/dashboard/data/dynamic/model/nelson/jbbardstown/latest.json",
    "Jim Beam Clermont": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/bullitt/jbclermont/latest.json",
    "Jim Beam Old Crow": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/franklin-ky/jboldcrow/latest.json",
    "Jim Beam Grand Dad": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/franklin-ky/jbgranddad/latest.json",
    "Woodford Courthouse": "https://cdn/weatherstem.com/dashboard/data/dynamic/model/woodford/courthouse/latest.json",
    "Adair County High School": "https://cdn/weatherstem.com/dashboard/data/dynamic/model/adair/achs/latest.json",
    "Clinton County High School": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/clinton/clintonhs/latest.json",
    "Novelis Guthrie": "https://cdn/weatherstem.com/dashboard/data/dynamic/model/todd/novelis/latest.json"
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
ws_list = []
for site, url in urls.items():
    try:
        resp = requests.get(url, timeout=10).json()
        rec = resp.get("records", [])
        time_obs = resp.get("time", "N/A")

        ws_list.append({
            "name": site,
            "wbgt_f": extract_value(rec, "Wet Bulb Globe Temperature"),
            "temp_f": extract_value(rec, "Thermometer"),
            "dew_f": extract_value(rec, "Dewpoint"),
            "wind_mph": extract_value(rec, "Anemometer"),
            "observation_time": time_obs,
            "source": "WeatherSTEM",
            "latitude": None,
            "longitude": None
        })

    except:
        ws_list.append({
            "name": site,
            "wbgt_f": None,
            "temp_f": None,
            "dew_f": None,
            "wind_mph": None,
            "observation_time": "Error",
            "source": "WeatherSTEM",
            "latitude": None,
            "longitude": None
        })

df_ws = pd.DataFrame(ws_list)

# ---------------------------------------------------------------
# Hard-coded Mesonet coords
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
for line in station_coords_text.strip().splitlines():
    sid, la, lo = line.split(',')
    station_coords[sid] = (float(la), float(lo))

# ---------------------------------------------------------------
# USGS IV Kentucky
# ---------------------------------------------------------------
usgs_url = "https://waterservices.usgs.gov/nwis/iv/?format=json&stateCd=ky&siteStatus=active"
try:
    usgs_json = requests.get(usgs_url).json()
except:
    usgs_json = {"value":{"timeSeries":[]}}

usgs_list = []
for ts in usgs_json["value"]["timeSeries"]:
    si = ts["sourceInfo"]
    lat = si["geoLocation"]["geogLocation"]["latitude"]
    lon = si["geoLocation"]["geogLocation"]["longitude"]
    name = si["siteName"]
    pid = si["siteCode"][0]["value"]

    var = ts["variable"]["variableName"]
    val = None
    if ts["values"] and ts["values"][0]["value"]:
        val = ts["values"][0]["value"][-1]["value"]

    usgs_list.append({
        "site_id": pid,
        "site_name": name,
        "lat": lat,
        "lon": lon,
        "param": var,
        "val": val
    })

df_usgs = pd.DataFrame(usgs_list)

# ---------------------------------------------------------------
# WBGT + Mesonet calculation (your logic kept verbatim)
# ---------------------------------------------------------------
def farenheit_to_celsius(t): return (t-32)*5/9
def celsius_to_farenheit(t): return t*9/5+32
def dbdp2wb(tc, dc, p): return (tc+dc)/2

def wbgt_calc(tf, mph, rad, bar, dpf):
    tc = farenheit_to_celsius(tf)
    dp = farenheit_to_celsius(dpf)
    mps = mph*0.44704
    tg = tc + (rad-30)/(0.0252*rad + 10.5*mps + 22.5) if rad else np.nan
    w = dbdp2wb(tc, dp, bar*3.38639)
    if np.isnan(tg) or np.isnan(w):
        return np.nan
    return celsius_to_farenheit(0.7*w + 0.2*tg + 0.1*tc)

year="2025"

def process_station_data(sid):
    try:
        murl=f"https://d266k7wxhw6o23.cloudfront.net/data/{sid}/{year}/manifest.json"
        mani=requests.get(murl).json()
        day=max(mani.keys())
        key=mani[day]["key"]
        data=requests.get(f"https://d266k7wxhw6o23.cloudfront.net/{key}").json()

        df=pd.DataFrame(data["rows"], columns=data["columns"])
        tair=df["TAIR"].dropna().iloc[-1]
        dw=df["DWPT"].dropna().iloc[-1]
        wsp=df["WSPD"].dropna().iloc[-1]
        sr=df["SRAD"].dropna().iloc[-1]
        pr=df["PRES"].dropna().iloc[-1]
        ts=df["UTCTimestampCollected"].dropna().iloc[-1]

        wb=wbgt_calc(celsius_to_farenheit(tair), wsp*2.23694, sr, pr*0.02953, celsius_to_farenheit(dw))

        lat,lon=station_coords.get(sid,(None,None))

        return {"name":sid,"latitude":lat,"longitude":lon,"wbgt_f":wb,"observation_time":ts,"source":"Mesonet"}
    except:
        return None

mesonet_list=[]
for sid in station_coords.keys():
    res=process_station_data(sid)
    if res:
        mesonet_list.append(res)

df_mesonet=pd.DataFrame(mesonet_list)

# ---------------------------------------------------------------
# Combine
# ---------------------------------------------------------------
df_ws = df_ws.rename(columns={"temp_f":"Temperature (°F)","dew_f":"Dewpoint (°F)","wind_mph":"Wind Speed (mph)"})
combined = pd.concat([df_mesonet, df_ws], ignore_index=True)

# ---------------------------------------------------------------
# Render MAP ONLY
# ---------------------------------------------------------------
st.subheader("WBGT Weather Map")

m = folium.Map(location=[37.4,-85.5], zoom_start=7, control_scale=True)

meso_layer = folium.FeatureGroup("Mesonet").add_to(m)
ws_layer = folium.FeatureGroup("WeatherSTEM").add_to(m)
usgs_layer = folium.FeatureGroup("USGS").add_to(m)

def wbgt_color(w):
    if w is None or pd.isna(w): return "gray"
    if w < 60: return "green"
    if w < 70: return "orange"
    if w < 80: return "red"
    return "black"

# WeatherSTEM + Mesonet
for _,r in combined.iterrows():
    lat=r["latitude"]; lon=r["longitude"]
    if pd.isna(lat) or pd.isna(lon): continue

    popup=f"{r['name']}<br>WBGT: {r['wbgt_f']}"
    marker = folium.CircleMarker(
        location=[lat,lon],
        radius=7,
        color=wbgt_color(r["wbgt_f"]),
        fill=True,
        fill_color=wbgt_color(r["wbgt_f"]),
        popup=popup
    )
    if r["source"]=="Mesonet":
        marker.add_to(meso_layer)
    else:
        marker.add_to(ws_layer)

# USGS markers
for _,r in df_usgs.iterrows():
    folium.Marker(
        location=[r["lat"],r["lon"]],
        popup=f"{r['site_name']}<br>{r['param']}: {r['val']}",
        icon=folium.Icon(color="blue")
    ).add_to(usgs_layer)

folium.LayerControl().add_to(m)

st_folium(m, width=1200, height=700)
