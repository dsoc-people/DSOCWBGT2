# app.py
# -*- coding: utf-8 -*-
# Kentucky WBGT + WeatherSTEM + Mesonet + USGS (no county view, source toggles)

import requests
import pandas as pd
import numpy as np
import folium
from streamlit_folium import st_folium
import streamlit as st
import inspect
from streamlit_autorefresh import st_autorefresh
import branca.colormap as cm

# -------------------------------------------------------------
#  STATION COORDINATES (YOUR ORIGINAL LIST)
# -------------------------------------------------------------
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
Makers Mark Warehouse,37.6333457845,-85.4075842212
Makers Mark St Mary,37.5707524233,-85.3743790708
Makers Mark Lebanon,37.5758692691,-85.2736659636
Maker's Mark Innovation Garden,37.64686,-85.34895
JimBeam Booker Noe,37.8127589004,-85.6849316392
JimBeam Bardstown,37.8344634433,-85.4711423977
JimBeam Clermont,37.9317945798,-85.6520369416
JimBeam Old Crow,38.1463823354,-84.8415031586
JimBeam Grand Dad,38.215725282,-84.8093261477
Woodford County Courthouse,38.052717,-84.73067
Clinton County High School,36.708211,-85.131276
Novelis Guthrie,36.6025431022,-87.7186136559
"""

# -------------------------------------------------------------
# NAME NORMALIZATION (Fixes Mesonet vs WeatherSTEM conflicts)
# -------------------------------------------------------------
name_variants = {
    "WKUChaos": "WKU Chaos",
    "WKUCHAOS": "WKU Chaos",
    "Etown": "E'town",
    "WKUIMFields": "WKU IM Fields",
    "Owensboro": "Owensboro",
    "Glasgow": "Glasgow",
    "WKU": "WKU",
}

name_variants.update({
    "Makers Mark Warehouse": "Maker's Mark Warehouse",
    "Makers Mark St Mary": "Maker's Mark St Mary",
    "Makers Mark Lebanon": "Maker's Mark Lebanon",

    "JimBeam Grand Dad": "Jim Beam Grand Dad",
    "JimBeam Old Crow": "Jim Beam Old Crow",
    "JimBeam Booker Noe": "Jim Beam Booker Noe",
    "JimBeam Bardstown": "Jim Beam Bardstown",
    "JimBeam Clermont": "Jim Beam Clermont",
})

# -------------------------------------------------------------
# WEATHERSTEM URL TABLE
# -------------------------------------------------------------
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

    "Woodford County Courthouse": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/woodford/courthouse/latest.json",
    "Clinton County High School": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/clinton/clintonhs/latest.json",
    "Novelis Guthrie": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/todd/novelis/latest.json",
}

# -------------------------------------------------------------
# STREAMLIT SETUP
# -------------------------------------------------------------
st.set_page_config(page_title="Kentucky WBGT Monitor", layout="wide")
st.title("🌡️ Kentucky WBGT / Weather Map Dashboard")

sig = inspect.signature(st_autorefresh)
if "rerun" in sig.parameters:
    st_autorefresh(interval=300000, key="wbgt_refresh", rerun=False)
else:
    st_autorefresh(interval=300000, key="wbgt_refresh")

year = "2025"

# -------------------------------------------------------------
# TOGGLES FOR DATA SOURCES
# -------------------------------------------------------------
with st.sidebar:
    st.header("Data Sources")
    show_weatherstem = st.checkbox("WeatherSTEM / White Squirrel Weather", True)
    show_mesonet    = st.checkbox("Mesonet", True)
    show_usgs       = st.checkbox("USGS River Gauges", True)

selected_var = st.sidebar.selectbox(
    "Variable to Display:",
    ["WBGT (°F)"]
)

# -------------------------------------------------------------
# LOAD STATION COORDINATES
# -------------------------------------------------------------
def load_station_coords():
    station_data = []
    for line in station_coords_text.strip().split("\n"):
        parts = line.split(",")
        if len(parts) == 3:
            station_data.append({
                "abbrev": parts[0],
                "lat": float(parts[1]),
                "lon": float(parts[2])
            })
    df = pd.DataFrame(station_data)
    coords = {r["abbrev"]: (r["lat"], r["lon"]) for _, r in df.iterrows()}
    return df, coords

stations_df, station_coords = load_station_coords()

# -------------------------------------------------------------
# WEATHERSTEM FETCHER
# -------------------------------------------------------------
def extract_value(records, target):
    for r in records:
        if target.lower() in r.get("sensor_name", "").lower():
            return r.get("value")
    return None

@st.cache_data(ttl=300)
def fetch_weatherstem_data():
    out = []
    for name, url in urls.items():
        try:
            j = requests.get(url, timeout=10).json()
            records = j.get("records", [])
            out.append({
                "name": name,
                "observation_time": j.get("time", "N/A"),
                "WBGT (°F)": extract_value(records, "Wet Bulb Globe Temperature"),
                "Temperature (°F)": extract_value(records, "Thermometer"),
                "Dewpoint (°F)": extract_value(records, "Dewpoint"),
                "Wind Speed (mph)": extract_value(records, "Anemometer"),
                "source": "WeatherSTEM"
            })
        except:
            out.append({
                "name": name,
                "observation_time": "Error",
                "WBGT (°F)": None,
                "Temperature (°F)": None,
                "Dewpoint (°F)": None,
                "Wind Speed (mph)": None,
                "source": "WeatherSTEM",
            })
    return pd.DataFrame(out)

# -------------------------------------------------------------
# MESONET FETCHING
# -------------------------------------------------------------
def farenheit_to_celsius(f): return (f - 32) * 5 / 9
def celsius_to_farenheit(c): return c * 9 / 5 + 32

def dbdp2wb(tc, dc, p): return (tc + dc) / 2

def wbgt_calc(tempF, mph, rad, bar, dpF):
    tempC = farenheit_to_celsius(tempF)
    mps = mph * 0.44704
    tempK = tempC + 273.15
    tempG = (
        np.nan if rad is None else tempK + (rad - 30) / (0.0252 * rad + 10.5 * mps + 22.5 + 1e-9) - 273.15
    )
    dpC = farenheit_to_celsius(dpF)
    p = bar * 3.38639
    wbc = dbdp2wb(tempC, dpC, p)
    wbgt_c = 0.7*wbc + 0.2*tempG + 0.1*tempC
    return celsius_to_farenheit(wbgt_c)

@st.cache_data(ttl=300)
def process_mesonet_station(st_id, coords):
    lat, lon = coords.get(st_id, (None, None))
    try:
        murl = f"https://d266k7wxhw6o23.cloudfront.net/data/{st_id}/{year}/manifest.json"
        manifest = requests.get(murl, timeout=10).json()
        latest_day = max(manifest.keys())
        key = manifest[latest_day]["key"]
        data = requests.get(f"https://d266k7wxhw6o23.cloudfront.net/{key}", timeout=10).json()
        df = pd.DataFrame(data["rows"], columns=data["columns"])
        tair, dew, wspd, srad, pres = df[["TAIR","DWPT","WSPD","SRAD","PRES"]].dropna().iloc[-1]
        pres_inhg = pres * 0.02953
        wb = wbgt_calc(celsius_to_farenheit(tair), wspd*2.23694, srad, pres_inhg, celsius_to_farenheit(dew))
        return {
            "name": st_id,
            "latitude": lat,
            "longitude": lon,
            "WBGT (°F)": wb,
            "Temperature (°F)": celsius_to_farenheit(tair),
            "Dewpoint (°F)": celsius_to_farenheit(dew),
            "Wind Speed (mph)": wspd*2.23694,
            "observation_time": df["UTCTimestampCollected"].dropna().iloc[-1],
            "source": "Mesonet"
        }
    except:
        return {
            "name": st_id,
            "latitude": lat,
            "longitude": lon,
            "WBGT (°F)": None,
            "Temperature (°F)": None,
            "Dewpoint (°F)": None,
            "Wind Speed (mph)": None,
            "observation_time": "Error",
            "source": "Mesonet"
        }

@st.cache_data(ttl=300)
def fetch_mesonet_data():
    mesonet_ids = []
    for abbr in stations_df["abbrev"]:
        normalized = name_variants.get(abbr, abbr)
        if normalized not in urls:     # SEPARATION FIX
            mesonet_ids.append(abbr)

    data = [process_mesonet_station(st, station_coords) for st in mesonet_ids]
    return pd.DataFrame(data)

# -------------------------------------------------------------
# USGS FETCHER (STATEWIDE)
# -------------------------------------------------------------
@st.cache_data(ttl=300)
def fetch_usgs_data():
    url = "https://waterservices.usgs.gov/nwis/iv/?format=json&stateCd=ky&siteStatus=active"
    try:
        js = requests.get(url, timeout=20).json()
        ts = js["value"]["timeSeries"]
        rows = []
        for s in ts:
            info = s["sourceInfo"]
            var = s["variable"]
            site_id = info["siteCode"][0]["value"]
            site_name = info["siteName"]
            lat = info["geoLocation"]["geogLocation"]["latitude"]
            lon = info["geoLocation"]["geogLocation"]["longitude"]
            latest = s["values"][0]["value"][-1]
            rows.append({
                "name": site_name,
                "Site ID": site_id,
                "latitude": lat,
                "longitude": lon,
                "Parameter Name": var["variableName"],
                "Latest Value": latest["value"],
                "Value Unit": var["unit"]["unitCode"],
                "source": "USGS"
            })
        return pd.DataFrame(rows)
    except:
        return pd.DataFrame(columns=[
            "name","Site ID","latitude","longitude",
            "Parameter Name","Latest Value","Value Unit","source"
        ])

# -------------------------------------------------------------
# FETCH SOURCES
# -------------------------------------------------------------
ws_df  = fetch_weatherstem_data() if show_weatherstem else pd.DataFrame()
mes_df = fetch_mesonet_data()     if show_mesonet    else pd.DataFrame()
usgs_df = fetch_usgs_data()       if show_usgs       else pd.DataFrame()

# -------------------------------------------------------------
# APPLY NAME NORMALIZATION + COORDINATE FIXES (WeatherSTEM)
# -------------------------------------------------------------
# Build lookup from station list
coord_map = {}
for line in station_coords_text.strip().split("\n"):
    parts = line.split(",")
    if len(parts)==3:
        name = name_variants.get(parts[0], parts[0])
        coord_map[name] = (float(parts[1]), float(parts[2]))

if not ws_df.empty:
    for i, row in ws_df.iterrows():
        raw = row["name"]
        normalized = name_variants.get(raw, raw)
        ws_df.loc[i, "name"] = normalized
        lat, lon = coord_map.get(normalized, (None, None))
        ws_df.loc[i, "latitude"] = lat
        ws_df.loc[i, "longitude"] = lon

# -------------------------------------------------------------
# COMBINE ALL SOURCES
# -------------------------------------------------------------
combined = pd.concat([df for df in [ws_df, mes_df] if not df.empty], ignore_index=True)

# -------------------------------------------------------------
# CREATE MAP
# -------------------------------------------------------------
# Determine map center
all_lats = pd.to_numeric(combined.get("latitude", pd.Series()), errors="coerce")
all_lons = pd.to_numeric(combined.get("longitude", pd.Series()), errors="coerce")
if show_usgs and not usgs_df.empty:
    all_lats = pd.concat([all_lats, pd.to_numeric(usgs_df["latitude"], errors="coerce")])
    all_lons = pd.concat([all_lons, pd.to_numeric(usgs_df["longitude"], errors="coerce")])

center_lat = all_lats.mean() if len(all_lats.dropna()) else 37.5
center_lon = all_lons.mean() if len(all_lons.dropna()) else -85.3

m = folium.Map(location=[center_lat, center_lon], zoom_start=7, control_scale=True)

# Coloring for WBGT / Temp / Dewpoint / Wind
def variable_color(val, var):
    if pd.isna(val): return "#808080"
    if var == "WBGT (°F)":
        if val < 66: return "#008000"
        elif val < 74: return "#FEF200"
        elif val < 83: return "#FF0000"
        else: return "#000000"
    if var in ["Temperature (°F)", "Dewpoint (°F)"]:
        return cm.LinearColormap(["#0000FF", "#00FF00", "#FF0000"], 30, 100)(val)
    if var == "Wind Speed (mph)":
        return cm.LinearColormap(["#FFFFFF", "#00FFFF", "#0000FF"], 0, 20)(val)
    return "#808080"

# -------------------------------------------------------------
# WeatherSTEM + Mesonet markers
# -------------------------------------------------------------
for _, row in combined.iterrows():
    lat = row.get("latitude")
    lon = row.get("longitude")
    if pd.isna(lat) or pd.isna(lon): 
        continue
    val = row.get(selected_var)
    color = variable_color(val, selected_var)
    popup = f"<b>{row['name']} ({row['source']})</b><br>{selected_var}: {val}"
    folium.CircleMarker(
        location=[lat, lon],
        radius=7,
        color="black" if row["source"]=="WeatherSTEM" else color,
        fill=True,
        fill_color=color,
        fill_opacity=0.85,
        popup=popup
    ).add_to(m)

# -------------------------------------------------------------
# USGS markers
# -------------------------------------------------------------
if show_usgs and not usgs_df.empty:
    for _, r in usgs_df.iterrows():
        lat, lon = r["latitude"], r["longitude"]
        if pd.isna(lat) or pd.isna(lon): continue
        popup = (
            f"<b>{r['name']}</b><br>"
            f"Site ID: {r['Site ID']}<br>"
            f"{r['Parameter Name']}: {r['Latest Value']} {r['Value Unit']}"
        )
        folium.Marker(
            location=[lat, lon],
            popup=popup,
            icon=folium.Icon(color="blue", icon="tint", prefix="fa")
        ).add_to(m)

# -------------------------------------------------------------
# DISPLAY MAP
# -------------------------------------------------------------
st_folium(m, width=1000, height=650)
