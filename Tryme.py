# app.py — Kentucky Mesonet + WeatherSTEM + USGS Integrated WBGT Map
# -*- coding: utf-8 -*-

import requests
import pandas as pd
import numpy as np
import folium
from shapely.geometry import Point
from streamlit_folium import st_folium
import streamlit as st
from streamlit_autorefresh import st_autorefresh
import geopandas as gpd
import branca.colormap as cm
import inspect
from datetime import datetime

# -------------------------------------------------------------------
# STREAMLIT SETUP
# -------------------------------------------------------------------
st.set_page_config(page_title="Kentucky WBGT Monitor", layout="wide")
st.title("🌡️ Kentucky WBGT / Weather Map Dashboard")

# Version-safe autorefresh
sig = inspect.signature(st_autorefresh)
if "rerun" in sig.parameters:
    st_autorefresh(interval=5 * 60 * 1000, limit=None, key="refresh", rerun=False)
else:
    st_autorefresh(interval=5 * 60 * 1000, limit=None, key="refresh")

if "last_map" not in st.session_state:
    st.session_state["last_map"] = None

year = str(datetime.utcnow().year)

# -------------------------------------------------------------------
# SIDEBAR
# -------------------------------------------------------------------
st.sidebar.header("Map Controls")
selected_var = st.sidebar.selectbox(
    "Variable to Display:",
    ["WBGT (°F)", "Temperature (°F)", "Dewpoint (°F)", "Wind Speed (mph)"]
)

show_ws = st.sidebar.checkbox("Show WeatherSTEM / WSW", value=True)
show_mesonet = st.sidebar.checkbox("Show Mesonet", value=True)
show_usgs = st.sidebar.checkbox("Show USGS", value=True)

# -------------------------------------------------------------------
# WBGT SUPPORT
# -------------------------------------------------------------------
def farenheit_to_celsius(f):
    return (f - 32) * 5 / 9

def celsius_to_farenheit(c):
    return c * 9 / 5 + 32

def dbdp2wb(tc, dc, p):
    return (tc + dc) / 2

def wbgt(tempF, mph, rad, bar, dpF):
    tempC = farenheit_to_celsius(tempF)
    dpC = farenheit_to_celsius(dpF)
    mps = mph * 0.44704
    tempK = tempC + 273.15
    p = bar * 3.38639

    tempG = np.nan if rad is None else (
        tempK + (rad - 30) / (0.0252 * rad + 10.5 * mps + 22.5 + 1e-9) - 273.15
    )
    wbc = np.nan if (np.isnan(tempC) or np.isnan(dpC) or np.isnan(p)) else dbdp2wb(tempC, dpC, p)

    wbgt_c = (
        0.7 * wbc + 0.2 * tempG + 0.1 * tempC
        if not (np.isnan(wbc) or np.isnan(tempG) or np.isnan(tempC))
        else np.nan
    )
    return celsius_to_farenheit(wbgt_c)

def extract_value(records, target):
    for r in records:
        if target.lower() in r.get("sensor_name", "").lower():
            return r.get("value")
    return None

# -------------------------------------------------------------------
# WEATHERSTEM URLs (DEFINES WEATHERSTEM STATIONS)
# -------------------------------------------------------------------
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
    "Adair County High School": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/adair/achs/latest.json",
    "Clinton County High School": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/clinton/clintonhs/latest.json",
    "Novelis Guthrie": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/todd/novelis/latest.json"
}

# -------------------------------------------------------------------
# STATION COORDINATES (ALL MESONET)
# -------------------------------------------------------------------
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

@st.cache_data
def load_station_coords():
    rows = []
    for line in station_coords_text.strip().split("\n"):
        parts = line.split(",")
        if len(parts) == 3:
            rows.append({"abbrev": parts[0], "lat": float(parts[1]), "lon": float(parts[2])})
    df = pd.DataFrame(rows)
    coords = {r["abbrev"]: (r["lat"], r["lon"]) for _, r in df.iterrows()}
    return df, coords

stations_df, station_coords = load_station_coords()

# -------------------------------------------------------------------
# FETCH WEATHERSTEM DATA
# -------------------------------------------------------------------
@st.cache_data(ttl=300)
def fetch_weatherstem():
    out = []
    for name, url in urls.items():
        try:
            j = requests.get(url, timeout=10).json()
            rec = j.get("records", [])
            out.append({
                "name": name,
                "WBGT (°F)": extract_value(rec, "Wet Bulb"),
                "Temperature (°F)": extract_value(rec, "Thermometer"),
                "Dewpoint (°F)": extract_value(rec, "Dewpoint"),
                "Wind Speed (mph)": extract_value(rec, "Anemometer"),
                "observation_time": j.get("time", None),
                "source": "WeatherSTEM"
            })
        except:
            out.append({
                "name": name,
                "WBGT (°F)": None,
                "Temperature (°F)": None,
                "Dewpoint (°F)": None,
                "Wind Speed (mph)": None,
                "observation_time": "Error",
                "source": "WeatherSTEM"
            })
    return pd.DataFrame(out)

# -------------------------------------------------------------------
# FETCH MESONET
# -------------------------------------------------------------------
@st.cache_data(ttl=300)
def process_mesonet_station(station_id):
    lat, lon = station_coords.get(station_id, (None, None))

    try:
        manifest_url = f"https://d266k7wxhw6o23.cloudfront.net/data/{station_id}/{year}/manifest.json"
        manifest = requests.get(manifest_url, timeout=10).json()
        latest_day = max(manifest.keys())
        key = manifest[latest_day]["key"]

        j = requests.get(f"https://d266k7wxhw6o23.cloudfront.net/{key}", timeout=10).json()
        df = pd.DataFrame(j["rows"], columns=j["columns"])

        tair = df["TAIR"].dropna().iloc[-1]
        dwpt = df["DWPT"].dropna().iloc[-1]
        wspd = df["WSPD"].dropna().iloc[-1]
        rad = df["SRAD"].dropna().iloc[-1]
        pres = df["PRES"].dropna().iloc[-1] * 0.02953
        obs = df["UTCTimestampCollected"].dropna().iloc[-1]

        wbgt_val = wbgt(celsius_to_farenheit(tair),
                        wspd * 2.23694,
                        rad,
                        pres,
                        celsius_to_farenheit(dwpt))

        return {
            "name": station_id,
            "latitude": lat,
            "longitude": lon,
            "WBGT (°F)": wbgt_val,
            "Temperature (°F)": celsius_to_farenheit(tair),
            "Dewpoint (°F)": celsius_to_farenheit(dwpt),
            "Wind Speed (mph)": wspd * 2.23694,
            "observation_time": obs,
            "source": "Mesonet"
        }
    except:
        return {
            "name": station_id,
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
def fetch_mesonet():
    rows = []
    for sid in stations_df["abbrev"]:
        if sid not in urls:  # WeatherSTEM stations are excluded
            rows.append(process_mesonet_station(sid))
    return pd.DataFrame(rows)

# -------------------------------------------------------------------
# FETCH USGS — ALL GAUGES IN KENTUCKY
# -------------------------------------------------------------------
@st.cache_data(ttl=300)
def fetch_usgs():
    url = (
        "https://waterservices.usgs.gov/nwis/iv/"
        "?format=json&stateCd=KY&parameterCd=00010,00020,00045,00060,00065"
    )

    try:
        j = requests.get(url, timeout=10).json()
        sites = j["value"]["timeSeries"]

        rows = []
        for s in sites:
            props = s["sourceInfo"]
            vals = s["values"][0]["value"]
            if not vals:
                continue

            latest = vals[-1]
            value = float(latest["value"]) if latest["value"] not in ["", None] else None
            timestamp = latest.get("dateTime", "")

            sid = props["siteCode"][0]["value"]
            lat = props["geoLocation"]["geogLocation"]["latitude"]
            lon = props["geoLocation"]["geogLocation"]["longitude"]
            variable = s["variable"]["variableName"]

            rows.append({
                "name": sid,
                "latitude": lat,
                "longitude": lon,
                "variable": variable,
                "value": value,
                "observation_time": timestamp,
                "source": "USGS"
            })
        return pd.DataFrame(rows)
    except:
        return pd.DataFrame(columns=["name","latitude","longitude","variable","value","observation_time","source"])

# -------------------------------------------------------------------
# LOAD DATA
# -------------------------------------------------------------------
df_ws = fetch_weatherstem() if show_ws else pd.DataFrame()
df_mes = fetch_mesonet() if show_mesonet else pd.DataFrame()
df_usgs = fetch_usgs() if show_usgs else pd.DataFrame()

# Attach coordinates to WeatherSTEM from input list
ws_coords = {}
for line in station_coords_text.strip().split("\n"):
    p = line.split(",")
    if len(p) == 3:
        ws_coords[p[0]] = (float(p[1]), float(p[2]))

df_ws["latitude"] = df_ws["name"].map(lambda x: ws_coords.get(x, (None,None))[0])
df_ws["longitude"] = df_ws["name"].map(lambda x: ws_coords.get(x, (None,None))[1])

# -------------------------------------------------------------------
# COMBINE MESONET + WEATHERSTEM
# -------------------------------------------------------------------
combined = pd.concat([df_ws, df_mes], ignore_index=True)

# -------------------------------------------------------------------
# COLOR SCALING
# -------------------------------------------------------------------
def variable_color(val, var):
    if pd.isna(val):
        return "#808080"
    if var == "WBGT (°F)":
        if val < 66: return "#008000"
        if val < 74: return "#FEF200"
        if val < 83: return "#FF0000"
        return "#000000"
    if var in ["Temperature (°F)", "Dewpoint (°F)"]:
        return cm.LinearColormap(["blue","lime","red"], vmin=30, vmax=100)(val)
    if var == "Wind Speed (mph)":
        return cm.LinearColormap(["white","cyan","blue"], vmin=0, vmax=20)(val)
    return "#808080"

# -------------------------------------------------------------------
# BASE MAP
# -------------------------------------------------------------------
lat_center = combined["latitude"].dropna().mean()
lon_center = combined["longitude"].dropna().mean()

m = folium.Map(location=[lat_center, lon_center], zoom_start=7, control_scale=True)

weather_layer = folium.FeatureGroup("Weather Stations")
usgs_layer = folium.FeatureGroup("USGS")

weather_layer.add_to(m)
usgs_layer.add_to(m)

# WEATHER / MESONET MARKERS
for _, row in combined.iterrows():
    lat, lon = row["latitude"], row["longitude"]
    if pd.isna(lat) or pd.isna(lon): continue
    val = row[selected_var]
    color = variable_color(val, selected_var)
    popup = f"<b>{row['name']} ({row['source']})</b><br>{selected_var}: {val}"
    folium.CircleMarker(
        location=[lat, lon],
        radius=7,
        color="black" if row["source"] == "WeatherSTEM" else color,
        weight=2 if row["source"] == "WeatherSTEM" else 1,
        fill=True,
        fill_color=color,
        fill_opacity=0.85,
        popup=popup
    ).add_to(weather_layer)

# USGS MARKERS
for _, row in df_usgs.iterrows():
    popup = f"""
    <b>USGS {row['name']}</b><br>
    {row['variable']}: {row['value']}<br>
    Observed: {row['observation_time']}
    """
    folium.CircleMarker(
        location=[row['latitude'], row['longitude']],
        radius=6,
        color="blue",
        weight=2,
        fill=True,
        fill_color="lightblue",
        fill_opacity=0.9,
        popup=popup
    ).add_to(usgs_layer)

st_folium(m, width=1000, height=650)

# -------------------------------------------------------------------
# COUNTY VIEW
# -------------------------------------------------------------------
counties_gdf = gpd.read_file("https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json")
counties_gdf = counties_gdf[counties_gdf["STATE"] == "21"]
counties_gdf["NAME"] = counties_gdf["NAME"].str.title()
county_list = sorted(counties_gdf["NAME"])

selected_county = st.sidebar.selectbox("Select a Kentucky County:", county_list)
geom = counties_gdf[counties_gdf["NAME"] == selected_county].geometry.iloc[0]

bounds = geom.bounds
county_map = folium.Map(
    location=[(bounds[1]+bounds[3])/2, (bounds[0]+bounds[2])/2],
    zoom_start=9
)

folium.GeoJson(
    geom.__geo_interface__,
    style_function=lambda x: {
        "fillColor": "#ff7800",
        "color": "black",
        "weight": 2,
        "fillOpacity": 0.25
    }
).add_to(county_map)

# Points inside county
pts = [Point(lon, lat) for lon, lat in zip(combined["longitude"], combined["latitude"])]
combined["in_county"] = [geom.contains(p) for p in pts]

subset = combined[combined["in_county"]]

for _, row in subset.iterrows():
    color = variable_color(row[selected_var], selected_var)
    folium.CircleMarker(
        location=[row.latitude, row.longitude],
        radius=8,
        fill=True,
        fill_color=color,
        color="black" if row["source"] == "WeatherSTEM" else color,
        weight=2,
        fill_opacity=0.9,
        popup=f"{row['name']}<br>{selected_var}: {row[selected_var]}"
    ).add_to(county_map)

st.markdown("### 🧭 County Focus View")
st_folium(county_map, width=850, height=450)
