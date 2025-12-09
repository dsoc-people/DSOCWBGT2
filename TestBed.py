# app.py — Kentucky WBGT / Multi-Source Weather + USGS Map Dashboard
# Clean version with:
#   ✔ WeatherSTEM toggle
#   ✔ Mesonet toggle
#   ✔ USGS river gauges toggle
#   ✔ NO county focus view
#   ✔ Blue USGS water-drop icons
#   ✔ Full USGS parameter popups (Style B)

import requests
import pandas as pd
import numpy as np
import folium
from streamlit_folium import st_folium
import streamlit as st
import inspect
from streamlit_autorefresh import st_autorefresh
from shapely.geometry import Point
import branca.colormap as cm


# ---------------- Station Coordinates (Mesonet + WeatherSTEM) ----------------
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


# ---------------- Streamlit Setup ----------------
st.set_page_config(page_title="Kentucky WBGT Monitor", layout="wide")
st.title("🌡️ Kentucky WBGT / Weather Map Dashboard")

# Safe autorefresh
sig = inspect.signature(st_autorefresh)
if "rerun" in sig.parameters:
    st_autorefresh(interval=300000, limit=None, key="wbgt_refresh", rerun=False)
else:
    st_autorefresh(interval=300000, limit=None, key="wbgt_refresh")

year = "2025"


# ---------------- Sidebar ----------------
selected_var = st.sidebar.selectbox(
    "Variable to Display:",
    ["WBGT (°F)", "Temperature (°F)", "Dewpoint (°F)", "Wind Speed (mph)"]
)

st.sidebar.subheader("Data Layers")
show_mesonet = st.sidebar.checkbox("Show Mesonet", value=True)
show_weatherstem = st.sidebar.checkbox("Show WeatherSTEM (White Squirrel)", value=True)
show_usgs = st.sidebar.checkbox("Show USGS River Gauges", value=True)


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
    tempG = np.nan if rad is None else (
        tempK + (rad - 30) / (0.0252 * rad + 10.5*mps + 22.5) - 273.15
    )
    p = bar * 3.38639
    dpC = farenheit_to_celsius(dpF)
    wbc = dbdp2wb(tempC, dpC, p)
    wbgt_c = 0.7*wbc + 0.2*tempG + 0.1*tempC
    return celsius_to_farenheit(wbgt_c)


# ---------------- Parse Coordinates ----------------
@st.cache_data
def load_station_coords():
    rows = []
    for line in station_coords_text.strip().split("\n"):
        parts = line.split(",")
        if len(parts) == 3:
            rows.append({
                "name": parts[0],
                "lat": float(parts[1]),
                "lon": float(parts[2])
            })
    df = pd.DataFrame(rows)
    return df, {r["name"]: (r["lat"], r["lon"]) for _, r in df.iterrows()}

stations_df, station_coords = load_station_coords()


# ---------------- WeatherSTEM URLs ----------------
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
    "Jim Beam Clermont": "https://cdn/weatherstem.com/dashboard/data/dynamic/model/bullitt/jbclermont/latest.json",
    "Jim Beam Old Crow": "https://cdn/weatherstem.com/dashboard/data/dynamic/model/franklin-ky/jboldcrow/latest.json",
    "Jim Beam Grand Dad": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/franklin-ky/jbgranddad/latest.json",
    "Woodford County Courthouse": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/woodford/courthouse/latest.json",
    "Clinton County High School": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/clinton/clintonhs/latest.json",
    "Novelis Guthrie": "https://cdn.weatherstem.com/dashboard/data/dynamic/model/todd/novelis/latest.json",
}


# ---------------- WeatherSTEM Fetch ----------------
@st.cache_data(ttl=300)
def fetch_weatherstem_data():
    rows = []
    for site, url in urls.items():
        try:
            j = requests.get(url, timeout=10).json()
            r = j.get("records", [])
            rows.append({
                "name": site,
                "observation_time": j.get("time", "N/A"),
                "WBGT (°F)": extract_value(r, "Wet Bulb Globe Temperature"),
                "Temperature (°F)": extract_value(r, "Thermometer"),
                "Dewpoint (°F)": extract_value(r, "Dewpoint"),
                "Wind Speed (mph)": extract_value(r, "Anemometer"),
                "source": "WeatherSTEM"
            })
        except Exception:
            rows.append({
                "name": site,
                "observation_time": "Error",
                "WBGT (°F)": None,
                "Temperature (°F)": None,
                "Dewpoint (°F)": None,
                "Wind Speed (mph)": None,
                "source": "WeatherSTEM"
            })
    return pd.DataFrame(rows)


# ---------------- Mesonet Fetch ----------------
@st.cache_data(ttl=300)
def fetch_mesonet_data(station_ids):
    rows = []
    for sid in station_ids:
        rows.append(process_mesonet_station(sid))
    return pd.DataFrame(rows)

def process_mesonet_station(sid):
    lat, lon = station_coords.get(sid, (None, None))
    try:
        manifest_url = f"https://d266k7wxhw6o23.cloudfront.net/data/{sid}/{year}/manifest.json"
        manifest = requests.get(manifest_url, timeout=10).json()
        day = max(manifest.keys())
        key = manifest[day]["key"]
        data = requests.get(f"https://d266k7wxhw6o23.cloudfront.net/{key}", timeout=10).json()

        df = pd.DataFrame(data["rows"], columns=data["columns"])
        tair_c = df["TAIR"].dropna().iloc[-1]
        dwpt_c = df["DWPT"].dropna().iloc[-1]
        wspd_mps = df["WSPD"].dropna().iloc[-1]
        srad = df["SRAD"].dropna().iloc[-1]
        pres_hpa = df["PRES"].dropna().iloc[-1]
        obs = df["UTCTimestampCollected"].dropna().iloc[-1]

        wbgt_f = wbgt(
            celsius_to_farenheit(tair_c),
            wspd_mps * 2.23694,
            srad,
            pres_hpa * 0.02953,
            celsius_to_farenheit(dwpt_c),
        )

        return {
            "name": sid,
            "latitude": lat,
            "longitude": lon,
            "WBGT (°F)": wbgt_f,
            "Temperature (°F)": celsius_to_farenheit(tair_c),
            "Dewpoint (°F)": celsius_to_farenheit(dwpt_c),
            "Wind Speed (mph)": wspd_mps * 2.23694,
            "observation_time": obs,
            "source": "Mesonet",
        }

    except Exception:
        return {
            "name": sid,
            "latitude": lat,
            "longitude": lon,
            "WBGT (°F)": None,
            "Temperature (°F)": None,
            "Dewpoint (°F)": None,
            "Wind Speed (mph)": None,
            "observation_time": "Error",
            "source": "Mesonet",
        }


# ---------------- USGS Fetch (Full Data Option 1) ----------------
@st.cache_data(ttl=300)
def fetch_usgs_data():
    url = "https://waterservices.usgs.gov/nwis/iv/?format=json&stateCd=ky&siteStatus=active"
    try:
        j = requests.get(url, timeout=15).json()
        series = j["value"]["timeSeries"]

        rows = []
        for s in series:
            info = s["sourceInfo"]
            var = s["variable"]

            lat = info["geoLocation"]["geogLocation"]["latitude"]
            lon = info["geoLocation"]["geogLocation"]["longitude"]

            name = info["siteName"]
            sid = info["siteCode"][0]["value"]
            var_name = var["variableName"]
            unit = var["unit"].get("unitAbbreviation", "")

            # latest reading
            val = None
            if s["values"]:
                vlist = s["values"][0]["value"]
                if vlist:
                    val = vlist[-1]["value"]

            rows.append({
                "Site ID": sid,
                "Site Name": name,
                "Latitude": lat,
                "Longitude": lon,
                "Parameter Name": var_name,
                "Latest Value": val,
                "Value Unit": unit,
            })

        return pd.DataFrame(rows)

    except Exception:
        return pd.DataFrame(columns=[
            "Site ID", "Site Name", "Latitude", "Longitude",
            "Parameter Name", "Latest Value", "Value Unit"
        ])


# ---------------- Fetch All Data ----------------
ws_df = fetch_weatherstem_data() if show_weatherstem else pd.DataFrame()
mesonet_ids = [n for n in station_coords.keys() if n not in urls]
mesonet_df = fetch_mesonet_data(mesonet_ids) if show_mesonet else pd.DataFrame()
usgs_df = fetch_usgs_data() if show_usgs else pd.DataFrame()

# Add lat/lon to WeatherSTEM
for i, row in ws_df.iterrows():
    lat, lon = station_coords.get(row["name"], (None, None))
    ws_df.loc[i, "latitude"] = lat
    ws_df.loc[i, "longitude"] = lon

combined = pd.concat([mesonet_df, ws_df], ignore_index=True)


# ---------------- Color Mapping ----------------
def variable_color(val, var):
    if pd.isna(val):
        return "#808080"
    if var in ["Temperature (°F)", "Dewpoint (°F)"]:
        cmap = cm.LinearColormap(["#0000FF", "#00FF00", "#FF0000"], vmin=30, vmax=100)
        return cmap(val)
    if var == "WBGT (°F)":
        if val < 66: return "#008000"
        if val < 74: return "#FEF200"
        if val < 83: return "#FF0000"
        return "#000000"
    if var == "Wind Speed (mph)":
        cmap = cm.LinearColormap(["#FFFFFF", "#00FFFF", "#0000FF"], vmin=0, vmax=20)
        return cmap(val)
    return "#808080"


# ---------------- MAIN MAP ----------------
all_lats = combined["latitude"].dropna().tolist()
if show_usgs:
    all_lats += usgs_df["Latitude"].dropna().tolist()

all_lons = combined["longitude"].dropna().tolist()
if show_usgs:
    all_lons += usgs_df["Longitude"].dropna().tolist()

center_lat = np.mean(all_lats) if all_lats else 37.5
center_lon = np.mean(all_lons) if all_lons else -85.5

m = folium.Map(location=[center_lat, center_lon], zoom_start=7, control_scale=True)


# --- Mesonet + WeatherSTEM markers ---
for _, row in combined.iterrows():
    src = row["source"]
    if src == "Mesonet" and not show_mesonet:
        continue
    if src == "WeatherSTEM" and not show_weatherstem:
        continue

    lat, lon = row["latitude"], row["longitude"]
    if pd.isna(lat) or pd.isna(lon):
        continue

    val = row.get(selected_var)
    popup = (
        f"<b>{row['name']} ({row['source']})</b><br>"
        f"{selected_var}: {val if pd.notna(val) else 'N/A'}<br>"
        f"Observed: {row['observation_time']}"
    )

    folium.CircleMarker(
        location=[lat, lon],
        radius=7,
        color="black" if row["source"] == "WeatherSTEM" else variable_color(val, selected_var),
        weight=2 if row["source"] == "WeatherSTEM" else 1,
        fill=True,
        fill_color=variable_color(val, selected_var),
        fill_opacity=0.85,
        popup=popup,
    ).add_to(m)


# --- USGS markers ---
if show_usgs and not usgs_df.empty:
    for sid in usgs_df["Site ID"].unique():
        subset = usgs_df[usgs_df["Site ID"] == sid]
        first = subset.iloc[0]

        lat = first["Latitude"]
        lon = first["Longitude"]

        popup_html = (
            f"<b>{first['Site Name']}</b><br>"
            f"Site ID: {sid}<br>"
        )

        # append parameters for this site
        for _, r in subset.iterrows():
            if r["Latest Value"] is not None:
                popup_html += f"{r['Parameter Name']}: {r['Latest Value']} {r['Value Unit']}<br>"

        folium.Marker(
            location=[lat, lon],
            icon=folium.Icon(color="blue", icon="tint", prefix="fa"),
            tooltip=f"{first['Site Name']} ({sid})",
            popup=popup_html
        ).add_to(m)


# ---------------- RENDER MAP ----------------
st_folium(m, width=1100, height=700)
