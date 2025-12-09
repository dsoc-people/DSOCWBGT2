# -*- coding: utf-8 -*-
"""
Multi-Source Weather Station Data Collection and Visualization
Streamlit App for WBGT (Wet Bulb Globe Temperature) Monitoring
"""

import streamlit as st
import requests
import pandas as pd
import numpy as np
import folium
from geopy.geocoders import Nominatim
from geopy.extra.rate_limiter import RateLimiter
from streamlit_folium import st_folium
from streamlit_autorefresh import st_autorefresh
import base64
import os
from datetime import datetime

# -----------------------------
# Page Configuration
# -----------------------------
st.set_page_config(
    page_title="Multi-Source Weather Station WBGT Monitor",
    page_icon="🌡️",
    layout="wide"
)

# -----------------------------
# Unit Conversions
# -----------------------------
def farenheit_to_celsius(temp_f):
    return (temp_f - 32) * 5 / 9

def celsius_to_farenheit(temp_c):
    return temp_c * 9 / 5 + 32

# -----------------------------
# Wet Bulb Temperature Calculation
# -----------------------------
def dbdp2wb(tempC, dpC, p):
    return (tempC + dpC) / 2

# -----------------------------
# WBGT Calculation
# -----------------------------
def wbgt(tempF, mph, rad, bar, dpF):
    if tempF is None:
        tempF = np.nan
    if mph is None:
        mph = np.nan
    if rad is None:
        rad = np.nan
    if bar is None:
        bar = np.nan
    if dpF is None:
        dpF = np.nan

    tempC = farenheit_to_celsius(tempF)
    mps = mph * 0.44704
    tempK = tempC + 273.15

    if rad is None or np.isnan(rad):
        tempG = np.nan
    else:
        tempG = tempK + (rad - 30) / (0.0252 * rad + 10.5 * mps + 22.5 + 1e-9)
        tempG = tempG - 273.15

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

    wbgt_f = celsius_to_farenheit(wbgt_c)
    return wbgt_f

# -----------------------------
# Extract a sensor value from WeatherSTEM
# -----------------------------
def extract_value(records, target):
    for r in records:
        if target.lower() in r.get("sensor_name", "").lower():
            return r.get("value")
    return None

# -----------------------------
# WeatherSTEM (White Squirrel) Fetch
# -----------------------------
@st.cache_data(ttl=300)
def fetch_weatherstem_data():

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

    data = []

    for site, url in urls.items():
        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            site_json = response.json()

            records = site_json.get("records", [])
            time_obs = site_json.get("time", "N/A")

            wbgt_val = extract_value(records, "Wet Bulb Globe Temperature")

            temp = extract_value(records, "Thermometer")
            if temp is None:
                temp = extract_value(records, "Temperature")
            if temp is None:
                temp = extract_value(records, "Air Temperature")

            dewpt = extract_value(records, "Dewpoint")
            if dewpt is None:
                dewpt = extract_value(records, "Dew Point")

            wind = extract_value(records, "Anemometer")
            if wind is None:
                wind = extract_value(records, "Wind")

            data.append({
                "Site": site,
                "Observation Time": time_obs,
                "WBGT (°F)": wbgt_val,
                "Temperature (°F)": temp,
                "Dewpoint (°F)": dewpt,
                "Wind Speed (mph)": wind,
                "source": "White Squirrel Weather"
            })

        except Exception:
            data.append({
                "Site": site,
                "Observation Time": "Error",
                "WBGT (°F)": None,
                "Temperature (°F)": None,
                "Dewpoint (°F)": None,
                "Wind Speed (mph)": None,
                "source": "White Squirrel Weather"
            })

    return pd.DataFrame(data)

# -----------------------------
# USGS Fetch
# -----------------------------
@st.cache_data(ttl=300)
def fetch_usgs_data():

    usgs_url = 'https://waterservices.usgs.gov/nwis/iv/?format=json&stateCd=ky&siteStatus=active'

    try:
        response = requests.get(usgs_url, timeout=10)
        response.raise_for_status()
        json_data = response.json()
    except Exception:
        return pd.DataFrame()

    all_rows = []

    if json_data and 'value' in json_data and 'timeSeries' in json_data['value']:
        for series in json_data['value']['timeSeries']:
            site_info = series['sourceInfo']
            variable_info = series['variable']

            site_id = site_info['siteCode'][0]['value']
            site_name = site_info['siteName']
            lat = site_info['geoLocation']['geogLocation']['latitude']
            lon = site_info['geoLocation']['geogLocation']['longitude']

            latest = "N/A"
            if series['values'] and series['values'][0]['value']:
                latest = series['values'][0]['value'][-1]['value']

            unit = variable_info['unit'].get('unitAbbreviation', 'N/A')

            all_rows.append({
                'Site ID': site_id,
                'Site Name': site_name,
                'Latitude': lat,
                'Longitude': lon,
                'Parameter Name': variable_info['variableName'],
                'Latest Value': latest,
                'Value Unit': unit
            })

    return pd.DataFrame(all_rows)

# -----------------------------
# Parse Mesonet Station Coordinates
# -----------------------------
def get_station_coordinates():

    text = """FARM,36.93,-86.47
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

    coords = {}
    abbrevs = []

    for line in text.strip().splitlines():
        p = line.split(',')
        if len(p) == 3:
            try:
                abbrev = p[0]
                coords[abbrev] = (float(p[1]), float(p[2]))
                abbrevs.append(abbrev)
            except:
                pass

    return coords, abbrevs

# -----------------------------
# Mesonet Station Processor
# -----------------------------
@st.cache_data(ttl=300)
def process_mesonet_station(station_id, year, coordmap):

    lat, lon = coordmap.get(station_id, (None, None))

    try:
        manifest_url = f"https://d266k7wxhw6o23.cloudfront.net/data/{station_id}/{year}/manifest.json"
        manifest = requests.get(manifest_url, timeout=10).json()

        latest_day = max(manifest.keys())
        key = manifest[latest_day]["key"]

        json_data = requests.get(
            f"https://d266k7wxhw6o23.cloudfront.net/{key}",
            timeout=10
        ).json()

        df = pd.DataFrame(json_data["rows"], columns=json_data["columns"])

        required = ["TAIR", "DWPT", "WSPD", "SRAD", "PRES", "UTCTimestampCollected"]
        if not all(c in df.columns for c in required):
            raise ValueError("Missing Mesonet columns")

        tair = df["TAIR"].dropna().iloc[-1]
        dew = df["DWPT"].dropna().iloc[-1]
        wspd = df["WSPD"].dropna().iloc[-1]
        srad = df["SRAD"].dropna().iloc[-1]
        pres_hpa = df["PRES"].dropna().iloc[-1]
        timestamp = df["UTCTimestampCollected"].dropna().iloc[-1]

        temp_f = celsius_to_farenheit(tair)
        dew_f = celsius_to_farenheit(dew)
        wspd_mph = wspd * 2.23694
        pres_inhg = pres_hpa * 0.02953

        wbgt_val = wbgt(
            temp_f,
            wspd_mph,
            srad,
            pres_inhg,
            dew_f
        )

        return {
            "name": station_id,
            "latitude": lat,
            "longitude": lon,
            "wbgt_f": wbgt_val,
            "Temperature (°F)": temp_f,
            "Dewpoint (°F)": dew_f,
            "Wind Speed (mph)": wspd_mph,
            "observation_time": timestamp,
            "source": "Mesonet"
        }

    except Exception as e:
        return {
            "name": station_id,
            "latitude": lat,
            "longitude": lon,
            "wbgt_f": None,
            "Temperature (°F)": None,
            "Dewpoint (°F)": None,
            "Wind Speed (mph)": None,
            "observation_time": "Error",
            "source": "Mesonet"
        }

# -----------------------------
# Mesonet Fetch
# -----------------------------
@st.cache_data(ttl=300)
def fetch_mesonet_data(year, coordmap, abbrevs):

    rows = []

    for station in abbrevs:
        row = process_mesonet_station(station, year, coordmap)
        rows.append(row)

    df = pd.DataFrame(rows)

    # ------------ FIX: FORCE MESONET TEMP + DEW TO NUMERIC ------------
    df["Temperature (°F)"] = pd.to_numeric(df["Temperature (°F)"], errors="coerce")
    df["Dewpoint (°F)"] = pd.to_numeric(df["Dewpoint (°F)"], errors="coerce")

    return df

# -----------------------------
# Geocode WeatherSTEM stations
# -----------------------------
@st.cache_data(ttl=3600)
def geocode_stations(df):

    known = {
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

    df = df.copy()
    df["latitude"] = None
    df["longitude"] = None

    geolocator = Nominatim(user_agent="streamlit-wbgt-app")
    geocode = RateLimiter(geolocator.geocode, min_delay_seconds=1)

    def lookup(site):
        if site in known:
            return known[site]
        q = f"{site}, Kentucky"
        loc = geocode(q)
        if loc:
            return (loc.latitude, loc.longitude)
        return (None, None)

    for i, r in df.iterrows():
        lat, lon = lookup(r["Site"])
        df.loc[i, "latitude"] = lat
        df.loc[i, "longitude"] = lon

    return df

# -----------------------------
# WBGT Color Scale
# -----------------------------
def wbgt_color(v):
    if v is None or pd.isna(v):
        return "#808080"
    if v < 60:
        return "#2ca02c"
    if v < 70:
        return "#ff7f0e"
    if v < 80:
        return "#d62728"
    return "#800026"

def marker_radius():
    return 10

# -----------------------------
# Folium Map Builder
# -----------------------------
def create_map(combined_df, usgs_df, selection):

    # -----------------------------
    # Determine Center of Map
    # -----------------------------
    all_lat = []
    all_lon = []

    if not combined_df.empty:
        df = combined_df.dropna(subset=["latitude", "longitude"])
        all_lat += df["latitude"].tolist()
        all_lon += df["longitude"].tolist()

    if not usgs_df.empty:
        df = usgs_df.dropna(subset=["Latitude", "Longitude"])
        all_lat += df["Latitude"].tolist()
        all_lon += df["Longitude"].tolist()

    if all_lat:
        center = [np.mean(all_lat), np.mean(all_lon)]
    else:
        center = [37.2, -85.7]

    m = folium.Map(location=center, zoom_start=8)

    folium.TileLayer("cartodbpositron", name="Light").add_to(m)
    folium.TileLayer("cartodbdark_matter", name="Dark").add_to(m)
    folium.TileLayer("openstreetmap", name="OSM").add_to(m)

    weather_layer = folium.FeatureGroup(name="Weather Stations").add_to(m)
    usgs_layer = folium.FeatureGroup(name="USGS Gauges").add_to(m)

    # -----------------------------
    # Weather Stations
    # -----------------------------
    if not combined_df.empty:
        for _, row in combined_df.iterrows():

            site = row.get("name")
            lat = row.get("latitude")
            lon = row.get("longitude")
            wb = row.get("wbgt_f")
            temp = row.get("Temperature (°F)")
            dew = row.get("Dewpoint (°F)")
            time = row.get("observation_time", "N/A")
            src = row.get("source")

            if pd.isna(lat) or pd.isna(lon):
                continue

            if selection == "WBGT":
                value = wb
                label = "WBGT"
            elif selection == "Temperature":
                value = temp
                label = "Temperature"
            else:
                value = dew
                label = "Dewpoint"

            col = wbgt_color(wb)

            popup = f"""
            <div>
                <h4>{site} ({src})</h4>
                <b>Observed:</b> {time}<br>
                <b>WBGT:</b> {wb if pd.notna(wb) else "N/A"} °F<br>
                <b>Temperature:</b> {temp if pd.notna(temp) else "N/A"} °F<br>
                <b>Dewpoint:</b> {dew if pd.notna(dew) else "N/A"} °F<br>
            </div>
            """

            tooltip = f"{site}: {label} {value if pd.notna(value) else 'N/A'}"

            folium.CircleMarker(
                location=[lat, lon],
                radius=10,
                color=col,
                fill=True,
                fill_color=col,
                fill_opacity=0.85,
                popup=popup,
                tooltip=tooltip,
            ).add_to(weather_layer)

    # -----------------------------
    # USGS
    # -----------------------------
    if not usgs_df.empty:
        for sid in usgs_df["Site ID"].unique():

            rows = usgs_df[usgs_df["Site ID"] == sid]
            r0 = rows.iloc[0]

            lat = r0["Latitude"]
            lon = r0["Longitude"]
            name = r0["Site Name"]

            popup = f"<div><h4>{name}</h4>Site ID: {sid}<br>"

            for _, rr in rows.iterrows():
                vn = rr["Parameter Name"]
                v = rr["Latest Value"]
                u = rr["Value Unit"]
                popup += f"<b>{vn}:</b> {v} {u}<br>"

            popup += "</div>"

            folium.Marker(
                location=[lat, lon],
                popup=popup,
                tooltip=f"{name} ({sid})",
                icon=folium.Icon(color="blue", icon="tint", prefix="fa")
            ).add_to(usgs_layer)

    # -----------------------------
    # Legend
    # -----------------------------
    legend = """
    <div style="
        position: fixed; 
        bottom: 30px; left: 30px; 
        z-index: 9999; background: white; 
        padding: 10px; border: 1px solid #ccc;">
        <b>WBGT Scale (°F)</b><br>
        <span style="background:#2ca02c;width:12px;height:12px;display:inline-block;"></span> <60<br>
        <span style="background:#ff7f0e;width:12px;height:12px;display:inline-block;"></span> 60–69.9<br>
        <span style="background:#d62728;width:12px;height:12px;display:inline-block;"></span> 70–79.9<br>
        <span style="background:#800026;width:12px;height:12px;display:inline-block;"></span> ≥80<br>
        <span style="background:#808080;width:12px;height:12px;display:inline-block;"></span> N/A<br>
    </div>
    """

    m.get_root().html.add_child(folium.Element(legend))
    return m

# -----------------------------
# GitHub Upload Helper
# -----------------------------
def upload_map_to_github(filename):

    if not os.path.exists(filename):
        return {"success": False, "error": "File not found"}

    token = None
    try:
        token = st.secrets["GITHUB_TOKEN"]
    except:
        token = os.getenv("GITHUB_TOKEN")

    if not token:
        return {"success": False, "error": "No GitHub token provided"}

    repo_owner = "dsoc-people"
    repo_name = "DSOCWBGT2"
    path = "Image storage/wbgt_map.html"

    url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/contents/{path}"
    headers = {"Authorization": f"token {token}"}

    # Read mapping
    content = open(filename, "r", encoding="utf-8").read()
    encoded = base64.b64encode(content.encode()).decode()

    r = requests.get(url, headers=headers, timeout=10)

    sha = r.json().get("sha") if r.status_code == 200 else None

    payload = {
        "message": "Update WBGT map",
        "content": encoded,
        "branch": "main"
    }
    if sha:
        payload["sha"] = sha

    r = requests.put(url, headers=headers, json=payload, timeout=10)

    if r.status_code in (200, 201):
        return {
            "success": True,
            "url": f"https://github.com/{repo_owner}/{repo_name}/blob/main/{path}"
        }

    return {"success": False, "error": r.text}

# -----------------------------
# MAIN APP
# -----------------------------
def main():

    st.title("🌡️ Multi-Source Weather Station WBGT Monitor")

    if "load_mesonet" not in st.session_state:
        st.session_state["load_mesonet"] = False

    # --------------------
    # SIDEBAR
    # --------------------
    with st.sidebar:
        st.header("Data Sources")

        show_white = st.checkbox("White Squirrel Weather", True)
        show_mesonet = st.checkbox("Mesonet (slow)", True)
        show_usgs = st.checkbox("USGS River Gauges", False)
        upload = st.checkbox("Upload map to GitHub", False)

        if show_mesonet:
            if not st.session_state["load_mesonet"]:
                if st.button("Load Mesonet Data"):
                    st.session_state["load_mesonet"] = True
                    st.rerun()
            else:
                st.caption("Mesonet enabled (cached).")

        st_autorefresh(interval=300000, key="refresh")

    # --------------------
    # Selection
    # --------------------
    selection = st.radio(
        "Select display metric",
        ["WBGT", "Temperature", "Dewpoint"],
        horizontal=True
    )

    # --------------------
    # Fetch Data
    # --------------------
    df_white = pd.DataFrame()
    df_mesonet = pd.DataFrame()
    df_usgs = pd.DataFrame()

    if show_white:
        df_white = fetch_weatherstem_data()
        df_white = geocode_stations(df_white)
        df_white = df_white.rename(columns={"Site": "name", "Observation Time": "observation_time"})
        df_white["source"] = "White Squirrel Weather"

    if show_usgs:
        df_usgs = fetch_usgs_data()

    if show_mesonet and st.session_state["load_mesonet"]:
        coords, abbrevs = get_station_coordinates()
        year = datetime.utcnow().year

        with st.spinner("Loading Mesonet data…"):
            df_mesonet = fetch_mesonet_data(year, coords, abbrevs)

    # --------------------
    # Combine Weather Data
    # --------------------
    dfs = []

    if not df_white.empty:
        dfs.append(df_white)

    if not df_mesonet.empty:
        dfs.append(df_mesonet)

    if dfs:
        combined = pd.concat(dfs, ignore_index=True)
    else:
        combined = pd.DataFrame()

    # --------------------
    # Summary
    # --------------------
    c1, c2, c3 = st.columns(3)
    c1.metric("White Squirrel Stations", len(df_white))
    c2.metric("Mesonet Stations", len(df_mesonet))
    c3.metric("USGS Gauges", df_usgs["Site ID"].nunique() if not df_usgs.empty else 0)

    # --------------------
    # Tables
    # --------------------
    if st.checkbox("Show Data Tables"):
        t1, t2, t3 = st.tabs(["White Squirrel", "Mesonet", "USGS"])

        with t1:
            st.dataframe(df_white)

        with t2:
            st.dataframe(df_mesonet)

        with t3:
            st.dataframe(df_usgs)

    # --------------------
    # Map
    # --------------------
    if not combined.empty or not df_usgs.empty:
        st.subheader("Interactive Map")

        m = create_map(combined, df_usgs, selection)

        filename = "wbgt_map.html"
        m.save(filename)

        if upload:
            status = upload_map_to_github(filename)
            if status["success"]:
                st.success(f"Map uploaded: {status['url']}")
            else:
                st.warning(f"Upload failed: {status['error']}")

        st_folium(m, height=600, width=None)

    else:
        st.warning("No data available.")

if __name__ == "__main__":
    main()
