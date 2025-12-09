# -*- coding: utf-8 -*-
"""
Multi-Source Weather Station Data Collection and Visualization
Streamlit App for WBGT (Wet Bulb Globe Temperature) Monitoring
(Patched Version – Mesonet Temp/Dewpoint Fix + WeatherSTEM Black Border)
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

# Page configuration
st.set_page_config(
    page_title="Multi-Source Weather Station WBGT Monitor",
    page_icon="🌡️",
    layout="wide"
)

# --------------------
# Unit conversion
# --------------------
def farenheit_to_celsius(f):
    return (f - 32) * 5 / 9

def celsius_to_farenheit(c):
    return c * 9 / 5 + 32

# --------------------
# WBGT Helpers
# --------------------
def dbdp2wb(tempC, dpC, p):
    return (tempC + dpC) / 2

def wbgt(tempF, mph, rad, bar, dpF):
    # Convert inputs
    tempC = farenheit_to_celsius(tempF) if tempF is not None else np.nan
    mps = mph * 0.44704 if mph is not None else np.nan
    tempK = tempC + 273.15

    # Globe temp
    if rad is None or np.isnan(rad):
        tempG = np.nan
    else:
        tempG = tempK + (rad - 30) / (0.0252 * rad + 10.5*mps + 22.5 + 1e-9)
        tempG -= 273.15

    p = bar * 3.38639 if bar is not None else np.nan
    dpC = farenheit_to_celsius(dpF) if dpF is not None else np.nan

    if np.isnan(tempC) or np.isnan(dpC) or np.isnan(p):
        wbc = np.nan
    else:
        wbc = dbdp2wb(tempC, dpC, p)

    if np.isnan(wbc) or np.isnan(tempG) or np.isnan(tempC):
        wbg_c = np.nan
    else:
        wbg_c = 0.7*wbc + 0.2*tempG + 0.1*tempC

    return celsius_to_farenheit(wbg_c)

# --------------------
# Extract value from WeatherSTEM JSON
# --------------------
def extract_value(records, target):
    for r in records:
        if target.lower() in r.get("sensor_name", "").lower():
            return r.get("value")
    return None

# --------------------
# Fetch WeatherSTEM / White Squirrel Weather
# (UNCHANGED as required)
# --------------------
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
            j = requests.get(url, timeout=10).json()
            records = j.get("records", [])
            t = j.get("time", "N/A")

            wbgt_val = extract_value(records, "Wet Bulb Globe Temperature")
            temp = extract_value(records, "Thermometer") or \
                   extract_value(records, "Temperature") or \
                   extract_value(records, "Air Temperature")
            dew = extract_value(records, "Dewpoint") or extract_value(records, "Dew Point")
            wind = extract_value(records, "Anemometer") or extract_value(records, "Wind")

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
            data.append({
                "name": site,
                "observation_time": "Error",
                "wbgt_f": None,
                "Temperature (°F)": None,
                "Dewpoint (°F)": None,
                "Wind Speed (mph)": None,
                "source": "White Squirrel Weather"
            })

    return pd.DataFrame(data)

# --------------------
# Mesonet station coordinate parser
# --------------------
def get_station_coordinates():
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
BEDD,38.63,-85.32"""
    station_coords = {}
    station_ids = []

    for line in station_coords_text.splitlines():
        parts = line.split(",")
        if len(parts) == 3:
            sid = parts[0]
            lat, lon = float(parts[1]), float(parts[2])
            station_coords[sid] = (lat, lon)
            station_ids.append(sid)

    return station_coords, station_ids

# --------------------
# Mesonet process logic (PATCHED)
# --------------------
@st.cache_data(ttl=300)
def process_mesonet_station(station_id, year, station_coords):
    lat, lon = station_coords.get(station_id, (None, None))

    try:
        manifest_url = f"https://d266k7wxhw6o23.cloudfront.net/data/{station_id}/{year}/manifest.json"
        manifest = requests.get(manifest_url, timeout=15).json()

        if not manifest:
            raise ValueError("Manifest empty")

        latest_day = max(manifest.keys())
        key = manifest[latest_day]["key"]

        j = requests.get(f"https://d266k7wxhw6o23.cloudfront.net/{key}", timeout=15).json()
        df = pd.DataFrame(j["rows"], columns=j["columns"])

        cols = ["TAIR", "DWPT", "WSPD", "SRAD", "PRES", "UTCTimestampCollected"]
        if not all(c in df.columns for c in cols):
            raise ValueError("Mesonet missing required data fields")

        tair_c = df["TAIR"].dropna().iloc[-1]
        dwpt_c = df["DWPT"].dropna().iloc[-1]
        wspd_mps = df["WSPD"].dropna().iloc[-1]
        srad = df["SRAD"].dropna().iloc[-1]
        pres_hpa = df["PRES"].dropna().iloc[-1]
        obs_time = df["UTCTimestampCollected"].dropna().iloc[-1]

        pres_inhg = pres_hpa * 0.02953

        wbgt_val = wbgt(
            celsius_to_farenheit(tair_c),
            wspd_mps * 2.23694,
            srad,
            pres_inhg,
            celsius_to_farenheit(dwpt_c)
        )

        return {
            "name": station_id,
            "latitude": lat,
            "longitude": lon,
            "wbgt_f": wbgt_val,
            "Temperature (°F)": celsius_to_farenheit(tair_c),
            "Dewpoint (°F)": celsius_to_farenheit(dwpt_c),
            "Wind Speed (mph)": wspd_mps * 2.23694,
            "observation_time": obs_time,
            "source": "Mesonet",
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
            "source": "Mesonet",
        }

# --------------------
# Fetch all Mesonet stations
# --------------------
@st.cache_data(ttl=300)
def fetch_mesonet_data(year, station_coords, station_ids):
    rows = []
    for sid in station_ids:
        rows.append(process_mesonet_station(sid, year, station_coords))
    return pd.DataFrame(rows)

# --------------------
# WeatherSTEM geocoding (UNCHANGED)
# --------------------
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
        "Novelis Guthrie": (36.6025431022, -87.7186136559),
    }

    df = df.copy()
    df["latitude"] = None
    df["longitude"] = None

    geolocator = Nominatim(user_agent="streamlit-wbgt")
    geocode = RateLimiter(geolocator.geocode, min_delay_seconds=1)

    for i, row in df.iterrows():
        site = row["name"]
        if site in known:
            df.loc[i, ["latitude", "longitude"]] = known[site]
        else:
            df.loc[i, ["latitude", "longitude"]] = (None, None)

    return df

# --------------------
# Marker color by WBGT
# --------------------
def wbgt_color(w):
    if w is None or pd.isna(w): return "#808080"
    if w < 60: return "#2ca02c"
    if w < 70: return "#ff7f0e"
    if w < 80: return "#d62728"
    return "#800026"

# --------------------
# Main map creation
# --------------------
def marker_radius():
    return 10

def create_map(combined_df, ky_sites_df, selected_measurement):
    all_lats, all_lons = [], []

    dfc = combined_df.dropna(subset=["latitude", "longitude"])
    all_lats.extend(dfc["latitude"].tolist())
    all_lons.extend(dfc["longitude"].tolist())

    dfu = ky_sites_df.dropna(subset=["Latitude", "Longitude"])
    all_lats.extend(dfu["Latitude"].tolist())
    all_lons.extend(dfu["Longitude"].tolist())

    center = [np.mean(all_lats), np.mean(all_lons)] if all_lats else [37.1, -85.9]
    m = folium.Map(location=center, zoom_start=8, control_scale=True)

    # Base layers
    folium.TileLayer(
        tiles="https://{s}.basemaps.cartocdn.com/rastertiles/voyager/{z}/{x}/{y}{r}.png",
        attr='CartoDB', name="CartoDB Voyager",
    ).add_to(m)

    folium.TileLayer(
        tiles="https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png",
        attr='CartoDB', name="Dark Matter"
    ).add_to(m)

    folium.TileLayer(
        tiles="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
        attr='Esri', name="Esri World Imagery"
    ).add_to(m)

    weather_layer = folium.FeatureGroup(name="Weather Stations").add_to(m)
    usgs_layer = folium.FeatureGroup(name="USGS River Gauges").add_to(m)

    # Weather markers
    for _, row in combined_df.iterrows():
        lat = row.get("latitude")
        lon = row.get("longitude")
        if pd.isna(lat) or pd.isna(lon): continue

        site = row.get("name")
        source = row.get("source")
        obs = row.get("observation_time")

        wbgt_val = row.get("wbgt_f")
        temp = row.get("Temperature (°F)")
        dew = row.get("Dewpoint (°F)")

        if selected_measurement == "WBGT":
            val = wbgt_val
            color = wbgt_color(wbgt_val)
            label = "WBGT"
        elif selected_measurement == "Temperature":
            val = temp
            color = wbgt_color(wbgt_val) if pd.notna(wbgt_val) else "#808080"
            label = "Temperature"
        else:
            val = dew
            color = wbgt_color(wbgt_val) if pd.notna(wbgt_val) else "#808080"
            label = "Dewpoint"

        popup_html = f"""
        <div style="font-family:system-ui; min-width:220px">
            <h4>{site} ({source})</h4>
            <b>Observed:</b> {obs}<br>
            <b>WBGT:</b> {wbgt_val if pd.notna(wbgt_val) else 'N/A'} °F<br>
            <b>Temperature:</b> {temp if pd.notna(temp) else 'N/A'} °F<br>
            <b>Dewpoint:</b> {dew if pd.notna(dew) else 'N/A'} °F
        </div>
        """

        # WEATHERSTEM gets black border
        edge_color = "black" if source == "White Squirrel Weather" else color
        edge_weight = 3 if source == "White Squirrel Weather" else 1

        folium.CircleMarker(
            location=[lat, lon],
            radius=marker_radius(),
            color=edge_color,
            weight=edge_weight,
            fill=True,
            fill_color=color,
            fill_opacity=0.85,
            popup=folium.Popup(popup_html, max_width=260),
        ).add_to(weather_layer)

    # USGS markers (unchanged)
    for site_id in ky_sites_df['Site ID'].unique():
        rows = ky_sites_df[ky_sites_df['Site ID'] == site_id]
        r0 = rows.iloc[0]
        lat, lon = r0['Latitude'], r0['Longitude']
        site_name = r0['Site Name']

        popup = f"<b>{site_name}</b><br>Site ID: {site_id}<br>"
        for _, r in rows.iterrows():
            popup += f"{r['Parameter Name']}: {r['Latest Value']} {r['Value Unit']}<br>"

        folium.Marker(
            location=[lat, lon],
            popup=popup,
            icon=folium.Icon(color="blue", icon="tint", prefix="fa")
        ).add_to(usgs_layer)

    return m

# --------------------
# GitHub Upload (unchanged)
# --------------------
def upload_map_to_github(map_filename):
    try:
        if not os.path.exists(map_filename):
            return {"success": False, "error": f"{map_filename} missing"}

        with open(map_filename, "r", encoding="utf-8") as f:
            content = f.read()

        repo_owner = "dsoc-people"
        repo_name = "DSOCWBGT2"
        file_path = "Image storage/wbgt_map.html"

        token = st.secrets.get("GITHUB_TOKEN", None) or os.getenv("GITHUB_TOKEN")
        if not token:
            return {"success": False, "error": "No GitHub token found"}

        url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/contents/{file_path}"
        headers = {"Authorization": f"token {token}"}

        get_r = requests.get(url, headers=headers)
        sha = get_r.json().get("sha") if get_r.status_code == 200 else None

        encoded = base64.b64encode(content.encode()).decode()
        data = {"message": "Update WBGT map", "content": encoded, "branch": "main"}
        if sha:
            data["sha"] = sha

        put_r = requests.put(url, json=data, headers=headers)
        if put_r.status_code in (200, 201):
            return {"success": True, "url": f"https://github.com/{repo_owner}/{repo_name}/blob/main/{file_path}"}

        return {"success": False, "error": put_r.text}

    except Exception as e:
        return {"success": False, "error": str(e)}

# --------------------
# MAIN APP
# --------------------
def main():

    st.title("🌡️ Multi-Source Weather Station WBGT Monitor")

    # Autorefresh
    st_autorefresh(interval=300000, key="refresh")

    # Sidebar
    with st.sidebar:
        show_ws = st.checkbox("White Squirrel Weather", True)
        show_mesonet = st.checkbox("Mesonet (slower)", True)
        show_usgs = st.checkbox("USGS River Gauges", False)
        upload_to_github = st.checkbox("Upload map to GitHub", False)

    selected_measurement = st.radio(
        "Value to display:", ["WBGT", "Temperature", "Dewpoint"], horizontal=True
    )

    # Fetch sources
    if show_ws:
        df_ws = fetch_weatherstem_data()
        df_ws = geocode_stations(df_ws)
    else:
        df_ws = pd.DataFrame()

    year = datetime.utcnow().year
    station_coords, station_ids = get_station_coordinates()

    df_mesonet = (
        fetch_mesonet_data(year, station_coords, station_ids) if show_mesonet else pd.DataFrame()
    )

    df_usgs = fetch_usgs_data() if show_usgs else pd.DataFrame()

    # Combine
    dfs = []
    if not df_ws.empty: dfs.append(df_ws)
    if not df_mesonet.empty: dfs.append(df_mesonet)
    combined = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

    # Summary
    c1, c2, c3 = st.columns(3)
    c1.metric("White Squirrel", len(df_ws))
    c2.metric("Mesonet", len(df_mesonet))
    c3.metric("USGS Gauges", df_usgs['Site ID'].nunique() if show_usgs else 0)

    # Data tables
    if st.checkbox("Show Data Tables"):
        tab1, tab2, tab3 = st.tabs(["White Squirrel", "Mesonet", "USGS"])
        with tab1: st.dataframe(df_ws)
        with tab2: st.dataframe(df_mesonet)
        with tab3: st.dataframe(df_usgs)

    # Map
    if not combined.empty or not df_usgs.empty:
        m = create_map(combined, df_usgs, selected_measurement)
        map_file = "wbgt_map.html"
        m.save(map_file)

        if upload_to_github:
            result = upload_map_to_github(map_file)
            if result["success"]:
                st.success(result["url"])
            else:
                st.warning(result["error"])

        st_folium(m, height=600)

    else:
        st.warning("No data available.")

# --------------------
if __name__ == "__main__":
    main()
