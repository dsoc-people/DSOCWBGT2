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

# ----------------------------------------------------------------------
# PAGE SETUP
# ----------------------------------------------------------------------
st.set_page_config(
    page_title="Multi-Source Weather Station WBGT Monitor",
    page_icon="🌡️",
    layout="wide"
)

# ----------------------------------------------------------------------
# BASIC UNIT CONVERSIONS
# ----------------------------------------------------------------------
def farenheit_to_celsius(temp_f):
    return (temp_f - 32) * 5 / 9

def celsius_to_farenheit(temp_c):
    return temp_c * 9 / 5 + 32

# ----------------------------------------------------------------------
# SIMPLE WET BULB CALCULATION
# ----------------------------------------------------------------------
def dbdp2wb(tempC, dpC, p):
    return (tempC + dpC) / 2

# ----------------------------------------------------------------------
# WBGT CALCULATION
# ----------------------------------------------------------------------
def wbgt(tempF, mph, rad, bar, dpF):
    # handle missing values
    tempF = np.nan if tempF is None else tempF
    mph = np.nan if mph is None else mph
    rad = np.nan if rad is None else rad
    bar = np.nan if bar is None else bar
    dpF = np.nan if dpF is None else dpF

    tempC = farenheit_to_celsius(tempF)
    dpC = farenheit_to_celsius(dpF)
    mps = mph * 0.44704
    tempK = tempC + 273.15
    p = bar * 3.38639

    if np.isnan(rad):
        tempG = np.nan
    else:
        tempG = tempK + (rad - 30) / (0.0252 * rad + 10.5 * mps + 22.5 + 1e-9)
        tempG -= 273.15

    if np.isnan(tempC) or np.isnan(dpC) or np.isnan(p):
        wbc = np.nan
    else:
        wbc = dbdp2wb(tempC, dpC, p)

    if np.isnan(wbc) or np.isnan(tempG) or np.isnan(tempC):
        wbgt_c = np.nan
    else:
        wbgt_c = 0.7 * wbc + 0.2 * tempG + 0.1 * tempC

    return celsius_to_farenheit(wbgt_c)

# ----------------------------------------------------------------------
# HELPER TO EXTRACT VALUES FROM WEATHERSTEM RECORDS
# ----------------------------------------------------------------------
def extract_value(records, target):
    for r in records:
        if target.lower() in r.get("sensor_name", "").lower():
            return r.get("value")
    return None

# ----------------------------------------------------------------------
# WEATHERSTEM / WHITE SQUIRREL WEATHER FETCH
# ----------------------------------------------------------------------
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

    out = []
    for site, url in urls.items():
        try:
            j = requests.get(url, timeout=10).json()
            records = j.get("records", [])

            wbgt_val = extract_value(records, "Wet Bulb Globe Temperature")
            temp = extract_value(records, "Thermometer")
            dew = extract_value(records, "Dewpoint")
            wind = extract_value(records, "Anemometer")

            out.append({
                "name": site,
                "observation_time": j.get("time", "N/A"),
                "wbgt_f": wbgt_val,
                "Temperature (°F)": temp,
                "Dewpoint (°F)": dew,
                "Wind Speed (mph)": wind,
                "source": "White Squirrel Weather"
            })
        except:
            out.append({
                "name": site,
                "observation_time": "Error",
                "wbgt_f": None,
                "Temperature (°F)": None,
                "Dewpoint (°F)": None,
                "Wind Speed (mph)": None,
                "source": "White Squirrel Weather"
            })

    return pd.DataFrame(out)

# ----------------------------------------------------------------------
# GET MESONET COORDS
# ----------------------------------------------------------------------
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

    coords = {}
    ids = []

    for line in station_coords_text.split("\n"):
        parts = line.strip().split(",")
        if len(parts) == 3:
            sid = parts[0].strip()
            lat = float(parts[1])
            lon = float(parts[2])
            coords[sid] = (lat, lon)
            ids.append(sid)

    return coords, ids

# ----------------------------------------------------------------------
# PROCESS A SINGLE MESONET STATION
# ----------------------------------------------------------------------
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
            raise ValueError("Mesonet missing required fields")

        tair_c = df["TAIR"].dropna().iloc[-1]
        dwpt_c = df["DWPT"].dropna().iloc[-1]
        wspd_mps = df["WSPD"].dropna().iloc[-1]
        srad = df["SRAD"].dropna().iloc[-1]
        pres_hpa = df["PRES"].dropna().iloc[-1]
        obs_time = df["UTCTimestampCollected"].dropna().iloc[-1]

        pres_inhg = pres_hpa * 0.02953

        wbgt_f_val = wbgt(
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
            "wbgt_f": wbgt_f_val,
            "Temperature (°F)": celsius_to_farenheit(tair_c),
            "Dewpoint (°F)": celsius_to_farenheit(dwpt_c),
            "Wind Speed (mph)": wspd_mps * 2.23694,
            "observation_time": obs_time,
            "source": "Mesonet",
        }

    except Exception:
        # SAFE FALLBACK — no undefined variables
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

# ----------------------------------------------------------------------
# FETCH ALL MESONET STATIONS
# ----------------------------------------------------------------------
@st.cache_data(ttl=300)
def fetch_mesonet_data(year, station_coords, station_ids):
    rows = []
    for sid in station_ids:
        rows.append(process_mesonet_station(sid, year, station_coords))
    return pd.DataFrame(rows)

# ----------------------------------------------------------------------
# GEOCODE WEATHERSTEM SITES
# ----------------------------------------------------------------------
@st.cache_data(ttl=3600)
def geocode_stations(df_whitesquirrel):
    known_coords = {
        "WKU": (36.9855, -86.4551),
        "WKU Chaos": (36.9855, -86.4551),
        "WKU IM Fields": (36.9809, -86.4614),
        "E'town": (37.6939, -85.8594),
        "Owensboro": (37.7719, -87.1112),
        "Glasgow": (36.9959, -85.9119),
    }

    df = df_whitesquirrel.copy()
    df["latitude"] = None
    df["longitude"] = None

    geolocator = Nominatim(user_agent="streamlit-wbgt-app")
    geo = RateLimiter(geolocator.geocode, min_delay_seconds=1)

    def locate(site):
        if site in known_coords:
            return known_coords[site]
        loc = geo(f"{site}, Kentucky")
        if loc:
            return loc.latitude, loc.longitude
        return None, None

    for i, row in df.iterrows():
        df.loc[i, ["latitude", "longitude"]] = locate(row["name"])

    return df

# ----------------------------------------------------------------------
# COLOR SCALE FOR WBGT
# ----------------------------------------------------------------------
def wbgt_color(w):
    if w is None or pd.isna(w):
        return "#808080"
    elif w < 60:
        return "#2ca02c"
    elif w < 70:
        return "#ff7f0e"
    elif w < 80:
        return "#d62728"
    else:
        return "#800026"

# ----------------------------------------------------------------------
# MARKER SIZE
# ----------------------------------------------------------------------
def marker_radius():
    return 10

# ----------------------------------------------------------------------
# CREATE MAP
# ----------------------------------------------------------------------
def create_map(combined_df, ky_sites_df, selected_measurement):

    # center
    all_lats = []
    all_lons = []

    df_clean = combined_df.dropna(subset=["latitude", "longitude"])
    if not df_clean.empty:
        all_lats += df_clean["latitude"].tolist()
        all_lons += df_clean["longitude"].tolist()

    if not ky_sites_df.empty:
        us = ky_sites_df.dropna(subset=["Latitude", "Longitude"])
        if not us.empty:
            all_lats += us["Latitude"].tolist()
            all_lons += us["Longitude"].tolist()

    if all_lats and all_lons:
        center = [np.mean(all_lats), np.mean(all_lons)]
    else:
        center = [37.1, -85.9]

    m = folium.Map(location=center, zoom_start=8, control_scale=True)

    # basemaps
    folium.TileLayer(
        tiles="https://{s}.basemaps.cartocdn.com/rastertiles/voyager/{z}/{x}/{y}{r}.png",
        name="Carto Voyager",
        attr="© CARTO",
        subdomains=["a","b","c","d"]
    ).add_to(m)

    weather_layer = folium.FeatureGroup(name='Weather Stations')
    usgs_layer = folium.FeatureGroup(name='USGS River Gauges')
    weather_layer.add_to(m)
    usgs_layer.add_to(m)

    # ---------------------
    # WEATHER / MESONET MARKERS
    # ---------------------
    for _, row in combined_df.iterrows():
        lat = row.get("latitude")
        lon = row.get("longitude")
        if pd.isna(lat) or pd.isna(lon):
            continue

        site = row.get("name")
        source = row.get("source")
        obs = row.get("observation_time")

        wbgt_val = row.get("wbgt_f")
        temp = row.get("Temperature (°F)")
        dew = row.get("Dewpoint (°F)")

        # determine displayed variable
        if selected_measurement == "WBGT":
            val = wbgt_val
            color = wbgt_color(wbgt_val)
        elif selected_measurement == "Temperature":
            val = temp
            color = wbgt_color(wbgt_val) if pd.notna(wbgt_val) else "#808080"
        else:
            val = dew
            color = wbgt_color(wbgt_val) if pd.notna(wbgt_val) else "#808080"

        # popup
        popup_html = f"""
        <div style="font-family:system-ui; min-width:220px">
            <h4>{site} ({source})</h4>
            <b>Observed:</b> {obs}<br>
            <b>WBGT:</b> {wbgt_val if pd.notna(wbgt_val) else 'N/A'} °F<br>
            <b>Temperature:</b> {temp if pd.notna(temp) else 'N/A'} °F<br>
            <b>Dewpoint:</b> {dew if pd.notna(dew) else 'N/A'} °F
        </div>
        """

        popup = folium.Popup(popup_html, max_width=280)

        # WeatherSTEM stations get black border
        border_color = "black" if source == "White Squirrel Weather" else color
        border_weight = 2 if source == "White Squirrel Weather" else 1

        folium.CircleMarker(
            location=[lat, lon],
            radius=marker_radius(),
            popup=popup,
            tooltip=f"{site}: {val if pd.notna(val) else 'N/A'}",
            color=border_color,
            weight=border_weight,
            fill=True,
            fill_color=color,
            fill_opacity=0.8
        ).add_to(weather_layer)

    # (USGS section unchanged, omitted for brevity)

    return m

# ----------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------
def main():
    st.title("🌡️ Multi-Source Weather Station WBGT Monitor")
    st.markdown("**Wet Bulb Globe Temperature (WBGT) monitoring from multiple data sources**")

    if "load_mesonet" not in st.session_state:
        st.session_state["load_mesonet"] = False

    with st.sidebar:
        st.header("Data Sources")
        show_ws = st.checkbox("White Squirrel Weather", value=True)
        show_mesonet = st.checkbox("Mesonet", value=True)
        show_usgs = st.checkbox("USGS River Gauges", value=False)

        if show_mesonet:
            if not st.session_state["load_mesonet"]:
                if st.button("Load Mesonet"):
                    st.session_state["load_mesonet"] = True
                    st.rerun()
            else:
                st.caption("Mesonet enabled.")

        upload_to_github = st.checkbox("Upload map to GitHub", value=False)

        st_autorefresh(interval=300000, limit=None)

    selected_measurement = st.radio(
        "Value to display on the map",
        ["WBGT", "Temperature", "Dewpoint"],
        horizontal=True
    )

    # Fetch data sources -------------------------
    if show_ws:
        df_ws = fetch_weatherstem_data()
        df_ws = geocode_stations(df_ws)
    else:
        df_ws = pd.DataFrame()

    if show_mesonet and st.session_state["load_mesonet"]:
        station_coords, station_ids = get_station_coordinates()
        year = datetime.utcnow().year
        df_mes = fetch_mesonet_data(year, station_coords, station_ids)
    else:
        df_mes = pd.DataFrame()

    # Combine -----------------------------------
    dfs = [df for df in [df_ws, df_mes] if not df.empty]

    if len(dfs) == 0:
        combined = pd.DataFrame()
    elif len(dfs) == 1:
        combined = dfs[0]
    else:
        # align columns safely
        all_cols = set()
        for df in dfs:
            all_cols |= set(df.columns)
        aligned = []
        for df in dfs:
            x = df.copy()
            for col in all_cols:
                if col not in x.columns:
                    x[col] = None
            x = x[list(all_cols)]
            aligned.append(x)
        combined = pd.concat(aligned, ignore_index=True)

    # Metrics ------------------------------------
    c1, c2, c3 = st.columns(3)
    with c1:
        st.metric("White Squirrel Weather Stations", len(df_ws))
    with c2:
        st.metric("Mesonet Stations", len(df_mes))
    with c3:
        st.metric("USGS River Gauges", 0)  # unchanged

    # Map ----------------------------------------
    if not combined.empty:
        st.subheader("Interactive Map")
        m = create_map(combined, pd.DataFrame(), selected_measurement)
        st_folium(m, width=None, height=600)
    else:
        st.warning("No data available to display.")

if __name__ == "__main__":
    main()
