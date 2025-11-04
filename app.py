# app.py
# -*- coding: utf-8 -*-
import requests
import pandas as pd
import numpy as np
import folium
from streamlit_folium import st_folium
import streamlit as st

# --- Configuration ---
st.set_page_config(page_title="Kentucky WBGT Monitor", layout="wide")
year = "2025"

# --- Unit Conversions ---
def farenheit_to_celsius(temp_f):
    return (temp_f - 32) * 5 / 9

def celsius_to_farenheit(temp_c):
    return temp_c * 9 / 5 + 32

# --- Approximate Wet Bulb Calculation (placeholder) ---
def dbdp2wb(tempC, dpC, p):
    # Replace with psychrometric calc later if needed
    return (tempC + dpC) / 2

# --- WBGT Function ---
def wbgt(tempF, mph, rad, bar, dpF):
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

    return celsius_to_farenheit(wbgt_c)

# --- Station Coordinates ---
@st.cache_data
def load_station_coords():
    url = "https://d266k7wxhw6o23.cloudfront.net/metadata/stations_468eb55962c18d1fc333160925381b9d6fb5eb86cd6fbbfbfc285b1d6fcfe7a0.json"
    data = requests.get(url).json()
    df = pd.DataFrame(data)
    coords = {row["abbrev"]: (row["lat"], row["lon"]) for _, row in df.iterrows()}
    return df, coords

stations_df, station_coords = load_station_coords()
station_abbrevs = stations_df["abbrev"].tolist()

# --- Fetch + Process WBGT Data ---
@st.cache_data(ttl=900)
def process_station_data(station_id):
    try:
        manifest_url = f"https://d266k7wxhw6o23.cloudfront.net/data/{station_id}/{year}/manifest.json"
        manifest = requests.get(manifest_url).json()
        latest_day = max(manifest.keys())
        data_key = manifest[latest_day]["key"]
        data_url = f"https://d266k7wxhw6o23.cloudfront.net/{data_key}"
        data = requests.get(data_url).json()
        df = pd.DataFrame(data['rows'], columns=data['columns'])

        cols = ["TAIR", "DWPT", "WSPD", "SRAD", "PRES", "UTCTimestampCollected"]
        if not all(c in df.columns for c in cols):
            return None

        tair_c = df["TAIR"].dropna().iloc[-1]
        dwpt_c = df["DWPT"].dropna().iloc[-1]
        wspd_mps = df["WSPD"].dropna().iloc[-1]
        srad = df["SRAD"].dropna().iloc[-1]
        pres_hpa = df["PRES"].dropna().iloc[-1]
        pres_inhg = pres_hpa * 0.02953
        obs_time = df["UTCTimestampCollected"].dropna().iloc[-1]

        tair_f = celsius_to_farenheit(tair_c)
        dwpt_f = celsius_to_farenheit(dwpt_c)
        wspd_mph = wspd_mps * 2.23694

        wbgt_f = wbgt(tair_f, wspd_mph, srad, pres_inhg, dwpt_f)
        lat, lon = station_coords.get(station_id, (None, None))

        return {
            "Station": station_id,
            "WBGT (°F)": wbgt_f,
            "Latitude": lat,
            "Longitude": lon,
            "Time": obs_time
        }
    except Exception as e:
        return None

# --- UI Layout ---
st.title("🌡️ Kentucky Mesonet WBGT Map")
st.markdown("Displays current **Wet Bulb Globe Temperature (°F)** across Mesonet stations.")

selected_stations = st.multiselect(
    "Select stations to process",
    station_abbrevs,
    default=["FARM", "RSVL", "MRRY"]
)

if st.button("Fetch WBGT Data"):
    with st.spinner("Processing station data..."):
        results = [process_station_data(s) for s in selected_stations]
        results = [r for r in results if r]
        if not results:
            st.error("No valid data found.")
        else:
            wbgt_df = pd.DataFrame(results)

            st.dataframe(
                wbgt_df.style.format({"WBGT (°F)": "{:.1f}"}),
                use_container_width=True
            )

            # --- Folium Map ---
            avg_lat, avg_lon = wbgt_df["Latitude"].mean(), wbgt_df["Longitude"].mean()
            fmap = folium.Map(location=[avg_lat, avg_lon], zoom_start=7)

            for _, row in wbgt_df.iterrows():
                folium.CircleMarker(
                    location=[row["Latitude"], row["Longitude"]],
                    radius=7,
                    color="red",
                    fill=True,
                    fill_opacity=0.7,
                    popup=f"<b>{row['Station']}</b><br>WBGT: {row['WBGT (°F)']:.1f} °F<br>{row['Time']}"
                ).add_to(fmap)

            st_data = st_folium(fmap, width=900, height=600)

else:
    st.info("Select stations and click **Fetch WBGT Data** to begin.")
