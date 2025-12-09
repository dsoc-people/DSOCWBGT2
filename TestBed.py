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
from folium.plugins import FloatImage
from geopy.geocoders import Nominatim
from geopy.extra.rate_limiter import RateLimiter
import math
from streamlit_folium import st_folium
from streamlit_autorefresh import st_autorefresh
import base64
import os

# Page configuration
st.set_page_config(
    page_title="Multi-Source Weather Station WBGT Monitor",
    page_icon="🌡️",
    layout="wide"
)

# --- Helper functions for unit conversion ---
def farenheit_to_celsius(temp_f):
    return (temp_f - 32) * 5 / 9

def celsius_to_farenheit(temp_c):
    return temp_c * 9 / 5 + 32

# --- Wet Bulb Temperature Calculation ---
def dbdp2wb(tempC, dpC, p):
    return (tempC + dpC) / 2

# --- Wet Bulb Globe Temperature Calculation ---
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

    wbgt_f = celsius_to_farenheit(wbgt_c)
    return wbgt_f

# --- Function to safely extract parameter value ---
def extract_value(records, target):
    for r in records:
        if target.lower() in r.get("sensor_name", "").lower():
            return r.get("value")
    return None

@st.cache_data(ttl=300)  # Cache for 5 minutes
def fetch_weatherstem_data():
    """Fetch data from WeatherSTEM stations (White Squirrel Weather)"""
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

    target_sensors = [
        "Wet Bulb Globe Temperature",
        "Thermometer",
        "Dewpoint",
        "Anemometer"
    ]

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
                "source": "White Squirrel Weather"  # Explicitly set source
            })

        except Exception as e:
            data.append({
                "Site": site,
                "Observation Time": "Error",
                "WBGT (°F)": None,
                "Temperature (°F)": None,
                "Dewpoint (°F)": None,
                "Wind Speed (mph)": None,
                "source": "White Squirrel Weather"  # Explicitly set source even on error
            })
            st.warning(f"Error fetching data from {site}: {e}")

    df = pd.DataFrame(data)
    return df[["Site", "Observation Time", "WBGT (°F)", "Temperature (°F)", "Dewpoint (°F)", "Wind Speed (mph)", "source"]]

@st.cache_data(ttl=300)
def fetch_usgs_data():
    """Fetch data from USGS river gauges"""
    usgs_iv_ky_url = 'https://waterservices.usgs.gov/nwis/iv/?format=json&stateCd=ky&siteStatus=active'
    
    try:
        response = requests.get(usgs_iv_ky_url, timeout=30)
        response.raise_for_status()
        ky_sites_iv_json = response.json()

        site_data = []
        if ky_sites_iv_json and 'value' in ky_sites_iv_json and 'timeSeries' in ky_sites_iv_json['value']:
            for series in ky_sites_iv_json['value']['timeSeries']:
                site_info = series['sourceInfo']
                variable_info = series['variable']

                site_id = site_info['siteCode'][0]['value']
                site_name = site_info['siteName']
                latitude = site_info['geoLocation']['geogLocation']['latitude']
                longitude = site_info['geoLocation']['geogLocation']['longitude']
                parameter_code = variable_info['variableCode'][0]['value']
                parameter_name = variable_info['variableName']

                latest_value = 'N/A'
                if series['values'] and series['values'][0]['value']:
                    latest_value = series['values'][0]['value'][-1]['value']

                value_unit = 'N/A'
                if 'unit' in variable_info and variable_info['unit']:
                    if 'unitAbbreviation' in variable_info['unit']:
                        value_unit = variable_info['unit']['unitAbbreviation']
                    elif 'unitCode' in variable_info['unit']:
                        value_unit = variable_info['unit']['unitCode']

                site_data.append({
                    'Site ID': site_id,
                    'Site Name': site_name,
                    'Latitude': latitude,
                    'Longitude': longitude,
                    'Parameter Code': parameter_code,
                    'Parameter Name': parameter_name,
                    'Latest Value': latest_value,
                    'Value Unit': value_unit
                })

        return pd.DataFrame(site_data)
    except Exception as e:
        st.error(f"Error fetching USGS data: {e}")
        return pd.DataFrame()

@st.cache_data(ttl=300)
def process_mesonet_station(station_id, year, station_coords):
    """Process a single Mesonet station and calculate WBGT - always returns a row"""
    lat, lon = station_coords.get(station_id, (None, None))
    
    try:
        manifest_url = f"https://d266k7wxhw6o23.cloudfront.net/data/{station_id}/{year}/manifest.json"
        manifest = requests.get(manifest_url, timeout=15).json()
        if not manifest:
            raise ValueError("Empty manifest")

        latest_day = max(manifest.keys())
        key = manifest[latest_day]["key"]
        data = requests.get(
            f"https://d266k7wxhw6o23.cloudfront.net/{key}",
            timeout=15
        ).json()
        df = pd.DataFrame(data["rows"], columns=data["columns"])
        cols = ["TAIR", "DWPT", "WSPD", "SRAD", "PRES", "UTCTimestampCollected"]
        if not all(c in df.columns for c in cols):
            raise ValueError("Missing required columns")

        # Use the exact approach from the working code
        tair_c, dwpt_c, wspd_mps, srad, pres_hpa = [
            df[c].dropna().iloc[-1] for c in cols[:-1]
        ]
        pres_inhg = pres_hpa * 0.02953
        obs_time = df["UTCTimestampCollected"].dropna().iloc[-1]

        tair_f = celsius_to_farenheit(tair_c)
        dwpt_f = celsius_to_farenheit(dwpt_c)
        wspd_mph = wspd_mps * 2.23694

        wbgt_f = wbgt(tair_f, wspd_mph, srad, pres_inhg, dwpt_f)

        # Explicitly return with all required columns - ensure values are not numpy types
        result = {
            "name": str(station_id),
            "latitude": float(lat) if lat is not None else None,
            "longitude": float(lon) if lon is not None else None,
            "wbgt_f": float(wbgt_f) if not pd.isna(wbgt_f) else None,
            "Temperature (°F)": float(tair_f) if not pd.isna(tair_f) else None,
            "Dewpoint (°F)": float(dwpt_f) if not pd.isna(dwpt_f) else None,
            "Wind Speed (mph)": float(wspd_mph) if not pd.isna(wspd_mph) else None,
            "observation_time": str(obs_time) if obs_time is not None else "N/A",
            "source": "Mesonet",
        }
        return result

    except Exception as e:
        # Always return a row even on error, matching the provided code
        # Log the error for debugging but don't break the app
        import traceback
        error_msg = str(e)
        # Only show first error to avoid spam
        if not hasattr(process_mesonet_station, '_error_shown'):
            st.warning(f"Error processing Mesonet station {station_id}: {error_msg}")
            process_mesonet_station._error_shown = True
        
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

@st.cache_data(ttl=300)
def fetch_mesonet_data(year, station_coords, station_abbreviations):
    """Fetch and process all Mesonet stations - always returns rows even on error"""
    mesonet_data_rows = []
    for station_id in station_abbreviations:
        station_data = process_mesonet_station(station_id, year, station_coords)
        # Always append - function now always returns a row
        if station_data is not None:
            mesonet_data_rows.append(station_data)
    
    if not mesonet_data_rows:
        # Return empty dataframe with correct columns
        return pd.DataFrame(columns=[
            "name", "latitude", "longitude", "wbgt_f", 
            "Temperature (°F)", "Dewpoint (°F)", "Wind Speed (mph)", 
            "observation_time", "source"
        ])
    
    # Create DataFrame and ensure all expected columns exist
    df = pd.DataFrame(mesonet_data_rows)
    
    # Explicitly ensure required columns exist (fill with None if missing)
    required_columns = [
        "name", "latitude", "longitude", "wbgt_f", 
        "Temperature (°F)", "Dewpoint (°F)", "Wind Speed (mph)", 
        "observation_time", "source"
    ]
    for col in required_columns:
        if col not in df.columns:
            df[col] = None
    
    return df

def get_station_coordinates():
    """Parse station coordinates from text"""
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
WKUCHAOS,36.98582726072027,-86.44967208166477
WKUChaos,36.98582632147347,-86.44968800031974
Etown,37.69563805082102,-85.88387790284976
Glasgow,36.9774781561,-85.916651431
WKUIMFields,36.9774781561,-85.9166514315
Owensboro,36.9774781561,-85.9166514315
WKU,36.9774781561,-85.9166514315
"""
    station_coords = {}
    station_abbreviations = []
    for line in station_coords_text.strip().splitlines():
        parts = line.split(',')
        if len(parts) == 3:
            abbrev = parts[0].strip()
            try:
                lat = float(parts[1].strip())
                lon = float(parts[2].strip())
                station_coords[abbrev] = (lat, lon)
                station_abbreviations.append(abbrev)
            except ValueError:
                pass
    return station_coords, station_abbreviations

@st.cache_data(ttl=3600)  # Cache geocoding for 1 hour
def geocode_stations(df_whitesquirrel):
    """Geocode station locations"""
    known_coords = {
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

    df_whitesquirrel = df_whitesquirrel.copy()
    df_whitesquirrel["latitude"] = None
    df_whitesquirrel["longitude"] = None

    geolocator = Nominatim(user_agent="streamlit-wbgt-app")
    geocode = RateLimiter(geolocator.geocode, min_delay_seconds=1, swallow_exceptions=True)

    def locate(site: str):
        if site in known_coords:
            return known_coords[site]
        query = site
        if site.startswith("WKU"):
            query = f"{site}, Western Kentucky University, Bowling Green, KY"
        elif site in {"E'town", "Etown", "Elizabethtown"}:
            query = "Elizabethtown, KY"
        else:
            query = f"{site}, Kentucky"
        loc = geocode(query)
        if loc:
            return (loc.latitude, loc.longitude)
        return (None, None)

    for i, row in df_whitesquirrel.iterrows():
        latlon = locate(row["Site"])
        df_whitesquirrel.loc[i, "latitude"], df_whitesquirrel.loc[i, "longitude"] = latlon

    return df_whitesquirrel

def wbgt_color(w):
    """Color function by WBGT (°F)"""
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

def marker_radius():
    """Marker size function - fixed size for all markers"""
    return 10  # Fixed size for all markers

def create_map(combined_df, ky_sites_df, selected_measurement):
    """Create the interactive Folium map"""
    # Calculate center
    all_lats = []
    all_lons = []

    if not combined_df.empty:
        df_for_center = combined_df.dropna(subset=["latitude", "longitude"])
        if not df_for_center.empty:
            all_lats.extend(df_for_center["latitude"].tolist())
            all_lons.extend(df_for_center["longitude"].tolist())

    if not ky_sites_df.empty:
        df_usgs_center = ky_sites_df.dropna(subset=["Latitude", "Longitude"])
        if not df_usgs_center.empty:
            all_lats.extend(df_usgs_center["Latitude"].tolist())
            all_lons.extend(df_usgs_center["Longitude"].tolist())

    if all_lats and all_lons:
        center = [np.mean(all_lats), np.mean(all_lons)]
    else:
        center = [37.1, -85.9]

    m = folium.Map(location=center, zoom_start=8, control_scale=True)

    # Add tile layers
    folium.TileLayer(
        tiles="https://{s}.basemaps.cartocdn.com/rastertiles/voyager/{z}/{x}/{y}{r}.png",
        attr='&copy; <a href="https://carto.com/attributions">CARTO</a> &copy; OpenStreetMap contributors',
        name="CartoDB Voyager",
        subdomains=["a","b","c","d"],
        overlay=False,
        control=True
    ).add_to(m)

    folium.TileLayer(
        tiles="https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png",
        attr='&copy; <a href="https://carto.com/attributions">CARTO</a> &copy; OpenStreetMap contributors',
        name="CartoDB Dark Matter",
        subdomains=["a","b","c","d"],
        overlay=False,
        control=True
    ).add_to(m)

    folium.TileLayer(
        tiles="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
        attr='Tiles &copy; Esri — Source: Esri, Maxar, Earthstar Geographics, and the GIS User Community',
        name="Esri World Imagery",
        overlay=False,
        control=True
    ).add_to(m)

    folium.TileLayer(
        tiles="https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png",
        attr='Map data: &copy; OpenStreetMap contributors, SRTM | Map style: &copy; OpenTopoMap (CC-BY-SA)',
        name="OpenTopoMap",
        subdomains=["a","b","c"],
        overlay=False,
        control=True
    ).add_to(m)

    folium.TileLayer(
        tiles="https://{s}.basemaps.cartocdn.com/rastertiles/light_all/{z}/{x}/{y}{r}.png",
        attr='&copy; <a href="https://carto.com/attributions">CARTO</a> &copy; OpenStreetMap contributors',
        name="CartoDB Positron",
        subdomains=["a","b","c","d"],
        overlay=False,
        control=True
    ).add_to(m)

    # Create a single layer group (no layer control - using sidebar checkboxes instead)
    weather_layer = folium.FeatureGroup(name='Weather Stations')
    usgs_layer = folium.FeatureGroup(name='USGS River Gauges')
    weather_layer.add_to(m)
    usgs_layer.add_to(m)

    # Add markers for weather data
    if not combined_df.empty:
        for index, row in combined_df.iterrows():
            wbgt = row.get("wbgt_f")
            temp = row.get("Temperature (°F)")
            dew = row.get("Dewpoint (°F)")
            site = row.get("name") if "name" in row else row.get("Site")
            obs_time = row.get("observation_time") if "observation_time" in row else row.get("Observation Time")
            source = row.get("source")
            latitude = row.get("latitude")
            longitude = row.get("longitude")

            if pd.isna(latitude) or pd.isna(longitude):
                continue

            # Determine color and value based on selected measurement
            if selected_measurement == "WBGT":
                value = wbgt
                color = wbgt_color(wbgt)
                value_label = "WBGT"
                value_unit = "°F"
            elif selected_measurement == "Temperature":
                value = temp
                color = wbgt_color(wbgt) if pd.notna(wbgt) else "#808080"  # Still use WBGT for color coding
                value_label = "Temperature"
                value_unit = "°F"
            elif selected_measurement == "Dewpoint":
                value = dew
                color = wbgt_color(wbgt) if pd.notna(wbgt) else "#808080"  # Still use WBGT for color coding
                value_label = "Dewpoint"
                value_unit = "°F"
            else:
                value = wbgt
                color = wbgt_color(wbgt)
                value_label = "WBGT"
                value_unit = "°F"

            radius = marker_radius()

            popup_html = f"""
            <div style="font-family:system-ui;min-width:220px">
              <h4 style="margin:0 0 6px 0">{site} ({source})</h4>
              <div><b>Observed:</b> {obs_time if pd.notna(obs_time) else 'N/A'}</div>
              <div><b>WBGT:</b> {f'{wbgt:.1f} °F' if pd.notna(wbgt) else 'N/A'}</div>
              <div><b>Temperature:</b> {f'{temp:.1f} °F' if pd.notna(temp) else 'N/A'}</div>
              <div><b>Dewpoint:</b> {f'{dew:.1f} °F' if pd.notna(dew) else 'N/A'}</div>
            </div>
            """

            tooltip_text = f"{site}: {value_label} {value:.1f}{value_unit}" if pd.notna(value) else f"{site}: {value_label} N/A"
            
            marker = folium.CircleMarker(
                location=[latitude, longitude],
                radius=radius,
                popup=folium.Popup(popup_html, max_width=280),
                tooltip=tooltip_text,
                color=color,
                fill=True,
                fill_color=color,
                fill_opacity=0.8,
                weight=1,
            )

            marker.add_to(weather_layer)

    # Add USGS river gauge markers
    if not ky_sites_df.empty:
        for site_id in ky_sites_df['Site ID'].unique():
            site_rows = ky_sites_df[ky_sites_df['Site ID'] == site_id]
            first_row = site_rows.iloc[0]

            latitude = first_row['Latitude']
            longitude = first_row['Longitude']
            site_name = first_row['Site Name']

            popup_html = f"<div style='font-family:system-ui;min-width:220px'><h4 style='margin:0 0 6px 0'>{site_name}</h4>"
            popup_html += f"<div><b>Site ID:</b> {site_id}</div>"

            for idx, row in site_rows.iterrows():
                param_name = row['Parameter Name']
                latest_value = row['Latest Value']
                value_unit = row['Value Unit']
                if latest_value != 'N/A':
                    popup_html += f"<div><b>{param_name}:</b> {latest_value} {value_unit}</div>"

            popup_html += "</div>"

            folium.Marker(
                location=[latitude, longitude],
                popup=folium.Popup(popup_html, max_width=280),
                tooltip=f"{site_name} ({site_id})",
                icon=folium.Icon(color='blue', icon='tint', prefix='fa')
            ).add_to(usgs_layer)

    # No layer control - using sidebar checkboxes instead

    # Add legend
    legend_html = f"""
    <div style="
        position: fixed;
        bottom: 30px; left: 30px; z-index: 9999;
        background: white; padding: 10px 12px;
        border: 1px solid #ccc; border-radius: 8px;
        font-family: system-ui; font-size: 12px;">
      <div style="font-weight:600; margin-bottom:6px;">WBGT (°F) Color Scale</div>
      <div><span style="display:inline-block;width:12px;height:12px;background:#2ca02c;margin-right:6px;border:1px solid #999;"></span><60</div>
      <div><span style="display:inline-block;width:12px;height:12px;background:#ff7f0e;margin-right:6px;border:1px solid #999;"></span>60–69.9</div>
      <div><span style="display:inline-block;width:12px;height:12px;background:#d62728;margin-right:6px;border:1px solid #999;"></span>70–79.9</div>
      <div><span style="display:inline-block;width:12px;height:12px;background:#800026;margin-right:6px;border:1px solid #999;"></span>≥80</div>
      <div><span style="display:inline-block;width:12px;height:12px;background:#808080;margin-right:6px;border:1px solid #999;"></span>N/A</div>
      <div style="margin-top:6px;"><em>Displaying: {selected_measurement}</em></div>
    </div>
    """
    m.get_root().html.add_child(folium.Element(legend_html))

    return m

def upload_map_to_github(map_filename):
    """Upload map HTML file to GitHub repository"""
    try:
        # Read the map file
        if not os.path.exists(map_filename):
            return {"success": False, "error": f"Map file {map_filename} not found"}
        
        with open(map_filename, 'r', encoding='utf-8') as f:
            map_content = f.read()
        
        # GitHub repository details
        repo_owner = "dsoc-people"
        repo_name = "DSOCWBGT2"
        file_path = "Image storage/wbgt_map.html"  # Path in GitHub repo
        
        # Try to get token from Streamlit secrets first, then environment variable
        try:
            github_token = st.secrets.get("GITHUB_TOKEN", None)
        except:
            github_token = None
        
        if not github_token:
            github_token = os.getenv("GITHUB_TOKEN")
        
        if not github_token:
            # Try to use git command as fallback
            import subprocess
            try:
                # Check if we're in a git repo
                result = subprocess.run(
                    ["git", "remote", "get-url", "origin"],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                if result.returncode == 0:
                    # Use git to commit and push
                    subprocess.run(["git", "add", map_filename], check=False, timeout=5)
                    commit_result = subprocess.run(
                        ["git", "commit", "-m", "Update WBGT map", map_filename],
                        capture_output=True,
                        text=True,
                        timeout=5
                    )
                    push_result = subprocess.run(
                        ["git", "push"],
                        capture_output=True,
                        text=True,
                        timeout=10
                    )
                    if push_result.returncode == 0:
                        return {
                            "success": True,
                            "method": "git",
                            "url": f"https://github.com/{repo_owner}/{repo_name}/blob/main/{file_path}"
                        }
            except Exception as git_error:
                return {
                    "success": False,
                    "error": f"GITHUB_TOKEN not set and git push failed: {str(git_error)}"
                }
        
        # Use GitHub API
        import requests
        
        # Get the file SHA if it exists
        url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/contents/{file_path}"
        headers = {
            "Authorization": f"token {github_token}",
            "Accept": "application/vnd.github.v3+json"
        }
        
        # Check if file exists
        response = requests.get(url, headers=headers, timeout=10)
        sha = None
        if response.status_code == 200:
            sha = response.json().get("sha")
        elif response.status_code == 404:
            # File doesn't exist yet, that's okay
            pass
        else:
            return {
                "success": False,
                "error": f"GitHub API error checking file: {response.status_code} - {response.text}"
            }
        
        # Encode content
        content_encoded = base64.b64encode(map_content.encode('utf-8')).decode('utf-8')
        
        # Prepare data
        data = {
            "message": "Update WBGT map",
            "content": content_encoded,
            "branch": "main"
        }
        if sha:
            data["sha"] = sha
        
        # Upload file
        response = requests.put(url, json=data, headers=headers, timeout=10)
        
        if response.status_code in [200, 201]:
            return {
                "success": True,
                "method": "api",
                "url": f"https://github.com/{repo_owner}/{repo_name}/blob/main/{file_path}"
            }
        else:
            return {
                "success": False,
                "error": f"GitHub API error: {response.status_code} - {response.text}"
            }
        
    except Exception as e:
        return {
            "success": False,
            "error": f"Upload failed: {str(e)}"
        }

# Main Streamlit app
def main():
    st.title("🌡️ Multi-Source Weather Station WBGT Monitor")
    st.markdown("**Wet Bulb Globe Temperature (WBGT) monitoring from multiple data sources**")

    # Sidebar controls
    with st.sidebar:
        st.header("Data Sources")
        show_whitesquirrel = st.checkbox("White Squirrel Weather", value=True)
        show_mesonet = st.checkbox("Mesonet", value=True)
        show_usgs = st.checkbox("USGS River Gauges", value=True)
        
        st.header("Settings")
        year = st.selectbox("Year", ["2025", "2024", "2023"], index=0)
        selected_measurement = st.selectbox(
            "Display Measurement",
            ["WBGT", "Temperature", "Dewpoint"],
            index=0,
            help="Select which measurement to display on the map markers"
        )
        
        st.info("🔄 Auto-refreshing every 5 minutes")
        # Always run auto-refresh every 5 minutes (300 seconds = 300000 milliseconds)
        count = st_autorefresh(interval=300000, limit=None, key="weather_refresh")
        if count > 0:
            st.cache_data.clear()
            st.rerun()

    # Fetch data
    with st.spinner("Fetching weather data..."):
        # Fetch White Squirrel Weather data
        if show_whitesquirrel:
            df_whitesquirrel = fetch_weatherstem_data()
            df_whitesquirrel = geocode_stations(df_whitesquirrel)
            # Rename columns to match Mesonet structure
            df_whitesquirrel = df_whitesquirrel.rename(columns={
                "WBGT (°F)": "wbgt_f", 
                "Site": "name",
                "Observation Time": "observation_time"
            })
            # Ensure source is explicitly set
            df_whitesquirrel['source'] = 'White Squirrel Weather'
        else:
            df_whitesquirrel = pd.DataFrame()

        # Fetch Mesonet data
        if show_mesonet:
            station_coords, station_abbreviations = get_station_coordinates()
            df_mesonet = fetch_mesonet_data(year, station_coords, station_abbreviations)
            # Source is already set in process_mesonet_station function
            # Verify Temperature and Dewpoint columns exist
            if not df_mesonet.empty:
                if "Temperature (°F)" not in df_mesonet.columns:
                    st.error(f"❌ Mesonet dataframe missing 'Temperature (°F)' column! Columns: {list(df_mesonet.columns)}")
                if "Dewpoint (°F)" not in df_mesonet.columns:
                    st.error(f"❌ Mesonet dataframe missing 'Dewpoint (°F)' column! Columns: {list(df_mesonet.columns)}")
                # Debug: Show first row structure
                if len(df_mesonet) > 0:
                    first_row = df_mesonet.iloc[0].to_dict()
                    st.write("🔍 First Mesonet row keys:", list(first_row.keys()))
        else:
            df_mesonet = pd.DataFrame()

        # Fetch USGS data
        if show_usgs:
            ky_sites_df = fetch_usgs_data()
        else:
            ky_sites_df = pd.DataFrame()

    # Combine dataframes based on selected sources
    dfs_to_combine = []
    if show_whitesquirrel and not df_whitesquirrel.empty:
        dfs_to_combine.append(df_whitesquirrel.copy())
    if show_mesonet and not df_mesonet.empty:
        dfs_to_combine.append(df_mesonet.copy())
    
    if len(dfs_to_combine) > 1:
        # Align columns
        all_columns = list(set([col for df in dfs_to_combine for col in df.columns.tolist()]))
        aligned_dfs = []
        for df in dfs_to_combine:
            df_aligned = df.copy()
            for col in all_columns:
                if col not in df_aligned.columns:
                    df_aligned[col] = None
            df_aligned = df_aligned[all_columns]
            aligned_dfs.append(df_aligned)
        combined_df = pd.concat(aligned_dfs, ignore_index=True)
    elif len(dfs_to_combine) == 1:
        combined_df = dfs_to_combine[0].copy()
    else:
        combined_df = pd.DataFrame()
    
    # Filter USGS data based on checkbox
    if not show_usgs:
        ky_sites_df = pd.DataFrame()

    # Display data summary
    col1, col2, col3 = st.columns(3)
    with col1:
        if not df_whitesquirrel.empty:
            st.metric("White Squirrel Weather Stations", len(df_whitesquirrel))
        else:
            st.metric("White Squirrel Weather Stations", 0)
    with col2:
        if not df_mesonet.empty:
            st.metric("Mesonet Stations", len(df_mesonet))
        else:
            st.metric("Mesonet Stations", 0)
    with col3:
        if not ky_sites_df.empty:
            st.metric("USGS River Gauges", ky_sites_df['Site ID'].nunique())
        else:
            st.metric("USGS River Gauges", 0)

    # Display data tables
    if st.checkbox("Show Data Tables"):
        tab1, tab2, tab3 = st.tabs(["White Squirrel Weather", "Mesonet", "USGS"])
        
        with tab1:
            if not df_whitesquirrel.empty:
                st.dataframe(df_whitesquirrel, use_container_width=True)
            else:
                st.info("No White Squirrel Weather data available.")
        
        with tab2:
            if not df_mesonet.empty:
                st.dataframe(df_mesonet, use_container_width=True)
            else:
                st.info("No Mesonet data available.")
        
        with tab3:
            if not ky_sites_df.empty:
                st.dataframe(ky_sites_df, use_container_width=True)
            else:
                st.info("No USGS data available.")

    # Create and display map
    if not combined_df.empty or not ky_sites_df.empty:
        st.subheader("Interactive Map")
        m = create_map(combined_df, ky_sites_df, selected_measurement)
        
        # Save map to local file
        map_filename = "wbgt_map.html"
        m.save(map_filename)
        
        # Upload to GitHub
        upload_status = upload_map_to_github(map_filename)
        if upload_status.get("success"):
            st.success(f"✅ Map uploaded to GitHub: {upload_status.get('url', '')}")
        else:
            st.warning(f"⚠️ Could not upload map to GitHub: {upload_status.get('error', 'Unknown error')}")
            with st.expander("Setup Instructions"):
                st.markdown("""
                **To enable GitHub upload:**
                
                1. **Create a GitHub Personal Access Token:**
                   - Go to https://github.com/settings/tokens
                   - Click "Generate new token (classic)"
                   - Give it a name like "WBGT Map Upload"
                   - Select scope: `repo` (full control of private repositories)
                   - Click "Generate token"
                   - Copy the token immediately (you won't see it again!)
                
                2. **Set the token as an environment variable:**
                   ```bash
                   export GITHUB_TOKEN=your_token_here
                   ```
                
                3. **Or set it in your Streamlit secrets:**
                   Create `.streamlit/secrets.toml`:
                   ```toml
                   GITHUB_TOKEN = "your_token_here"
                   ```
                
                4. **Restart your Streamlit app**
                """)
        
        st_folium(m, width=None, height=600)
    else:
        st.warning("No data available to display on map. Please enable at least one data source.")

    # Debug info (can be removed in production)
    if st.checkbox("Show Debug Info"):
        st.subheader("Source Verification")
        if not combined_df.empty:
            st.write("Sources in combined dataframe:")
            st.write(combined_df['source'].value_counts())
            st.write("**All columns in combined dataframe:**")
            st.write(list(combined_df.columns))
            st.write("**Sample rows from combined dataframe:**")
            st.dataframe(combined_df.head(10))
        
        st.subheader("Mesonet Data Debug")
        if not df_mesonet.empty:
            st.write("**Mesonet dataframe columns:**")
            st.write(list(df_mesonet.columns))
            st.write("**Mesonet dataframe shape:**")
            st.write(df_mesonet.shape)
            st.write("**Sample Mesonet rows:**")
            st.dataframe(df_mesonet.head(10))
            st.write("**Checking for Temperature and Dewpoint columns:**")
            st.write(f"Has 'Temperature (°F)': {'Temperature (°F)' in df_mesonet.columns}")
            st.write(f"Has 'Dewpoint (°F)': {'Dewpoint (°F)' in df_mesonet.columns}")
            if 'Temperature (°F)' in df_mesonet.columns:
                st.write("**Temperature values:**")
                st.write(df_mesonet['Temperature (°F)'].head(10))
            if 'Dewpoint (°F)' in df_mesonet.columns:
                st.write("**Dewpoint values:**")
                st.write(df_mesonet['Dewpoint (°F)'].head(10))
        else:
            st.write("Mesonet dataframe is empty")
        
        st.subheader("White Squirrel Weather Data Debug")
        if not df_whitesquirrel.empty:
            st.write("**White Squirrel dataframe columns:**")
            st.write(list(df_whitesquirrel.columns))
            st.write("**Sample White Squirrel rows:**")
            st.dataframe(df_whitesquirrel.head(5))

if __name__ == "__main__":
    main()

