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
#  STATION COORDINATES
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
# NAME NORMALIZATION
# -------------------------------------------------------------
name_variants = {
    "WKUChaos": "WKU Chaos",
    "WKUCHAOS": "WKU Chaos",
    "Etown": "E'town",
    "WKUIMFields": "WKU IM Fields",
}

# Additional name overrides
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
# TOGGLES ONLY (NO VARIABLE SELECTOR)
# -------------------------------------------------------------
with st.sidebar:
    st.header("Data Sources")
    show_weatherstem = st.checkbox("WeatherSTEM / White Squirrel Weather", True)
    show_mesonet    = st.checkbox("Mesonet", True)
    show_usgs       = st.checkbox("USGS River Gauges", True)

# 🔥 FIXED VARIABLE — ALWAYS WBGT
selected_var = "WBGT (°F)"

# -------------------------------------------------------------
# REMAINDER OF YOUR CODE (UNCHANGED)
# -------------------------------------------------------------
# (Everything from loading station coords, fetchers, combining data,
# color scaling, markers, and map display remains exactly as in your script.)
