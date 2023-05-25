#!/usr/bin/env python3

# ----------------------------------------------------------------- #
# EXTRACT NOAA WEATHER DATA (ISD) TO MERGE WITH PANEL DATA          #
# ----------------------------------------------------------------- #
# Hans Elliott
# ----------------------------------------------------------------- #
# 1. Extract list of all US weather stations from ISD metadata:
#   - https://www.ncei.noaa.gov/pub/data/noaa/isd-history.txt
# 2. Extract raw datasets from: https://www.ncei.noaa.gov/data/global-hourly/
#   - download full tar.gz archives
#   - extract the list of filenames and then uncompress only the US stations
# 3. For each year, aggregate into a station-day longitudinal dataset
#   - summarize each var. For ex, average cloud cover for station i on day t
# 4. Save longitudinal dataset for each year.
#
# Variables:
# We mostly want "sky-cover-layer", ie GA1, and the "coverage code" specifically
#   - variables are stored as comma separated lists, ex: GA1='02,1,30000,00,1'
#   - the coverage code is the 1st element:
#     "the fraction of the total celestial dome covered", in oktas
#   - the quality code is the 2nd element: 
#     it describes the quality of this observation
#   - note on the unit of measurement, okta: https://en.wikipedia.org/wiki/Okta
#   - GA1 is not reported in all station-year datasets, so we skip csvs which do not have it
# We also extract temperature and dewpoint which are guaranteed to be in every CSV
#   - again, we are interested in the 1st and 2nd fields, the value and quality code respectively
# Variable descriptions: 
#   - https://www.ncei.noaa.gov/data/global-hourly/doc/isd-format-document.pdf
#   - sky-cover-layer (GA1) is on pg. 54
#   - temperature (TMP) is on pg. 10
#   - dewpoint temp (DEW) is on pg. 11
import os
import re
import time
import shutil
import tarfile
# install
import requests #http requests
import numpy as np #nd arrays
import pandas as pd #dataframes
from pyHere import Here #filepath management
# options
verbose = True # print lots of info? 
start_year = 2003 # first year to extract data for
end_year = 2023 # last year to extract data for
here = Here("nielsen-research") # project root


class Timer:
    """Simple timer to help comply with web etiquette
    """
    def __init__(self) -> None:
        self.t0 = time.time()
    def restart(self):
        self.t0 = time.time()
    def et(self):
        return time.time() - self.t0
    def sleep_til(self, elaps_sec):
        et = self.et()
        if et < elaps_sec:
            time.sleep(elaps_sec - et)


# ============================================================================ #
# GENERATE LIST OF US WEATHER STATIONS FROM ISD METADATA
# ============================================================================ #

def _get_next_idx(line):
    llen = len(line)
    regex = re.search(r"\s[+-]\d", line)
    if regex:
        return regex.start(0)
    else:
        return -1
    
def parse_lat_lon(line):
    idx = _get_next_idx(line)
    coords = []
    while idx != -1:
        coords.append( line[idx:].split()[0] )
        line = ' ' + ' '.join(line[idx:].split()[1:])
        idx = _get_next_idx(line)
    
    keys = ["lat", "lon", "elev_m"]
    out = {keys[i] : [coords[i]] for i in range(len(coords))}
    
    return out

def get_us_stations():
    """Extract ISD Station Metadata. Generate US Station Dataset.
    Saves CSV containing the US stations and returns the pandas dataframe.
    Dataset: station (identifier), longitude, latitude, elevation (meters)
    """
    us_stations_path = here.here("weather/data/us_stations.csv")
    if os.path.exists(us_stations_path):
        return pd.read_csv(us_stations_path, dtype=str)
    
    headers = {"user-agent": 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/112.0.0.0 Safari/537.36'}
    url = "https://www.ncei.noaa.gov/pub/data/noaa/isd-history.txt"
    r = requests.get(url, headers=headers)
    with open(here.here("weather/data", "isd-history.txt"), "w") as f:
        f.write(r.text) #save raw txt as backup
        
    # Parse station id and coordinates from text file
    dfs = []
    for line in r.text.splitlines():
        if " US " in line and line.strip() != "ST = State for US stations":
            # extract the usaf station id and the ncdc wban id
            station = ''.join(line.split()[0:2])
            # extract the lat and lon values
            coords = parse_lat_lon(line)
            coords["station"] = [station]
            dfs.append(pd.DataFrame(coords))
    us_stations = pd.concat(dfs).reset_index(drop=True)
    us_stations.to_csv(us_stations_path, index=False)
    
    return us_stations



# ============================================================================ #
#  EXTRACT & AGGREGATE NOAA WEATHER DATA FOR RELEVANT STATIONS
# ============================================================================ #

def clean_tmp(path):
    if os.path.exists(path):
        shutil.rmtree(path)

def download_raw_isd(tmp_dir, year, stations, verbose=False):
    """Download compressed tar.gz file from NOAA containing all stations for the
    given 'year'. Extract csv files from the tar file into 'tmp_dir' for all 
    stations that are also in the 'stations' iterable.
    """
    # Download compressed data
    url = f"https://www.ncei.noaa.gov/data/global-hourly/archive/csv/{year}.tar.gz"
    name = url.split("/")[-1] 
    if verbose: print(f"- Downloading {name} from NOAA archive.")
    tar_path = here.here("weather", name)
    r = requests.get(url, stream=True)
    if r.ok:
        # stream compressed data into local tarfile
        with open(tar_path, "wb") as f:
            f.write(r.raw.read())
    
    # Extract CSVs, but only for the stations we care about
    if verbose: print("- Extracting members.")
    clean_tmp(tmp_dir)
    os.makedirs(tmp_dir, exist_ok=False)
    with tarfile.open(tar_path, mode="r:gz") as tarf:
        members = tarf.getmembers()
        extract = set(m.name.split(".")[0] for m in members).intersection(
            set(stations))
        extract = [f+".csv" for f in extract]
        tarf.extractall(tmp_dir, 
                        members=[m for m in members if m.name in extract])
    # Remove tar file
    os.remove(tar_path)


def extract_station_data(csv):
    """Given a path to a CSV which contains the data for one weather station for 
    one year, processes the data and returns a dataset that is summarized at the 
    daily level. 
    Returns 'None' if the station's data is missing the key variable (GA1).
    See NOAA's data documentation for variable details: 
    https://www.ncei.noaa.gov/data/global-hourly/doc/isd-format-document.pdf
    - GA1/sky-cover-layer: pg. 54
    """ 
    d = pd.read_csv(csv, dtype=str)
    d["day"] = d["DATE"].str.split("T").str[0]
    
    if "GA1" not in d.columns:
        # We really care about sky cover, so if not available don't waste time.
        # The other vars are guaranteed to be in all files. 
        return None

    # parse sky cloud cover (oktas)
    d["sky_cover_layer"] = d["GA1"].str.split(",").str[0].apply(
        pd.to_numeric, errors="coerce").apply(
            lambda x: np.nan if x > 8 else x)
    d["quality_code"] = d["GA1"].str.split(",").str[1].apply(
        pd.to_numeric, errors="coerce")
    d.loc[d.quality_code.isin((3., 7.)), "sky_cover_layer"] = np.nan
    
    # parse temperature (celsius)
    d["temp_c"] = d["TMP"].str.split(",").str[0].apply(
        pd.to_numeric, errors="coerce")
    d["quality_code"] = d["TMP"].str.split(",").str[1].apply(
        pd.to_numeric, errors="coerce")
    d.loc[d.quality_code.isin((3., 7.)), "temp_c"] = np.nan

    # parse dew point temp (celsius)
    d["dewpoint_c"] = d["DEW"].str.split(",").str[0].apply(
        pd.to_numeric, errors="coerce")
    d["quality_code"] = d["DEW"].str.split(",").str[1].apply(
        pd.to_numeric, errors="coerce")
    d.loc[d.quality_code.isin((3., 7.)), "dewpoint_c"] = np.nan

    # summarize at daily level
    d = d.groupby(["STATION", "LONGITUDE", "LATITUDE", "ELEVATION", "NAME", "day"], as_index=False).agg({
        "sky_cover_layer" : "mean",
        "temp_c" : "mean",
        "dewpoint_c" : "mean",
        "DATE" : "count" #use 'DATE' to count the n of hrs reported per day
    })
    d.columns = [c.lower() for c in d.columns]
    d.rename(columns={"date" : "n_hourly_obs"}, inplace=True)
    return d 


def collect_isd_year(csv_dir, station_set, year, verbose=False):
    """Download and process one year of global-hourly-isd data.
    """
    assert len(str(year)) == 4, "Provide a 4-digit year, like: year=2010"
    # download raw data files
    download_raw_isd(tmp_dir=csv_dir,
                     year=year, stations=station_set, verbose=verbose)
    file_list = os.listdir(csv_dir)
    n_files = len(file_list)
    if verbose: print(f"Found {n_files :,} relevant weather station files.")
    # iterate through downloaded files and process them into one dataset
    # - each row contains weather features for one station-date
    dfs = []
    t = Timer()
    for i, file in enumerate(file_list):
        if verbose: print(f"[{i+1}/{n_files}] {year} {file} (et={t.et() :.2f}s)")
        fp = os.path.join(csv_dir, file)
        d = extract_station_data(fp)
        if d is not None:
            dfs.append(d)

    isd = pd.concat(dfs, axis=0).reset_index(drop=True)
    isd.to_csv(here.here(f"weather/data/isd_{year}.csv"),
               index=False)
    clean_tmp(csv_dir) # rm tmp files


# MAIN
def main(start_year, end_year, verbose=False):
    os.makedirs(here.here("weather/data"), exist_ok=True) 
    # US Stations: list of weather stations in the ISD data
    us_stations = get_us_stations()
    
    # Scrape NOAA ISD for all of US the weather station data
    station_set = set(us_stations.station)

    tmp_dir = here.here("weather/tmp")
    for year in range(start_year, end_year+1):
        if verbose: print(year)
        collect_isd_year(csv_dir=tmp_dir, station_set=station_set, year=year,
                         verbose=verbose)

if __name__=="__main__":
    main(start_year, end_year, verbose)

