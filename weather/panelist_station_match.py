#!/usr/bin/env python
import os
import math
import time
import sqlite3
import warnings
from tqdm import trange
import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree 
from pyHere import Here
# ============================================================================
# CONSTS
# ============================================================================
N_NEIGHBORS = 300 # no. of nearest stations to search for, per panelist zip code
START_YEAR = 2004 # first year of panelist data to match with
END_YEAR = 2020 # last year of panelist data to match with
VERBOSE = True # print info
here = Here("nielsen-research")
# INPUTS
## data matching panelist code to zip code
panelist_zipcodes_path = here.here("panelists/data/panelist_zips.tsv")
## data matching zip code to (lon,lat) centroid
zipcode_centroids_path = here.here("geographic/crosswalks/data",
                                   "zip_centroids.csv")
## path to df w/ all us weather stations and their coords
us_stations_path = here.here("weather/data/us_stations.csv")
## db containing scraped noaa isd data
isd_db_path = here.here("weather/isd.db")
## dir containing panelist cross-sections, with all possible zip code-dates
panelist_cs_dir = "/gpfs/projects/lcb/shared/cross-section/" 

# OUTPUTS
## file where zip-station nearest-neighbors match will be stored
zip_station_match = here.here("weather/data",
                                        f"zipcode_station_{N_NEIGHBORS}nn.csv")
## directory where annual zip-day-with-weather files will be stored
zip_date_station_match_dir = here.here("weather/data/zip_day_match")

# ============================================================================
# HELPERS
# ============================================================================
def _pad_id(id, length):
    # add leading 0s to identifier to reach str length of 'length'
    return "0"*(length - len(str(id))) + str(id)

def query_isd(year):
    """Query the NOAA ISD database for weather station data.
    Returns year-1, year, and year+1's data, with some preprocessing.
    - Each panelist cross-section will contain dates that include all of 'year'
      and likely some of the end of 'year-1' and beginning of 'year+1' (at most,
      a few days on either end).
    """
    conn = sqlite3.connect(isd_db_path)
    if conn:
        isd = []
        for yr in [year-1, year, year+1]:
            d = pd.read_sql_query(
                sql="SELECT * FROM isd_" + str(yr),
                con=conn 
            )
            if yr == year-1:
                d["month"] = pd.to_datetime(d["day"]).dt.month
                d = d.loc[d.month > 10, :].reset_index(drop=True)
                d.drop(labels=["month"], axis=1, inplace=True)
            elif yr == year+1:
                d["month"] = pd.to_datetime(d["day"]).dt.month
                d = d.loc[d.month < 3, :].reset_index(drop=True)
                d.drop(labels=["month"], axis=1, inplace=True)
            isd.append(d)
        # ---
        isd = pd.concat(isd, axis=0).drop_duplicates().reset_index(drop=True)
        # temp vars need to be rescaled
        isd["temp_c"] = isd["temp_c"].astype(float) / 10
        isd["dewpoint_c"] = isd["dewpoint_c"].astype(float) / 10
        isd.rename(columns={
            "longitude" : "station_lon",
            "latitude" : "station_lat",
            "elevation" : "station_elev_m",
        }, inplace=True)
        conn.close()
    else:
        raise SystemExit("Database connection failed.")
        
    return isd

# ============================================================================
# MATCH ZIP CODES TO THEIR K NEAREST WEATHER STATIONS
# ============================================================================
# To simplify downstream searches, we can match our relatively small set of zip
# codes to a large number (K) of weather stations, using a nearest-neighbor
# algorithm.
# I use the panelist masterfile to determine the set of all zip codes which
# appear in the data, and then match these zips to their closest stations.
# This is much less computationally demanding compared to trying to match
# zip-dates to closest station-dates. 

def balltree_haversine_nn(A, B, k=1):
    """For each point in A, query B for the k closest neighbors in B, measured 
    by haversine distance.
    A and B must be shape (n_samples, 2) - each point should be a 
    (latitude, longitude) pair in radians, for compatibility with 
    sklearn.metrics.parwise.haversine_distances
    """
    # ball-trees are slower than kd-trees in low (2) dimensions, but we can only 
    # use haversine as the metric w/ sklearn's ball-tree, and the performance is 
    # still decent.
    msg = "Inputs must have shape: (n_samples, 2)"
    assert A.shape[1] == 2, msg
    assert B.shape[1] == 2, msg
    b_tree = BallTree(B, leaf_size=16, metric="haversine")
    dist, idx = b_tree.query(A, k=k) 
    return dist, idx


def match_zip_to_station(us_stations, k, verb):
    """Uses a K-nearest-neighbor search to find the K closest weather stations
    to each zip code (centroid), based on haversine distance.
    """
    if os.path.exists(zip_station_match):
        if verb: print(f"Loading zip-station knn match for k={k}.")
        return pd.read_csv(zip_station_match, dtype=str)
    
    # get list of all zip codes which appear in panelist data
    panel_zips = pd.read_csv(panelist_zipcodes_path,
                             usecols=["zip"],
                             sep="\t", dtype=str
                             ).drop_duplicates().reset_index(drop=True)
    # import zip code to centroid coordinate mapping & merge w/ panelist zips
    zips = pd.read_csv(zipcode_centroids_path, dtype=str)
    zips["zip"] = zips.zip.apply(lambda x: _pad_id(x, length=5))
    panel_zips["zip"] = panel_zips.zip.apply(lambda x: _pad_id(x, length=5))
    panel_zips = panel_zips.merge(zips, on="zip", how="left")
    # 7 panelist zip codes do not match to zips in the zip-to-coordinate df
    # we can try and impute those later, but for now we just drop them    
    panel_zips.dropna(axis=0, how="any", subset=["latitude", "longitude"], inplace=True)
    panel_zips.reset_index(drop=True, inplace=True)

    # prep stations (drop stations w/out lat lon coords)
    us_stations.dropna(axis=0, how="any", subset=["lon", "lat"], inplace=True)
    us_stations.reset_index(drop=True, inplace=True)

    # Create arrays for nearest neighbor search
    panel_arr = panel_zips.loc[:, ["latitude", "longitude"]].astype(float)
    panel_arr["latitude"] = panel_arr["latitude"].apply(math.radians)
    panel_arr["longitude"] = panel_arr["longitude"].apply(math.radians)
    panel_arr = panel_arr.to_numpy()

    station_arr = us_stations.loc[:, ["lat", "lon"]].astype(float)
    station_arr["lat"] = station_arr["lat"].apply(math.radians)
    station_arr["lon"] = station_arr["lon"].apply(math.radians)
    station_arr = station_arr.to_numpy()

    # NN Search - query stations for the k NNs to each panelist, using
    # haversine distance as the distance metric
    if verb: print(f"Computing k nearest-neighbor search for k={k}")
    t0 = time.time()
    dist, idx = balltree_haversine_nn(A=panel_arr, B=station_arr, k=k)
    if verb: print(f"elapsed time={time.time() - t0 :.3f}s")

    station_ids = us_stations.loc[:, "station"].to_numpy()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=pd.errors.PerformanceWarning)
        for nn in range(k):
            panel_zips.loc[:,
                "nearest_station_" + str(nn + 1)] = station_ids[idx[:, nn]]
            panel_zips.loc[:,
                "distance_km_" + str(nn + 1)] = dist[:, nn] * 6371000/1000

    panel_zips.to_csv(zip_station_match, index=False)
    return panel_zips


# ============================================================================
# MATCH ZIP CODE - DATES TO CLOSEST STATION WITH AVAILABLE DATA
# ============================================================================
# For each zip-date, we need to know the closest of the zip's K closest stations
# which actually has weather data recorded for the given day.
#
def zip_day_closest(zip_station, year, max_neighbors, verb=False):
    """Using the zip code-to-weather station KNN match as a list of stations to
    search for each zip code, determine the closest station which has data 
    available for each zip code-date that appears in our panelist data. 
    """
    fp = os.path.join(panelist_cs_dir, f"cross_section_{year}.tsv")
    # prepare cross-section which contains the zipcode-date "keys"
    cs = pd.read_csv(fp, sep="\t", dtype=str,
                     usecols=["Panelist_ZipCd", "purchase_date"]
                     ).drop_duplicates()
    cs.rename(columns={"purchase_date" : "day", "Panelist_ZipCd" : "zip"},
              inplace=True)
    cs.dropna(axis=0, how="any", inplace=True)
    cs.reset_index(drop=True, inplace=True)
    n_rows = cs.shape[0]

    # load in weather data (noaa isd) and isolate unique station-day pairs
    # (loads data for year and year-1, since cs includes dates in both)
    isd = query_isd(year)
    station_day = isd.loc[:, ["station", "day"]]

    # iteratively match to closest available station
    if verb: print("- Finding zip-day matches.")
    matches = []
    for i in trange(0, max_neighbors):
        near_col = "nearest_station_" + str(i + 1)
        dist_col = "distance_km_" + str(i + 1)
        
        # extract id of next nearest station for each zip code, merge to cs
        nearest_i = zip_station.loc[:, ["zip", near_col]].rename(
            columns={near_col : "station"})
        cs = cs.merge(nearest_i, how="left", on="zip")
        
        # merge on station date indicator to determine which have avail. data
        station_day["nearest_idx"] = i + 1
        cs = cs.merge(station_day, how="left", on=["station", "day"])

        # extract successful matches & add other info about station        
        cs_s = cs.dropna(axis=0, subset=["nearest_idx"]).reset_index(drop=True)
        if not cs_s.empty:
            cs_s = cs_s.merge(
                zip_station.loc[:, ["zip", dist_col]].rename(
                    columns={dist_col : "distance_km"}), 
                how="left", on="zip"
            )
            matches.append(cs_s)

        # isolate unsussesful matches for next iteration
        cs = cs.loc[cs.nearest_idx.isna(), :].drop(
            ["station", "nearest_idx"], axis=1).reset_index(drop=True)

        if cs.empty:
            if verb: print("- Found all matches by index", i + 1)
            break
    # --- 
    # concat the matched dfs back together
    matches = pd.concat(matches, axis=0)

    # determine if any zip-days were left unmatched    
    # (if there are any, add them to the df so we know that no match was found)
    n_unmatched = int(cs.shape[0])
    if n_unmatched:
        matches = pd.concat((matches, cs), axis=0)
        if verb: 
            print(f"- Rows left unmatched: {n_unmatched} ({cs.zip.nunique()} zips)")

    matches.reset_index(drop=True, inplace=True)
    assert matches.shape[0] == n_rows, "Missing rows after matching"
    
    # left-join weather variables for the station which matched to each zip-day
    matches = matches.merge(isd, how="left", on=["station", "day"])
    # save
    out = f"zd_weather_{year}.csv"
    matches.to_csv(os.path.join(zip_date_station_match_dir, out), index=False)
    if verb: print(f"- Saved zip-date weather as {out}")

# ============================================================================
# MAIN
# ============================================================================
def main(start_year, end_year, verbose=True):
    # ensure inputs exist
    for inp, src in {panelist_zipcodes_path : "panelists/panelist_zips.sh", 
                     zipcode_centroids_path : "geographic/crosswalks",
                     us_stations_path : "weather/extract_noaa_isd.py",
                     isd_db_path : "weather/isd_db.py"}.items():
        assert os.path.exists(inp), \
            f"Cannot find {inp}. Ensure prerequisites have been executed: {src}"

    # prep directory for final outputs
    os.makedirs(zip_date_station_match_dir, exist_ok=True) 
    
    # US Stations: list of weather stations in the ISD data
    us_stations = pd.read_csv(us_stations_path, dtype=str)

    # Match zip code centroids to closest weather station
    ## This loads the combined annual panelist files which contains the panelist
    ## id and the panelist zip code for every year.
    ## We match each zip code centroid to the k closest weather stations.
    zip_station = match_zip_to_station(us_stations, k=N_NEIGHBORS, verb=verbose)
    
    # make sure that zip code will merge correctly with other datasets
    zip_station["zip"] = zip_station.zip.apply(lambda x: _pad_id(x, 5))

    # for each year:
    # load panelist cross-section to isolate all valid zip code-date pairings
    # for each zip-date, find the closest weather station *with avail. data*
    # save final datasets matching zip-dates to weather vars at closest station
    print(f"Saving zip-date-weather to: {zip_date_station_match_dir}")
    for yr in range(start_year, end_year + 1):
        print(yr)
        zip_day_closest(zip_station, year=yr, max_neighbors=N_NEIGHBORS,
                        verb=verbose)


if __name__=="__main__":
    main(start_year=START_YEAR, end_year=END_YEAR, verbose=VERBOSE)
