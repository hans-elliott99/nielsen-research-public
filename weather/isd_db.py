#!/usr/bin/env python3

# Created SQLite Database for ISD Data
import re
import os
import sqlite3
import pandas as pd
from pyHere import Here
here = Here("nielsen-research")


ISD_DIR = here.here("weather/data")
DB_PATH = here.here("weather/isd.db")

def make_table_for_year(file, conn:sqlite3.Connection) -> None:
    """create table in the database for given year"""
    year = re.findall(r"\d+", file)[0]
    dat = pd.read_csv(file, dtype=str)
    tbl_name = f"isd_{year}" 
    dat.to_sql(tbl_name, conn, if_exists="replace", index=False)
    return tbl_name

if __name__ == "__main__":
    try:
        conn = sqlite3.connect(DB_PATH)
    except Exception as e:
        print(f"[ERROR] Failed to initialize or connect to {DB_PATH}")
        print(e)
    
    print("Successful db connection.")
    for file in os.listdir(ISD_DIR):
        if "isd" in file and file.endswith(".csv"):
            fp = os.path.join(ISD_DIR, file)
            try:
                t = make_table_for_year(fp, conn)
                print("added table:", t)
            except Exception as e:
                print("Failed for file:", file)
                print(e)    
    
    conn.close()
