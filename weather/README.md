# Weather

The goal of this sub-project was to acquire weather variables and match them to our panelists.  
We accomplished this by extrtacting all of the US weather station data from NOAA's Integrated Surface Database (ISD) and then matching panelists to their closest weather station.  
- We extracted all of the weather data from US stations for the years in our analysis.  
- We 'geocoded' the panelists using their (5-digit) ZIP code by merging on population weighted ZIP code centroids.  
- We then computed a nearest-neighbor match between the panelist centroids and K different station coordinates (using haversine distance as our metric).  
- Finally, for each panelist-day in our dataset, we determine the closest of the K stations which has data available for this day.  

