#!/usr/bin/env Rscript

# ============================================================================
# Match panelist zip codes to their closest weather station and generate maps
#
# Hans Elliott
# ============================================================================
if (!require("pacman")) install.packages("pacman") # package-manager
pacman::p_load(here, tools,
               # database
               DBI, RSQLite,
               # plotting
               USAboundaries, ggtext, ggrepel, ggplot2,
               # data containers & geospatial ops
               data.table, sf)

# CONSTS
verbose <- TRUE
n_neighbors <- 100
## path to panelist-zipcode match
panelist_path_ <- here::here("panelists/data/panelist_zips.tsv")
## path to zipcode-weather station match
zip_station_path_ <- here::here(paste0("weather/data/zipcode_station_",
                                                        n_neighbors, "nn.csv"))
## path to noaa isd database w/ weather station data
isd_db_path_ <- here::here("weather/isd.db")
## state bounding boxes, for plotting
state_bbox_path_ <- here::here("geographic/crosswalks/data/",
                                                    "state_bounding_boxes.csv")
state_bbox_ <- data.table::fread(state_bbox_path_)
wgs84_crs_ <- 4326

# =============== #
#     HELPERS     #
# =============== #
if (verbose) {
    v <- \(...) { message(paste0(...)) }
} else {
    v <- \(...) {}
}

# get basemaps
get_states <- function() {
    suppressMessages({ USAboundaries::us_states() })
}

get_cities <- function() {
    suppressMessages({ USAboundaries::us_cities() })
}

# transform sf geometries to specific coord. ref. systems
transform_wgs <- function(sf_df) {
    # world geodetic system 84 crs
    return(sf::st_transform(sf_df, crs = 4326))
}

transform_albers <- function(sf_df) {
    # albers equal area conic crs
    return(sf::st_transform(sf_df, crs = 5070))
}

# wrapper for using ggrepel labeling with geom_sf
geom_sf_text_repel <- function(mapping,
                                data,
                                fun.geometry = sf::st_point_on_surface,
                                geom_fun = ggrepel::geom_text_repel,
                                ...) {
    if (is.null(mapping$geometry)) {
        g <- attr(data, "sf_column")
        if (is.null(g)) g <- "geometry"
        mapping$geometry <- as.name(g)
    }
    geom_fun(
        mapping = mapping,
        data = data,
        stat = ggplot2::StatSfCoordinates,
        fun.geometry = fun.geometry,
        ...)
}


# =========================================== #
#     MATCH PANELISTS TO WEATHER STATIONS     #
# =========================================== #
# Each row in the "panelist" data is a panelist(id)-year. Its cols include the
# N (100) closest weather stations to that panelist's centroid, as computed in
# panelist_station_match.py (distance measured via haversine formula).
# Of these N stations, we want to determine the closest weather station that has
# data for the given year and assign it to the panelist.

find_closest_station <- function(zip_station, station_list, neighbors) {
    # Finds the closest station *with data* to each zip-code.
    # Args:
    # zip_station: data.table, data containing zip code to station matches
    #   for one year.
    # station_list: character vector, list of distinct station ids with data for
    #   the given year.
    # neighbors: integer, max number of nearest neighbors to check.
    zip_station[, bool_mask := FALSE]
    zip_station[, closest_station := NA_character_]
    zip_station[, closest_dist_km := NA_real_]
    zip_station[, closest_idx := NA_integer_]

    for (idx in 1:neighbors) {
        # Determine which rows still need to be matched to a weather station:
        # If the next closest station has available weather data and we haven't
        # already found a match for this zip, then the next closest station is
        # our match.
        zip_station[, bool_mask :=
                    (get(paste0("nearest_station_", idx)) %in% station_list)
                    & (is.na(closest_idx))
            ]
        if (sum(is.na(zip_station$closest_idx)) == 0) {
            v("- Found all zip-station matches by index ", idx - 1)
            # all zips will have matched to one of their closest stations
            break
        }
        zip_station[, closest_station := data.table::fifelse(
                    bool_mask,
                    get(paste0("nearest_station_", idx)),
                    closest_station)]
        zip_station[, closest_dist_km := data.table::fifelse(
                    bool_mask,
                    get(paste0("distance_km_", idx)),
                    closest_dist_km)]
        zip_station[, closest_idx := data.table::fifelse(
                    bool_mask,
                    idx,
                    closest_idx)]
    }
    no_match <- sum(is.na(zip_station$closest_idx))
    if (no_match) {
        v("- No. of zips not matched to a station: ", no_match)
    }

    # return necessary cols
    return(zip_station[, .(zip, latitude, longitude,
                           closest_station, closest_dist_km, closest_idx)])
}


# ============================================== #
#     PLOT MAPS OF PANELIST-STATION NETWORKS     #
# ============================================== #
# Each map displays panelists (as points, based on centroid) and a line segment
# connecting them to their closest weather station (also displayed as a point).
# The line segment is colored based on the value of the haversine distance (km)
# that we calculated between the 2 points.

get_state_bbox <- function(state) {
    # returns a ggplot2::coord_sf object to bound the given state
    stopifnot(is.character(state))
    if (nchar(state) == 2) {
        state <- toupper(state)
        bb <- as.numeric(state_bbox_[STUSPS == state,
                                     .(xmin, xmax, ymin, ymax)])
    } else {
        state <- tools::toTitleCase(state)
        bb <- as.numeric(state_bbox_[NAME == state, .(xmin, xmax, ymin, ymax)])
    }
    coord <- ggplot2::coord_sf(xlim = bb[1:2], ylim = bb[3:4],
                               default_crs = sf::st_crs(wgs84_crs_))

    return(coord)
}


plot_network <- function(map_sf, panel_sf, station_sf, loc, title = "") {
    # Plot map of location (loc) that shows panelist zip code centroids and
    # a connection to their closest available weather station.
    stopifnot(nchar(loc) == 2)
    loc <- tolower(loc)
    is_state <- TRUE
    if (loc == "us") is_state <- FALSE

    # prep basemap(s), bounding boxes, and other attributes
    basemap <- get_states()
    if (is_state) {
        transform_crs <- transform_wgs
        cities <- get_cities()
        # label the 5 most populated cities in the state
        cities <- cities[cities$state_abbr == toupper(loc), ]
        cities <- head(cities[order(-cities$population), ], n = 5)
        cities <- transform_crs(cities)
        city_labels <- geom_sf_text_repel(
            data = cities, mapping = aes(label = city),
            box.padding = 0.5, max.overlaps = Inf,
            size = 4, fontface = "bold")
        bbox <- get_state_bbox(loc)
        zip_pt_size <- 1

    } else {
        transform_crs <- transform_albers
        city_labels <- NULL
        bbox <- ggplot2::coord_sf(xlim = c(-125, -65), ylim = c(24, 50),
                                  default_crs = sf::st_crs(wgs84_crs_))
        zip_pt_size <- 0.5
    }
    # transform the crs of each layer
    basemap <- transform_crs(basemap)
    map_sf <- transform_crs(map_sf)
    panel_sf <- transform_crs(panel_sf)
    station_sf <- transform_crs(station_sf)

    # make the plot
    subtitle <- "<span style='color:blue;'>Panelists (blue)</span> and closest
     <span style='color:red;'>station (red)</span>."
    p <- ggplot2::ggplot() +
        # state polygon base map
        geom_sf(data = basemap, fill = "gray95") +
        # panelist to station linestring
        geom_sf(data = map_sf, mapping = aes(color = closest_dist_km),
                alpha = 0.6) +
        # panelist centroids
        geom_sf(data = panel_sf, size = zip_pt_size,
                color = "blue", alpha = 0.4) +
        # station centroids
        geom_sf(data = station_sf, mapping = aes(size = n_panelists),
                color = "red", alpha = 0.5) +
        scale_size_continuous(range = c(1, 3)) +
        scale_color_viridis_c() +
        theme_minimal() +
        theme(
            plot.subtitle = ggtext::element_markdown(),
            panel.grid.minor = element_blank(),
            axis.text = element_blank()
        ) +
        labs(title = title,
             x = "", y = "",
             subtitle = subtitle,
             caption = "Distance estimated by the haversine formula.", 
             color = "Distance (km)",
             size = "Number of Matched Panelists") +
        guides(
            color = guide_colorbar(order = 1),
            size = guide_legend(order = 2)
        ) +
        bbox +
        city_labels

    return(p)
}

# ========= #
#   MAIN    #
# ========= #
# Load data for the given year, create plot inputs, generate and save plots
make_year <- function(year) {
    # load panelist zip code data for the year
    # - matches an individual panelist to the zip code where they live
    panel <- data.table::fread(panelist_path_,
                               sep = "\t", keepLeadingZeros = TRUE,
                               )[panel_year == year, ]
    # group by zip since panelists in the same zip have the same coords
    panel <- panel[, .(n_panelists = .N), by = .(zip)]

    # load weather data for the year
    # - contains all stations - and their coords - which have data for this year
    db <- DBI::dbConnect(RSQLite::SQLite(), isd_db_path_)
    ## (can't use parameterized query to include table name, have to paste)
    isd_uniq <- data.table(DBI::dbGetQuery(db, paste0(
        "SELECT DISTINCT station, longitude, latitude FROM isd_", year),
    ))
    DBI::dbDisconnect(db)

    # load the zip code to weather station match
    # - matches a zip code to the k closest weather stations
    # reduce dims by filtering out zip codes that are not in this panel year
    zip_station <- data.table::fread(zip_station_path_, keepLeadingZeros = TRUE)
    zip_station <- zip_station[zip %in% unique(panel$zip), ]

    # For each zip, find the closest weather station with available data this yr
    zip_station <- find_closest_station(zip_station = zip_station,
                                        station_list = isd_uniq$station,
                                        neighbors = n_neighbors)

    # left-join the closest station data onto panelist ids
    panel <- merge(panel, zip_station, by = "zip", all.x = TRUE)
    ## ideally, all zips in this yr match to a station, but some might not, and
    ## we are missing some zip centroids, so drop na
    panel <- panel[!is.na(panel$closest_station), ]

    # left-join the closest station's lon-lat coordinates to panelist data
    panel <- merge(panel,
                   isd_uniq[, .(closest_station = station,
                                station_lon = longitude,
                                station_lat = latitude)],
                   by = "closest_station",
                   all.x = TRUE)

    # Create 2 individual point geometries, then combine into one df so a
    # linestring can be created connecting panelist zip codes to stations
    panel_sf <- sf::st_as_sf(panel[, .(zip, longitude, latitude)],
                             coords = c("longitude", "latitude"),
                             crs = wgs84_crs_)

    station_sf <- sf::st_as_sf(panel[, .(closest_station, closest_dist_km,
                                        closest_idx, station_lon, station_lat,
                                        n_panelists)],
                               coords = c("station_lon", "station_lat"),
                               crs = wgs84_crs_)
    map_sf <- cbind(panel_sf, station_sf)

    # Union the geometries (lon,lat points) and cast them as a linestring
    v("- Connecting panelist and station coords")
    network <- sf::st_sfc(
        mapply(
            \(a, b) { sf::st_cast(sf::st_union(a, b), "LINESTRING") },
            map_sf$geometry, map_sf$geometry.1,
            SIMPLIFY = FALSE
        )
    )
    sf::st_geometry(map_sf) <- network

    # rm extra geometry & manually assign a crs to avoid transformation errors
    map_sf$geometry.1 <- NULL
    sf::st_crs(map_sf) <- sf::st_crs(wgs84_crs_)

    # group station_sf so that each row is a unique station & count the number
    # of connected panelists
    # (can't group by geometry so have to convert to DT, then back to SF)
    data.table::setDT(station_sf)
    station_uniq <- station_sf[, .(n_panelists = sum(n_panelists)),
                                by = .(closest_station)]
    station_sf <- sf::st_as_sf(
        merge(station_uniq,
              station_sf[, .(closest_station, geometry)],
              by = "closest_station")
    )

    # Create plots
    for (loc in c("us", "or", "wa", "ca", "ny")) {
        v("- Plotting: ", toupper(loc))
        p <- plot_network(map_sf, panel_sf, station_sf,
                          loc = loc,
                          title = paste0(toupper(loc), ", ", year))

        ggplot2::ggsave(
            filename = here::here(paste0(
               "weather/analysis/figures/", loc, "_", year, "_panel_station.png"
              )),
            plot = p,
            bg = "white",
            dpi = 400, width = 7, height = 7, units = "in"
        )
    }
}

# ------------------------------------
for (yr in seq(2004, 2020, by = 4)) {
    v("Making year: ", yr)
    make_year(yr)
}
