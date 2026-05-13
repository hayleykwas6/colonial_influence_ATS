################################################################################
# WATER MASS CLASSIFICATION PIPELINE  v2 — Ross Sea Extension
#
# Classifies fishing capture locations into Southern Ocean water masses using
# three independent methods, then gap-fills with the manual method if agreement
# is sufficiently high (>80%) among co-located observations.
#
# Methods
#   1. Manual      – Latitude/depth/longitude zone boundaries (no external data)
#   2. Argo        – Nearest real-time float profile (year-matched to capture)
#   3. WOA         – World Ocean Atlas 2018 seasonal climatology (season-matched)
#
# Water masses — open ocean (sigma0-based)
#   AASW  Antarctic Surface Water          sigma0 < 26.8
#   SAMW  Subantarctic Mode Water    26.8 ≤ sigma0 < 27.1
#   AAIW  Antarctic Intermediate W.  27.1 ≤ sigma0 < 27.4
#   CDW   Circumpolar Deep Water     27.4 ≤ sigma0 < 27.9
#   AABW  Antarctic Bottom Water          sigma0 ≥ 27.9
#
# Water masses — Ross Sea continental shelf (T/S-based)
#   AASW     Antarctic Surface Water        depth < 150 m
#   MCDW     Modified Circumpolar Deep Water T > -1.5°C, depth ≥ 150 m
#   SW_salty Dense Shelf Water (saline)     T < -1.7°C, S ≥ 34.62
#   SW_fresh Dense Shelf Water (fresh)      T < -1.7°C, 34.3 ≤ S < 34.62
#   ISW      Ice Shelf Water                T < -1.9°C, S < 34.52
#   AABW     Exported Antarctic Bottom Water depth > 700 m (shelf break)
#
# Key references for Ross Sea thresholds:
#   Orsi & Wiederwohl (2009) Deep Sea Res II 56:778-795
#   Jacobs et al. (1985) Antarctic Research Series 43:59-85
#   Ainley & Jacobs (1981) Deep Sea Res 28:1173-1185
#   Smith et al. (2012) Oceanography 25(3):90-103
#
# Expected columns in edge_df
#   Lat       – decimal degrees (negative = South)
#   Long360   – longitude 0–360°
#   Depth     – capture depth in metres
#   Year      – integer capture year  (e.g. 2019)
#   Month     – integer capture month (1–12); if absent, derived from Date
#   Date      – optional Date/POSIXct; used only if Month is absent
################################################################################


# ==============================================================================
# SECTION 0 · LIBRARIES
# ==============================================================================

library(argoFloats)   # Argo float data access
library(gsw)          # TEOS-10 thermodynamic functions
library(oce)          # geodDist()
library(ncdf4)        # WOA NetCDF files
library(dplyr)        # Data manipulation
library(ggplot2)      # Visualisation
library(patchwork)    # Multi-panel plots
library(knitr)        # Summary tables
library(purrr)        # pmap_chr in Li depth analysis


# ==============================================================================
# SECTION 1 · SHARED WATER MASS CLASSIFICATION FUNCTIONS
# ==============================================================================

# ---------------------------------------------------------------------------- #
# 1a. Sigma0-based classifier (open ocean; retained for backward compatibility)
# ---------------------------------------------------------------------------- #
# CDW combines the former UCDW and LCDW categories for consistency across
# the manual, Argo, and WOA methods in non-Ross-Sea regions.

classify_water_mass <- function(sigma0) {
  dplyr::case_when(
    sigma0 < 26.8 ~ "AASW",
    sigma0 < 27.1 ~ "SAMW",
    sigma0 < 27.4 ~ "AAIW",
    sigma0 < 27.9 ~ "CDW",
    TRUE          ~ "AABW"
  )
}

# Helper: compute sigma0 from practical salinity, in-situ temperature, depth
compute_sigma0 <- function(salinity, temperature, depth, lon, lat) {
  SA <- gsw::gsw_SA_from_SP(salinity, depth, lon, lat)
  CT <- gsw::gsw_CT_from_t(SA, temperature, depth)
  gsw::gsw_sigma0(SA, CT)
}

# ---------------------------------------------------------------------------- #
# 1b. Ross Sea geographic detector
# ---------------------------------------------------------------------------- #
# Returns TRUE for points on the Ross Sea continental shelf.
# Bounds follow Orsi & Wiederwohl (2009) and Smith et al. (2012) Box 1:
#   Latitude  : south of -65° (captures full shelf + outer gyre influence)
#   Longitude : 155–220° in Long360 (≈ 155°E to 140°W)
#   Shelf break defined at ~700 m (Smith et al. 2012 Box 1)

is_ross_sea_shelf <- function(lat, lon360, depth = NULL,
                              lat_cut  = -65,
                              lon_min  = 155,
                              lon_max  = 220,
                              shelf_m  = 700) {
  in_lat   <- !is.na(lat)    & lat    <= lat_cut
  in_lon   <- !is.na(lon360) & lon360 >= lon_min & lon360 <= lon_max
  on_shelf <- if (!is.null(depth)) is.na(depth) | depth <= shelf_m else TRUE
  in_lat & in_lon & on_shelf
}

# ---------------------------------------------------------------------------- #
# 1c. T/S-based classifier (used by Argo and WOA passes)
# ---------------------------------------------------------------------------- #
# On the Ross Sea continental shelf, temperature tests take priority over sigma0
# because MCDW, SW, and ISW overlap heavily in density space but are thermally
# distinct.  Off-shelf, falls back to the original sigma0 thresholds.
#
# Thresholds:
#   ISW      T < -1.9°C, S < 34.52   Jacobs et al. (1985); sub-ice-shelf outflow
#            T below surface freezing point, relatively fresh
#   SW_salty T < -1.7°C, S ≥ 34.62  Brine-rejected dense water (Orsi & Wiederwohl 2009)
#   SW_fresh T < -1.7°C, 34.3 ≤ S < 34.62  Eastern/seaward shelf variety
#   MCDW     T > -1.5°C, depth ≥ 150 Modified CDW intruding onto the shelf
#            (Core MCDW: T > 0.5°C; Ainley & Jacobs 1981)
#   AASW     depth < 150 OR sigma0 < 27.0

classify_water_mass_ts <- function(temp, sal, sigma0, depth, lat, lon360) {
  
  on_rs <- is_ross_sea_shelf(lat, lon360, depth)
  
  if (on_rs) {
    
    # Surface layer — always AASW regardless of T/S
    if (!is.na(depth) && depth < 150) return("AASW")
    
    # ISW: sub-surface, colder than surface freezing, relatively fresh
    # Emerges at mid-depths from beneath the Ross Ice Shelf
    if (!is.na(temp) && !is.na(sal) && temp < -1.9 && sal < 34.52)
      return("ISW")
    
    # Salty Shelf Water: brine-enriched, near-freezing, dense
    # Fills continental shelf bottom; primary AABW precursor
    if (!is.na(temp) && !is.na(sal) && temp < -1.7 && sal >= 34.62)
      return("SW_salty")
    
    # Fresh Shelf Water: near-freezing, lower salinity (eastern shelf)
    if (!is.na(temp) && !is.na(sal) && temp < -1.7 && sal >= 34.3 && sal < 34.62)
      return("SW_fresh")
    
    # MCDW: subsurface warm intrusion from the ACC via the shelf break
    # T > -1.5°C is the broad definition; T > 0.5°C marks the core
    if (!is.na(temp) && !is.na(depth) && temp > -1.5 && depth >= 150)
      return("MCDW")
    
    # Fallback: sigma0-based for remaining Ross Sea shelf points
    # sigma0 < 27.0 is NOT assigned AASW here — any point reaching this
    # fallback is at depth > 150 m and is better described as SW_fresh
    # (cold, intermediate density, meltwater-influenced)
    if (!is.na(sigma0)) {
      if (sigma0 >= 27.88) return("AABW")  # rare exported type via troughs
      return("SW_fresh")                    # cold, intermediate density
    }
    
  } else {
    
    # ── Open ocean / off-shelf: sigma0 thresholds ────────────────────────────
    # AASW: surface only — depth guard prevents anomalously fresh deep Argo
    #        profiles from being miscategorised as surface water
    # SAMW/AAIW: Subantarctic and Polar Frontal Zone only (lat >= -60°) —
    #        these water masses do not exist in the Antarctic Zone and any
    #        sigma0 in their range south of -60° is meltwater-modified CDW
    if (!is.na(sigma0)) {
      return(dplyr::case_when(
        sigma0 < 26.8 & !is.na(depth) & depth < 300       ~ "AASW",
        sigma0 < 27.1 & !is.na(lat)   & lat   >= -60      ~ "SAMW",
        sigma0 < 27.4 & !is.na(lat)   & lat   >= -60      ~ "AAIW",
        sigma0 < 27.9                                      ~ "CDW",
        TRUE                                               ~ "AABW"
      ))
    }
  }
  
  NA_character_
}

# ==============================================================================
# SECTION 2 · MANUAL CLASSIFICATION (latitude + depth + longitude zones)
# ==============================================================================
# Ross Sea continental shelf receives a dedicated branch using depth zones
# calibrated to the regional water mass structure (Smith et al. 2012 Fig 3).
# The remaining branches follow the classical Southern Ocean front/zone scheme.

classify_water_mass_manual <- function(lat, depth, lon360 = NA) {
  
  on_rs <- !is.na(lon360) && is_ross_sea_shelf(lat, lon360, depth)
  
  # ── Ross Sea continental shelf ────────────────────────────────────────────
  # Boundaries from Orsi & Wiederwohl (2009); Smith et al. (2012)
  if (on_rs) {
    return(dplyr::case_when(
      
      # Surface mixed layer — always AASW
      depth < 150                                     ~ "AASW",
      
      # ISW: emerges at mid-depths near the Ross Ice Shelf front
      # Sub-ice-shelf outflow appears at ~200–500 m (Jacobs & Giulivi 1998)
      depth >= 150 & depth < 500 & lat < -77          ~ "ISW",
      
      # MCDW: intrudes at intermediate depths along western bank flanks
      # Most conspicuous 150–500 m range (Orsi & Wiederwohl 2009)
      depth >= 150 & depth < 500 & lat >= -77         ~ "MCDW",
      
      # Below ~500 m on shelf: Shelf Water dominates the bottom layer
      # (salty variety fills most of the shelf bottom; Smith et al. 2012)
      depth >= 500 & depth <= 700                     ~ "SW_salty",
      
      # Below shelf break (~700 m): rare AABW exported through troughs
      # (Drygalski, Joides, Glomar Challenger; Orsi et al. 1999)
      depth > 700                                     ~ "AABW",
      
      TRUE                                            ~ NA_character_
    ))
  }
  
  # ── Antarctic Zone south of -60° (non-Ross Sea) ──────────────────────────
  # NOTE: original code incorrectly labelled 200–800 m at lat < -70 as AABW;
  # corrected here — that depth/lat range is CDW in the open Antarctic Zone.
  if (lat < -60) {
    return(dplyr::case_when(
      depth < 200                  ~ "AASW",
      depth >= 200  & depth < 1200 ~ "CDW",
      depth >= 1200 & depth < 3500 ~ "CDW",
      depth >= 3500                ~ "AABW",
      TRUE                         ~ NA_character_
    ))
  }
  
  # ── Polar Frontal Zone -60° to -52° ─────────────────────────────────────
  if (lat >= -60 & lat < -52) {
    return(dplyr::case_when(
      depth < 200                  ~ "AASW",
      depth >= 200  & depth < 600  ~ "SAMW",
      depth >= 600  & depth < 1500 ~ "AAIW",
      depth >= 1500 & depth < 4000 ~ "CDW",
      depth >= 4000                ~ "AABW",
      TRUE                         ~ NA_character_
    ))
  }
  
  # ── Subantarctic Zone north of -52° ─────────────────────────────────────
  dplyr::case_when(
    depth < 200                  ~ "AASW",
    depth >= 200  & depth < 700  ~ "SAMW",
    depth >= 700  & depth < 1500 ~ "AAIW",
    depth >= 1500 & depth < 4000 ~ "CDW",
    depth >= 4000                ~ "AABW",
    TRUE                         ~ NA_character_
  )
}

run_manual_classification <- function(df) {
  cat("=== MANUAL CLASSIFICATION (lat + depth + lon zones) ===\n")
  
  df$water_mass_manual <- NA_character_
  
  # Long360 is now required for the Ross Sea branch
  has_lon <- "Long360" %in% names(df)
  if (!has_lon) warning("Long360 column not found — Ross Sea branch will not activate.")
  
  valid <- !is.na(df$Lat) & !is.na(df$Depth)
  
  for (i in which(valid)) {
    lon_i <- if (has_lon) df$Long360[i] else NA_real_
    df$water_mass_manual[i] <- classify_water_mass_manual(
      df$Lat[i], df$Depth[i], lon_i
    )
  }
  
  cat(sprintf("Classified %d / %d locations\n", sum(!is.na(df$water_mass_manual)), nrow(df)))
  cat("Distribution:\n")
  print(table(df$water_mass_manual))
  
  # Latitude zone label (useful for diagnostics)
  df$lat_zone <- dplyr::case_when(
    is_ross_sea_shelf(df$Lat, df$Long360)        ~ "Ross Sea Shelf",
    df$Lat < -60                                 ~ "Antarctic Zone (<-60°)",
    df$Lat >= -60 & df$Lat < -52                 ~ "Polar Frontal Zone (-60 to -52°)",
    df$Lat >= -52                                ~ "Subantarctic Zone (>-52°)",
    TRUE                                         ~ NA_character_
  )
  
  df
}


# ==============================================================================
# SECTION 3 · ARGO CLASSIFICATION (year-matched to capture year)
# ==============================================================================

# ---------------------------------------------------------------------------- #
# 3a. Download Argo index and cache profiles for the capture year range
# ---------------------------------------------------------------------------- #

get_argo_profiles_for_years <- function(df,
                                        lat_col    = "Lat",
                                        year_col   = "Year",
                                        buffer_deg = 2,
                                        cache_file = "argo_profiles_cache.rds") {
  
  if (file.exists(cache_file)) {
    cat("Loading cached Argo profiles from disk...\n")
    profiles <- readRDS(cache_file)
    cat(sprintf("Loaded %d profiles from cache.\n", length(profiles)))
    return(profiles)
  }
  
  capture_years <- df[[year_col]][!is.na(df[[year_col]])]
  if (length(capture_years) == 0) stop("No valid capture years found in '", year_col, "' column.")
  
  date_from <- as.Date(paste0(min(capture_years) - 1, "-01-01"))
  date_to   <- as.Date(paste0(max(capture_years) + 1, "-12-31"))
  lat_range <- range(df[[lat_col]], na.rm = TRUE)
  
  cat("=== DOWNLOADING ARGO DATA ===\n")
  cat(sprintf("  Capture years: %d to %d\n", min(capture_years), max(capture_years)))
  cat(sprintf("  Download window: %s to %s\n", date_from, date_to))
  cat(sprintf("  Latitude band: %.1f° to %.1f° (full longitude)\n",
              lat_range[1] - buffer_deg, lat_range[2] + buffer_deg))
  cat("  Expected time: 20–40 minutes\n\n")
  
  ai <- argoFloats::getIndex()
  
  ai_sub <- argoFloats::subset(ai, rectangle = list(
    longitude = c(-180, 180),
    latitude  = c(lat_range[1] - buffer_deg, lat_range[2] + buffer_deg)
  ))
  cat(sprintf("Profiles in latitude band: %d\n", length(ai_sub[["longitude"]])))
  
  ai_sub <- argoFloats::subset(ai_sub, time = list(from = date_from, to = date_to))
  cat(sprintf("Profiles after time filter: %d\n", length(ai_sub[["longitude"]])))
  
  if (length(ai_sub[["longitude"]]) == 0) {
    warning("No Argo profiles found for the specified region and time window.")
    return(NULL)
  }
  
  cat("Reading profiles (may take 10–20 min)...\n")
  profiles <- argoFloats::readProfiles(argoFloats::getProfiles(ai_sub))
  
  cat("Caching profiles to disk...\n")
  saveRDS(profiles, cache_file)
  cat(sprintf("Saved to: %s\n", cache_file))
  
  profiles
}

# ---------------------------------------------------------------------------- #
# 3b. Extract T/S from the nearest temporally-matched profile
# ---------------------------------------------------------------------------- #

get_ts_from_profiles <- function(profiles, lat, lon, depth,
                                 capture_year = NULL,
                                 max_dist_km  = 300,
                                 year_window  = 1) {
  
  if (is.null(profiles) || length(profiles) == 0) {
    return(list(temp = NA, sal = NA, source = "none"))
  }
  
  best_temp <- NA
  best_sal  <- NA
  min_dist  <- Inf
  
  profile_list <- profiles@data$argos
  n_profiles   <- length(profile_list)
  
  for (i in seq_len(n_profiles)) {
    prof <- profile_list[[i]]
    
    if (is.null(prof@data$pressure)    ||
        is.null(prof@data$temperature) ||
        is.null(prof@data$salinity))   next
    
    prof_lon <- prof@data$longitude[1]
    prof_lat <- prof@data$latitude[1]
    
    if (is.na(prof_lon) || is.na(prof_lat)) next
    
    if (!is.null(capture_year) && !is.na(capture_year)) {
      prof_time <- tryCatch(
        as.numeric(format(as.POSIXct(prof@metadata$time, origin = "1970-01-01",
                                     tz = "UTC"), "%Y")),
        error = function(e) NA_real_
      )
      if (!is.na(prof_time) &&
          abs(prof_time - capture_year) > year_window) next
    }
    
    dist <- oce::geodDist(lon, lat, prof_lon, prof_lat)
    if (is.na(dist) || is.infinite(dist) || dist > max_dist_km) next
    
    pressure    <- prof@data$pressure
    temperature <- prof@data$temperature
    salinity    <- prof@data$salinity
    
    valid <- !is.na(pressure) & !is.na(temperature) & !is.na(salinity)
    pressure    <- pressure[valid]
    temperature <- temperature[valid]
    salinity    <- salinity[valid]
    
    if (length(pressure) < 2) next
    
    max_p <- max(pressure, na.rm = TRUE)
    if (is.na(max_p) || is.infinite(max_p) || max_p < depth) next
    
    if (dist < min_dist) {
      temp_i <- stats::approx(pressure, temperature, xout = depth)$y
      sal_i  <- stats::approx(pressure, salinity,    xout = depth)$y
      
      if (!is.na(temp_i) && !is.na(sal_i)) {
        best_temp <- temp_i
        best_sal  <- sal_i
        min_dist  <- dist
      }
    }
  }
  
  list(
    temp   = best_temp,
    sal    = best_sal,
    source = ifelse(is.na(best_temp), "none", "argo")
  )
}

# ---------------------------------------------------------------------------- #
# 3c. Orchestrate the full Argo classification pass
# ---------------------------------------------------------------------------- #

run_argo_classification <- function(df,
                                    profiles,
                                    year_col    = "Year",
                                    max_dist_km = 300,
                                    year_window = 1) {
  
  cat("=== ARGO CLASSIFICATION (year-matched, T/S-aware) ===\n")
  
  has_year <- year_col %in% names(df)
  if (!has_year) warning("Column '", year_col, "' not found; temporal matching disabled.")
  
  df$temperature_argo <- NA_real_
  df$salinity_argo    <- NA_real_
  df$density_argo     <- NA_real_
  df$water_mass_argo  <- NA_character_
  df$source_argo      <- NA_character_
  
  valid_rows <- !is.na(df$Lat) & !is.na(df$Long360) & !is.na(df$Depth)
  cat(sprintf("Processing %d valid locations...\n", sum(valid_rows)))
  
  pb <- utils::txtProgressBar(min = 0, max = nrow(df), style = 3)
  
  for (i in seq_len(nrow(df))) {
    utils::setTxtProgressBar(pb, i)
    if (!valid_rows[i]) next
    
    cap_year <- if (has_year) df[[year_col]][i] else NULL
    
    ts <- get_ts_from_profiles(
      profiles, df$Lat[i], df$Long360[i], df$Depth[i],
      capture_year = cap_year,
      max_dist_km  = max_dist_km,
      year_window  = year_window
    )
    
    if (!is.na(ts$temp)) {
      sigma0 <- tryCatch(
        compute_sigma0(ts$sal, ts$temp, df$Depth[i], df$Long360[i], df$Lat[i]),
        error = function(e) NA_real_
      )
      df$temperature_argo[i] <- ts$temp
      df$salinity_argo[i]    <- ts$sal
      df$density_argo[i]     <- sigma0
      
      # ── UPDATED: T/S-aware classifier (replaces sigma0-only) ──────────────
      df$water_mass_argo[i] <- classify_water_mass_ts(
        ts$temp, ts$sal, sigma0,
        df$Depth[i], df$Lat[i], df$Long360[i]
      )
      df$source_argo[i] <- "argo"
    }
  }
  
  close(pb)
  
  n_class <- sum(!is.na(df$water_mass_argo))
  n_valid <- sum(valid_rows)
  cat(sprintf("\nClassified %d / %d (%.1f%%)\n", n_class, n_valid, 100 * n_class / n_valid))
  cat("Water mass distribution (Argo):\n")
  print(table(df$water_mass_argo))
  
  df
}


# ==============================================================================
# SECTION 4 · WOA CLASSIFICATION (season-matched to capture month)
# ==============================================================================

month_to_woa_season <- function(month) {
  dplyr::case_when(
    month %in% c(12, 1, 2)  ~ "13",  # Dec–Feb (austral summer)
    month %in% c(3, 4, 5)   ~ "14",  # Mar–May
    month %in% c(6, 7, 8)   ~ "15",  # Jun–Aug
    month %in% c(9, 10, 11) ~ "16",  # Sep–Nov
    TRUE                     ~ "00"   # fallback: annual
  )
}

# ---------------------------------------------------------------------------- #
# 4a. Download the required WOA season files
# ---------------------------------------------------------------------------- #

download_woa_seasonal <- function(season_codes, data_dir = "woa_data") {
  
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
  
  base_t <- "https://www.ncei.noaa.gov/thredds-ocean/fileServer/ncei/woa/temperature/decav/1.00/"
  base_s <- "https://www.ncei.noaa.gov/thredds-ocean/fileServer/ncei/woa/salinity/decav/1.00/"
  
  files <- list()
  
  for (sc in unique(season_codes)) {
    t_name <- sprintf("woa18_decav_t%s_01.nc", sc)
    s_name <- sprintf("woa18_decav_s%s_01.nc", sc)
    t_path <- file.path(data_dir, t_name)
    s_path <- file.path(data_dir, s_name)
    
    if (!file.exists(t_path)) {
      cat(sprintf("Downloading WOA temperature (season %s)...\n", sc))
      utils::download.file(paste0(base_t, t_name), t_path, mode = "wb", quiet = FALSE)
    } else {
      cat(sprintf("WOA temperature (season %s) already cached.\n", sc))
    }
    
    if (!file.exists(s_path)) {
      cat(sprintf("Downloading WOA salinity (season %s)...\n", sc))
      utils::download.file(paste0(base_s, s_name), s_path, mode = "wb", quiet = FALSE)
    } else {
      cat(sprintf("WOA salinity (season %s) already cached.\n", sc))
    }
    
    files[[sc]] <- list(temp = t_path, sal = s_path)
  }
  
  files
}

# ---------------------------------------------------------------------------- #
# 4b. Orchestrate the full WOA classification pass
# ---------------------------------------------------------------------------- #

run_woa_classification <- function(df,
                                   month_col = "Month",
                                   data_dir  = "woa_data") {
  
  cat("=== WOA CLASSIFICATION (season-matched, T/S-aware) ===\n")
  
  has_month <- month_col %in% names(df)
  if (has_month) {
    df$woa_season_code <- month_to_woa_season(df[[month_col]])
    cat("Using capture Month column for seasonal WOA matching.\n")
  } else {
    df$woa_season_code <- "00"
    cat("No Month column found; using WOA annual climatology for all records.\n")
  }
  
  season_codes <- unique(df$woa_season_code[!is.na(df$woa_season_code)])
  cat(sprintf("Unique season codes to download: %s\n\n", paste(season_codes, collapse = ", ")))
  
  woa_files <- download_woa_seasonal(season_codes, data_dir)
  
  df$temperature_woa <- NA_real_
  df$salinity_woa    <- NA_real_
  df$density_woa     <- NA_real_
  df$water_mass_woa  <- NA_character_
  
  valid_rows <- !is.na(df$Lat) & !is.na(df$Long360) & !is.na(df$Depth)
  cat(sprintf("\nProcessing %d valid locations...\n", sum(valid_rows)))
  
  nc_handles <- list()
  open_nc <- function(sc) {
    if (is.null(nc_handles[[sc]])) {
      nc_handles[[sc]] <<- list(
        temp = ncdf4::nc_open(woa_files[[sc]]$temp),
        sal  = ncdf4::nc_open(woa_files[[sc]]$sal)
      )
    }
    nc_handles[[sc]]
  }
  
  nc0 <- ncdf4::nc_open(woa_files[[season_codes[1]]]$temp)
  woa_lon   <- ncdf4::ncvar_get(nc0, "lon")
  woa_lat   <- ncdf4::ncvar_get(nc0, "lat")
  woa_depth <- ncdf4::ncvar_get(nc0, "depth")
  ncdf4::nc_close(nc0)
  
  pb <- utils::txtProgressBar(min = 0, max = nrow(df), style = 3)
  
  for (i in seq_len(nrow(df))) {
    utils::setTxtProgressBar(pb, i)
    if (!valid_rows[i]) next
    
    lat   <- df$Lat[i]
    lon   <- df$Long360[i]
    depth <- df$Depth[i]
    
    if (!is.na(lon) && lon < 0) lon <- lon + 360
    
    sc  <- df$woa_season_code[i]
    nch <- open_nc(sc)
    
    lon_idx   <- which.min(abs(woa_lon   - lon))
    lat_idx   <- which.min(abs(woa_lat   - lat))
    depth_idx <- which.min(abs(woa_depth - depth))
    
    temp <- tryCatch(
      ncdf4::ncvar_get(nch$temp, "t_an",
                       start = c(lon_idx, lat_idx, depth_idx, 1),
                       count = c(1, 1, 1, 1)),
      error = function(e) NA_real_
    )
    sal <- tryCatch(
      ncdf4::ncvar_get(nch$sal, "s_an",
                       start = c(lon_idx, lat_idx, depth_idx, 1),
                       count = c(1, 1, 1, 1)),
      error = function(e) NA_real_
    )
    
    if (!is.na(temp) && !is.na(sal) &&
        temp > -10 && temp < 40 &&
        sal  >   0 && sal  < 50) {
      
      sigma0 <- tryCatch(
        compute_sigma0(sal, temp, depth, df$Long360[i], lat),
        error = function(e) NA_real_
      )
      df$temperature_woa[i] <- temp
      df$salinity_woa[i]    <- sal
      df$density_woa[i]     <- sigma0
      
      # ── UPDATED: T/S-aware classifier (replaces sigma0-only) ──────────────
      df$water_mass_woa[i] <- classify_water_mass_ts(
        temp, sal, sigma0,
        depth, lat, df$Long360[i]
      )
    }
  }
  
  close(pb)
  
  for (sc in names(nc_handles)) {
    ncdf4::nc_close(nc_handles[[sc]]$temp)
    ncdf4::nc_close(nc_handles[[sc]]$sal)
  }
  
  n_class <- sum(!is.na(df$water_mass_woa))
  n_valid <- sum(valid_rows)
  cat(sprintf("\nClassified %d / %d (%.1f%%)\n", n_class, n_valid, 100 * n_class / n_valid))
  cat("Water mass distribution (WOA):\n")
  print(table(df$water_mass_woa))
  
  df
}


# ==============================================================================
# SECTION 5 · AGREEMENT ANALYSIS & MANUAL GAP-FILL
# ==============================================================================

run_agreement_and_gap_fill <- function(df, agreement_threshold = 0.80) {
  
  cat("\n=== AGREEMENT ANALYSIS ===\n")
  
  covered        <- !is.na(df$water_mass_argo) | !is.na(df$water_mass_woa)
  uncovered_argo <- is.na(df$water_mass_argo) & !is.na(df$Lat)
  uncovered_woa  <- is.na(df$water_mass_woa)  & !is.na(df$Lat)
  
  cat(sprintf("Records with Argo classification:  %d / %d\n",
              sum(!is.na(df$water_mass_argo)), nrow(df)))
  cat(sprintf("Records with WOA classification:   %d / %d\n",
              sum(!is.na(df$water_mass_woa)), nrow(df)))
  cat(sprintf("Records missing Argo (temporal gap or no nearby float): %d\n",
              sum(uncovered_argo)))
  cat(sprintf("Records missing WOA:               %d\n\n",
              sum(uncovered_woa)))
  
  complete_rows <- !is.na(df$water_mass_manual) &
    !is.na(df$water_mass_argo)   &
    !is.na(df$water_mass_woa)
  
  n_complete <- sum(complete_rows)
  cat(sprintf("Records with all 3 methods (used for agreement): %d\n\n", n_complete))
  
  agreement_argo <- NA_real_
  agreement_woa  <- NA_real_
  
  if (n_complete > 0) {
    sub <- df[complete_rows, ]
    
    agree_manual_argo <- mean(sub$water_mass_manual == sub$water_mass_argo, na.rm = TRUE)
    agree_manual_woa  <- mean(sub$water_mass_manual == sub$water_mass_woa,  na.rm = TRUE)
    agree_argo_woa    <- mean(sub$water_mass_argo   == sub$water_mass_woa,  na.rm = TRUE)
    
    agreement_argo <- agree_manual_argo
    agreement_woa  <- agree_manual_woa
    
    cat(sprintf("Manual vs Argo agreement: %.1f%%\n", 100 * agree_manual_argo))
    cat(sprintf("Manual vs WOA  agreement: %.1f%%\n", 100 * agree_manual_woa))
    cat(sprintf("Argo   vs WOA  agreement: %.1f%%\n\n", 100 * agree_argo_woa))
    
    cat("Confusion matrix – Manual vs Argo:\n")
    print(table(Manual = sub$water_mass_manual, Argo = sub$water_mass_argo))
    cat("\nConfusion matrix – Manual vs WOA:\n")
    print(table(Manual = sub$water_mass_manual, WOA  = sub$water_mass_woa))
  } else {
    cat("Warning: No records with all three methods – cannot compute agreement.\n")
    cat("         Manual gap-fill will NOT be applied.\n")
  }
  
  df$water_mass_argo_filled <- df$water_mass_argo
  df$water_mass_woa_filled  <- df$water_mass_woa
  df$argo_fill_source       <- ifelse(!is.na(df$water_mass_argo), "argo", NA_character_)
  df$woa_fill_source        <- ifelse(!is.na(df$water_mass_woa),  "woa",  NA_character_)
  
  fill_argo <- !is.na(agreement_argo) &&
    agreement_argo >= agreement_threshold &&
    sum(uncovered_argo) > 0
  
  fill_woa <- !is.na(agreement_woa) &&
    agreement_woa  >= agreement_threshold &&
    sum(uncovered_woa) > 0
  
  if (fill_argo) {
    cat(sprintf(
      "\nManual-vs-Argo agreement (%.1f%%) exceeds threshold (%.0f%%).\n",
      100 * agreement_argo, 100 * agreement_threshold))
    cat(sprintf("Applying manual classification to %d Argo-uncovered records.\n",
                sum(uncovered_argo)))
    df$water_mass_argo_filled[uncovered_argo] <- df$water_mass_manual[uncovered_argo]
    df$argo_fill_source[uncovered_argo]       <- "manual_fill"
  } else if (sum(uncovered_argo) > 0) {
    cat(sprintf(
      "\nManual-vs-Argo agreement (%.1f%%) is BELOW threshold (%.0f%%).\n",
      100 * agreement_argo, 100 * agreement_threshold))
    cat("Argo-uncovered records remain NA.\n")
  }
  
  if (fill_woa) {
    cat(sprintf(
      "\nManual-vs-WOA agreement (%.1f%%) exceeds threshold (%.0f%%).\n",
      100 * agreement_woa, 100 * agreement_threshold))
    cat(sprintf("Applying manual classification to %d WOA-uncovered records.\n",
                sum(uncovered_woa)))
    df$water_mass_woa_filled[uncovered_woa] <- df$water_mass_manual[uncovered_woa]
    df$woa_fill_source[uncovered_woa]       <- "manual_fill"
  } else if (sum(uncovered_woa) > 0) {
    cat(sprintf(
      "\nManual-vs-WOA agreement (%.1f%%) is BELOW threshold (%.0f%%).\n",
      100 * agreement_woa, 100 * agreement_threshold))
    cat("WOA-uncovered records remain NA.\n")
  }
  
  df
}


# ==============================================================================
# SECTION 6 · CONSENSUS CLASSIFICATION
# ==============================================================================

compute_consensus <- function(df) {
  
  cat("\n=== CONSENSUS CLASSIFICATION ===\n")
  
  df$water_mass_consensus <- NA_character_
  
  for (i in seq_len(nrow(df))) {
    
    # Depth override: fish at < 160 m in the Antarctic Zone are in AASW
    # regardless of T/S method output — otolith chemistry validates this
    if (!is.na(df$Depth[i]) && !is.na(df$Lat[i]) &&
        df$Depth[i] < 160 && df$Lat[i] < -52 &&
        !is.na(df$water_mass_manual[i]) &&
        df$water_mass_manual[i] == "AASW") {
      df$water_mass_consensus[i] <- "AASW"
      next  # skip majority vote for this record
    }
    
    # Standard majority vote for everything else
    methods <- c(
      df$water_mass_manual[i],
      df$water_mass_argo_filled[i],
      df$water_mass_woa_filled[i]
    )
    methods <- methods[!is.na(methods)]
    
    if (length(methods) >= 2) {
      tbl <- table(methods)
      if (max(tbl) >= 2) {
        df$water_mass_consensus[i] <- names(tbl)[which.max(tbl)]
      }
    } else if (length(methods) == 1) {
      df$water_mass_consensus[i] <- methods[1]
    }
  }  # ← single closing brace for the one loop
  
  n_consensus <- sum(!is.na(df$water_mass_consensus))
  cat(sprintf("Consensus assigned: %d / %d (%.1f%%)\n",
              n_consensus, nrow(df), 100 * n_consensus / nrow(df)))
  cat("Consensus distribution:\n")
  print(table(df$water_mass_consensus))
  
  df
}


# ==============================================================================
# SECTION 7 · SUMMARY TABLES
# ==============================================================================

print_summary_tables <- function(df) {
  
  cat("\n\n=== SUMMARY TABLE 1: Dataset Coverage ===\n")
  cov <- data.frame(
    Method = c("Manual (lat/depth/lon)", "Argo (year-matched)",
               "Argo (after gap-fill)", "WOA (season-matched)",
               "WOA (after gap-fill)", "Consensus (≥2 agree)"),
    n = c(
      sum(!is.na(df$water_mass_manual)),
      sum(!is.na(df$water_mass_argo)),
      sum(!is.na(df$water_mass_argo_filled)),
      sum(!is.na(df$water_mass_woa)),
      sum(!is.na(df$water_mass_woa_filled)),
      sum(!is.na(df$water_mass_consensus))
    )
  )
  cov$Pct <- round(100 * cov$n / nrow(df), 1)
  print(knitr::kable(cov, col.names = c("Method", "n", "%"), format = "simple"))
  
  cat("\n=== SUMMARY TABLE 2: Agreement Rates (records with all 3 methods) ===\n")
  complete <- !is.na(df$water_mass_manual) &
    !is.na(df$water_mass_argo)   &
    !is.na(df$water_mass_woa)
  if (sum(complete) > 0) {
    sub <- df[complete, ]
    agr <- data.frame(
      Comparison = c("Manual vs Argo", "Manual vs WOA", "Argo vs WOA"),
      Agree = c(
        sum(sub$water_mass_manual == sub$water_mass_argo),
        sum(sub$water_mass_manual == sub$water_mass_woa),
        sum(sub$water_mass_argo   == sub$water_mass_woa)
      ),
      Total = rep(nrow(sub), 3)
    )
    agr$Pct <- round(100 * agr$Agree / agr$Total, 1)
    print(knitr::kable(agr, col.names = c("Comparison", "Agree", "Total", "%"),
                       format = "simple"))
  }
  
  # ── UPDATED: expanded water mass levels including Ross Sea categories ───────
  cat("\n=== SUMMARY TABLE 3: Water Mass Distribution by Method ===\n")
  wm_levels <- c("AASW", "SAMW", "AAIW", "MCDW", "CDW",
                 "SW_fresh", "SW_salty", "ISW", "AABW")
  wm_tbl <- data.frame(Water_Mass = wm_levels)
  for (col in c("water_mass_manual", "water_mass_argo_filled",
                "water_mass_woa_filled", "water_mass_consensus")) {
    wm_tbl[[gsub("water_mass_", "", col)]] <-
      sapply(wm_levels, function(wm) sum(df[[col]] == wm, na.rm = TRUE))
  }
  print(knitr::kable(wm_tbl, format = "simple"))
  
  cat("\n=== SUMMARY TABLE 4: Consensus Water Mass Characteristics ===\n")
  if (sum(!is.na(df$water_mass_consensus)) > 0) {
    stats_tbl <- df %>%
      dplyr::filter(!is.na(water_mass_consensus)) %>%
      dplyr::group_by(water_mass_consensus) %>%
      dplyr::summarise(
        n          = dplyr::n(),
        Mean_Depth = round(mean(Depth, na.rm = TRUE), 0),
        SD_Depth   = round(sd(Depth,   na.rm = TRUE), 0),
        Mean_Lat   = round(mean(Lat,    na.rm = TRUE), 2),
        Lat_Range  = paste0(round(min(Lat, na.rm = TRUE), 1),
                            " to ",
                            round(max(Lat, na.rm = TRUE), 1)),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(n))
    print(knitr::kable(stats_tbl, format = "simple",
                       col.names = c("Water Mass", "n", "Mean Depth (m)",
                                     "SD Depth (m)", "Mean Lat (°)", "Lat Range (°)")))
  }
  
  # ── Ross Sea shelf breakdown ───────────────────────────────────────────────
  cat("\n=== SUMMARY TABLE 5: Ross Sea Shelf — Water Mass × Period ===\n")
  rs_cols <- c("AASW", "MCDW", "SW_fresh", "SW_salty", "ISW", "AABW")
  rs_wm <- df %>%
    dplyr::filter(is_ross_sea_shelf(Lat, Long360, Depth)) %>%
    dplyr::mutate(
      wm = factor(water_mass_consensus, levels = rs_cols)
    ) %>%
    dplyr::count(wm, .drop = FALSE) %>%
    dplyr::mutate(Pct = round(100 * n / sum(n), 1))
  print(knitr::kable(rs_wm, col.names = c("Water Mass", "n", "%"), format = "simple"))
}


# ==============================================================================
# SECTION 8 · COLOUR PALETTE & VISUALISATION
# ==============================================================================

# ── UPDATED: expanded palette including Ross Sea water masses ─────────────────
WM_COLORS <- c(
  "AABW"    = "#E8977A",   # salmon        (open ocean bottom water)
  "AAIW"    = "#A4C969",   # green         (intermediate water)
  "AASW"    = "#4CB5BD",   # teal          (surface water)
  "CDW"     = "#9D7BC7",   # purple        (open ocean deep water)
  "SAMW"    = "#F5C842",   # yellow        (mode water)
  # ── Ross Sea ────────────────────────────────────────────────────────────────
  "MCDW"    = "#E8612A",   # burnt orange  (modified CDW on shelf; Smith et al. Fig 3)
  "ISW"     = "#7FB3D3",   # ice blue      (ice shelf water)
  "SW_salty"= "#5B2D8E",   # deep purple   (dense shelf water; Fig 3)
  "SW_fresh"= "#B39DDB"    # lavender      (fresher shelf water variant)
)

# Ordered factor for plots (surface → deep, light → dense)
WM_ORDER <- c("AASW", "SAMW", "AAIW", "MCDW", "CDW",
              "SW_fresh", "SW_salty", "ISW", "AABW")

plot_water_mass_panels <- function(df, save_path = "water_mass_comparison.pdf") {
  
  make_panel <- function(data, col, title, show_legend = FALSE) {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = Lat, y = Depth,
                                            color = .data[[col]])) +
      ggplot2::geom_point(alpha = 0.7, size = 1.8) +
      ggplot2::scale_y_reverse() +
      ggplot2::scale_color_manual(
        values   = WM_COLORS,
        name     = "Water Mass",
        na.value = "#AAAAAA",
        limits   = WM_ORDER,
        drop     = FALSE
      ) +
      ggplot2::labs(x = "Latitude (°)", y = "Depth (m)", title = title) +
      ggplot2::theme_classic() +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    
    if (!show_legend) p <- p + ggplot2::theme(legend.position = "none")
    p
  }
  
  p_manual    <- make_panel(df, "water_mass_manual",      "Manual (lat/depth/lon)")
  p_argo      <- make_panel(df, "water_mass_argo_filled", "Argo (year-matched, gap-filled)")
  p_woa       <- make_panel(df, "water_mass_woa_filled",  "WOA (season-matched, gap-filled)")
  p_consensus <- make_panel(df, "water_mass_consensus",   "Consensus (≥2 agree)",
                            show_legend = TRUE)   # single legend on bottom-right panel
  
  combined <- (p_manual | p_argo) / (p_woa | p_consensus) +
    patchwork::plot_annotation(
      title = "Water Mass Classification – All Methods",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    )
  
  print(combined)
  ggplot2::ggsave(save_path, combined, width = 14, height = 10)
  cat(sprintf("\nPlot saved to: %s\n", save_path))
  
  combined
}


# ==============================================================================
# SECTION 9 · GEOTRACES IRON DATA (supplementary; run independently if needed)
# ==============================================================================

run_geotraces_iron <- function(nc_path = "Downloads/geotraces.nc") {
  
  library(ncdf4)
  library(tidyverse)
  
  nc_data <- ncdf4::nc_open(nc_path)
  
  lat     <- ncdf4::ncvar_get(nc_data, "latitude")
  lon     <- ncdf4::ncvar_get(nc_data, "longitude")
  iron    <- ncdf4::ncvar_get(nc_data, "Fe_D_CONC_BOTTLE")
  iron_qc <- ncdf4::ncvar_get(nc_data, "Fe_D_CONC_BOTTLE_qc")
  depth   <- ncdf4::ncvar_get(nc_data, "Bot_Depth")
  
  n_depth   <- dim(iron)[1]
  n_station <- dim(iron)[2]
  
  df_fe <- expand.grid(sample_id = seq_len(n_depth),
                       station_id = seq_len(n_station)) %>%
    dplyr::mutate(
      lat     = lat[station_id],
      lon     = lon[station_id],
      iron    = as.vector(iron),
      iron_qc = as.vector(iron_qc)
    ) %>%
    dplyr::filter(iron > -1e9, iron_qc == 49) %>%
    dplyr::mutate(lon_360 = ifelse(lon < 0, lon + 360, lon))
  
  assign_region <- function(lon360) {
    dplyr::case_when(
      lon360 >= 40  & lon360 < 150  ~ "East Antarctic",
      lon360 >= 160 & lon360 < 200  ~ "Ross",
      (lon360 >= 220 & lon360 <= 360) | (lon360 >= 0 & lon360 < 20) ~ "Weddell-Bell-Amundsen",
      TRUE ~ "Other"
    )
  }
  
  df_fe <- df_fe %>%
    dplyr::mutate(region = assign_region(lon_360)) %>%
    dplyr::filter(region != "Other")
  
  summary_stats <- df_fe %>%
    dplyr::group_by(region) %>%
    dplyr::summarise(
      n          = dplyr::n(),
      n_stations = dplyr::n_distinct(station_id),
      mean_Fe    = mean(iron, na.rm = TRUE),
      median_Fe  = median(iron, na.rm = TRUE),
      sd_Fe      = sd(iron, na.rm = TRUE),
      min_Fe     = min(iron, na.rm = TRUE),
      max_Fe     = max(iron, na.rm = TRUE),
      .groups    = "drop"
    )
  
  print(summary_stats)
  
  p_violin <- ggplot2::ggplot(df_fe, ggplot2::aes(x = region, y = iron, fill = region)) +
    ggplot2::geom_violin() +
    ggplot2::geom_boxplot(width = 0.1, alpha = 0.5) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(title = "Dissolved Iron Distribution by Region",
                  x = "Region", y = "Dissolved Fe (nmol/kg)") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "none")
  
  print(p_violin)
  list(data = df_fe, summary = summary_stats)
}


# ==============================================================================
# SECTION 10 · MAIN PIPELINE  ── run from here ──
# ==============================================================================

cat("===========================================================\n")
cat(" WATER MASS CLASSIFICATION PIPELINE  v2\n")
cat("===========================================================\n\n")

edge_file <- "~/Documents/Data/OtoData/TE/edge_FINAL.csv"
edge_df   <- readr::read_csv(edge_file, show_col_types = FALSE)

# ── 0. Derive Month if absent ────────────────────────────────────────────────
if (!"Month" %in% names(edge_df) && "Date" %in% names(edge_df)) {
  edge_df$Month <- as.integer(format(as.Date(edge_df$Date), "%m"))
  cat("Derived Month from Date column.\n")
} else if (!"Month" %in% names(edge_df)) {
  cat("Note: No Month or Date column found – WOA will use annual climatology.\n")
}

# ── 1. Manual classification ─────────────────────────────────────────────────
edge_df <- run_manual_classification(edge_df)

# ── 2. Argo classification (year-matched) ────────────────────────────────────
# Delete stale cache if you want a fresh download (uncomment below):
# if (file.exists("argo_profiles_cache.rds")) file.remove("argo_profiles_cache.rds")

profiles <- get_argo_profiles_for_years(
  df         = edge_df,
  lat_col    = "Lat",
  year_col   = "Year",
  buffer_deg = 2,
  cache_file = "argo_profiles_cache.rds"
)

edge_df <- run_argo_classification(
  df          = edge_df,
  profiles    = profiles,
  year_col    = "Year",
  max_dist_km = 300,
  year_window = 1
)

# ── 3. WOA classification (season-matched) ───────────────────────────────────
edge_df <- run_woa_classification(
  df        = edge_df,
  month_col = "Month",
  data_dir  = "woa_data"
)

# ── 4. Agreement analysis & manual gap-fill ──────────────────────────────────
edge_df <- run_agreement_and_gap_fill(
  df                  = edge_df,
  agreement_threshold = 0.80
)

# ── 5. Consensus classification ──────────────────────────────────────────────
edge_df <- compute_consensus(edge_df)

# ── 6. Summary tables ────────────────────────────────────────────────────────
print_summary_tables(edge_df)

# ── 7. Save outputs ──────────────────────────────────────────────────────────
write.csv(edge_df, "edge_rs_wm_FINAL.csv", row.names = FALSE)
cat("\n✓ Final classified dataset saved to: edge_water_mass_classified_FINAL.csv\n")

# ── 8. Plots ─────────────────────────────────────────────────────────────────
plot_water_mass_panels(edge_df, save_path = "water_mass_comparison.pdf")

cat("\n===========================================================\n")
cat(" PIPELINE COMPLETE\n")
cat("===========================================================\n")
cat("\nColumn guide:\n")
cat("  water_mass_manual       : Lat/depth/lon zones (always available)\n")
cat("  water_mass_argo         : Argo float T/S (year-matched; NAs = no float nearby)\n")
cat("  water_mass_argo_filled  : Argo + manual gap-fill where agreement ≥80%\n")
cat("  argo_fill_source        : 'argo' | 'manual_fill' | NA\n")
cat("  water_mass_woa          : WOA season T/S (NAs = missing climatology)\n")
cat("  water_mass_woa_filled   : WOA + manual gap-fill where agreement ≥80%\n")
cat("  woa_fill_source         : 'woa' | 'manual_fill' | NA\n")
cat("  water_mass_consensus    : Majority vote across all three methods\n")
cat("  lat_zone                : Geographic zone label (incl. Ross Sea Shelf)\n")
cat("\nRoss Sea categories (Long360 155-220°, Lat ≤ -65°, Depth ≤ 700m):\n")
cat("  AASW     : depth < 150 m\n")
cat("  MCDW     : T > -1.5°C at depth ≥ 150 m (Modified CDW intrusion)\n")
cat("  SW_fresh : T < -1.7°C, 34.3 ≤ S < 34.62 (fresh shelf water)\n")
cat("  SW_salty : T < -1.7°C, S ≥ 34.62 (dense brine-rejected shelf water)\n")
cat("  ISW      : T < -1.9°C, S < 34.52 (ice shelf water outflow)\n")
cat("  AABW     : depth > 700 m (exported via Drygalski/Joides/GC Troughs)\n")


# ==============================================================================
# SECTION 8b · ELEMENT:Ca BOXPLOTS WITH SIGNIFICANCE STARS
# ==============================================================================

library(rstatix)
library(ggpubr)
library(tidyr)

plot_element_ca_by_watermass <- function(
    df,
    element_cols   = c("raw_ratio_X135Ba", "raw_ratio_X23Na", "raw_ratio_X25Mg",
                       "raw_ratio_X55Mn",  "raw_ratio_X57Fe", "raw_ratio_X88Sr"),
    element_labels = c("Ba", "Na", "Mg", "Mn", "Fe", "Sr"),
    wm_col         = "water_mass_manual",
    ref_group      = "AABW",
    wm_order       = c("AABW", "CDW", "MCDW", "SW_salty", "SW_fresh",
                       "ISW", "AAIW", "AASW", "SAMW"),
    save_path      = "element_ca_watermass.pdf"
) {
  
  # ── UPDATED: expanded palette matching Section 8 ───────────────────────────
  LOCAL_WM_COLORS <- c(
    "AABW"    = "#E8977A",
    "CDW"     = "#9D7BC7",
    "MCDW"    = "#E8612A",
    "SW_salty"= "#5B2D8E",
    "SW_fresh"= "#B39DDB",
    "ISW"     = "#7FB3D3",
    "AAIW"    = "#A4C969",
    "AASW"    = "#4CB5BD",
    "SAMW"    = "#F5C842"
  )
  
  stopifnot(length(element_cols) == length(element_labels))
  missing_cols <- setdiff(c(element_cols, wm_col), names(df))
  if (length(missing_cols) > 0)
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  
  label_map <- setNames(element_labels, element_cols)
  
  plot_df <- df %>%
    dplyr::select(dplyr::all_of(c(wm_col, element_cols))) %>%
    dplyr::rename(WaterMass = dplyr::all_of(wm_col)) %>%
    dplyr::filter(!is.na(WaterMass)) %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(element_cols),
      names_to  = "Element",
      values_to = "Ratio"
    ) %>%
    dplyr::filter(!is.na(Ratio), Ratio > 0) %>%
    dplyr::mutate(
      log_Ratio = log10(Ratio),
      Element   = factor(label_map[Element], levels = element_labels),
      WaterMass = factor(WaterMass, levels = intersect(wm_order, unique(WaterMass)))
    )
  
  y_maxima <- plot_df %>%
    dplyr::group_by(Element, WaterMass) %>%
    dplyr::summarise(q95 = quantile(log_Ratio, 0.95, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(Element) %>%
    dplyr::summarise(y_top = max(q95, na.rm = TRUE), .groups = "drop")
  
  stat_results <- plot_df %>%
    dplyr::group_by(Element) %>%
    rstatix::wilcox_test(log_Ratio ~ WaterMass, p.adjust.method = "BH") %>%
    rstatix::add_significance(
      p.col     = "p.adj",
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols   = c("***", "**", "*", "ns")
    ) %>%
    dplyr::filter(p.adj.signif != "ns") %>%
    rstatix::add_xy_position(x = "WaterMass", fun = "max") %>%
    dplyr::left_join(y_maxima, by = "Element") %>%
    dplyr::group_by(Element) %>%
    dplyr::mutate(
      bracket_rank = dplyr::row_number(),
      y.position   = y_top + 0.12 * bracket_rank,
      xmin         = as.numeric(factor(group1, levels = levels(plot_df$WaterMass))),
      xmax         = as.numeric(factor(group2, levels = levels(plot_df$WaterMass)))
    ) %>%
    dplyr::ungroup()
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = WaterMass, y = log_Ratio, fill = WaterMass)
  ) +
    ggplot2::geom_boxplot(
      outlier.shape  = 16,
      outlier.size   = 0.9,
      outlier.alpha  = 0.45,
      outlier.colour = "grey40",
      linewidth      = 0.45,
      width          = 0.55,
      colour         = "grey20"
    ) +
    ggpubr::stat_pvalue_manual(
      stat_results,
      label          = "p.adj.signif",
      tip.length     = 0.005,
      bracket.size   = 0.45,
      size           = 3.8,
      hide.ns        = TRUE,
      remove.bracket = FALSE,
      inherit.aes    = FALSE
    ) +
    ggplot2::scale_fill_manual(values = LOCAL_WM_COLORS, guide = "none") +
    ggplot2::facet_wrap(~ Element, scales = "free_y", ncol = 3) +
    ggplot2::labs(
      x = "Water Mass",
      y = expression(log[10]~"(Element:Ca ratio)")
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      panel.background   = ggplot2::element_rect(fill = "white", colour = NA),
      panel.border       = ggplot2::element_rect(fill = NA, colour = "grey30",
                                                 linewidth = 0.5),
      panel.grid.major.y = ggplot2::element_line(colour = "grey93", linewidth = 0.3),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.spacing      = ggplot2::unit(0.8, "lines"),
      strip.text         = ggplot2::element_text(face = "bold", size = 12,
                                                 margin = ggplot2::margin(b = 4)),
      strip.background   = ggplot2::element_blank(),
      strip.placement    = "outside",
      axis.text.x        = ggplot2::element_text(size = 9, colour = "grey20",
                                                 angle = 35, hjust = 1),
      axis.text.y        = ggplot2::element_text(size = 9,  colour = "grey20"),
      axis.title         = ggplot2::element_text(size = 11),
      axis.line          = ggplot2::element_line(colour = "grey30", linewidth = 0.4),
      axis.ticks         = ggplot2::element_line(colour = "grey30", linewidth = 0.4),
      plot.background    = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin        = ggplot2::margin(10, 15, 10, 10)
    )
  
  ggplot2::ggsave(save_path, p, width = 13, height = 7.5, dpi = 300)
  cat(sprintf("Plot saved to: %s\n", save_path))
  
  cat("\n=== PAIRWISE SIGNIFICANCE (Wilcoxon, BH-adjusted) ===\n")
  print(
    stat_results %>%
      dplyr::select(Element, group1, group2, p.adj, p.adj.signif) %>%
      dplyr::arrange(Element, group2),
    n = Inf
  )
  
  invisible(p)
}

# ── Call ──────────────────────────────────────────────────────────────────────
plot_element_ca_by_watermass(
  df             = edge_df,
  element_cols   = c("raw_ratio_X135Ba", "raw_ratio_X23Na", "raw_ratio_X25Mg",
                     "raw_ratio_X55Mn",  "raw_ratio_X57Fe", "raw_ratio_X88Sr"),
  element_labels = c("Ba", "Na", "Mg", "Mn", "Fe", "Sr"),
  wm_col         = "water_mass_consensus",
  ref_group      = "AABW",
  save_path      = "element_ca_watermass_consensus.pdf"
)


# ==============================================================================
# SECTION 11 · HIGH-Li DEPTH × WATER MASS ANALYSIS (Ross Sea)
# ==============================================================================
# Connects the Li meltwater signal (from the glacial influx analysis) to the
# water mass domain, testing whether the 2000–2004 (B-15) and 2019–2024 (WAIS)
# events occupy distinct water masses as predicted by their spatial patterns.

plot_ross_li_depth_watermass <- function(
    df,
    li_col        = "conc_ppm_X7Li",
    li_threshold  = 0.84,
    save_path     = "ross_li_depth_watermass.pdf"
) {
  
  ross_highli <- df %>%
    dplyr::filter(
      Long360 >= 160 & Long360 <= 220,
      Lat     <= -70,
      !is.na(.data[[li_col]]),
      .data[[li_col]] > li_threshold
    ) %>%
    dplyr::mutate(
      period = dplyr::case_when(
        CatchYr >= 2000 & CatchYr <= 2004 ~ "2000–2004 (B-15)",
        CatchYr >= 2019 & CatchYr <= 2024 ~ "2019–2024 (WAIS)",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(period))
  
  # ── Summary table ──────────────────────────────────────────────────────────
  cat("=== HIGH-Li ROSS SEA: DEPTH × WATER MASS × PERIOD ===\n")
  summary_tbl <- ross_highli %>%
    dplyr::group_by(period, water_mass_consensus) %>%      # use pre-computed column
    dplyr::summarise(
      n          = dplyr::n(),
      mean_depth = round(mean(Depth, na.rm = TRUE), 0),
      sd_depth   = round(sd(Depth,   na.rm = TRUE), 0),
      mean_Lat   = round(mean(Lat, na.rm = TRUE), 0),
      sd_Lat     = round(sd(Lat,   na.rm = TRUE), 0),
      mean_Long360 = round(mean(Long360, na.rm = TRUE), 0),
      sd_Long360   = round(sd(Long360,   na.rm = TRUE), 0),
      mean_Li    = round(mean(.data[[li_col]], na.rm = TRUE), 3),
      .groups    = "drop"
    ) %>%
    dplyr::arrange(period, water_mass_consensus)           # fixed: was wm_manual
  print(summary_tbl, n = Inf)
  
  # ── Plot ───────────────────────────────────────────────────────────────────
  p <- ross_highli %>%
    dplyr::filter(!is.na(water_mass_consensus)) %>%        # fixed: was wm_manual
    ggplot2::ggplot(ggplot2::aes(
      x      = .data[[li_col]],
      y      = Depth,
      colour = water_mass_consensus,                       # fixed: was wm_manual
      shape  = period
    )) +
    ggplot2::geom_point(alpha = 0.75, size = 2.8) +
    ggplot2::scale_y_reverse(name = "Capture depth (m)") +
    ggplot2::scale_x_continuous(name = "Otolith Li (ppm)") +
    ggplot2::scale_colour_manual(
      values = WM_COLORS,
      name   = "Water mass",
      limits = WM_ORDER,                                   # consistent legend
      drop   = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = c("2000–2004 (B-15)" = 16, "2019–2024 (WAIS)" = 17),
      name   = "Period"
    ) +
    ggplot2::geom_vline(
      xintercept = li_threshold,
      linetype = "dashed", colour = "grey40", linewidth = 0.4
    ) +
    ggplot2::facet_wrap(~ period, ncol = 2) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(title = "High-Li Ross Sea fish: depth × water mass × period")
  
  ggplot2::ggsave(save_path, p, width = 10, height = 5.5, dpi = 300)
  cat(sprintf("\nPlot saved to: %s\n", save_path))
  
  invisible(list(plot = p, data = ross_highli, summary = summary_tbl))
}

plot_ross_li_depth_watermass(
  df           = edge_df,
  li_col       = "conc_ppm_X7Li",
  li_threshold = 0.84,
  save_path    = "ross_li_depth_watermass.pdf"
) %>%
  print(width = Inf)
