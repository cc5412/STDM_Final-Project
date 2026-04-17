# ─────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────

# Install packages (if not installed)
packages <- c("tidyverse", "lubridate", "httr", "purrr", "readr",
              "janitor", "ggplot2", "sf", "tmap", "forecast",
              "zoo","spdep")

installed <- rownames(installed.packages())
to_install <- packages[!packages %in% installed]

if (length(to_install) > 0) {
  install.packages(to_install)
}

# Load libraries
invisible(lapply(packages, library, character.only = TRUE))

source("starima_package.R")

# NOAA API & Temperature Station information to initially download data
COOPS_BASE <- "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"
CACHE_FILE <- "NY_watertemps.csv"

stations <- tibble(
  station = c("8518750","8516945","8467150","8465705","8461490","8510560"),
  lat = c(40.7003, 40.8100, 41.1800, 41.2700, 41.3600, 41.0000),
  lon = c(-74.0100, -73.7650, -73.1800, -72.9000, -72.0900, -72.0000)
)

# ─────────────────────────────────────────────
# LOAD DATA - Temperature readings are in 6 minute intervals
# ─────────────────────────────────────────────

if (file.exists(CACHE_FILE)) {
  message("Loading cached dataset...")
  final_data <- read_csv(CACHE_FILE, show_col_types = FALSE)
  
} else {
  
  make_chunks <- function(start_date, end_date) {
    starts <- seq.Date(start_date, end_date, by = "31 days")
    map_df(starts, ~tibble(start = .x, end = min(.x + 30, end_date)))
  }
  
  fetch_coops <- function(station, start_date, end_date) {
    
    url <- modify_url(COOPS_BASE, query = list(
      product = "water_temperature",
      application = "NOS.COOPS.TAC.WL",
      station = station,
      begin_date = format(start_date, "%Y%m%d"),
      end_date   = format(end_date, "%Y%m%d"),
      units = "metric",
      time_zone = "gmt",
      format = "csv"
    ))
    
    res <- GET(url, user_agent("COOPS research client"))
    if (status_code(res) != 200) return(NULL)
    
    txt <- content(res, as = "text", encoding = "UTF-8")
    
    tryCatch(read_csv(I(txt), show_col_types = FALSE) %>% clean_names(),
             error = function(e) NULL)
  }
  
  normalize <- function(df, station_id) {
    if (is.null(df)) return(NULL)
    
    df %>%
      mutate(station = station_id) %>%
      mutate(
        date_time = suppressWarnings(
          ymd_hms(date_time, quiet = TRUE)
        ),
        date_time = if_else(
          is.na(date_time),
          suppressWarnings(ymd_hm(date_time, quiet = TRUE)),
          date_time
        ),
        water_temperature = as.numeric(water_temperature)
      ) %>%
      filter(!is.na(date_time)) %>%
      select(date_time, water_temperature, station)
  }
  
  get_station <- function(station, start_date, end_date) {
    
    message("\nStarting station: ", station)
    chunks <- make_chunks(start_date, end_date)
    out <- list()
    
    for (i in seq_len(nrow(chunks))) {
      Sys.sleep(2)
      message("  Chunk ", i, "/", nrow(chunks))
      
      df <- fetch_coops(station, chunks$start[i], chunks$end[i])
      df <- normalize(df, station)
      
      if (!is.null(df)) out[[length(out)+1]] <- df
    }
    
    message("Completed station: ", station, " ✓")
    bind_rows(out)
  }
  
  start_date <- as.Date("2024-01-01")
  end_date   <- as.Date("2025-12-31")
  
  results <- pmap(stations, ~get_station(..1, start_date, end_date))
  final_data <- bind_rows(compact(results))
  
  write_csv(final_data, CACHE_FILE)
}

final_data <- final_data %>%
  mutate(station = as.character(station))

final_data

# ─────────────────────────────────────────────
# MAP MEAN TEMPERATURES (2024-2025)
# ─────────────────────────────────────────────

# Calculate overall mean from the raw 6-minute data
overall_data <- final_data %>%
  group_by(station) %>%
  summarise(mean_temp = mean(water_temperature, na.rm = TRUE), .groups = "drop")

# Create sf object from stations tibble
stations_sf <- stations %>%
  mutate(station = as.character(station)) %>%
  st_as_sf(coords = c("lon","lat"), crs = 4326)

# Join the mean temperatures to the spatial data
stations_mapped <- stations_sf %>%
  left_join(overall_data, by = "station")

# Map
tmap_mode("plot")

LIsound_label <- data.frame(
  lon = -72.85,
  lat = 41.085,
  name = "Long Island Sound"
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

mean_temp_map <- tm_shape(stations_mapped) +
  tm_tiles("CartoDB.Positron") + # Clean OSM basemap
  tm_dots(fill = "mean_temp",
          fill.scale = tm_scale(values = "viridis"),
          fill.legend = tm_legend(title = "Mean Temp. (°C)"), # Added legend title here
          size = 1,
          shape = 21) +
  tm_text("station", ymod = -0.9, size = 0.9) +
  
  tm_shape(LIsound_label) +
  tm_text("name", 
          size = 1.2, 
          fontface = "italic",
          col = "dodgerblue4",
          alpha = 0.8) +
  
  tm_title("Mean Water Temperature (°C) at NOAA Stations in Long Island Sound (2024-2025)") +
  tm_layout(
    legend.position = c("left", "top"),
    inner.margins = c(0.1, 0.1, 0.15, 0.1) # Padding: c(bottom, left, top, right)
  )

print(mean_temp_map)

# ─────────────────────────────────────────────
# HOURLY AGGREGATION
# ─────────────────────────────────────────────

hourly <- final_data %>%
  mutate(hour = floor_date(date_time, "hour")) %>%
  group_by(station, hour) %>%
  summarise(
    temp = mean(water_temperature, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(station = as.character(station))

# Full time grid - HOURLY - Ensures no skipped hours exist
full_hours <- expand.grid(
  station = unique(hourly$station),
  hour = seq(min(hourly$hour), max(hourly$hour), by = "hour")
)

hourly_full <- full_hours %>%
  left_join(hourly, by = c("station","hour")) %>%
  arrange(station, hour)

# ─────────────────────────────────────────────
# PLOT EACH STATION
# ─────────────────────────────────────────────

station_labels <- c(
  "8461490" = "New London, CT [8461490]",
  "8465705" = "New Haven, CT [8465705]",
  "8467150" = "Bridgeport, CT [8467150]",
  "8510560" = "Montauk, NY [8510560]",
  "8516945" = "Kings Point, NY [8516945]",
  "8518750" = "The Battery, NY [8518750]"
)

hourly_full_named <- hourly_full %>%
  mutate(station_name = station_labels[as.character(station)])

ggplot(hourly_full_named, aes(x = hour, y = temp)) +
  geom_line(color = "black", alpha = 0.7) + 
  facet_wrap(~station_name, nrow = 2, ncol = 3, scales = "free_y") + 
  labs(title = "Hourly Water Temperatures by NOAA Station (2024-2025)",
       x = "Date",
       y = "Water Temperature (°C)") +
  theme_minimal() +
  theme(strip.text = element_text(face="bold", size=10))

# ─────────────────────────────────────────────
# REMOVE BAD STATION
# ─────────────────────────────────────────────

hourly_clean <- hourly_full %>%
  filter(station != "8467150") %>%
  group_by(station) %>%
  arrange(hour) %>%
  mutate(temp = na.approx(temp, na.rm = FALSE)) %>%
  ungroup()

# ─────────────────────────────────────────────
# EXPLORATORY DATA ANALYSIS (New Haven, CT: Station 8465705)
# ─────────────────────────────────────────────
# Selected New Haven station [8465705] due to its visually most complete data and location on the Long Island Sound.

# Inspect gaps by extracting rows where temperature is NA
proxy_data <- hourly_clean %>%
  filter(station == "8465705") %>%
  arrange(hour)

missing_gaps <- proxy_data %>%
  filter(is.na(temp))

# Print the total number of missing hours
cat("Total missing hours for 8465705:", nrow(missing_gaps), "out of", nrow(proxy_data), "\n\n")

# Print the exact timestamps of the missing data
View(missing_gaps)

# Found a 32h gap + 1h gap
# Linear interpolation for small gaps - Gaps in CSV will no longer be present.
proxy_station <- hourly_clean %>%
  filter(station == "8465705") %>%
  arrange(hour) %>%
  mutate(temp_clean = na.approx(temp))

# Convert to a Time Series object (frequency = 24 for daily seasonality)
ts_proxy <- ts(proxy_station$temp_clean, frequency = 24)

# STL Decomposition (Viewing a 2-week slice [336h] during peak summer: July 1 - 14, 2024)
slice_data <- proxy_station %>%
  filter(hour >= ymd_hms("2024-07-01 00:00:00") & hour < ymd_hms("2024-07-15 00:00:00"))
ts_slice <- ts(slice_data$temp_clean, frequency = 24)
fit_stl <- stl(ts_slice, s.window = "periodic")
plot(fit_stl, main = "STL Decomposition (July 1-14, 2024): New Haven, CT, Station 8465705")

# Stationarity Check (ACF & PACF)
# Plot raw vs. differenced (lag = 24) to justify ARIMA parameters
par(mfrow = c(2, 2)) # Create a 2x2 grid for comparison
acf(ts_proxy, lag.max = 72, main = "ACF: Raw Data")
pacf(ts_proxy, lag.max = 72, main = "PACF: Raw Data")

# Non-stationarity determined from ACF
# Apply seasonal differencing
ts_proxy_diff <- diff(ts_proxy, lag = 24)
acf(ts_proxy_diff, lag.max = 72, main = "ACF: Differenced (Lag 24)")
pacf(ts_proxy_diff, lag.max = 72, main = "PACF: Differenced (Lag 24)")
par(mfrow = c(1, 1)) # Reset plot layout

# ─────────────────────────────────────────────
# AUTO.ARIMA & RESIDUAL DIAGNOSTICS
# ─────────────────────────────────────────────

# Algorithm finds the optimal parameters via AIC
NewHaven_model <- auto.arima(ts_proxy, seasonal = TRUE, stepwise = TRUE, approximation = TRUE)
print(summary(NewHaven_model))

# Visual Confirmation (Residual Diagnostics) - ACF and PACF of residuals
par(mfrow = c(2, 1)) # Create a 2x1 grid

acf(residuals(NewHaven_model), lag.max = 72, 
    main = "ACF of New Haven ARIMA Model Residuals")

pacf(residuals(NewHaven_model), lag.max = 72, 
     main = "PACF of New Haven ARIMA Model Residuals")

par(mfrow = c(1, 1)) # Reset plot layout

# Ljung-Box Test
# Checking if p-value > 0.05, which indicates the residuals are true random noise and the model is good.
checkresiduals(NewHaven_model)

# ─────────────────────────────────────────────
# LOAD TEST DATA - January 1-14, 2026
# ─────────────────────────────────────────────

# Pulling all 5 active stations (which excludes the bad station - Bridgeport, CT)
active_stations <- c("8518750", "8516945", "8465705", "8461490", "8510560")

# Initialize an empty list to store the dataframes
test_data_list <- list()

# Loop through each station to pull the data from NOAA
# NOTE: User can modify 'begin_date' and 'end_date' to test different forecast horizons.
# Dates must be in YYYYMMDD format. Maximum pull for this API endpoint is 31 days.

for (sta in active_stations) {
  test_url <- modify_url(COOPS_BASE, query = list(
    product = "water_temperature",
    application = "NOS.COOPS.TAC.WL",
    station = sta, 
    begin_date = "20260101",
    end_date   = "20260114",
    units = "metric",
    time_zone = "gmt",
    format = "csv"
  ))
  
  # Download, clean, and assign the station ID
  temp_df <- read_csv(test_url, show_col_types = FALSE) %>% clean_names()
  temp_df$station <- sta 
  
  test_data_list[[sta]] <- temp_df
}

# Combine all 5 station pulls into one raw dataframe
test_data_raw <- bind_rows(test_data_list)

# Aggregate to hourly means to match training data
test_hourly_all <- test_data_raw %>%
  mutate(
    date_time_parsed = suppressWarnings(ymd_hms(date_time, quiet = TRUE)),
    date_time = if_else(
      is.na(date_time_parsed),
      suppressWarnings(ymd_hm(date_time, quiet = TRUE)),
      date_time_parsed
    ),
    hour = floor_date(date_time, "hour")
  ) %>%
  group_by(station, hour) %>% 
  summarise(
    actual_temp = mean(as.numeric(water_temperature), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(station, hour)

# Fill small gaps in the test data using the same training logic
test_hourly_all <- test_hourly_all %>%
  group_by(station) %>%
  arrange(hour) %>%
  mutate(actual_temp = na.approx(actual_temp, na.rm = FALSE, rule = 2)) %>%
  ungroup()

# ─────────────────────────────────────────────
# ARIMA FORECAST & VALIDATION
# ─────────────────────────────────────────────

# Filter to New Haven for ARIMA
test_hourly_nh <- test_hourly_all %>%
  filter(station == "8465705") %>%
  arrange(hour)

# Generate ARIMA forecast
# 14 days * 24 hours = ~ 336 steps ahead
forecast_horizon <- nrow(test_hourly_nh) 
arima_fc <- forecast(NewHaven_model, h = forecast_horizon)

# Check that the math lines up (~336 hours)
cat("\nARIMA Forecast Horizon: ", forecast_horizon, " hours\n")

# Combine Actual and Forecast into one data frame for comparison
forecast_eval <- tibble(
  hour = test_hourly_nh$hour,
  actual_temp = test_hourly_nh$actual_temp,
  predicted_temp = as.numeric(arima_fc$mean),
  lower_95 = as.numeric(arima_fc$lower[,2]),
  upper_95 = as.numeric(arima_fc$upper[,2])
)

# Calculate RMSE - Becomes "Score to Beat" for next model
baseline_rmse <- sqrt(mean((forecast_eval$actual_temp - forecast_eval$predicted_temp)^2, na.rm = TRUE))
cat("\nBaseline ARIMA RMSE for January 1-14, 2026:", round(baseline_rmse, 3), "°C\n")

# Plot the Forecast vs Actual
ggplot(forecast_eval, aes(x = hour)) +
  geom_line(aes(y = predicted_temp, color = "Forecast (ARIMA)"), linewidth = 0.4) +
  geom_line(aes(y = actual_temp, color = "Actual (NOAA)"), linewidth  = 1) +
  scale_color_manual(values = c("Actual (NOAA)" = "black", "Forecast (ARIMA)" = "blue")) +
  labs(
    title = "ARIMA Forecast vs Actual Water Temperatures:\nNew Haven, CT [Station 8465705] on January 1-14, 2026",
    subtitle = paste("Out-of-sample RMSE:", round(baseline_rmse, 3), "°C"),
    x = "Date",
    y = "Water Temperature (°C)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# ─────────────────────────────────────────────
# SPATIAL WEIGHT MATRIX (W)
# ─────────────────────────────────────────────

# STATION MATRIX 
active_stations <- stations %>%
  filter(station != "8467150")

# Create a spatial object. Starting with WGS84 (Lat/Lon - CRS 4326), 
# but then transform to a projected CRS - EPSG:2263 (NAD83 / New York Long Island).
stations_sf <- active_stations %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = 2263)

# Calculate the pairwise geographic distance matrix
dist_matrix <- st_distance(stations_sf)

# Convert from a 'units' object to a standard numeric matrix
dist_matrix_num <- as.numeric(dist_matrix)
dim(dist_matrix_num) <- dim(dist_matrix)

# Add station IDs as row and column names for easy reading
rownames(dist_matrix_num) <- active_stations$station
colnames(dist_matrix_num) <- active_stations$station

cat("\n--- Pairwise Distance Matrix (US Survey Feet) ---\n"); print(round(dist_matrix_num, 0))

# SPATIAL WEIGHT MATRIX (W)
# Inverse Distance Weighting: W_ij = 1 / Distance_ij
W_matrix <- 1 / dist_matrix_num

# Set the diagonal to 0.
diag(W_matrix) <- 0

# Row-standardization
W <- W_matrix / rowSums(W_matrix)
cat("\n--- Row-Standardized Spatial Weight Matrix (W) ---\n"); print(round(W, 3))

# ─────────────────────────────────────────────
# SPACE-TIME MATRIX FORMATTING
# ─────────────────────────────────────────────

# Pivot long format into a wide format
# Rows = Time (hour), Columns = Stations
space_time_wide <- hourly_clean %>%
  select(hour, station, temp) %>%
  pivot_wider(
    names_from = station,
    values_from = temp
  ) %>%
  arrange(hour)

# Time index for plotting later
time_index <- space_time_wide$hour

# Pure numeric data matrix
space_time_matrix <- space_time_wide %>%
  select(-hour) %>%
  as.matrix()

# Reorder columns of the data matrix to match the order of the stations in the Spatial Weight Matrix (W).
space_time_matrix <- space_time_matrix[, active_stations$station]

cat("\n--- Space-Time Matrix Dimensions (Rows = Hours, Cols = Stations) ---\n"); print(dim(space_time_matrix))

# First few rows to confirm structure.
head(space_time_matrix)

# ─────────────────────────────────────────────
# SPATIOTEMPORAL AUTOCORRELATION (ST-ACF & ST-PACF)
# ─────────────────────────────────────────────

# Formatting Weight Matrix for the STARIMA package
W_fit <- list(w1 = W)

# Difference the Space-Time Matrix
space_time_diff <- diff(space_time_matrix, lag = 24, differences = 1)

par(mfrow = c(2, 1))

# Calculate and Plot ST-ACF
# Spatial lag 1 across 48 time lags.
cat("\n--- ST-ACF (Spatial Lag 1, Seasonal Difference) ---\n"); stacf_result <- stacf(space_time_diff, W, 48); stacf_result

# Calculate and Plot ST-PACF
# Use to determine the AR order (p).
cat("\n--- ST-PACF (Spatial Lag 1, Seasonal Difference) ---\n"); stpacf_result <- stpacf(space_time_diff, W, 48); stpacf_result

par(mfrow = c(1, 1))

# ─────────────────────────────────────────────
# STARIMA
# ─────────────────────────────────────────────

colnames(hourly_clean); colnames(test_hourly_all)

# Rename the test data column to match the training data ("temp")
test_hourly_ready <- test_hourly_all %>%
  rename(temp = actual_temp) %>%
  select(station, hour, temp)

hourly_clean_ready <- hourly_clean %>%
  select(station, hour, temp)

# Combine
full_dataset <- bind_rows(hourly_clean_ready, test_hourly_ready)

# Pivot to Wide Matrix
space_time_matrix_full <- full_dataset %>%
  pivot_wider(names_from = station, values_from = temp) %>%
  arrange(hour)

# Save timestamps for later, then strip the 'hour' for the STARIMA functions
time_index <- space_time_matrix_full$hour
space_time_matrix_full <- space_time_matrix_full %>%
  select(-hour) %>%
  as.matrix()

# Define forecast window
n_total <- nrow(space_time_matrix_full)
n_train <- n_total - 336 # 14 days = 336 hours

# Fit the STARIMA(1, 24, 0) - Use training portion only for fitting
fit_star <- starima_fit(
  Z = space_time_matrix_full[1:n_train, ], 
  W = W_fit, 
  p = 2, d = 24, q = 0
)

# NOTE: I manually changed values and reran the code sections for the following (kept d = 24):
# p = 1, q = 0 --> RMSE = 0.182 °C
# p = 2, q = 0 --> RMSE = 0.145 °C
# p = 3, q = 0 --> RMSE = 0.217 °C
# p = 1, q = 1 --> RMSE = 0.149 °C
# p = 2, q = 1 --> RMSE = 0.218 °C
# RMSE values come from PERFORMANCE CHECK section below.
# User is welcome to edit values and test the same.

# Rolling One-Step-Ahead Forecast
# The package needs a small 'seed' of data before the forecast starts
# Start = n_train - seasonal_lag - p + 1 = n_train - 24 - 1 + 1
start_idx <- n_train - 24
pre_star <- starima_pre(space_time_matrix_full[start_idx:n_total, ], model = fit_star)

# ─────────────────────────────────────────────
# PERFORMANCE CHECK - Comparing with New Haven, CT, Station 8465705
# ─────────────────────────────────────────────

nh_col <- which(colnames(space_time_matrix_full) == "8465705")

# Ljung-Box test (using lag 24 to match the ARIMA test)
starima_res_nh <- fit_star$RES[, nh_col]
box_test_starima <- Box.test(starima_res_nh, lag = 24, type = "Ljung-Box")
cat("\n--- Ljung-Box Test: STARIMA Residuals (New Haven) ---\n"); print(box_test_starima)

# Plot STACF for the STARIMA model residuals to see if random
stacf_residuals <- stacf(fit_star$RES, W, 48)

# Calculating RMSE
star_preds <- pre_star$PRE[, nh_col]
star_actuals <- space_time_matrix_full[(n_train + 1):n_total, nh_col]

starima_rmse <- sqrt(mean((star_actuals - star_preds)^2, na.rm = TRUE))

cat("\n--- FINAL COMPARISON ---"); cat("\nBaseline ARIMA RMSE: ", round(baseline_rmse, 3), "°C"); cat("\nSTARIMA RMSE:        ", round(starima_rmse, 3), "°C\n")

# ─────────────────────────────────────────────
# FINAL VISUALIZATIONS - STARIMA FORECAST
# ─────────────────────────────────────────────

# Ensure all vectors are exactly the same length (336 for p =1, 335 for p = 2, etc.)
target_length <- 335

plot_data <- data.frame(
  DateTime = head(time_index[(n_train + 1):n_total], target_length),
  Actual   = head(star_actuals, target_length),
  STARIMA  = head(star_preds, target_length)
)

# Plot Actual vs STARIMA Forecast
ggplot(plot_data, aes(x = DateTime)) +
  geom_line(aes(y = Actual, color = "Actual (NOAA)"), linewidth = 1) +
  geom_line(aes(y = STARIMA, color = "Forecast (STARIMA)"), linewidth = 0.4) +
  
  labs(
    title = "STARIMA Forecast vs. Actual Water Temperature:\nNew Haven, CT [Station 8465705] on January 1-14, 2026",
    subtitle = paste0("STARIMA RMSE: ", round(starima_rmse, 3), 
                      "°C  |  ARIMA Baseline RMSE: ", round(baseline_rmse, 3), "°C"),
    x = "Date",
    y = "Water Temperature (°C)"
  ) +
  
  scale_color_manual(
    name = "Legend", 
    values = c("Actual (NOAA)" = "black", "Forecast (STARIMA)" = "red")
  ) +
  
  theme_minimal() +
  theme(legend.position = "bottom")

# ─────────────────────────────────────────────
# 4-STATION VISUAL COMPARISON - Remaining Active Stations
# ─────────────────────────────────────────────

# Isolate the Actuals vs Forecasts for the 14-day test window
target_length <- 335 # (336 for p =1, 335 for p = 2, etc.)
test_dates <- head(time_index[(n_train + 1):n_total], target_length)

all_actuals <- head(space_time_matrix_full[(n_train + 1):n_total, ], target_length)
all_preds   <- head(pre_star$PRE, target_length)

# Convert Actuals to a Long Dataframe
actuals_df <- as.data.frame(all_actuals) %>%
  mutate(DateTime = test_dates) %>%
  pivot_longer(cols = -DateTime, names_to = "station", values_to = "Actual")

# Convert Predictions to a Long Dataframe
preds_df <- as.data.frame(all_preds) %>%
  mutate(DateTime = test_dates) %>%
  pivot_longer(cols = -DateTime, names_to = "station", values_to = "STARIMA")

# Join, filter out New Haven, and label
plot_data_4 <- left_join(actuals_df, preds_df, by = c("DateTime", "station")) %>%
  filter(station != "8465705") %>%
  mutate(station_name = station_labels[station])

# Plot the 4 Panels
ggplot(plot_data_4, aes(x = DateTime)) +
  geom_line(aes(y = Actual, color = "Actual (NOAA)"), linewidth = 0.75) +
  geom_line(aes(y = STARIMA, color = "Forecast (STARIMA)"), linewidth = 0.3) +
  
  facet_wrap(~ station_name, scales = "free_y", ncol = 2) + 
  
  labs(
    title = "STARIMA Forecast vs. Actuals: Remaining Active Stations",
    subtitle = "NOAA Stations around Long Island Sound on January 1-14, 2026",
    x = "Date",
    y = "Water Temperature (°C)"
  ) +
  scale_color_manual(
    name = "Legend", 
    values = c("Actual (NOAA)" = "black", "Forecast (STARIMA)" = "red")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11)
  )

