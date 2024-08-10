library(tidynab)

df_raw <- tidynab::load_data()
df <- tidynab::parse_data(df_raw)

set.seed(1)
df_station <- tidynab::geolocate_stations() %>%
  drop_na(lat, lon) %>%
  rename(station = id) %>%
  sample_n(10)

df_raw <- tidynab::load_data()
df <- tidynab::parse_data(df_raw) %>%
  rename(station = stationid) %>%
  filter(station %in% df_station$station) %>%
  filter(taxa %in% c("Acer", "Quercus", "Populus", "Betula", "Ulmus", "Fraxinus", "Morus")) %>%
  mutate(date = lubridate::ymd(date)) %>%
  spread(key = "taxa", value = "count")


# download daymet data for all lat lon coordinates
ls_df_daymet <- vector("list", length = nrow(df_station))
for (i in 1:nrow(df_station)) {
  date_range <- df %>%
    filter(station == df_station$station[i]) %>%
    pull(date) %>%
    range() %>%
    lubridate::year()
  ls_df_daymet[[i]] <- daymetr::download_daymet(
    lat = df_station$lat[i],
    lon = df_station$lon[i],
    start = date_range[1],
    end = date_range[2]
  )$data %>%
    mutate(date = lubridate::ymd(str_c(year, "-01-01")) + yday - 1) %>%
    select(-year, -yday) %>%
    mutate(station = df_station$station[i]) %>%
    select(date, station, everything())
}
df_daymet <- bind_rows(ls_df_daymet)
