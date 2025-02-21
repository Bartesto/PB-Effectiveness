## Single AOI Fuel Age


# Packages ----------------------------------------------------------------

# NOTE FireHistory is an R package which can be installed from the dbca-wa 
# github repository (for more package details see https://github.com/dbca-wa/FireHistory).
# The script below uses some of its utility functions.

# To install FireHistory, uncomment and run the below:

# install.packages("devtools")
# devtools::install_github("dbca-wa/FireHistory")

library(FireHistory)
library(terra)
library(tidyterra)
library(sf)
library(dplyr)
library(readr)
library(glue)
library(tools)
library(fs)


# User Input Section ------------------------------------------------------

# Different pixel resolutions will affect results and processing times. Over time 
# the source data for mapping the fires in DBCA_060 has changed and of particular 
# note, mapping from satellite data has had 2 different pixel sizes relating to 
# the type of sensor employed. Landsat at 30m and Sentinel at 10 or 20m has been 
# used. Whilst the user has free range to choose whatever pixel resolution they 
# like for the conversion of vectors to raster, a good middle ground is 20m.

# Please input required pixel resolution:
resolution <- 20

# Please provide complete file path to the DBCA_060 data set
dbca_060_path <- "E:/pb_effectiveness/dbca_060/DBCA_Fire_History_DBCA_060.shp"

# Please input complete file path to your AOI shapefile
aoi_shp_path <- "E:/pb_effectiveness/vectors/blackwood_ten_diss.shp"

# Please provide file path to desired output location
out_path <- "E:/pb_effectiveness"


aoi_name <- tools::file_path_sans_ext(basename(aoi_shp_path))

# Assemble Data -----------------------------------------------------------

# takes users aoi and formats it for use (also projects it to EPSG:7844)
aoi <- FireHistory::user_aoi(aoi_shp_path, name = aoi_name)

# subsets to aoi extent on read in the very large DBCA_060 for dates FYfrom and FYto
data <- FireHistory::assemble_data(fire_path = dbca_060_path,
                                   FYfrom = 1900, FYto = 2024, aoi = aoi,
                                   accessed_on = "02/09/2024")# date that you downloaded DBCA_060


# Create Folder Structure  ------------------------------------------------

# Folder structure for outputs
cli::cli_progress_step("Creating Folder Structure")


tdir <- fs::path(out_path, paste0(aoi_name, "_res_",resolution, "/"))

by_fol <- fs::path(tdir, "01_by")
yob_fol <- fs::path(tdir, "02_yob")
tsf_fol <- fs::path(tdir, "03_tsf")
wfm_fol <- fs::path(tdir, "04_wfm")
wffa_fol <- fs::path(tdir, "05_wffa")
ofm_fol <- fs::path(tdir, "04_ofm")
offa_fol <- fs::path(tdir, "05_offa")

fols <- c(by_fol, yob_fol, tsf_fol, wfm_fol, wffa_fol, ofm_fol, offa_fol)
fs::dir_create(fols)


# Create Individual Burn Years - byYYYY -----------------------------------

cli::cli_progress_step("Creating Annual Burn Year Rasters")

# initial data only subset to extent - requires further cropping
f_vecs <- terra::vect(data[["fh_alb"]]) %>%
  terra::crop(terra::vect(data$aoi_alb)) %>%
  dplyr::arrange(fin_y)

# unique fire years
u_fyrs <- f_vecs %>%
  tidyterra::pull(var = fin_y) %>%
  unique()

# creates appropriate raster template
template <- terra::rast(data[["aoi_alb"]], res = resolution)
aoi_msk <- terra::rasterize(terra::vect(data[["aoi_alb"]]), template)

for(i in seq_along(u_fyrs)){
  by <- f_vecs %>%
    tidyterra::filter(fin_y == u_fyrs[i]) %>%
    terra::rasterize(template, field = "fin_y") %>%
    terra::crop(aoi_msk, mask = TRUE)
  names(by) <- u_fyrs[i]
  nom <- paste0("by",  u_fyrs[i], ".tif")
  by_nom <- fs::path(by_fol, nom)
  terra::writeRaster(by, by_nom)
}

# find non burn years and add an infill blank year
zero_rst <- template
terra::values(zero_rst) <- NaN #0
zero_yr <- zero_rst %>%
  terra::crop(aoi_msk, mask = TRUE)
minyr <- min(u_fyrs)
maxyr <- max(u_fyrs)
full_fyrs <- minyr:maxyr
missing <- full_fyrs[!(full_fyrs %in% u_fyrs)]

for(i in seq_along(missing)){
  nom <- glue::glue("by",  missing[i], ".tif")
  my_nom <- fs::path(by_fol, nom)
  terra::writeRaster(zero_yr, my_nom)
}


# Create Individual Year of Burn - yobYYYY --------------------------------

cli::cli_progress_step("Creating Year of Burn Rasters")

yob_noms <- fs::path(yob_fol, gsub("by", "yob", dir(by_fol)))

# start year same as first burn year
fs::file_copy(path = fs::dir_ls(by_fol)[1], new_path = yob_noms[1])

# 2nd yob
rst1 <- terra::rast(fs::dir_ls(by_fol)[1])
rst2 <- terra::rast(fs::dir_ls(by_fol)[2])
by_stk <- c(rst1, rst2)
yob <- terra::app(by_stk, fun ="max", na.rm = TRUE)
yob_nom <- yob_noms[2]
terra::writeRaster(yob, yob_nom)

# subsequent yob
yob_iter <- length(full_fyrs)-1
for(i in 2:yob_iter){
  if(i < yob_iter){
    rst1 <- terra::rast(yob_noms[i])
    rst2 <- terra::rast(fs::dir_ls(by_fol)[i+1])
    by_stk <- c(rst1, rst2)
    yob <- terra::app(by_stk, fun ="max", na.rm = TRUE)
    terra::writeRaster(yob, yob_noms[i+1])
  } else {
    rst1 <- terra::rast(yob_noms[i])
    rst2 <- terra::rast(fs::dir_ls(by_fol)[i+1])
    by_stk <- c(rst1, rst2)
    yob <- terra::app(by_stk, fun ="max", na.rm = TRUE)
    yob[yob == 0] <- NA
    terra::writeRaster(yob, yob_noms[i+1])
  }
}


# Create Individual Time Since Fire - tsfYYYY -----------------------------

cli::cli_progress_step("Creating Annual Time Since Fire Rasters")

tsf_noms <- fs::path(tsf_fol, gsub("yob", "tsf", dir(yob_fol)))

for(i in seq_along(tsf_noms)){
  yob <- terra::rast(fs::dir_ls(yob_fol)[i])
  tsf <- full_fyrs[i] - yob
  terra::writeRaster(tsf, tsf_noms[i])
}


# Create All Fires Fuel Age Matrix ----------------------------------------

cli::cli_progress_step("Calculating Fuel Age Matrix")

tsf_stk <- terra::rast(fs::dir_ls(tsf_fol))
names(tsf_stk) <- full_fyrs

pix_ha <- resolution^2/10000

# total area
tot_area <- dplyr::as_tibble(terra::freq(aoi_msk)) %>%
  dplyr::mutate(area_ha = count * pix_ha) %>%
  dplyr::pull(area_ha)
# calc fuel age area stats
fuel_df <- tibble::tibble()
for(i in seq_along(full_fyrs)){
  out <- dplyr::as_tibble(terra::freq(tsf_stk[[i]])) %>%
    dplyr::mutate(layer = full_fyrs[i],
                  value = paste0("fa", value),
                  area_ha = count * pix_ha) %>%
    dplyr::rename(fuel_age = value,
                  year = layer) %>%
    dplyr::select(-count)
  fuel_df <- dplyr::bind_rows(fuel_df, out)
}

# fuel age area matrix
fuel_mat <- fuel_df %>%
  tidyr::pivot_wider(names_from = fuel_age, values_from = area_ha) %>%
  replace(is.na(.), 0) %>%
  dplyr::mutate(total = rowSums(pick(where(is.numeric), -year)),
                unknown = tot_area - total,
                id = dplyr::row_number()) %>%
  dplyr::arrange(desc(id)) %>%
  dplyr::select(-id)

fa_mat_nom <- fs::path(tdir, paste0(basename(tdir), "_all_fire_fuel_age_area_matrix_", minyr, "-", maxyr, ".csv"))

readr::write_csv(fuel_mat, file = fa_mat_nom)



# Create Wild Fires Fuel Age Matrix ---------------------------------------

cli::cli_progress_step("Creating Annual Wild Fire Masks")

wf_vecs <- terra::vect(data[["fh_alb"]]) %>%
  tidyterra::filter(fih_fire_t == "WF") %>%
  terra::crop(terra::vect(data$aoi_alb)) %>%
  dplyr::arrange(fin_y)

# unique fire years
u_wfyrs <- wf_vecs %>%
  tidyterra::pull(var = fin_y) %>%
  unique()


wfm_noms <- fs::path(wfm_fol, paste0("wfm", u_wfyrs, ".tif"))

for(i in seq_along(u_wfyrs)){
  wf1 <- wf_vecs %>%
    tidyterra::filter(fin_y == u_wfyrs[i]) %>%
    terra::rasterize(template) %>%
    terra::crop(aoi_msk, mask = TRUE)
  names(wf1) <- u_wfyrs[i]
  terra::writeRaster(wf1, wfm_noms[i])
}


# fuel age burnt by wild fire ---------------------------------------------
cli::cli_progress_step("Calculating Wild Fire Fuel Age Matrix")

wffa_noms <- fs::path(wffa_fol, gsub("wfm", "wffa", dir(wfm_fol)))

# can't calculate for wf if present in first year
if(u_wfyrs[1] == u_fyrs[1]){
  wf_iter <- u_wfyrs[-1]
} else {
  wf_iter <- u_wfyrs
}

# get df of total hectares of burn that year for "wf"
wf_totha_df <- tibble::tibble()

for(yr in wf_iter){
  wfm <- terra::rast(fs::path(wfm_fol, paste0("wfm", yr, ".tif")))
  out_df <- wfm %>%
    terra::freq() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(year = paste0("wf", yr),
                  mask_area_ha = count * pix_ha) %>%
    dplyr::select(year, mask_area_ha)
  wf_totha_df <- dplyr::bind_rows(wf_totha_df, out_df)
  wf_burnt <- terra::rast(fs::path(tsf_fol, paste0("tsf", yr-1, ".tif"))) %>%
    terra::crop(wfm, mask = TRUE)
  terra::writeRaster(wf_burnt, fs::path(wffa_fol, paste0("wffa", yr, ".tif")))
}

wffa_stk <- terra::rast(fs::dir_ls(wffa_fol))
names(wffa_stk) <- wf_iter

# calc fuel age area stats
wffa_df <- tibble::tibble()
for(i in seq_along(wf_iter)){
  out <- dplyr::as_tibble(terra::freq(wffa_stk[[i]])) %>%
    dplyr::mutate(layer = paste0("wf", wf_iter[i]),
                  value = paste0("fa", value),
                  area_ha = count * pix_ha) %>%
    dplyr::rename(fuel_age = value,
                  year = layer) %>%
    dplyr::select(-count)
  wffa_df <- dplyr::bind_rows(wffa_df, out)
}

## re-order fuel ages when out of sync
levs <- paste0("fa", 0:100)


wffa_df$fuel_age <- factor(wffa_df$fuel_age, levels = levs)

# fuel age area matrix
wffa_mat <- wffa_df %>%
  tidyr::pivot_wider(names_from = fuel_age, values_from = area_ha,
                     names_sort = TRUE) %>%
  replace(is.na(.), 0) %>%
  dplyr::mutate(fa_burned = rowSums(pick(where(is.numeric), -year)),
                id = dplyr::row_number()) %>%
  dplyr::arrange(desc(id)) %>%
  dplyr::left_join(wf_totha_df, by = "year") %>%
  dplyr::mutate(fa_unknown = mask_area_ha - fa_burned) %>%
  dplyr::select(-id, -mask_area_ha)


wf_mat_nom <- fs::path(tdir, paste0(basename(tdir), "_wild_fire_fuel_age_area_matrix_", minyr, "-", maxyr, ".csv"))

readr::write_csv(wffa_mat, file = wf_mat_nom)



# Create Other Fires Fuel Age Matrix --------------------------------------

cli::cli_progress_step("Creating Other Fire Masks")

of_vecs <- terra::vect(data[["fh_alb"]]) %>%
  tidyterra::filter(fih_fire_t != "WF") %>%
  terra::crop(terra::vect(data$aoi_alb)) %>%
  dplyr::arrange(fin_y)

# unique fire years
u_ofyrs <- of_vecs %>%
  tidyterra::pull(var = fin_y) %>%
  unique()

ofm_noms <- fs::path(ofm_fol, paste0("ofm", u_ofyrs, ".tif"))

for(i in seq_along(u_ofyrs)){
  of1 <- of_vecs %>%
    tidyterra::filter(fin_y == u_ofyrs[i]) %>%
    terra::rasterize(template) %>%
    terra::crop(aoi_msk, mask = TRUE)
  names(of1) <- u_ofyrs[i]
  terra::writeRaster(of1, ofm_noms[i])
}


# fuel age burnt by other fire --------------------------------------------
cli::cli_progress_step("Calculating Other Fire Fuel Age Matrix")

offa_noms <- fs::path(offa_fol, gsub("ofm", "offa", dir(ofm_fol)))

# can't calculate for of if present in first year
if(u_ofyrs[1] == u_fyrs[1]){
  of_iter <- u_ofyrs[-1]
} else {
  of_iter <- u_ofyrs
}

# get df of total hectares of burn that year for "other"
of_totha_df <- tibble::tibble()

for(yr in of_iter){
  ofm <- terra::rast(fs::path(ofm_fol, paste0("ofm", yr, ".tif")))
  out_df <- ofm %>%
    terra::freq() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(year = paste0("other", yr),
                  mask_area_ha = count * pix_ha) %>%
    dplyr::select(year, mask_area_ha)
  of_totha_df <- dplyr::bind_rows(of_totha_df, out_df)
  of_burnt <- terra::rast(fs::path(tsf_fol, paste0("tsf", yr-1, ".tif"))) %>%
    terra::crop(ofm, mask = TRUE)
  terra::writeRaster(of_burnt, fs::path(offa_fol, paste0("offa", yr, ".tif")))
}

offa_stk <- terra::rast(fs::dir_ls(offa_fol))
names(offa_stk) <- of_iter

# calc fuel age area stats
offa_df <- tibble::tibble()
for(i in seq_along(of_iter)){
  out <- dplyr::as_tibble(terra::freq(offa_stk[[i]])) %>%
    dplyr::mutate(layer = paste0("other", of_iter[i]),
                  value = paste0("fa", value),
                  area_ha = count * pix_ha) %>%
    dplyr::rename(fuel_age = value,
                  year = layer) %>%
    dplyr::select(-count)
  offa_df <- dplyr::bind_rows(offa_df, out)
}

offa_df$fuel_age <- factor(offa_df$fuel_age, levels = levs)

# fuel age area matrix
offa_mat <- offa_df %>%
  tidyr::pivot_wider(names_from = fuel_age, values_from = area_ha,
                     names_sort = TRUE) %>%
  replace(is.na(.), 0) %>%
  dplyr::mutate(fa_burned = rowSums(pick(where(is.numeric), -year)),
                id = dplyr::row_number()) %>%
  dplyr::arrange(desc(id)) %>%
  dplyr::left_join(of_totha_df, by = "year") %>%
  dplyr::mutate(fa_unknown = mask_area_ha - fa_burned) %>%
  dplyr::select(-id, -mask_area_ha)


of_mat_nom <- fs::path(tdir, paste0(basename(tdir), "_other_fire_fuel_age_area_matrix_", minyr, "-", maxyr, ".csv"))

readr::write_csv(offa_mat, file = of_mat_nom)

