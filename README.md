# PB-Effectiveness

This repository contains code that constructs fuel age matrices from the 
Department of Biodiversity, Conservation and Attractions Fire History dataset (DBCA_060). 
DBCA_060 is freely available from [data WA
link](https://catalogue.data.wa.gov.au/dataset/dbca-fire-history).

### Requirements

- All code is written in R. You will need R (version > 4.0) to run the code.

- You will need to have downloaded the DBCA_060 locally to your computer. **NOTE** 
The data is periodically updated so it always pays to refresh your copy of DBCA_060.

- A shapefile of your area of interest (AOI).

### The AOI Shapefile

The AOI shapefile used in the development of this code was a DBCA region, further 
subset to DBCA tenure. If using something similar ensure that if it has single part 
geometry it has been converted to multipart and dissolved. 

To further focus on specific vegetation types, the above shapefile should then been 
further subset to reflect this.

The shapefile can have any CRS as it will be transformed to EPSG:9473 
(GDA2020 Australian Albers) to allow for calculations in meters. Note output units 
of area are hectares.

### Notes on Processing

The R script in this repository, `single_aoi_fa.R`, has been written to take one 
aoi and produce three csv fuel age matrices (as csv files), all fires fuel age, 
wild fires fuel age and other fires fuel age.

The script can be readily adapted to loop over a number of aoi's if required but 
bear in mind that large rasters can be very memory intensive and slow to process.

The DBCA_060 data set is very large and unwieldy and covers the whole State. Using 
the spatial extents of your aoi, the data set is subset as it is read into memory. 
As it is reading in data based on extent, where possible avoid constructing aoi's 
that might contain very large areas of no interest between polygons as these areas 
outside of the aoi polygons will still be processed (but not included in the results).

If the above behaviour can't be avoided (perhaps the requirement is to process and 
report on the entire south west for example) there is a work around. The script 
`00_preprocess_firescar_vectors_to_raster_base_data.R` found in the [Remote Sensing 
and Spatial Analysis Program's github](https://github.com/RSPaW/fire_metrics/tree/main) 
can be modified to produce the three essential annual products, `byYYYY`, `yobYYYY` and
`tsfYYYY`. The pre-amble for the script contains 3 functions that will assess the user 
aoi size and if necessary, spilt it into chunks that can be looped over to process 
large areas. The `00_preprocess_firescar_vectors_to_raster_base_data.R` isn't completely 
interchangeable with this workflow but it should contain enough to get you started 
on modifications.

### The Basic Process

- The spatial extent of your AOI will be used to extract historical fire vectors 
from DBCA_060.

- The vector data is rasterised and various annual raster products are made. These 
annual (YYYY) products include:
  * `byYYYY`    - raster fire scars that occured in that year with cell values of that year.
  * `yobYYYY`   - raster fire scars that occured in that year and all preceding years, with 
  rasters overlying each other with most recent year on top.
  * `tsfYYYY`   - current year minus yobYYYY for current year producing raster with cell 
  values of time since fire (or year since last burn).
  * `wfmYYYY`   - wild fire raster mask for that year. Wild fires are subset from 
  DBCA_060 by vector attribute.
  * `ofmYYYY`   - other fire mask for that year. All fires not attributed as wild fire 
  from DBCA_060.
  * `wffaYYYY`  - tsfYYYY masked by wfmYYYY that produces the wild fire fuel age raster 
  for that year.
  * `offaYYYY`  - tsfYYYY masked by ofmYYYY that produces the other fire fuel age raster 
  for that year.

- Area stats are then calculated by multiplying the pixel count of unique fuel ages 
by the pixel resolution and exported as csv.

### The YOB to TSF to Fuel Age process in more detail

If you were to look at all of the fire scar data from the 1991 to the 1994 
in your AOI it might look like the below. To arrive at this we have stacked the 
burn year data (byYYYY) on top of each other. This is the year of burn (yobYYYY) 
data.

![**Year of burn for 1994**](./pics/yob_stack.png)  


If you then subtracted the year of burn data from 1994 you would get the below, which 
is the time since fire (tsfYYYY) or also known as year since last burn.

![**Time since fire for 1994**](./pics/tsf_stack.png)


Next you can construct a fire mask (all, wild fire or other) for 1995 and use it to 
crop the time since fire for 1994. This will result in a classified raster that 
has the fuel ages that the wild fire in 1995 burnt across. Something like below.

![**1995 wild fire fuel age**](./pics/fuel_age.png)

The last process is to extract the count of pixels for each of the remaining fuel 
ages and convert them to area.

![**1995 wild fire fuel age area stats**](./pics/stats.png)