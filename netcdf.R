
## -- Módulo 24: Dados espaciais e mapas



# 14: Como ler arquivos NetCDF com o pacote raster ------------------------

library(magrittr)
library(raster)
library(ncdf4)  # functions for reading and writing netCDF files

c("ncmeta","fields") %in% installed.packages() 
# Check if the ncmeta and fields packages are installed

library(ncmeta) # Provides a set of tools to obtain metadata from NetCDF files
library(fields) # Provides functions for spatial data analysis, geostatistics...
library(ggplot2)

setwd("D:/Daniel/courses/curso_R/modulo24-dados_espaciais_e_mapas/")

ar_temp_nc <- "input/raster/temp.nc"


# Obter informações do arquivo nc

raster("input/raster/temp.nc")  # Fazer leitura do arquivo

nc_vars(ar_temp_nc)   # Generate a table of all variables

nc_dims(ar_temp_nc)   # Information about the dimensions in a NetCDF source


# Leitura e manipulação do arquivo
ar_temp <- raster(ar_temp_nc)

plot(ar_temp) # A temperatura do ar esta em graus K

plot(ar_temp-273.15) # Agora a temperatura e plotada em graus celsius

plot(ar_temp, col= c("blue","green","yellow","orange","red"))
# Modificação das cores

plot(ar_temp, col = fields::tim.colors(n = 150))




###### air.mon.mean file  #######

air.mon.nc <- "input/raster/air.mon.mean.nc"

raster(air.mon.nc)
# Segondo as informações do anterior comando, o arquivo não tem CRS

air.mon <- raster(air.mon.nc)
raster::crs(air.mon)              # CRS arguments: NA 


# Atribuir um CRS para o arquivo.
#   Para isso obter a definição PROJ.4 desde o site epsg.io

raster::crs(air.mon) <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"
raster::crs(air.mon)


## Plotar o arquivo

plot(air.mon)  # Coordenadas x (longitude) de 0 a 360 (frequently used in global climate models)

air.mon.180 <- raster::rotate(air.mon)

plot(air.mon.180, col = fields::tim.colors(n = 150))

ncmeta::nc_vars(air.mon.nc)
ncmeta::nc_dims(air.mon.nc)


## Open the NetCDF file

library(ncdf4)
nc_air_mon <- ncdf4::nc_open(air.mon.nc)

nc_air_mon

## Read dimension data

nc_dims(air.mon.nc)  # Get a summary of the dimensions. The input is the file route

lat <- ncvar_get(nc_air_mon,"lat")
lat

lon <- ncvar_get(nc_air_mon, "lon")
lon

time <- ncvar_get(nc_air_mon, "time")
time

## Read variable data

nc_vars(air.mon.nc) # Get a summary of the vars. The input is the file route

air <- ncvar_get(nc_air_mon, "air")


# Get the dimensions information

lat_dim <- nc_dim(nc_air_mon, "lat")
lon_dim <- nc_dim(nc_air_mon, "lon")
time_dim <- nc_dim(nc_air_mon, "time")

nc_dims(nc_air_mon)  # This time the input is not the file route, but the opened nc file


## Create raster stack

raster_stack <- raster::brick(air, 
                              xmn = min(lon), xmx = max(lon),
                              ymn = min(lat), ymx = max(lat),
                              crs = sp::CRS("+proj=longlat +datum=WGS84"))

raster_stack

## Plot a single layer

plot(raster_stack[[1]])


## Plot time series

time_series <- raster::extract(raster_stack, cbind(lon, lat)) # results in a matrix

class(time_series)
dim(time_series)

time_series <- t(time_series)
class(time_series)

time_series_df <- data.frame(time = as.POSIXct(time, origin = "1970-01-01"),
                             matrix(time_series, ncol = ncol(time_series),
                                    byrow = TRUE))

## Plot time series

ggplot(time_series_df, aes(x = time, y = X2)) +
  geom_line()


## Close the NetCDF file

nc_close(nc_air_mon)


# ersst file

ersst.nc <- "input/raster/ersst.nc"

nc_vars(ersst.nc)
nc_dims(ersst.nc)

ersst_ssta <- raster(ersst.nc, varname = "ssta")
ersst_sst <- raster(ersst.nc, varname = "sst")

plot(ersst_ssta)
plot(ersst_sst)

# Swap 0 to 360 scale to -180 to 180 one
ssta_180 <- rotate(ersst_ssta)

plot(ssta_180)




# 15: Como ler arquivos NetCDF com o pacote raster (2) --------------------

setwd("D:/Daniel/courses/curso_R/modulo24-dados_espaciais_e_mapas/")

library(raster)
library(ncmeta)


era.nc <- "input/raster/era.nc"

nc_vars(era.nc)
nc_dims(era.nc)

era_brick <- raster::brick(era.nc)

class(era_brick)   # RasterBrick: a multi-layer raster object.







