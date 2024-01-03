
## -- Módulo 24: Dados espaciais e mapas

setwd("D:/Daniel/courses/curso_R/modulo24-dados_espaciais_e_mapas/")
  # Descarga de dados desde: https://srtm.csi.cgiar.org/


# 6: Dados raster ---------------------------------------------------------

library(rgdal)
library(raster)

darien <- raster::raster("input/raster/srtm_21_11.tif")
darien
class(darien)

plot(darien)

methods(class = 'RasterLayer') # Revisar os métodos aplicáveis para RasterLayer



# 7: Dados vetoriais ------------------------------------------------------

library(sf)

colo <- sf::st_read("input/vector/Limite Departamental.shp")
class(colo)  ## sf : spatial features
plot(colo)

col_geom <- colo$geometry #  --> extrair a geometria sem os atributos
class(col_geom) ## sfc : spatial features collection
plot(col_geom)


guajira <- col_geom[[12]]
class(guajira) ## sfg : somente um polígono
plot(guajira)




# 8: Sistema de coordenadas na prática ------------------------------------


  # Sistema de coordenadas
    # Formas de representação:
    # - WKT2: Well Known Text
    # - Proj4, Proj4string
    # - EPSG: European Petroleum Survey Group


library(sf)
library(raster)


sergipe <- st_read("input/vector/SE_Microrregioes_2022.shp")
sergipe

class(sergipe)
plot(sergipe)

    # Consultar o Sistema de Referencia de Coordenadas em dados vetoriais
crs1 <- st_crs(sergipe)
crs2 <- crs(sergipe)

crs1$wkt
crs1$epsg
crs1$proj4string

crs2  # Mesmo resultado de crs1$proj4string



  # Consultar o Sistema de Referencia de Coordenadas em dados raster

ndvi <- raster("input/raster/MOD_NDVI_M_2023-11-01_rgb_720x360.FLOAT.TIFF")
ndvi
dim(ndvi)

plot(ndvi) # -> Mapa não e preciso ja que tem dados de 99999 para os oceanos

ndvi[,1] # -> revisar os dados de 99999 para a coluna 1


  # Mudar os valores de 99999 para NA e plotar novamente o mapa global de ndvi

ndvi[ndvi == 99999.000] <- NA
ndvi[,1]

plot(ndvi)


  # Revisar o CRS do mapa raster de ndvi
raster::crs(ndvi)




# 9: Atribuição de Sistema de Coordenadas ---------------------------------


  # Como atribuir CRS aos dados que não posuim

    # Vector

lagos <- st_read("input/vector/wtrbodies_utm.shp")
lagos
plot(lagos)
plot(lagos$geometry)

st_crs(lagos) # Coordinate Reference System: NA
crs(lagos)

st_crs(lagos) <- "+proj=utm +zone=11 +datum=NAD27 +units=m +no_defs"
                # Extrair o CRS desde https://epsg.io/
st_crs(lagos)
st_crs(lagos)$proj4string
crs(lagos)


    # Raster

amazonas <- raster("input/raster/srtm_23_13.tif")
amazonas

plot(amazonas)

crs(amazonas) <- "+proj=longlat +ellps=GRS80 +no_defs" # Aviso que não tem datum 




# 10: Transformação de Sistema de Coordenadas -----------------------------


library(sf)
library(raster)


setwd("D:/Daniel/courses/curso_R/modulo24-dados_espaciais_e_mapas/")

ndvi <- raster("input/raster/MOD_NDVI_M_2023-11-01_rgb_720x360.FLOAT.TIFF")
ndvi # crs = +proj=longlat +datum=WGS84 +no_defs --> Coordenadas geodésicas

plot(ndvi)

ndvi[ndvi == 99999.0] <- NA

plot(ndvi)

    # Transformação do CRS para coordenadas planas

# Obter notações de projeções no site proj.org -> Coordinate operations -> 
# Projections

ndvi.proj <- projectRaster(ndvi, crs = '+proj=ortho +long_0=10')

plot(ndvi.proj)                


    # Transformação do CRS em dados vector

countries <- st_read("input/vector/ne_110m_admin_0_countries_lakes.shp")


plot(countries$geometry) # Dados espaciais em coordenadas geodesicas

countries_proj <- st_transform(countries, crs = "+proj=moll")

plot(countries_proj$geometry)


  # Plotar os dois mapas juntos, uno acima do outro

plot(ndvi.proj)
plot(countries_proj, add = TRUE)


  # Gerar outra projeção

ortho.proj <- st_transform(countries, crs = "+proj=ortho +lon_0=-85")
plot(ortho.proj$geometry)


ndvi.ortho <- projectRaster(ndvi, crs = '+proj=ortho +lon_0=-85')

plot(ndvi.ortho)
plot(ortho.proj$geometry, add= TRUE)




# 11: Como acessar dados do IBGE (pacote geobr) ---------------------------

library(geobr)
library(magrittr)

br <- read_country(year = 2015)
plot(br)

estados <- read_state(code_state = 'all')

estados %>% 
  filter(name_region=='Norte')


plot(estados$geom)


biomas <- read_biomes()
plot(biomas$geom)


amazonas <- read_amazon()
plot(amazonas$geom)


# 12: Como ler dados de GPS -----------------------------------------------

library(sf)
library(ggplot2)

setwd("D:/Daniel/courses/curso_R/modulo24-dados_espaciais_e_mapas/")


  # Listar as camadas do arquivo gpx
st_layers("input/vector/passeio_bike.gpx")


  # Leitura da camada "tracks" do arquivo gpx
tracks <- st_read("input/vector/passeio_bike.gpx", layer = "tracks")


  # Plotar o objeto tracks
plot(tracks$geometry)

ggplot(tracks) +
  geom_sf(col = 'red')



# 13: Como ler dados do google earth --------------------------------------

library(sf)
library(ggplot2)
library(magrittr)

setwd("D:/Daniel/courses/curso_R/modulo24-dados_espaciais_e_mapas/")


  # Carregar o arquivo kml

tota <- st_read("input/vector/tota.kml")
plot(tota$geometry)

ggplot(tota) +
  geom_sf(fill = 'white', alpha = 0.2 ) +
  theme_minimal()


  # Carregar o arquivo kmz

tota2 <- unzip("input/vector/tota.kmz") %>% 
             st_read()

plot(tota2$geometry)  

ggplot(tota2) +
  geom_sf(fill = 'dodgerblue1', alpha = 0.2 ) +
  theme_void()


  
  # Working with google kml sample data

download.file(url = "https://developers.google.com/kml/documentation/KML_Samples.kml",
              destfile = "input/vector/kml_google_sample.kml")

st_layers("input/vector/kml_google_sample.kml")

google_campus <- st_read("input/vector/kml_google_sample.kml", 
                         layer = "Google Campus")


ggplot(google_campus) +
  geom_sf()




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


