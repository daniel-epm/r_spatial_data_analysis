
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




