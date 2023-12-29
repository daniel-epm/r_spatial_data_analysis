
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



# Sistema de coordenadas na prática ---------------------------------------

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



# Atribuição de Sistema de Coordenadas ------------------------------------

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

