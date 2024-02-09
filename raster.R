
## -- Módulo 24: Dados espaciais e mapas



# 18: Operações com raster - juntar dados ---------------------------------


setwd("D:/Daniel/courses/curso_R/modulo24-dados_espaciais_e_mapas/")

library(raster)
library(ggplot2)


  # Juntar arquivos tif do estado de Sergipe


tifs <- list.files("input/raster/srtm_sergipe/", pattern = "*.tif")


raster_list <- list()

# Use a for loop to read each file and store in the list
for (i in seq_along(tifs)) {
  # Construct the variable name (r1, r2, ..., r6)
  variable_name <- paste0("r", i)
  
  # Read the TIFF file using the raster function
  raster_obj <- raster(file.path("input/raster/srtm_sergipe", tifs[i]))
  
  # Assign the raster object to the variable
  assign(variable_name, raster_obj)
  
  # Save the raster object in the list
  raster_list[[i]] <- raster_obj
}


  # In order to merge raster files, theu must have the same origin:
  ## In the following for loop we make sure of the same origin for every file

for (i in raster_list) {
  print(raster::origin(i))
}


    # Plotting each raster file
plot(r1)
plot(r2)
plot(r3)
plot(r4)
plot(r5)
plot(r6)


  # Merging the raster files
sergipe <- raster::merge(r1,r2, r3, r4, r5, r6)  # Pode-se usar a função raster::mosaic() também

  # Plotting the resultant merged file
plot(sergipe)


  # Plotting with ggplot and adding a mask for the sergipe state
brasil <- geobr::read_state(code_state = 'all')

brasil[brasil$name_state == 'Sergipe', ]

st_sergipe <- geobr::read_state(code_state = 28)


ggplot() +
  geom_raster(data = as.data.frame(sergipe$layer, xy = TRUE),
              aes(x = x, y = y, fill = layer)) +
  geom_sf(data = st_sergipe, fill = NA, lwd = 0.7, color= 'white')




abc <- c('a','b','c')


nums <- 1:3

tidyr::expand_grid(abc, nums)




# 19: Operações com raster (juntar todos os dados) ------------------------

library(raster)
library(fields)

setwd("D:/Daniel/courses/curso_R/modulo24-dados_espaciais_e_mapas/")


tifs_sergipe <- list.files(path = "input/raster/srtm_sergipe/" , 
                           pattern = "*.tif", recursive = T, full.names = T)

rasters_sergipe <- lapply(tifs_sergipe, raster)
  # lapply: applies a function to each element of a list X, returning a list
  # with same length as initial X


  # Join the raster files using merge function
merged_rasters <- do.call(what = merge, args = rasters_sergipe)

plot(merged_rasters)


  # Join the raster files using mosaic function

    # Specify the fun and tolerance arguments required in the mosaic function
rasters_sergipe$fun <- mean  # Function to use on the overlayed points
rasters_sergipe$fun(c(1,2,3,4,5)) # Using the incorporated mean function

      # Tolerance is required whenever there's difference in the origin values
      # from the raster files. We can check those values as follows:

for (i in rasters_sergipe) {
  print(origin(i))
}
      # Here the values are the same, but in case they are different we set the 
      # tolerance argument in the mosaic function. Here we include it in the list
      # to be passed to the mosaic function through the do.call function.
rasters_sergipe$tolerance <- 0.2


rasters_mosaic <- do.call(what = mosaic, args = rasters_sergipe)

plot(rasters_mosaic, col = tim.colors(n = 50))



# 20: Operações com raster (recortar e mascarar dados) --------------------

library(raster)
library(geobr)
library(sf)
library(magrittr)

setwd("D:/Daniel/courses/curso_R/modulo24-dados_espaciais_e_mapas/")

    # Extrair os dados específicos para a extensão territorial de sergipe
# Using the mosaic previously generated: rasters_mosaic



  # Filtrar dados de Sergipe
brasil <- geobr::read_state(code_state = "all")
brasil[brasil$name_state == "Sergipe",]

  # Criar um objeto com a extensão territorial do estado de Sergipe
sgp <- geobr::read_state(code_state = 28)

plot(sgp$geom)


  # Revisar o CRS das capas
sf::st_crs(sgp)$proj4string  # +proj=longlat +ellps=GRS80
raster::crs(rasters_mosaic)  # +proj=longlat +datum=WGS84

    # O CRS não é o mesmo nas duas capas, é preciso mudar o CRS de uma delas
sgp <- sgp %>% 
  sf::st_transform(crs = "+proj=longlat +datum=WGS84")

st_crs(sgp)$proj4string  # [1] "+proj=longlat +datum=WGS84 +no_defs"


  # Após de ter o mesmo CRS é possivel juntar o mapa de relevo com os limites 
  # territoriais do estado de Sergipe
plot(rasters_mosaic)
plot(sgp$geom, add = T)


  # Crop
m1_crop <- raster::crop(rasters_mosaic, sgp)
plot(m1_crop)
plot(sgp$geom, add= T)


  # Mask
m2_mask <- raster::mask(rasters_mosaic, sgp) # Argument inverse = T for ext area
plot(m2_mask)
plot(sgp$geom, add=T)

