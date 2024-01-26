
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










