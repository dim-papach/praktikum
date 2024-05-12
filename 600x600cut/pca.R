  rm(list = ls())
  options(digits = 16)
  require(compiler)
  require(beepr)
  require(rgl)
  require(lattice)
  require(stringr)
  require(INLA)
  require(sp)
  require(fields)
  require(spam)
  require(ggplot2)
  require(viridis)
  require(yaml)
  require(colorspace)
  require(spatstat)
  require(data.table)
  require(gtools)
  require(FITSio)
  require(IDPmisc)
  library(reticulate)
  corTesteA <- colorspace::diverge_hsv
  set.seed(134775813)
  #weighted.median com w = inla
  
  require(jjb)
  require(classInt)
  library(reticulate)
  library ( pls )
  library(ggfortify)
  library( lattice )
  library(gridExtra)
  library(grid)
  
  
  
  ################################################################### 
  ######################## USER INTERACTION #########################
  ###################################################################
  
  
  angle<-c('face','0.0')
  #angle<-c('edge','90')
  
  broj<-'3e08'

  fit<-'halo16_ph3e08_Inc=0.0_Azi=165.34_total.fits'

  #original input
  input_fits <- readFITS(fit)
  
  dims =  dim(input_fits$imDat)
  features <- 1:dims[3]#
  
  original_pre <- (input_fits$imDat[1201:1800,1201:1800, features])

  rm(input_fits)
  # 
  dim_cube<-dim(original_pre)
  input_dataset_size <- dim_cube[1]*dim_cube[2]
  # 
  original_post <- array_reshape(original_pre, dim = c(input_dataset_size,dim_cube[3]))
  
     
  #pca
      
      
        start_time <- Sys.time()
        pca<-prcomp(original_post, center = TRUE,scale. = TRUE)
        end_time <- Sys.time()
        time<-end_time - start_time
  
  
        saveRDS(pca, file = "pca.rds")
        
  