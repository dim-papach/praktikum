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
    library(NMF)
    
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
      
      
      nmf_comps<-10
      nmf_dir <-paste(nmf_comps)
      
      
      mainDir <- paste0('/work/msmole/auriga/')

      
      fit<-'halo16_ph3e08_Inc=0.0_Azi=165.34_total.fits'
      
      
      #original input
      input_fits <- readFITS(fit)
      
      dims =  dim(input_fits$imDat)
      features <- 1:dims[3]#
      
      original_pre <- (input_fits$imDat[1201:1800,1201:1800, features])
      #dims =  dim(original_pre)
      
      rm(input_fits)
      # 
      dim_cube<-dim(original_pre)
      input_dataset_size <- dim_cube[1]*dim_cube[2]
      # 
      original_post <- array_reshape(original_pre, dim = c(input_dataset_size,dim_cube[3]))
      
      # ###########################3
     
      #normalization by lambda
      
      norm_lambda<- array(data=NA, dim= dim(original_post))
      for (i in 1:dim_cube[3]){
        norm_lambda[,i] <- (original_post[,i]-min(original_post[,i]))/(max(original_post[,i])-min(original_post[,i]))
      }
      
      set.seed(24527)
      
      
      start_time <- Sys.time()
      nmf_sed<-nmf(norm_lambda, nmf_comps, 'brunet')
      end_time <- Sys.time()
      time<-(end_time - start_time)
      print(time) 
      
      write.csv(time,"time_nmf.csv", row.names = FALSE)
      
      saveRDS(nmf_sed, file = paste0('nmf_',nmf_comps,'.rds'))
      
     