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


library(INLA)
inla.setOption(pardiso.license = "/work/msmole/inla_skirt/pardiso.lic")
inla.pardiso.check()

    
      
    ################################################################### 
    ######################## USER INTERACTION #########################
    ###################################################################
    

    cutoff_set <-3
    cutOff<-cutoff_set
    print("Define mesh cutoff distance.")

    n_sample<-0.1

    angle<-c('face','0.0')
    #angle<-c('edge','90')
    
    broj<-'3e08'


    ###################################################################
    ###################### END OF USER INTERACTION ####################
    ###################################################################
    
    
    #######################
    # Loading input files #
    #######################

    
    #sed
    sed<-'halo16_ph3e08_Inc=0.0_Azi=165.34_sed.dat'
    dat<-read.table(sed)
    inputScale <- dat$V1

    #fits
    fit<-'halo16_ph3e08_Inc=0.0_Azi=165.34_total.fits'
    input_fits <- readFITS(fit)

    dims =  dim(input_fits$imDat)
    features <- 1:dims[3]#
    
    input_pre = (input_fits$imDat[1201:1800,1201:1800, features])
    dims =  dim(input_pre)
  

    ##################################################
    # some variables that will be handy further down #
    ##################################################
    feature_size <- length(features)                 #
    spatial_size <- dims[1] * dims[2]                #
    ##################################################
    
    ########################################################
    # Turning map (without unwanted wavelengths) into line #
    ########################################################

    
    input_post = array_reshape(input_pre,
                               dim = c(spatial_size, feature_size))

    
    # sampling spaxels
    N <- n_sample * spatial_size
    sampled_ind<-c(1:spatial_size)
    input_Net_indeces <- sample(sampled_ind, N)#
    
    
    #######################################
    # Setting up data to input to INLA #
    #######################################
    input_Net <- array(NaN, dim = c(spatial_size,feature_size)) 
    
    
    input_Net[input_Net_indeces, ] <- input_post[input_Net_indeces, ]
    
    
    input_Net<-  array_reshape(input_Net,
                                          dim = c(dims[2], dims[1], feature_size))
    
    #######################################
 
   
    ##################################
    # Variables to be used with INLA #
    ##################################
    nSize <- max(dims[1], dims[2], na.rm = TRUE)
    X <- rep(1:dims[1], times = dims[2])
    Y <- rep(1:dims[2], each = dims[1])
    ##################################
    
    ####################################################
    # Arrays to hold data for and from INLA procedures #
    ####################################################
    time_holder <- array(NaN, dim = feature_size)       #
    INPUT_PRE_arr <- array(NaN, dim = c(dims[2], dims[1], feature_size))
    INPUT_POST_arr <- array(NaN, dim = c(dims[2], dims[1], feature_size))
    MEAN_arr <- array(NaN, dim = c(dims[2], dims[1], feature_size))
    SD_arr <- array(NaN, dim =c(dims[2], dims[1], feature_size))
    ####################################################
    
    ##########################
    ##########################
    ## MAIN INLA WORK CYCLE ##
    ##########################
    ##########################
    for(s in 1:feature_size){##
    #for(s in 31:31){##
        
      ##  
      # DATA TRANSFORM      ##
      tempZ <- input_Net[, , s]
      INPUT_PRE_arr[,,s] <-   input_pre[,,s]
      Z <- as.vector(tempZ)      
      
      
      #normalization
      
      normZ <- max(Z,na.rm=T)
      #fits_df <- data.frame(x=X, y=Y, z=Z/normZ, z_err=Z_ERR/normZ)
      Z<-Z/normZ
      
      
      
      #DATA TRANSFORM      ##
      Z[which(is.infinite(Z), arr.ind = TRUE)] <- NaN
      ##
      maxZ <- max(Z, na.rm = TRUE)
      minZ <- min(Z, na.rm = TRUE)
      ##
      R <- maxZ - minZ      ##
      Z <- Z - minZ         ##
      ##
      minZ_2 <- min(Z[which(Z > 0, arr.ind = TRUE)], na.rm = TRUE)
      minZ_3 <- minZ_2 / R^2##
      Z <- Z + minZ_3       ##
      
      
    
      
      fits_df <- data.frame(x=X, y=Y, z=Z)
      
      INPUT_POST_arr[,,s] <- fits_df$z
      
      
      rm(tempZ)             ##
      ##
      # Setup INLA Params   ##
      zoom = 1              ##
      p_range = c(2, 0.2)   ##
      p_sigma = c(2, 0.2)   ##



      new_exmag <- na.omit(fits_df)
      exCoord <- cbind(new_exmag$x, new_exmag$y)
      ##
      # INLA Procedures     ##
      start_time <- Sys.time()
      
      mesh <- inla.mesh.2d(exCoord, max.n = nSize, 
                           cutoff = cutOff)
      end_time <- Sys.time()##
      mesh_time <- end_time - start_time
      print( paste0('INLA mesh.2d() took: ', mesh_time, ' seconds') )
      
    
      
      ##
      start_time <- Sys.time()
      A <- inla.spde.make.A(mesh = mesh, loc = exCoord)
      end_time <- Sys.time()##
      make_time <- end_time - start_time
      print( paste0('INLA spde.make.A() took: ', make_time, ' seconds'))
      ##  
      start_time <- Sys.time()
      projection <- inla.mesh.projector(mesh, xlim = c(1, dims[1]), 
                                        ylim = c(1, dims[2]), 
                                        dim = zoom*c(dims[1] + 1,
                                                     dims[2] + 1))
      end_time <- Sys.time()##
      proj_time <- end_time - start_time
      print( paste0('INLA mesh.projector() took: ', proj_time, 
                    ' seconds'))
      ##
      start_time <- Sys.time()
      spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, 
                                  prior.range = p_range, 
                                  prior.sigma = p_sigma)
      end_time <- Sys.time()##
      pcmatern_time <- end_time - start_time
      print( paste0('INLA spde2.pcmatern() took: ', pcmatern_time, 
                    ' seconds') )
      ##  
      start_time <- Sys.time()
      stk <- inla.stack(    ##
        data = list(par = new_exmag$z), A = list(A, 1),
        effects = list(i = 1:spde$n.spde, 
                       intercept = rep(1, length(new_exmag$x))), 
        tag = 'est')        ## 
      end_time <- Sys.time()##
      stack_time <- end_time - start_time
      print( paste0('INLA stack() took: ', stack_time, ' seconds') )
      ##
      start_time <- Sys.time()
      possibleError <- tryCatch(
        res <- inla(new_exmag$z ~ -1 + intercept + f(i, model = spde), 
                    data = inla.stack.data(stk), 
                    control.predictor = list(A = inla.stack.A(stk), 
                                             compute = TRUE), 
                    #scale = (1 / (new_exmag$z_err)),
                    verbose = F),
        error = function(e)
          e)                ##
      
     
      ##
      if (inherits(possibleError, "error")) {
        print(paste('INLA reconstruction error: '))
        next                ##
      }                     ##
      end_time <- Sys.time()##
      inla_time <- end_time - start_time
      print( paste0('INLA inla() took: ', inla_time, ' seconds') )
      ##
      start_time <- Sys.time()
      ##
      # OUTPUT MEAN VALUE   ##
      output_mean <- inla.mesh.project(
        inla.mesh.projector(mesh, xlim = c(1, dims[1]), 
                            ylim = c(1, dims[2]), 
                            dim = zoom * c(dims[1], dims[2])),
        res$summary.random$i$mean) +
        t(matrix(as.numeric(res$summary.fixed$mean[1]), 
                 nrow = zoom * (dims[2]), ncol = zoom * (dims[1])))
      ##
      # OUTPUT SD           ##
      output_sd <- inla.mesh.project(
        inla.mesh.projector(mesh, xlim = c(1, dims[1]), 
                            ylim = c(1, dims[2]), 
                            dim = c(dims[1], dims[2])),
        res$summary.random$i$sd)
      ##
      end_time <- Sys.time()##
      project_time <- end_time - start_time
      print( paste0('INLA projecting solution took: ', project_time, 
                    ' seconds') )
      total_inla_time <- project_time + inla_time + stack_time + pcmatern_time + proj_time + make_time + mesh_time
      print( paste0('INLA in total took: ', total_inla_time, ' seconds') )
    ##########################
    ##########################
      #nekad je za neki pixel rekonstrukcija negativna, menja u 0
      #output_mean[which(output_mean <= 0)] <- 0 #because inla is going beyond the imposed inferior limit
      
      Z_mean <- output_mean - minZ_3  + minZ
      Z_sd <- output_sd - minZ_3  + minZ
      
      Z_mean<-Z_mean*normZ
      Z_sd<-Z_sd*normZ
      
      ##################################################################
      ################## END OF DATA DISPLAY TRANSFORM #################
      ##################################################################
      time_holder[s] <- total_inla_time
      
      print(paste0("End of iteration ",s," out of ", feature_size))
      
      
      MEAN_arr[,,s] <- Z_mean
      SD_arr[,,s] <- Z_sd
      #RES_arr[,,s] <- residuals_ratio
      #INPUT_RES_arr[,,s] <- input_residuals_ratio
      
      # RES_MOD_arr[,,s] <- residuals_mod_ratio 
      # INPUT_RES_MOD_arr[,,s] <- input_residuals_mod_ratio
      # RES_MOD_arr[1,1,s] <- rMr_counter
      # INPUT_RES_MOD_arr[1,1,s] <- IrMr_counter
      
      }
    ###################################################################
    ###################### OUTPUT FILES CREATION ######################
    ################################################################### 
    
    #writeFITSim(INPUT_PRE_arr, file = paste0('INPUT_PRE_T.fits'))
    #writeFITSim(INPUT_POST_arr, file = paste0('INPUT_POST_T.fits'))
    writeFITSim(MEAN_arr, file = paste0('pure_inla_Mean.fits'))
    #writeFITSim(SD_arr, file = paste0('INLA_sd.fits'))
    
    ###################################################################
    ################### END OF OUTPUT FILES CREATION ##################
    ###################################################################
    
    time_holder=c(time_holder/60,sum(time_holder, na.rm =TRUE)/60)
    
    write.csv(time_holder,"time_pure_inla.csv", row.names = FALSE)
