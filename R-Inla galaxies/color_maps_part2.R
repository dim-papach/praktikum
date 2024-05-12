rm(list = ls())
options(digits = 16)
require(ggplot2)
require(gtools)
require(FITSio)
require(keras)


require(jjb)
require(classInt)
library(reticulate)
library ( pls )
library(ggfortify)
library( lattice )
library(gridExtra)
library(grid)

library(cowplot)
################################################################### 
######################## USER INTERACTION #########################
###################################################################
###################################################################


 n_sample<-0.25

angle<-c('face','0.0')
#angle<-c('edge','90')

num.comp<-10


fit<-'halo16_ph3e08_Inc=0.0_Azi=165.34_total.fits'
fit_ref<-'Au16_Inc=0.0_Azi=165.34_total.fits'


#original fits file , input_fits input_fits <- readFITS(input_path_fit)

input_fits <- readFITS(fit)

dims =  dim(input_fits$imDat)
features <- 1:dims[3]#

original_pre <- input_fits$imDat[1201:1800,1201:1800, features]
#dims =  dim(original_pre)

original_pre[which(original_pre < 0)] <- 0

dim_cube<-dim(original_pre)
input_dataset_size <- dim_cube[1]*dim_cube[2]

original_post <- array_reshape(original_pre, dim = c(input_dataset_size,dim_cube[3]))


rm(input_fits)


###########################3
#original 1e8 fits file 

input_fits <- readFITS(fit_ref)

dims <-  dim(input_fits$imDat)    

original_ref_pre <- input_fits$imDat[1201:1800,1201:1800, features]


rm(input_fits)
original_ref_post <- array_reshape(original_ref_pre, dim = c(input_dataset_size,dim_cube[3]))

original_ref_post[which(original_ref_post < 0)] <- 0



#sed
sed<-'halo16_ph3e08_Inc=0.0_Azi=165.34_sed.dat'
dat<-read.table(sed)
inputScale <- dat$V1



# pure INLA

latent_size<-dim_cube[3]


input_fits <- readFITS( 'pure_inla_Mean.fits')

pure_inla_pre <- input_fits$imDat[1:dim_cube[1],1:dim_cube[2],1:latent_size]

pure_inla_post <- array_reshape(pure_inla_pre, dim = c(input_dataset_size,dim_cube[3]))

#true zeroes zeroes at all lambda for reference 

true_zero_ref<-apply(original_ref_post, 1, sum )

true_zero_pixels_ref <- which(true_zero_ref == 0)

pure_inla_post[true_zero_pixels_ref,] <- 0
pure_inla_pre <- array_reshape(pure_inla_post, dim = dim_cube)


#PCA+INLA

latent_size<-num.comp
input_fits <- readFITS('pca+inla_Mean.fits')


pca_pre <- input_fits$imDat[1:dim_cube[1],1:dim_cube[2],1:dim_cube[3]]


pca_post <- array_reshape(pca_pre, dim = c(input_dataset_size,dim_cube[3]))

rm(input_fits)

pca_post[true_zero_pixels_ref,] <- 0


#  NMF+INLA

input_fits <- readFITS('nmf+inla_Mean.fits')



nmf_pre <- input_fits$imDat[1:dim_cube[1],1:dim_cube[2],1:dim_cube[3]]


nmf_post <- array_reshape(nmf_pre, dim = c(input_dataset_size,dim_cube[3]))

rm(input_fits)



nmf_post[true_zero_pixels_ref,] <- 0


  pure_inla_post[which(pure_inla_post < 0)] <- 0    


  pca_post[which(pca_post < 0)] <- 0    
  pca_post <- array_reshape(pca_post, dim = dim_cube)

  nmf_post[which(nmf_post < 0)] <- 0    
  nmf_post <- array_reshape(nmf_post, dim = dim_cube)



######################################

original_ref_pre_r <- original_ref_pre                       #
original_ref_pre_r[which(original_ref_pre_r == 0)] <- 1



########## color maps

# sampled original input
set.seed(134775813)
spatial_size<-dim_cube[1]*dim_cube[2]
N <- n_sample * spatial_size
sampled_ind<-c(1:spatial_size)
input_Net_indeces <- sample(sampled_ind, N)#






#######################################
input_Net <- array(NaN, dim = c(spatial_size,dim_cube[3])) 


input_Net[input_Net_indeces, ] <- original_post[input_Net_indeces, ]


input_Net<-  array_reshape(input_Net,
                           dim = dim_cube)

rm(original_post, original_ref_post)

lambda<-c(1:50) #index

#for (i in 1:1) {
for (i in 1:length(lambda)) {
  png(file=paste('color_map',lambda[i],'_part2.png',sep=''))
  
  ref_granice_pre<-log10(original_ref_pre[,,lambda[i]])
  ref_granice_pre[which(is.infinite(ref_granice_pre))] <- NA 
  
  ref_granice_pre[which(ref_granice_pre < ( mean(ref_granice_pre, na.rm = T) - 2.5* sd(ref_granice_pre, na.rm = T)))]<-mean(ref_granice_pre, na.rm = T) - 2.5* sd(ref_granice_pre, na.rm = T)
  ref_granice_pre[which(ref_granice_pre > ( mean(ref_granice_pre, na.rm = T) + 2.5* sd(ref_granice_pre, na.rm = T)))]<-mean(ref_granice_pre, na.rm = T) + 2.5* sd(ref_granice_pre, na.rm = T)
  
  granice_ref<-seq(min(ref_granice_pre, na.rm = TRUE), max(ref_granice_pre, na.rm = TRUE),length.out = 100)
  
  # residuals original input vs ref
  
  res_original_map<- abs(original_ref_pre[,,lambda[i]] - original_pre[,,lambda[i]]) / original_ref_pre_r[,,lambda[i]] * 100
  res_original_map[which(is.infinite(res_original_map))] <- NA
  
  # median of residuals (number on plot)
  res_median_original<-res_original_map
  res_median_original[which(res_median_original == 0)] <- NA
  res_median_original<-median(res_median_original, na.rm = TRUE)
  
  res_original_map<-log10(res_original_map)
  
  res_map_granice_original<-res_original_map    
  
  res_map_granice_original[which(res_original_map < ( mean(res_original_map,na.rm = T) - 2.5* sd(res_original_map,na.rm = T)))]<-mean(res_original_map,na.rm = T) - 2.5* sd(res_original_map,na.rm = T)
  res_map_granice_original[which(res_original_map > ( mean(res_original_map,na.rm = T) + 2.5* sd(res_original_map,na.rm = T)))]<-mean(res_original_map,na.rm = T) + 2.5* sd(res_original_map,na.rm = T)
  
  
  ################################ pure inla 
  ref_granice_inla<-log10(pure_inla_pre[,,lambda[i]])
  ref_granice_inla[which(is.infinite(ref_granice_inla))] <- NA 
  
  ref_granice_inla[which(ref_granice_inla < ( mean(ref_granice_pre, na.rm = T) - 2.5* sd(ref_granice_pre, na.rm = T)))]<-mean(ref_granice_pre, na.rm = T) - 2.5* sd(ref_granice_pre, na.rm = T)
  ref_granice_inla[which(ref_granice_inla > ( mean(ref_granice_pre, na.rm = T) + 2.5* sd(ref_granice_pre, na.rm = T)))]<-mean(ref_granice_pre, na.rm = T) + 2.5* sd(ref_granice_pre, na.rm = T)
  #
  
  
  #residuals
  res_map<- abs(original_ref_pre[,,lambda[i]] - pure_inla_pre[,,lambda[i]]) / original_ref_pre_r[,,lambda[i]] * 100
  res_map[which(is.infinite(res_map))] <- NA  
  
  # median of residuals (number on plot)
  res_median<-res_map
  res_median[which(res_median == 0)] <- NA
  res_median<-median(res_median, na.rm = TRUE)
  
  res_map<-log10(res_map)
  
  #  inla vs ref
  
  res_map_granice<-res_map  
  
  res_map_granice[which(res_map < ( mean(res_map,na.rm = T) - 2.5* sd(res_map,na.rm = T)))]<-mean(res_map,na.rm = T) - 2.5* sd(res_map,na.rm = T)
  res_map_granice[which(res_map > ( mean(res_map,na.rm = T) + 2.5* sd(res_map,na.rm = T)))]<-mean(res_map,na.rm = T) + 2.5* sd(res_map,na.rm = T)
  
  # color scale za residuals
  granice<-seq(min(res_map_granice, na.rm = TRUE), max(res_map_granice, na.rm = TRUE),length.out = 100)
  
  
  ###################################### PCA
  ref_granice_pca <-log10(pca_post[,,lambda[i]])
  ref_granice_pca[which(is.infinite(ref_granice_pca))] <- NA 
  
  ref_granice_pca[which(ref_granice_pca < ( mean(ref_granice_pre, na.rm = T) - 2.5* sd(ref_granice_pre, na.rm = T)))]<-mean(ref_granice_pre, na.rm = T) - 2.5* sd(ref_granice_pre, na.rm = T)
  ref_granice_pca[which(ref_granice_pca > ( mean(ref_granice_pre, na.rm = T) + 2.5* sd(ref_granice_pre, na.rm = T)))]<-mean(ref_granice_pre, na.rm = T) + 2.5* sd(ref_granice_pre, na.rm = T)
  #
    
  
  #residuals
  res_map_pca<- abs(original_ref_pre[,,lambda[i]] - pca_post[,,lambda[i]]) / original_ref_pre_r[,,lambda[i]] * 100
  
  res_map_pca[which(is.infinite(res_map_pca))] <- NA  
  
  # median of residuals (number on plot)
  res_median_pca<-res_map_pca
  res_median_pca[which(res_median_pca == 0)] <- NA
  res_median_pca<-median(res_median_pca, na.rm = TRUE)
  
  res_map_pca<-log10(res_map_pca)
  
  # vrednosti za plot residual pca vs ref
  
  res_map_granice_pca<-res_map_pca  
  res_map_granice_pca[which(is.infinite(res_map_granice_pca))] <- NA  
  
  
  res_map_granice_pca[which(res_map_pca < ( mean(res_map,na.rm = T) - 2.5* sd(res_map,na.rm = T)))]<-mean(res_map,na.rm = T) - 2.5* sd(res_map,na.rm = T)
  res_map_granice_pca[which(res_map_pca > ( mean(res_map,na.rm = T) + 2.5* sd(res_map,na.rm = T)))]<-mean(res_map,na.rm = T) + 2.5* sd(res_map,na.rm = T)
  

  ########################################
  
  

  ###################################### NMF
  ref_granice_nmf <-log10(nmf_post[,,lambda[i]])
  ref_granice_nmf[which(is.infinite(ref_granice_nmf))] <- NA 
  
  ref_granice_nmf[which(ref_granice_nmf < ( mean(ref_granice_pre, na.rm = T) - 2.5* sd(ref_granice_pre, na.rm = T)))]<-mean(ref_granice_pre, na.rm = T) - 2.5* sd(ref_granice_pre, na.rm = T)
  ref_granice_nmf[which(ref_granice_nmf > ( mean(ref_granice_pre, na.rm = T) + 2.5* sd(ref_granice_pre, na.rm = T)))]<-mean(ref_granice_pre, na.rm = T) + 2.5* sd(ref_granice_pre, na.rm = T)
  #
  
  
  #residuals
  res_map_nmf<- abs(original_ref_pre[,,lambda[i]] - nmf_post[,,lambda[i]]) / original_ref_pre_r[,,lambda[i]] * 100
  
  res_map_nmf[which(is.infinite(res_map_nmf))] <- NA  
  
  # median of residuals (number on plot)
  res_median_nmf<-res_map_nmf
  res_median_nmf[which(res_median_nmf == 0)] <- NA
  res_median_nmf<-median(res_median_nmf, na.rm = TRUE)
  
  res_map_nmf<-log10(res_map_nmf)
  
  # residual nmf vs ref
  
  res_map_granice_nmf<-res_map_nmf 
  res_map_granice_nmf[which(is.infinite(res_map_granice_nmf))] <- NA  
  
  
  res_map_granice_nmf[which(res_map_nmf < ( mean(res_map,na.rm = T) - 2.5* sd(res_map,na.rm = T)))]<-mean(res_map,na.rm = T) - 2.5* sd(res_map,na.rm = T)
  res_map_granice_nmf[which(res_map_nmf > ( mean(res_map,na.rm = T) + 2.5* sd(res_map,na.rm = T)))]<-mean(res_map,na.rm = T) + 2.5* sd(res_map,na.rm = T)
  
  
  ########################################

 
  ####################################

  
  Van_Gogh <- c("#263C8B","#4E74A6", "#BDBF78", "#BFA524", "#2E231F")
  colpal<-colorRampPalette(Van_Gogh)(100)
  cutColor <- 100
  
  
  my.padding <- list(layout.heights= list(top.padding=0.0, bottom.padding=0.0),
       layout.widths= list(left.padding=0.0, right.padding=0.0))
  

 
  p5<-levelplot(ref_granice_pca,col.regions=colpal,
                cut=cutColor,at=granice_ref, main=list(paste('PCA+INLA')), par.settings = my.padding, 
                scales=list(draw=FALSE),
                xlab='', ylab=''
                #xlab=expression(bold(paste(Delta, 't= 5h 40m'),))
                
               
                
  ) 
  
  p6<-levelplot(res_map_granice_pca,col.regions=colpal,
                cut=cutColor,at=granice, 
                main=list(paste('Res. [%] ',
                                'median=', trunc(res_median_pca)),cex=0.8),
                # main=list(paste('Res. [%] (INLA to HPN ref.)', 
                #    'median=', format(round(res_median, 2), nsmall = 2)),cex=0.8),
                
                par.settings = my.padding, scales=list(draw=FALSE),
                xlab='', ylab='')
  
  
  p7<-levelplot(ref_granice_nmf,col.regions=colpal,
                cut=cutColor,at=granice_ref, main=list(paste('NMF+INLA')), par.settings = my.padding, 
                scales=list(draw=FALSE),
                xlab='', ylab=''
                #xlab=expression(bold(paste(Delta, 't= 5h 40m'),))
                
                
                
  ) 
  
  p8<-levelplot(res_map_granice_nmf,col.regions=colpal,
                cut=cutColor,at=granice, 
                main=list(paste('Res. [%] ',
                                'median=', trunc(res_median_nmf)),cex=0.8),
                # main=list(paste('Res. [%] (INLA to HPN ref.)', 
                #    'median=', format(round(res_median, 2), nsmall = 2)),cex=0.8),
                
                par.settings = my.padding, scales=list(draw=FALSE),
                xlab='', ylab='')
  
  
  grid.arrange(grobs = list( p5, p7, p6, p8), ncol=2, common.legend = TRUE, legend="right")
  
  
  
  dev.off()
  
}




