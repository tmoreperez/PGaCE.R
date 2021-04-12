#This script will extract occurrence data from the GBIF database, retrieves climate
#data for each occurrence, and remove potentially erroneus occurrence data that is
#outside of a 95% confidence interval for selected climatic variables. For species
#with more than 15 occurrences, the mean climatic distribution is calculated. 

#The user will have to provide the file path to the folder where their climate data
#is stored. 

#Load required libraries
library(spocc)
library(raster)
library(sp)
library(foreach)


#Load global climate data. If you don't have climwte data, it can be obtained from 
#this website: https://www.worldclim.org/

#Note: that there are several different databases with cliamte data. The gridded 
#cliamte data that can be downloaded from each database also has different spatial 
#resolution. In this example, I use 2.5km grid cells. You should also check out 
#the documentation for this function. Individual coordinates are needed for for 
#grid cells with 0.5km resolution. 

#Also please note that this version of code does not deal with any potentially 
#duplciated entries. The user should decide what they think is the best way to 
#deal with duplicates and add the appropriate code. I usually just use the 
#"duplicated" function to identify dupliated data. 

climdat =getData("worldclim", var="bio", res=2.5, path="/your path to your favorite climate data/")
str(climdat)
wcd <- climdat[[c(1:19)]] #1:19 retrieves all of the available Worldclim climate variables

#For large datasets, multicore processing may be more effective, but I actually 
#found this to be fairly slow for the occurrences of ~150 species. I played 
#around with the apply family functions and the doParallel package, but they 
#curiosuly were slower than the code below.

#Provide the names of the species you'd like occurrence data and climate data for
splist=unique(a.vector.species.names) 
splist=as.character("Acer rubrum")
#Detect cores, make a cluster with one free core as to not overload your computer
registerDoParallel(makeCluster(detectCores()[1]-1))
sclimout=foreach(i=1:length(splist), .combine=rbind, .packages = c("spocc", "sp", "raster")) %dopar% {
  #Get geographic occurence data from the GBIF database. See the documentation 
  #"?occ" for more information about this function
  fwdist=as.data.frame(spocc::occ2df(occ(splist[1], from="gbif", limit=100000, gbifopts = list(hasCoordinate = T))))
  
  if( length(fwdist[,1])>1){ #Make sure that there is occurrence data for a species 
    fwdist$unscrubbed.name=rep(splist[i], length(fwdist[,1]))#Paste the original name you entered into the extracted data.frame (GBIF will return corrected/"scrubbed" names)
    if( length(which( complete.cases(fwdist[,1:3])==T) )>=15 ){ #Get complete cases and only use if there are 15 or more occurrences
      fwdist2=fwdist[complete.cases(fwdist[,1:3]),] #drop columns that aren't needed
      
      #Extract the climate data for each coordinate
      coords = data.frame(longitude=fwdist2$longitude, latitude=fwdist2$latitude)
      points = sp::SpatialPoints(coords, proj4string =wcd@crs)
      values= raster::extract(wcd,points)
      climdat2 = cbind.data.frame(raster::coordinates(points),values, 
                                  names= rep(splist[i], length(fwdist2[,1])),
                                   retrieved_species_binomial=fwdist2$name,
                                   date_collected= fwdist2$date)
      
      #Calculate the upper and lower confidence intervals the following climatic variables:
      #BIO1 = Annual Mean Temperature
      uqb1=quantile(na.omit(climdat2$bio1), 0.975)[[1]]
      lqb1=quantile(na.omit(climdat2$bio1), 0.025)[[1]]
      #BIO5 = Max Temperature of Warmest Month
      uqb5=quantile(na.omit(climdat2$bio5), 0.975)
      lqb5=quantile(na.omit(climdat2$bio5), 0.025)
      #BIO6 = Min Temperature of Coldest Month 
      uqb6=quantile(na.omit(climdat2$bio6), 0.975)
      lqb6=quantile(na.omit(climdat2$bio6), 0.025)
      #BIO12 = Annual Precipitation
      uqb12=quantile(na.omit(climdat2$bio12), 0.975)
      lqb12=quantile(na.omit(climdat2$bio12), 0.025)
      
      #Remove occurrences for a species if it falls outside of the 95% CI of 
      #climatic distribution
      cleaned.fwdist.clim=climdat2[-which(na.omit(climdat2$bio1)>=uqb1 | na.omit(climdat2$bio1)<=lqb1 | 
                                            na.omit(climdat2$bio5)>=uqb5 | na.omit(climdat2$bio5)<=lqb5|
                                            na.omit(climdat2$bio6)>=uqb6 | na.omit(climdat2$bio6)<=lqb6|
                                            na.omit(climdat2$bio12)>=uqb12 | na.omit(climdat2$bio12)<=lqb12),]
      
      if(length(cleaned.fwdist.clim[,1])>=15){ #Again, only use if there are 15 or more occurrences.
      
      #Calculate the mean of the climatic distribitions for variables of interest.
      sclimout[[i]]=data.frame(retrieved_species_binomial=unique(cleaned.fwdist.clim$retrieved_species_binomial),
                               unscrubbed.name=unique(cleaned.fwdist.clim$names),
                               bio1.mn=mean(na.omit(cleaned.fwdist.clim$bio1)),
                               bio2.mn=mean(na.omit(cleaned.fwdist.clim$bio2)),
                               bio3.mn=mean(na.omit(cleaned.fwdist.clim$bio3)),
                               bio4.mn=mean(na.omit(cleaned.fwdist.clim$bio4)),
                               bio5.mn=mean(na.omit(cleaned.fwdist.clim$bio5)),
                               bio6.mn=mean(na.omit(cleaned.fwdist.clim$bio6)),
                               bio7.mn=mean(na.omit(cleaned.fwdist.clim$bio7)),
                               bio12.mn=mean(na.omit(cleaned.fwdist.clim$bio12)))}
    }
  }
}
stopImplicitCluster() #stop parallel processing

#Here's the output:
str(sclimout)
