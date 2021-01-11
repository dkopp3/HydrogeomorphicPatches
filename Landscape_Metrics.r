# StreamCLimes stream Classification
# calculate landscape metrics for resolutions 

setwd("C:/Users/AllenLabWorkstation/Dropbox/Dissertation/StreamClimes/Manuscript/Data")
source("NHDHR_Functions_Example.R")

#networks outlets form StreamClimes 
#outlet needs to be formatted exactly as below
x <- read.csv("NEON_outlet_example.csv", colClasses ="character" )
x <- x[,c("site", "HUC4", "NHDPLUSID")]
names(x) <- c("Network_ID", "HUC4", "NHDPlusID")

# Extent
#####
for(i in x$Network_ID){
  T1 <- Sys.time()
  #i <- x$Network_ID[1]
  
  print(i)
  # data 
  #####
  file.dir <- grep(i, list.dirs(), value = T)
  
  #HGP data 
  HGP <- read.csv(paste0(file.dir, "/Network_HGP.csv"), 
                  colClasses = c("NHDPlID" = "character"))
  #impute data
  impute <- read.csv(paste0(file.dir, "/FPZ_impute.csv"), 
                     colClasses = c("split_fID" = "character"))
  
  network <- read_sf(paste0(file.dir,"/Reach_segments.shp"))
  
  patches <- merge(HGP[,c("split_fID","NHDPlID")], impute, by = "split_fID")
  network <- merge(network, patches, by = "split_fID")

  #"scaling_grain.csv" 
  ext <- read.csv(paste0(file.dir,"/scaling_extent.csv"),
                    colClasses = c("NHDPlusID" = "character", 
                                   "net_NHDPlusID" ="character" ))
  
  #resample scaling values - one network every 10m, where possible 
  #this will change for msl
  da <- unique(ext[,c("NHDPlusID", "fac_cat")])
  da <- da[da$fac_cat > 10, ]
  
  da$group <- cut(da$fac_cat, seq(0, max(da$fac_cat), by = 10), labels = FALSE)
  maxval <- da[which.max(da$fac_cat),]
  da <- do.call(rbind, lapply(split(da, f = da$group), function(x) x[sample(1:nrow(x), size = 1),]))
  
  da <- rbind(da, maxval)
  da <- da[!duplicated(da$fac_cat), ]
  
  ext <- ext[ext$NHDPlusID %in% da$NHDPlusID, ]
  ext$fac_cat <- round(ext$fac_cat) 
  #############
  
  # calculate 
  #####
  # identify subnetworks 
  q <- split(ext, ext$fac_cat)
  #length(q)
  #q[["244"]]
  #class_scale_metrics(network = network, subnetID = q[["244"]]$net_NHDPlusID)
  #T1 <- Sys.time()
  #for (test in 1:length(q)){
  z <- lapply(q, function(x) class_scale_metrics(network, subnetID = x$net_NHDPlusID))
  #}
  #T2 <- Sys.time()
  #T2-T1
  
  #compile landscape scale functions
  LS_metrics <- lapply(z, function(x) 
    data.frame(tlen.patch_m = sum(x$tlen.patch_m),
               SHDI = (-1 * (sum((x$tlen.patch_m/sum(x$tlen.patch_m))*
                                   (log(x$tlen.patch_m/sum(x$tlen.patch_m)))))), 
               No.Patch = sum(x$No.Patch), 
               mean.patch_m = mean(x$mean.patch_m),
               richness = nrow(x),
               DCI = mean(x$DCI, na.rm = T), 
               meandist_m = mean(x$meandist, na.rm = T),
               mediandist_m = median(x$medianDist, na.rm = T)))
  
  LS_metrics <- do.call(rbind, LS_metrics)
  LS_metrics <- data.frame(extent = as.numeric(rownames(LS_metrics)), LS_metrics)
  
  cL_metrics <- do.call(rbind, z)
  extent <- do.call(c, lapply(strsplit(row.names(cL_metrics),"\\."),"[[",1))
  cL_metrics <- data.frame(extent, cL_metrics)
  #############
  
  write.csv(cL_metrics, paste0(file.dir,"/cL_metrics.csv"), row.names = F)
  write.csv(LS_metrics, paste0(file.dir,"/LS_metrics.csv"), row.names = F)
  T2<-Sys.time()
  print(T2-T1)
}
#############
  
# GRAIN
#####
for(i in x$Network_ID){
  T1 <- Sys.time()
  #i <- x$Network_ID[3]
  #i<-"LECO"
  print(i)
  # data 
  #####
  file.dir <- grep(i, list.dirs(), value = T)
  
  #HGP data 
  HGP <- read.csv(paste0(file.dir, "/Network_HGP.csv"), 
                  colClasses = c("NHDPlID" = "character"))
  #impute data
  impute <- read.csv(paste0(file.dir, "/FPZ_impute.csv"), 
                     colClasses = c("split_fID" = "character"))
  
  network <- read_sf(paste0(file.dir, "/Reach_segments.shp"))
  
  patches <- merge(HGP[,c("split_fID","NHDPlID")], impute, by = "split_fID")
  network <- merge(network, patches, by = "split_fID")
  
  grain <- read.csv(paste0(file.dir, "/scaling_grain.csv"),
                    colClasses = c("NHDPlID" = "character"))
  names(grain)[2] <- "NHDPlusID"
  
  #unique(grain$mainstemlength)
  #scaleing values
  msl <- unique(grain[,c("NHDPlusID","mainstemlength")])
  #head(msl)
  ###############
  
  # calculate 
  #####
  # identify subnetworks 
  q <- split(grain, grain$mainstemlength)
  
  #q[[62]]
  #T1 <- Sys.time()
  #for (test in 1:length(q)){
  z <- lapply(q, function(x) class_scale_metrics(network, subnetID = x$NHDPlusID))  
   # print(test)
  #} 
  #T2 <- Sys.time()
  #T2-T1
  
  #compile landscape scale functions
  LS_metrics <- lapply(z, function(x) 
    data.frame(tlen.patch_m = sum(x$tlen.patch_m),
               SHDI = (-1 * (sum((x$tlen.patch_m/sum(x$tlen.patch_m))*
                                   (log(x$tlen.patch_m/sum(x$tlen.patch_m)))))), 
               No.Patch = sum(x$No.Patch), 
               mean.patch_m = mean(x$mean.patch_m),
               richness = nrow(x),
               DCI = mean(x$DCI, na.rm = T), 
               meandist_m = mean(x$meandist, na.rm = T),
               mediandist_m = median(x$medianDist, na.rm = T)))
  
  LS_metrics <- do.call(rbind, LS_metrics)
  LS_metrics <- data.frame(msl = as.numeric(rownames(LS_metrics)), LS_metrics)
  
  cL_metrics <- do.call(rbind, z)
  msl <- do.call(c, lapply(strsplit(row.names(cL_metrics),"\\."),"[[",1))
  cL_metrics <- data.frame(msl, cL_metrics)
  #############
  
  #write.csv(cL_metrics, paste0(file.dir,"/cL_metrics_grain.csv"), row.names = F)
  #write.csv(LS_metrics, paste0(file.dir,"/LS_metrics_grain.csv"), row.names = F)
  T2<-Sys.time()
  print(T2-T1)
}  
#############

# Thematic
#####
orderDA <- c("MAYF")

for(i in orderDA){# x$Network_ID){
  T1 <- Sys.time()
  #i <- x$Network_ID[3]
  print(i)
  # data 
  #####
  file.dir <- grep(i, list.dirs(), value = T)
  #list.files(file.dir)
  #HGP data 
  HGP <- read.csv(paste0(file.dir, "/Network_HGP.csv"), 
                  colClasses = c("NHDPlID" = "character"))
  #impute data
  impute <- read.csv(paste0(file.dir, "/FPZ_impute_all.csv"), 
                     colClasses = c("split_fID" = "character"))
  
  network <- read_sf(paste0(file.dir, "/Reach_segments.shp"))
  patches <- merge(HGP[,c("split_fID","NHDPlID")], impute, by = "split_fID")
  network <- merge(network, patches, by = "split_fID")
  
  #cluster scaling
  clusters <- grep("fpz|FPZ", names(impute), value = T)
  ###############

  # calculate 
  #####
  # identify subnetworks 
  
  #iterate through each column 
  z <- list()
  clusterindex <- unique(c(seq(1, which(clusters == "fpz_30"), by = 2), 
                           which(clusters == "fpz_30")))
  index <- 1
  for (j in clusterindex){
    #j<-2
    #rename j thematic level to FPZ.sil, so it lays nicely with the function
    FPZ.sil <- which(names(network) == clusters[j])
    network2 <- network[, c(1:4, FPZ.sil)]
    names(network2)[5] <- "FPZ.sil"
    ps <- class_scale_metrics(network2, subnetID = network2$NHDPlID)
    z[[index]] <- ps
    index <- index + 1
    Tmid <- Sys.time()
    print(paste("Finished", clusters[j]))
    #print(Tmid-T1)
  }

  names(z) <- clusters[clusterindex]

  #compile landscape scale functions
  LS_metrics <- lapply(z, function(x) 
    data.frame(tlen.patch_m = sum(x$tlen.patch_m),
               SHDI = (-1 * (sum((x$tlen.patch_m/sum(x$tlen.patch_m))*
                                   (log(x$tlen.patch_m/sum(x$tlen.patch_m)))))), 
               No.Patch = sum(x$No.Patch), 
               mean.patch_m = mean(x$mean.patch_m),
               richness = nrow(x),
               DCI = mean(x$DCI, na.rm = T), 
               meandist_m = mean(x$meandist, na.rm = T),
               mediandist_m = median(x$medianDist, na.rm = T)))
  
  LS_metrics <- do.call(rbind, LS_metrics)
  LS_metrics <- data.frame(clusters = as.numeric(substring(rownames(LS_metrics),5)), LS_metrics)
  
  cL_metrics <- do.call(rbind, z)
  cust <- do.call(c, lapply(strsplit(row.names(cL_metrics),"\\."),"[[",1))
  cL_metrics <- data.frame(clusters = cust, cL_metrics)
  ###############
  
  #write.csv(cL_metrics, paste0(file.dir,"/cL_metrics_thematic.csv"), row.names = F)
  #write.csv(LS_metrics, paste0(file.dir,"/LS_metrics_thematic.csv"), row.names = F)
  T2 <- Sys.time()
  print(T2-T1)
}
#############




