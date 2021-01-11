# StreamCLimes stream Classification
# Cluster module - identify hydrogeomorphic patches 
# must have run create and attribute modules 


# clustering routine
# writes hclust object 
#####

setwd("C:/Users/AllenLabWorkstation/Dropbox/Dissertation/StreamClimes/Manuscript/Data")
source("NHDHR_Functions_Example.R")

#networks outlets form StreamClimes 
#outlet needs to be formatted exactly as below
x <- read.csv("NEON_outlet_example.csv", colClasses = "character" )
x <- x[,c("site", "HUC4", "NHDPLUSID")]
names(x) <- c("Network_ID", "HUC4", "NHDPlusID")

#iterates thru each network
for (site in 1:nrow(x)) {
  #site <- 1
  print(x[site, ])
    
  # Initialize directory
  #####
  write.dir <- c(paste0(getwd(), "/", x[site, "Network_ID"]))
  dir.create(write.dir)
  
  setwd("./NHDPLUSHR_3_7_2019")
  raster.dir <- paste0("./NHDPLUS_H_", x[site, "HUC4"], "_HU4_RASTER/HRNHDPlusRasters", x[site, "HUC4"])
  ancillary.dir <- "./ancillary_data"
  
  ###################################
    
  # cluster - saves cluster results and dendrograms
  #####
  # classify all points within network
  # eight variables form Maasri 2019 + temperature and meander length
  HGP <- read.csv(paste0(write.dir, "/Network_HGP.csv"), 
                    colClasses = c("NHDPlID" = "character"))
  
  FPZ_vars <- c("ele", "PRISMppt30", "PRISMtemp30", "mean_slope", 
                "v_wid", "vf_wid", "v_vf_ratio", "ch_sinu", 
                "mean_meander_length", "v_slope_precent",
                "dtb_cm", "kwfact", "ph1t12_") 
    
  # need complete case for FPZ_vars otherwise unclassified sections
  if(!all(complete.cases(HGP[, FPZ_vars]))){
    #save sites with missing variables - impute classification 
      vars_NA <- HGP[!complete.cases(HGP[,FPZ_vars]), ]
      HGP <- HGP[complete.cases(HGP[, FPZ_vars]), ]
    }
  
  #reclassify erosion polygons
  #Medium textured soils, such as the silt loam soils, have a moderate K values, 
  #about 0.25 to 0.4, because they are moderately susceptible 
  #to detachment and they produce moderate runoff.
  HGP$kwfact <- cut(HGP$kwfact, breaks = c(0, 0.25, 0.4, 0.6), labels = c("high", "medium", "low"))
  
  #Most living organisms, especially aquatic life, 
  #function at the optimal pH range of 6.5 to 8.5. (USEPA NRSA Acidification)
  HGP$ph1t12_ <- cut(HGP$ph1t12_, breaks = c(0, 6.50, 8.50, 14.00), labels = c("acidic", "neutral", "basic"))
    
  #read NHDWaterbody Attribute table  
  WB.attr <- suppressWarnings(st_read(dsn = grep(paste0(x[site, "HUC4"], "_HU4_GDB.gdb"), list.dirs(), value = T),
                                      layer = "NHDWaterbody", 
                                      quiet = T,  stringsAsFactors = T, as_tibble=F))

  #only include those that were listed as reservoir (fcode = 436)
  WB.attr <- WB.attr[WB.attr$FType == 436, ]
  #remove waterbdies = already classified "waterbodies" no need for clustering
  if(any(HGP[!is.na(HGP$WBA_P_I), "WBA_P_I"] %in%  as.character(WB.attr$Permanent_Identifier))){
    WB <- HGP[!is.na(HGP$WBA_P_I), ]
    WB <- WB[WB$WBA_P_I %in% as.character(WB.attr$Permanent_Identifier),]
    HGP <- HGP[!HGP$WBA_P_I %in% WB$WBA_P_I, ]
    } else {
      WB <- data.frame()
    }
  
  #identify FPZ, nmax is the maximium iterations to find avg siloutee width  
  nmax <- ifelse(nrow(HGP) > 25, 25, nrow(HGP) - 1)
    
  fpz <- FPZ(x = HGP, FPZ_vars = FPZ_vars, 
             normalize = F, nmax = nmax, 
             saveplot = T, savecluster = T, 
             write.dir = write.dir)

  fpz <- apply(fpz, 2, function(x) as.character(x))
  HGP <- data.frame(HGP, FPZ.sil = as.character(fpz), stringsAsFactors = F)
    
  #add back in Waterbodies
  if (nrow(WB) > 0){
      fpzs <- fpz[1,]
      WB <- data.frame(WB, t(data.frame(fpzs)), stringsAsFactors = F)
      WB[,names(fpzs)] <- "Waterbody"
      HGP <- rbind(HGP, WB)
      
      WB <- WB[WB$WBA_P_I %in% "emptydataframe",]
    }

  #add back in missing data reaches
  if (nrow(vars_NA) > 0){
      fpzs <- fpz[1,]
      vars_NA <- data.frame(vars_NA, t(data.frame(fpzs)), stringsAsFactors = F)
      vars_NA[,names(fpzs)] <- "Unclassified"
      HGP <- suppressWarnings(rbind(HGP, vars_NA))
      vars_NA <- vars_NA[vars_NA$SegID %in%"emptydataframe", ]
      unique(HGP$FPZ.sil)
    }
  
  ###################################

  #print(levels(HGP$FPZ))
  write.csv(HGP[, c("split_fID", colnames(fpz))], 
            paste0(write.dir, "/FPZ_2.csv"), row.names = F)
} 

###########

# impute any missing fpz classes
#####
# tyically because missing values in sideslopes
# fit randomforest on groups with variables 
# that are present and use those to predict 
# unclassified. 
# warning OK. Initated at "mods" (replecated variable vector)
for (site in 1:nrow(x)) {
  #site <- 3
  print(x[site, ])
  
  
  #Impute
  #####
  #options(warn = 1)
  HGP <- read.csv(paste0(write.dir, "/Network_HGP.csv"), 
                  colClasses = c("NHDPlID" = "character"))
  
  FPZ_2 <- read.csv(paste0(write.dir, "/FPZ_2.csv"), 
                    colClasses = c("split_fID" = "character"))
  
  FPZ_vars <- c("ele", "PRISMppt30", "PRISMtemp30", "mean_slope", 
                "v_wid", "vf_wid", "v_vf_ratio", "ch_sinu", 
                "mean_meander_length", "v_slope_precent",
                "dtb_cm", "kwfact", "ph1t12_") 
  
  z <- merge(HGP, FPZ_2, by = "split_fID")
  
  #segments with missing vlaues
  Unclass_segs <- z[z$FPZ.sil == "Unclassified", ]
  
  #names of missing variables 
  mods <- apply(Unclass_segs[,FPZ_vars], 1, 
                function(x) names(x[is.na(x)]))
  
  #random forest classification
  if(class(mods) == "list"){
    mods <- do.call(rbind, mods)
    dropvars <- unique(mods)
    impute <- data.frame()
    
    for (i in 1:nrow(dropvars)){
      #i <- 1
      dv <- dropvars[i,]
      newvars <- FPZ_vars[!FPZ_vars %in% dv]
      
      #predict unclassified sites
      ucrows <- t(apply(mods, 1, function(x) all(x == dv)))
      predictdata <- Unclass_segs[which(ucrows), newvars]
      
      # Iterate through each of the responses 
      # for different clustering solutions...
      outz <- data.frame(split_fID = Unclass_segs[which(ucrows), "split_fID"])
      for (j in grep("FPZ|fpz", names(z), value = T)){ 
        #j <- "fpz_8" 
        
        fitdata <- z[z$FPZ.sil != "Unclassified" & 
                       z$FPZ.sil != "Waterbody", 
                     c(j, newvars)]
        
        #renames response variable to work with RF
        names(fitdata)[1] <- "FPZ.sil"
        
        #fit rf model 
        fitdata$FPZ.sil <- as.factor(as.character(fitdata$FPZ.sil)) 
        rf.mod <- randomForest(FPZ.sil ~ ., data = fitdata)
        
        temp <- data.frame(predict(rf.mod, newdata = predictdata))
        names(temp) <- j
        
        #temp <- data.frame(split_fID = Unclass_segs[which(ucrows), "split_fID"], 
        #                 FPZ.sil = predict(rf.mod, newdata = predictdata))
        outz <- data.frame(outz, temp)
      }
      
      #save resuts
      impute <- rbind(impute, outz)
    }
    
  } else {
    
    dropvars <- unique(mods)
    impute <- data.frame()
    
    for (i in 1:length(dropvars)){
      # i <- 1
      dv <- dropvars[i]
      newvars <- FPZ_vars[!FPZ_vars %in% dv]
      
      #predict unclassified sites
      ucrows <- mods == dv
      predictdata <- Unclass_segs[which(ucrows), newvars]
      
      # Iterate through each of the responses 
      # for different clustering solutions...
      
      outz <- data.frame(split_fID = Unclass_segs[which(ucrows), "split_fID"])
      for (j in grep("FPZ|fpz", names(z), value = T)){
        
        fitdata <- z[z$FPZ.sil != "Unclassified" & 
                       z$FPZ.sil != "Waterbody", 
                     c(j, newvars)]
        
        
        #renames response variable to work with RF
        names(fitdata)[1] <- "FPZ.sil"
        
        #fit rf model 
        fitdata$FPZ.sil <- as.factor(as.character(fitdata$FPZ.sil)) 
        rf.mod <- randomForest(FPZ.sil ~ ., data = fitdata)
        
        temp <- data.frame(predict(rf.mod, newdata = predictdata))
        names(temp) <- j
        
        #temp <- data.frame(split_fID = Unclass_segs[which(ucrows), "split_fID"], 
        #                 FPZ.sil = predict(rf.mod, newdata = predictdata))
        outz <- data.frame(outz, temp)
      }
      
      impute <- rbind(impute, outz)
    }
  }
  
  FPZ_2 <- FPZ_2[FPZ_2$FPZ.sil != "Unclassified", ]
  FPZ_2 <- rbind(FPZ_2, impute)
  ###################################  
  
  write.csv(FPZ_2, paste0(write.dir, "/FPZ_impute.csv"), row.names = F)
}

############


# manipulate scales of each network
# identify different  extent (drainage area),
# grain (minumum path length) and thematic resoluton
# these files are used to quantify spatial pattern at different scales 

# extent
#####
# StreamCLimes stream Classification
setwd("C:/Users/AllenLabWorkstation/Dropbox/Dissertation/StreamClimes/Manuscript/Data")
source("NHDHR_Functions_Example.R")

#networks outlets form StreamClimes 
#outlet needs to be formatted exactly as below
x <- read.csv("NEON_outlet_example.csv", colClasses ="character" )
x <- x[,c("site", "HUC4", "NHDPLUSID")]
names(x) <- c("Network_ID", "HUC4", "NHDPlusID")

# sample watersheds - scaling extent
# iterate through each each watershed
for (site in 1:nrow(x)) {
  #site <- 3
  print(x[site, ])
  
  # Initialize
  #####
  setwd("C:/Users/AllenLabWorkstation/Dropbox/Dissertation/StreamClimes/Manuscript/Data")
  write.dir <- c(paste0(getwd(), "/", x[site, "Network_ID"]))
  dir.create(write.dir)
  
  setwd("./NHDPLUSHR_3_7_2019")
  raster.dir <- paste0("./NHDPLUS_H_", x[site, "HUC4"], "_HU4_RASTER/HRNHDPlusRasters", x[site, "HUC4"])
  ancillary.dir <- "./ancillary_data"
  
  ###################################
  
  #### sample extent
  ###############
  HGP <- read.csv(paste0(write.dir, "/Network_HGP.csv"), 
                  colClasses = c("NHDPlID" = "character"))
  FPZs <- read.csv(paste0(write.dir, "/FPZ_impute.csv"))
  
  #the summary function requires FPZ to be asssociated with HGP
  HGP <- merge(HGP, FPZs, by = "split_fID")
  
  #Idenfiy comids along catchment area gradient
  HGP <- HGP[order(HGP$fac_cat), ]
  
  #stratified random sample with a 
  #random sample of 1% on the watersheds within a class
  sample_index <- HGP[HGP$fac_cat > 0 & 
                        HGP$FPZ != "Waterbody"&
                        HGP$FPZ != "Unclassified", 
                      c("split_fID", "NHDPlID","fac_cat")]
  
  #always include max value
  sample_index$group <- cut(sample_index$fac_cat,
                            seq(0, max(sample_index$fac_cat), by = 10), 
                            labels = FALSE)
  maxval <- sample_index[which.max(sample_index$fac_cat),]
  sample_index <- lapply(split(sample_index, f = sample_index$group), function(x)
    x[sample(1:nrow(x), size = ceiling(nrow(x)/100)),])
  sample_index <- do.call(rbind, sample_index)
  #test[[1]][sample(1:nrow(test[[1]]), size = ceiling(nrow(test[[1]])/100)),]
  
  sample_index <- rbind(sample_index, maxval)
  sample_index <- sample_index[!duplicated(sample_index$NHDPlID), ]
  
  # all flowlines above the root node
  z <- net_delin_HR(huc = x[site, "HUC4"], NHDPlusID = as.character(sample_index$NHDPlID))
  sample_index <- merge(z$network, sample_index, by.x = "NHDPlusID", by.y = "NHDPlID" )
  
  #scaling_FPZ <- net_summary(HPG_Sample = sample_index, HGP = HGP, show_progress = F)
  
  ###################################
  #to keep track of the upstream sampling points
  write.csv(sample_index[, c("NHDPlusID", "net_NHDPlusID", "huc", "split_fID", "fac_cat")], 
            paste0(write.dir,"/scaling_extent.csv"), row.names = F)
}
############

# grain
#####

# StreamCLimes stream Classification
setwd("C:/Users/AllenLabWorkstation/Dropbox/Dissertation/StreamClimes/Manuscript/Data")
source("NHDHR_Functions_Example.R")

#networks outlets form StreamClimes 
#outlet needs to be formatted exactly as below
x <- read.csv("NEON_outlet_example.csv", colClasses ="character" )
x <- x[,c("site", "HUC4", "NHDPLUSID")]
names(x) <- c("Network_ID", "HUC4", "NHDPlusID")


# iterate through each each watershed
for (site in 1:nrow(x)) {
  #site <- 3
  print(x[site, ])
  
  # Initialize
  #####
  setwd("C:/Users/AllenLabWorkstation/Dropbox/Dissertation/StreamClimes/Manuscript/Data")
  write.dir <- c(paste0(getwd(), "/", x[site, "Network_ID"]))
  dir.create(write.dir)
  
  setwd("./NHDPLUSHR_3_7_2019")
  raster.dir <- paste0("./NHDPLUS_H_", x[site, "HUC4"], "_HU4_RASTER/HRNHDPlusRasters", x[site, "HUC4"])
  ancillary.dir <- "./ancillary_data"
  
  ###################################
  
  #### sample/grain
  ###############
  
  network <- read_sf(paste0(write.dir,"/NHDHR_Network.shp"))
  huc = x[site,"HUC4"]
  
  fileGDB <- grep(paste0(huc, "_HU4_GDB.gdb"), list.dirs(), value = T)
  
  #suppress spatial warnings bc its just a table
  vaa.tbl <- suppressWarnings(st_read(dsn = fileGDB, layer = "NHDPlusFlowlineVAA", quiet = T))
  #vaa.tbl$NHDPlusID <- as.character(vaa.tbl$NHDPlusID)
  #drop divergent flowpaths
  vaa.tbl <- vaa.tbl[vaa.tbl[, "StreamOrde"] == vaa.tbl[, "StreamCalc"], ]
  vaa.tbl <- vaa.tbl[,c("NHDPlusID","FromNode", "ToNode", "HydroSeq", "LevelPathI")]
  vaa.tbl <- data.frame(apply(vaa.tbl, 2, as.character))
  
  network <- merge(network, vaa.tbl, by.x ="NHDPlID" , by.y ="NHDPlusID" )
  pathlengths <- lapply(split(network, as.character(network$LevelPathI)), 
                        function(x) data.frame(LevelPathI = unique(as.character(x$LevelPathI)), 
                                               length = as.numeric(sum(st_length(x)))))
  pathlengths <- do.call(rbind, pathlengths)
  pathlengths <- pathlengths[order(pathlengths$length, decreasing = T),]
  
  
  #always include max value - split by 500m, sample path lengths... 
  pathlengths$group <- cut(pathlengths$length,
                           seq(0, max(pathlengths$length), by = 500), 
                           labels = FALSE)
  maxval <- pathlengths[which.max(pathlengths$length),]
  minval <- pathlengths[which.min(pathlengths$length),]
  
  #sample lengths
  samplelengths <- lapply(split(pathlengths, f = pathlengths$group), function(x)
    x[sample(1:nrow(x), size = ceiling(nrow(x)/100)),])
  samplelengths <- do.call(rbind, samplelengths)
  
  samplelengths <- rbind(samplelengths, maxval, minval)
  samplelengths$length <- round(samplelengths$length) 
  
  out <- data.frame()
  for (q in sort(samplelengths$length, decreasing = T)){
    #q <- max(samplelengths$length)[1]
    path.inclu <- as.character(pathlengths [round(pathlengths$length) >= q, "LevelPathI"])
    net_NHDPlID <- st_set_geometry(network[network$LevelPathI %in%path.inclu, "NHDPlID"], NULL)
    temp <- data.frame(mainstemlength = q, net_NHDPlID)
    out <- rbind(out, temp)
    
    #plot(st_geometry(st_zm(network[network$LevelPathI %in% path.inclu,  ], drop = T , what ="ZM")))
    #plot(st_geometry(st_zm(network[network$LevelPathI %in%path.inclu,  ], drop = T , what ="ZM")), 
    #     add = T, col = "black")
  }
  
  ###################################
  #to keep track of the upstream smapling points
  write.csv(out, paste0(write.dir,"/scaling_grain.csv"), row.names = F)
  #write.csv(scaling_FPZ, paste0(write.dir, "/scaling_FPZ_2.csv"), row.names = F)
}
############

# thematic resolution
#####
# StreamCLimes stream Classification
setwd("C:/Users/AllenLabWorkstation/Dropbox/Dissertation/StreamClimes/Manuscript/Data")
source("NHDHR_Functions_Example.R")

#networks outlets form StreamClimes 
#outlet needs to be formatted exactly as below
x <- read.csv("NEON_outlet_example.csv", colClasses = "character" )
x <- x[,c("site", "HUC4", "NHDPLUSID")]
names(x) <- c("Network_ID", "HUC4", "NHDPlusID")

 
for (site in c(1:nrow(x))) {
  #site <- 2
  print(x[site, ])
  
  # Initialize
  #####
  setwd("C:/Users/AllenLabWorkstation/Dropbox/Dissertation/StreamClimes/Manuscript/Data")
  write.dir <- c(paste0(getwd(), "/", x[site, "Network_ID"]))
  dir.create(write.dir)
  
  setwd("./NHDPLUSHR_3_7_2019")
  raster.dir <- paste0("./NHDPLUS_H_", x[site, "HUC4"], "_HU4_RASTER/HRNHDPlusRasters", x[site, "HUC4"])
  ancillary.dir <- "./ancillary_data"
  
  ###################################
  
  #thematic scaling - identify 1:50 clusters
  #as the the number of clusters, the subordinate 
  #members become more similar. 
  #####
  load(paste0(write.dir, "/hclust"))
  fpz <- read.csv(paste0(write.dir, "/FPZ_2.csv"))
  
  #need to remppve waterbodies and unclseeified reasches apppended 
  #at the end of last script 
  to.add <- fpz[which(fpz$FPZ.sil == "Unclassified" | 
                     fpz$FPZ.sil == "Waterbody"), ]
  fpz <- fpz[which(fpz$FPZ.sil != "Unclassified" & 
                   fpz$FPZ.sil != "Waterbody"), ]
  
  sil.value <- max(as.numeric(as.character(fpz$FPZ.sil)), na.rm = T)
  
  #iterate from the sil value too the 50 
  #anything over 50patches might not be tractable 
  #from a management perspective
  
  for (k in (sil.value + 1):50){
    #k <- sil.value + 1
    fpz.new <- cutree(hcluster, k)
    fpz <- data.frame(fpz, fpz.new)
    names(fpz)[ncol(fpz)] <- paste0("fpz_", k)
  }
  
  #add back in the unclassified sections
  to.add <- data.frame(to.add$split_fID, 
                       matrix(rep(to.add$FPZ.sil, (length((sil.value):50))),
                              ncol = (length((sil.value):50))))
  names(to.add)<-names(fpz)
  fpz <- rbind(fpz, to.add)
  ######################
  write.csv(fpz, paste0(write.dir, "/FPZ_2.csv"), row.names = F)
}
###########

# impute missing fpz all thematic resolutons
#####
# tyically bc sideslope
# fit randomforest on groups with variables 
# that are present and use those to predict 
# unclassified. 
# warning OK. Initated at "mods" (replecated variable vector)

for (site in 1:nrow(x)) {
  #site <- 3
  print(x[site, ])
  
  #Impute
  #####
  #options(warn = 1)
  HGP <- read.csv(paste0(write.dir, "/Network_HGP.csv"), 
                  colClasses = c("NHDPlID" = "character"))
  
  FPZ_2 <- read.csv(paste0(write.dir, "/FPZ_2.csv"), 
                    colClasses = c("split_fID" = "character"))
  
  FPZ_vars <- c("ele", "PRISMppt30", "PRISMtemp30", "mean_slope", 
                "v_wid", "vf_wid", "v_vf_ratio", "ch_sinu", 
                "mean_meander_length", "v_slope_precent",
                "dtb_cm", "kwfact", "ph1t12_") 
  
  z <- merge(HGP, FPZ_2, by = "split_fID")
  
  Unclass_segs <- z[z$FPZ.sil == "Unclassified", ]
  
  mods <- apply(Unclass_segs[,FPZ_vars], 1, 
                function(x) names(x[is.na(x)]))
  
  if(class(mods) == "list"){
    mods <- do.call(rbind, mods)
    dropvars <- unique(mods)
    impute <- data.frame()
    
    for (i in 1:nrow(dropvars)){
      #i <- 1
      dv <- dropvars[i,]
      newvars <- FPZ_vars[!FPZ_vars %in% dv]
      
      #predict unclassified sites
      ucrows <- t(apply(mods, 1, function(x) all(x == dv)))
      predictdata <- Unclass_segs[which(ucrows), newvars]
      
      # Iterate through each of the responses 
      # for different clustering solutions...
      outz <- data.frame(split_fID = Unclass_segs[which(ucrows), "split_fID"])
      for (j in grep("FPZ|fpz", names(z), value = T)){ 
        #j <- "fpz_8" 
        
        fitdata <- z[z$FPZ.sil != "Unclassified" & 
                       z$FPZ.sil != "Waterbody", 
                     c(j, newvars)]
        
        #renames response variable to work with RF
        names(fitdata)[1] <- "FPZ.sil"
        
        #fit rf model 
        fitdata$FPZ.sil <- as.factor(as.character(fitdata$FPZ.sil)) 
        rf.mod <- randomForest(FPZ.sil ~ ., data = fitdata)
        
        temp <- data.frame(predict(rf.mod, newdata = predictdata))
        names(temp) <- j
        
        #temp <- data.frame(split_fID = Unclass_segs[which(ucrows), "split_fID"], 
        #                 FPZ.sil = predict(rf.mod, newdata = predictdata))
        outz <- data.frame(outz, temp)
      }
      
      #save resuts
      impute <- rbind(impute, outz)
    }
    
  } else {
    
    dropvars <- unique(mods)
    impute <- data.frame()
    
    for (i in 1:length(dropvars)){
      # i <- 1
      dv <- dropvars[i]
      newvars <- FPZ_vars[!FPZ_vars %in% dv]
      
      #predict unclassified sites
      ucrows <- mods == dv
      predictdata <- Unclass_segs[which(ucrows), newvars]
      
      # Iterate through each of the responses 
      # for different clustering solutions...
      
      outz <- data.frame(split_fID = Unclass_segs[which(ucrows), "split_fID"])
      for (j in grep("FPZ|fpz", names(z), value = T)){
        
        fitdata <- z[z$FPZ.sil != "Unclassified" & 
                       z$FPZ.sil != "Waterbody", 
                     c(j, newvars)]
        
        
        #renames response variable to work with RF
        names(fitdata)[1] <- "FPZ.sil"
        
        #fit rf model 
        fitdata$FPZ.sil <- as.factor(as.character(fitdata$FPZ.sil)) 
        rf.mod <- randomForest(FPZ.sil ~ ., data = fitdata)
        
        temp <- data.frame(predict(rf.mod, newdata = predictdata))
        names(temp) <- j
        
        #temp <- data.frame(split_fID = Unclass_segs[which(ucrows), "split_fID"], 
        #                 FPZ.sil = predict(rf.mod, newdata = predictdata))
        outz <- data.frame(outz, temp)
      }
      
      impute <- rbind(impute, outz)
    }
  }
  
  FPZ_2 <- FPZ_2[FPZ_2$FPZ.sil != "Unclassified", ]
  FPZ_2 <- rbind(FPZ_2, impute)
  ###################################  
 
  write.csv(FPZ_2, paste0(write.dir, "/FPZ_impute_all.csv"), row.names = F)
  
}

###########



