# StreamCLimes stream Classification

setwd("C:/Users/AllenLabWorkstation/Dropbox/Dissertation/StreamClimes/Manuscript/Data")
source("NHDHR_Functions_Example.R")

#networks outlets form StreamClimes 
#outlet needs to be formatted exactly as below
x <- read.csv("NEON_outlet_example.csv", colClasses ="character" )
x <- x[,c("site", "HUC4", "NHDPLUSID")]
names(x) <- c("Network_ID", "HUC4", "NHDPlusID")


#iterate through each each watershed
for (site in 2:nrow(x)) {
 # site<-1
  print(x[site, ])
  
  # Initialize
  #####
  write.dir <- c(paste0(getwd(), "/", x[site, "Network_ID"]))
  dir.create(write.dir)
  
  setwd("./NHDPLUSHR_3_7_2019")
  raster.dir <- paste0("./NHDPLUS_H_", x[site, "HUC4"], "_HU4_RASTER/HRNHDPlusRasters", x[site, "HUC4"])
  ancillary.dir <- "./ancillary_data"
  
  ###################################
    
  # Prepare Ancillary Data
  # only provide projected datasets with example
  # all logincal statements below are FALSE
  #####

  #these are national coverages so only need one
  #used for attributes below
  newCRS <- crs(raster(paste0(raster.dir, "/elev_cm.tif")))
    
  if (!"ppt_ras.tif" %in% list.files(ancillary.dir)){
    ras <- raster(paste0(ancillary.dir, "/PRISM_ppt_30yr_normal_800mM2_all_asc/PRISM_ppt_30yr_normal_800mM2_annual_asc.asc"))
    ras <- projectRaster(ras, crs = newCRS)
    writeRaster(ras, paste0(ancillary.dir, "/ppt_ras.tif"))
    } else {
      #print("ppt_ras.tif already exists")
    }
    
    if (!"temp_ras.tif" %in% list.files(ancillary.dir)){
      ras <- raster(paste0(ancillary.dir, "/PRISM_tmean_30yr_normal_800mM2_all_asc/PRISM_tmean_30yr_normal_800mM2_annual_asc.asc"))
      ras <- projectRaster(ras, crs = newCRS)
      writeRaster(ras, paste0(ancillary.dir,"/temp_ras.tif"))
      rm(ras)
    } else {
      #print("ppt_ras.tif already exists")
    }
    
    #lithology classes
    if(!"geo.shp" %in% list.files(ancillary.dir)){
      geo <- read_sf(paste0(ancillary.dir, "/USGS_SGMC_Shapefiles/USGS_SGMC_Shapefiles/SGMC_Geology.shp"))
      geo <- st_transform(geo, crs = 5070)
      write_sf(geo, paste0(ancillary.dir,"/geo.shp")) #throws warings on lenght of area field no a problem
      rm(geo)
    } else {
      #print("geo.shp already exists")
    }
    
    #SSURGO/STATSGO data. Whole soil erosion and pH
    if(!"soilmu.shp" %in% list.files(ancillary.dir)){
      soilmu <- read_sf(paste0(ancillary.dir, "/wss_gsmsoil_US_[2016-10-13]/wss_gsmsoil_US_[2016-10-13]/spatial/gsmsoilmu_a_us.shp"))
      soilmu <- st_transform(siolmu, crs = 5070)
    
      #component table - map units contain multiple components
      comp <- read.table(paste0(ancillary.dir,"/wss_gsmsoil_US_[2016-10-13]/wss_gsmsoil_US_[2016-10-13]/tabular/comp.txt"), sep = "|")    
      #2 comppct_r RV Integer Smallint No 
      #6 majcompflag Major Component Boolean Varchar No 3
      #31 geomdesc Geomorphic Description Narrative Text No
      #84 taxorder Order Choice Varchar No 254 
      #88 taxpartsize Particle Size Choice No 254 
      #108 mukey Mapunit Key String Varchar Yes 30
      #109 cokey Component Key String Varchar Yes 30
      comp <- comp[,c(2, 31, 84, 88, 108, 109)]
      names(comp) <- c("comppct_r", "geomdesc", "taxorder", "taxpartsize", "mukey", "cokey")
    
      #horizon table - has whole soil erosion factor and soil pH valuse for each horison in component 
      chorizon <- read.table(paste0(ancillary.dir, "/wss_gsmsoil_US_[2016-10-13]/wss_gsmsoil_US_[2016-10-13]/tabular/chorizon.txt"), sep = "|")
      #1 hzname Designation String Varchar No 12
      #112 kwfact Kw Choice Varchar No 254
      #136 ph1to1h2o_r RV Float Real No 1 1.8 11
      #170 cokey Component Key String Varchar Yes 30
      chorizon <- chorizon[,c(1, 112, 136, 170)] 
      names(chorizon) <- c("hzname", "kwfact", "ph1to1h2o_r", "cokey")
    
      #add zeros to NA values and mean horizons (soils are in 3 demensions)
      chorizon[is.na(chorizon)] <- 0
      chorizon <- aggregate(chorizon[, c("kwfact", "ph1to1h2o_r")], 
                            by = list(cokey = chorizon$cokey), function(x) round(mean(x), 2))

      
      soil_attr <- Reduce(function(x, y) merge(x, y, by = "cokey"), list(comp,chorizon))
      soil_attr <- split(soil_attr, soil_attr$mukey)
      
      #use majority rules for categorical variables and weighted averages for continuous
      #removed "geomdesc" was long text field 
      soil_attr <- lapply(soil_attr, function(x) 
        data.frame(mukey = unique(x$mukey), 
                   x[which.max(x$comppct_r), c("taxorder", "taxpartsize")], 
                   kwfact = round(weighted.mean(x$kwfact, x$comppct_r),2),
                   ph1to1h2o_r = round(weighted.mean(x$ph1to1h2o_r, x$comppct_r),2))) 
  
    soil_attr <- do.call(rbind, soil_attr)
    soilmu <- merge(soilmu, soil_attr, by.x = "MUKEY", by.y = "mukey")
    
    write_sf(soilmu, paste0(ancillary.dir, "/soilmu.shp"))
    } else {
      #print("soilmu.shp already exists")
    }
      
    #absolute depth to bedrock    
    if (!"BDTICM_ras.tif" %in% list.files(ancillary.dir)){
      #youll have to check the wd for raster
      ras <- raster(paste0(ancillary.dir, "/BDTICM_M_1km_ll/BDTICM_M_1km_ll.tif"))
      #US extent to clop DTB raster before projection. 
      #May speed things up a little
      US_extent <- extent(c(-124.9022, -66.23839, 24.51224, 50.63472))
      ras <- crop(ras, US_extent)
      
      ras <- projectRaster(ras, crs = newCRS)

      writeRaster(ras, paste0(ancillary.dir, "/BDTICM_ras.tif"))
      
    } else {
      #print("BDTICM_ras.tif already exists")
    }
    
    #identify the max fac within all catchments (gridcodes), used to calculate 
    #difference between sampling point and outlet drainage area
    if(!"max_fac_cat.csv" %in% list.files(raster.dir)){
      cat <- raster(paste0(raster.dir,"/cat.tif"))
      fac <- raster(paste0(raster.dir,"/fac.tif"))
      zone_max <- data.frame(zonal(fac, cat, 'max'))
      names(zone_max) <- c("gridcode", "fac")
      write.csv(zone_max, paste0(raster.dir, "/max_fac_cat.csv"), row.names = F)
    } else {
      print("max_fac_cat.csv already exists")
    }
    
  ###################################
    
  ##### Attribute points
  #########
    
    Reach_Points <- read_sf(paste0(write.dir, "/Reach_Points.shp"))
    
    nhdhr_net <- read_sf(paste0(write.dir,"/NHDHR_Network.shp"), as_tibble = F)
    # merege with vaa with nhdhr network
    vaa <- suppressWarnings(st_read(dsn = grep(paste0(x[site, "HUC4"], "_HU4_GDB.gdb"), list.dirs(), value = T),
                                    layer = "NHDPlusFlowlineVAA", 
                                    quiet = T,  stringsAsFactors = T, as_tibble=F))
    vaa$NHDPlusID <- as.character(vaa$NHDPlusID)
    
    
    #assign vaa to reconditioned network (reach points)
    attr_nhd <- attr_nhdplusHR(s_pts = Reach_Points,
                               nhdhr_net = nhdhr_net, 
                               vaa_tbl = vaa, 
                               vaa_attr = c("StreamLeve", "StreamOrde", 
                                            "StreamCalc", "TotDASqKm"))
    
    adjcat <- adj_cat_area(s_pts = Reach_Points, 
                           cat = raster(paste0(raster.dir,"/cat.tif")), 
                           fac = raster(paste0(raster.dir,"/fac.tif")),
                           max_cat = read.csv(paste0(raster.dir, "/max_fac_cat.csv")))
    
    attr_dtf <- attr_dtf_vaa(s_pts = Reach_Points, 
                             dtf = 4, 
                             dem = raster(paste0(raster.dir,"/elev_cm.tif")),
                             nhd_attr = attr_nhd,
                             cat_adjust = adjcat)
    
    #PRISM Temp and Climate
    PRISMppt30 <- extract(raster(paste0(ancillary.dir, "/ppt_ras.tif")), 
                          sf::st_coordinates(Reach_Points))
    PRISMppt30 <- data.frame(split_fID = Reach_Points$split_fID, PRISMppt30)
    
    PRISMtemp30 <- extract(raster(paste0(ancillary.dir, "/temp_ras.tif")), 
                           sf::st_coordinates(Reach_Points))
    PRISMtemp30 <- data.frame(split_fID = Reach_Points$split_fID, PRISMtemp30)
    
    #extract deth to bedrock 
    dtb_cm <- extract(raster(paste0(ancillary.dir, "/BDTICM_ras.tif")), 
                      sf::st_coordinates(Reach_Points))
    dtb_cm <- data.frame(split_fID = Reach_Points$split_fID, dtb_cm)
    
    #extract ssurgo/statsgo data
    soilmu <- read_sf(paste0(ancillary.dir, "/soilmu.shp"))
    st_crs(soilmu) <- "EPSG:5070"
    
    soilmu_index <- st_within(Reach_Points, soilmu)
    soilmu_index <- lapply(soilmu_index, function(x) ifelse(length(x) == 0, NA, x))
    soilmu_index <- do.call(rbind, soilmu_index)
    
    soilmu_index <- data.frame(split_fID = Reach_Points$split_fID, soilmu_index)
    
    #ssurgo stat is incomplete in some locations - needed to make NA in output
    if(any(!complete.cases(soilmu_index))){
      NoData <- soilmu_index[is.na(soilmu_index$soilmu_index), "split_fID"]
      
      Nodata_soilmu <- st_set_geometry(soilmu[1,], NULL)
      Nodata_soilmu[1,] <- NA
      NoData <- data.frame(split_fID = NoData, Nodata_soilmu)
      
      soilmu_index <- soilmu_index[!is.na(soilmu_index$soilmu_index), ]
      
      soilmu <- data.frame(split_fID = soilmu_index[,"split_fID"], 
                           st_set_geometry(soilmu[soilmu_index[,2], ], NULL))
      soilmu <- rbind(soilmu, NoData)
      
      } else {
        
        soilmu <- data.frame(split_fID = Reach_Points$split_fID,
                           st_set_geometry(soilmu[soilmu_index[,2], ], NULL))
      }
    
    # major geology
    geo <- read_sf(paste0(ancillary.dir, "/geo.shp"))
    reclass_geo <- read.csv(paste0(ancillary.dir, "/geology_reclass.csv"))
    #not reading in crs, need to set
    st_crs(geo) <- "EPSG:5070"
    major <- st_within(Reach_Points, geo)
    major <- do.call(rbind, major) #dealing with sgbc class (list)
    
    #some geology classes are NA, like water, lable with UNIT name
    MjrGeo <- st_set_geometry(geo[major[,1], c("MAJOR1", "UNIT_NAME")], NULL)
    MjrGeo <- merge(MjrGeo, reclass_geo, by = "MAJOR1")
    
    MjrGeo$Reclass <- as.character(MjrGeo$Reclass)
    MjrGeo[is.na(MjrGeo$MAJOR1), "MAJOR1"] <- MjrGeo[is.na(MjrGeo$MAJOR1), "UNIT_NAME"]
    MjrGeo[is.na(MjrGeo$Reclass), "Reclass"] <- MjrGeo[is.na(MjrGeo$Reclass), "UNIT_NAME"]
    
    MjrGeo <- data.frame(split_fID = Reach_Points$split_fID, MjrGeo)
    
    attr_points <- Reduce(function(x,y) merge(x,y, by = "split_fID"), 
                          list(attr_nhd, attr_dtf, PRISMppt30, 
                               PRISMtemp30, MjrGeo,  dtb_cm, soilmu))
    
    ###############################
    write.csv(attr_points, paste0(write.dir, "/attr_Points.csv"), row.names = F)
    
    ##### Attribute Valley Segments 
    ##########
    Valley_Segments <- read_sf(paste0(write.dir, "/Valley_Segments.shp"))
    attr_valley <- dwn_v_slope(v_seg = Valley_Segments, 
                               r = raster(paste0(raster.dir, "/elev_cm.tif")))
   
    ################################
    write.csv(attr_valley, paste0(write.dir, "/attr_valley.csv"), row.names = F)
    
    ##### Attribute Reach Transects
    ##################
    Reach_Transects <- read_sf(paste0(write.dir, "/Reach_Transects.shp"))
    v_side_slope <- valley_s_slope(vw = Reach_Transects,
                                   r = raster(paste0(raster.dir,"/elev_cm.tif")))
    
    v_floor_width <- valleyfloor_w(s_pts = read_sf(paste0(write.dir, "/Reach_Points.shp")), 
                                   dem = raster(paste0(raster.dir,"/elev_cm.tif")), 
                                   v_tran = Reach_Transects, 
                                   attr_Points = read.csv(paste0(write.dir, "/attr_Points.csv"))) 
    
    attr_transects <- merge(v_side_slope,v_floor_width, by = "split_fID")
    # not running channel belts - takes a long time 
    # channelBelts <- read_sf(paste0(write.dir, "channelBelts.shp"))
    # Reach_Transects <- belt_w(Reach_Transects = Reach_Transects, belts = channelBelts)
    
    ################################
    write.csv(attr_transects, paste0(write.dir, "/attr_transects.csv"), row.names = F)
    
    ##### Attribute reach segments
    ############
    Reach_segments <- read_sf(paste0(write.dir, "/Reach_segments.shp"))
    attr_reach <- channel_shp(reach_seg = Reach_segments)

    ################################
    write.csv(attr_reach, paste0(write.dir, "/attr_reach.csv"), row.names = F)
    
    ##### Merge
    ##############
    points <- read.csv(paste0(write.dir, "/attr_Points.csv"))
    transects <- read.csv(paste0(write.dir, "/attr_transects.csv"))
    reach <- read.csv(paste0(write.dir, "/attr_reach.csv"))
    
    out <- Reduce(function(x,y) merge(x,y, by = "split_fID"), list(points,transects,reach))    
    
    valley <- read.csv(paste0(write.dir, "/attr_valley.csv"))
    out <- merge(out, valley, by = "SegID")
    
    ################################    
    write.csv(out, paste0(write.dir, "/Network_HGP.csv"), row.names = F)
}
