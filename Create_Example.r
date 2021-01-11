# StreamCLimes Stream Classification
#create module

setwd("C:/Users/AllenLabWorkstation/Dropbox/Dissertation/StreamClimes/Manuscript/Data")
source("NHDHR_Functions_Example.R")

#networks outlets form StreamClimes 
#outlet needs to be formatted exactly as below
x <- read.csv("NEON_outlet_example.csv", colClasses ="character" )
x <- x[,c("site", "HUC4", "NHDPLUSID")]
names(x) <- c("Network_ID", "HUC4", "NHDPlusID")

#iterate through each each watershed
for (site in 1:nrow(x)) {
  #site <- 1
  print(x[site, ])
  
  # Initialize
  #####
  
  write.dir <- c(paste0(getwd(), "/", x[site, "Network_ID"]))
  dir.create(write.dir)
  
  setwd("./NHDPLUSHR_3_7_2019")
  raster.dir <- paste0("./NHDPLUS_H_", x[site, "HUC4"], "_HU4_RASTER/HRNHDPlusRasters", x[site, "HUC4"])
  ancillary.dir <- "./ancillary_data"
  
  ###################################

  # Create   
  #####

  # Delineate network
  net <- net_delin_HR(huc = as.character(x[site, "HUC4"]), 
                    NHDPlusID =	as.character(x[site, "NHDPlusID"]))
  write_sf(net$net.SF, paste0(write.dir, "/NHDHR_Network.shp"))

  # Segmentize Network
  net <- net_seg_HR(netdelin = net) 
  write_sf(net$net.SF, paste0(write.dir, "/Valley_Segments.shp"))
  write_sf(net$cat.SF, paste0(write.dir, "/Valley_Catchments.shp"))

  # Split Valley
  Reach_segments <- sample_reach(input_lines = net$net.SF, 
                                 max_length = 1000, id = "SegID")
  write_sf(Reach_segments, paste0(write.dir, "/Reach_segments.shp"))

  # Sample Network Points
  Reach_points <- sample_pts(reach = Reach_segments, t_dist = 5000)
  write_sf(Reach_points$Points, paste0(write.dir,"/Reach_Points.shp"))
  #write_sf(Reach_points$Transects, paste0(write.dir,"/Reach_Transects_long.shp"))

  # Valley TRANSECTS
  Reach_Transects <- lapply(split(Reach_points$Transects, 
                                  Reach_points$Transects$split_fID),
                            function(x) valley_w(x = x, 
                                      y = net$cat.SF[net$cat.SF$SegID == x$SegID, ], 
                                      pt = Reach_points$Points[Reach_points$Points$split_fID == x$split_fID, ], 
                                      inner = T))

  Reach_Transects <- do.call(rbind, Reach_Transects)
  write_sf(Reach_Transects, paste0(write.dir, "/Reach_Transects.shp"))

  ##### Not used in analysis (Be patient... 18,000 points took ~4hrs) 
  #belts <- channel_belt(reach_seg = Reach_segments)
  #write_sf(belts, paste0(write.dir, "channelBelts.shp"))
  ###################################
  
}

