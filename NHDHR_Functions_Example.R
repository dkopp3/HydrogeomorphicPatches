# functions for stream network classification

#required packages
library(sf)
library(rgdal)
library(dplyr)
library(raster)
library(gstat)
library(lwgeom)
library(smoothr)
library(units)
library(cluster)
library(dendextend)
library(tidygraph)
library(igraph)
library(randomForest)
library(sf)
library(MuMIn)


#' @title NHDHR network delineation 
#' @description Returns NHDHR comid's upstream of comid, flowlines and catchments form NHD HR
#' @param HUC 4-digit hydrologic unit code.
#' @param NHDPlusID comid of reach
#' @return Returns NHDHR stream network
#' @importFrom add sf and rgdal libraries here
#' @export
#'
net_delin_HR <- function(huc, NHDPlusID){
  
  if(!is.character(NHDPlusID)){
    stop("NHDPlusID is not character")
  }
  
  fileGDB <- grep(paste0(huc,"_HU4_GDB.gdb"), list.dirs(), value = T)
  #supress simple features warning with flow table in .gdb. 
  #I know there is no spatial information 
  flotbl <- suppressWarnings(st_read(dsn = fileGDB, layer = "NHDPlusFlow", quiet = T))
  flotbl[, "FromNHDPID"] <- as.character(flotbl[, "FromNHDPID"])
  flotbl[, "ToNHDPID"] <- as.character(flotbl[, "ToNHDPID"])
  
  network <- data.frame(group_NHDPlusID = character(), 
                        net_NHDPlusID = character(), 
                        huc = character())

  for (i in 1:length(NHDPlusID)) {
    net <- NHDPlusID[i]
    fcomid <- NHDPlusID[i]
    vpus <- as.character(flotbl[flotbl[,"ToNHDPID"] == NHDPlusID[i], "FromVPUID"])
  
    while (length(flotbl[flotbl[, "ToNHDPID"] %in% fcomid, "FromNHDPID"]) >= 
           1 & all(!is.na(net))) {
      fcomid <- flotbl[flotbl[, "ToNHDPID"] %in% fcomid, "FromNHDPID"]
      fcomid <- fcomid[fcomid != 0]
      vpu <- as.character(flotbl[flotbl[, "ToNHDPID"] %in% fcomid, "FromVPUID"])
    
      net <- append(fcomid, net, length(net))
      vpus <- append(vpu, vpus, length(net))
    
      #check to see if network jumps vpu
      if(length(unique(vpus[!is.na(vpus)]))>1){
        stop(print(unique(vpus)))
      }
    }
    
    group.comid <- NHDPlusID[i]
    net.comid <- unique(net[order(net)])
    #M <- ifelse(group.comid == net.comid, Ms[i], 1)
    network <- rbind(network, 
                     data.frame(NHDPlusID = NHDPlusID[i], 
                                net_NHDPlusID = as.character(net.comid), 
                                huc = huc))
  }
  
  floln <- st_read(dsn = fileGDB, layer = "NHDFlowline", quiet = T)
  floln$NHDPlusID <- as.character(floln$NHDPlusID)
  net.shp <- floln[as.character(floln$NHDPlusID) %in% as.character(network$net_NHDPlusID), ]
  GNIS_Name <- floln[as.character(floln$NHDPlusID) %in% as.character(NHDPlusID[i]), "GNIS_Name"]
  st_geometry(GNIS_Name) <- NULL
  
  cats <- st_read(dsn = fileGDB, layer = "NHDPlusCatchment", quiet = T)
  cats$NHDPlusID <- as.character(cats$NHDPlusID)
  cat.shp<-cats[as.character(cats$NHDPlusID) %in% as.character(network$net_NHDPlusID),]

  out<-list(name = GNIS_Name, 
            network = network, 
            net.SF = net.shp, 
            cat.SF = cat.shp)
  return(out)
}


#' @title Identify valley segments form NHD HR 
#' @description identify valley segments NHDHR form netdelin object 
#' @param netdelin netdelin object
#' @param NHDHR.Dir directory of NHDHR data
#' @return Returns NHDHR netdelin object with SegID field
#' @importFrom add sf and rgdal libraries here
#' @export
#'
net_seg_HR <- function (netdelin = z, NHDHR.Dir = getwd()){
  huc <- levels(netdelin$network[, "huc"])
  
  if(length(huc) > 1){
    stop("cross HUC boundary")
  }
  
  oldwd <- getwd() # save old directory
  setwd(NHDHR.Dir)
  fileGDB <- grep(paste0(huc, "_HU4_GDB.gdb"), list.dirs(), value = T)
  
  if(length(fileGDB) != 1){
    stop("check NHDHR.dir directory is correct")
  }
  
  #suppress spatial warnings bc its just a table
  vaa.tbl <- suppressWarnings(st_read(dsn = fileGDB, layer = "NHDPlusFlowlineVAA", quiet = T))
  #drop divergent flowpaths
  vaa.tbl <- vaa.tbl[vaa.tbl[, "StreamOrde"] == vaa.tbl[, "StreamCalc"], ]
  vaa.tbl <- vaa.tbl[,c("NHDPlusID","FromNode", "ToNode", "HydroSeq", "LevelPathI")]
  #change values to character for matching
  vaa.tbl <- apply(vaa.tbl, 2, as.character)
  
  netdelin$network <- apply(netdelin$network, 2, as.character)
  
  #merging reduces the number of VAA records
  vaa.tbl <- merge(netdelin$network, vaa.tbl, 
                   by.x = "net_NHDPlusID", by.y = "NHDPlusID")

  #split flowpaths by the LevelPathID
  net.path <- split(vaa.tbl, as.factor(as.character(vaa.tbl[,"LevelPathI"])))
  #net.path[1]
  #count the number of upstream segments 
  #match Fnode to Tonode and return the Fnodes
  #then group by toNodes and count the number of upstream

  #for each pathlevel count the number of nodes upstream of each comid node (>two indicate confluence) 
  node.count <- lapply(net.path, function (i) 
    vaa.tbl[vaa.tbl[,"ToNode"] %in% i[, c("FromNode")], c("ToNode", "FromNode")])
  single.path <- names(node.count[which(lapply(node.count, function(i) dim(i)[1]) == 0)])
  node.count <- node.count[which(lapply(node.count,function(i) dim(i)[1]) != 0)]
  node.count <- lapply(node.count, function(i) {
    aggregate(i[,"FromNode"], list(FromNode = i[,"ToNode"]), length)
  })

  #these are paths that have multiple comids
  multi.paths <- net.path[names(net.path) %in% names(node.count)]

  node.count.df <- do.call(rbind.data.frame, node.count)
  rownames(node.count.df) <- NULL
  
  func <- function(x, y){merge(x, y, by = "FromNode", all.x = T)}
  multi.paths <- lapply(multi.paths, func, node.count.df)
  
  #need to add headwaters that are part of a path
  new.x <- lapply(multi.paths, function(x) replace(x["x"], is.na(x["x"]), 2))
  multi.paths<-mapply(function(x,y) "[<-"(x, "x", value = y),
                        multi.paths, new.x, SIMPLIFY = FALSE)
  #sort dataframes within list 
  multi.paths <- lapply(multi.paths, function(x) x[order(x$HydroSeq), ])

  #add id's to multi.paths 
  IDs <- lapply(multi.paths, function(z) paste0(z[,"LevelPathI"], "_", seq(1, dim(z)[1])))
  names(IDs) <- NULL
  multi.paths <- mapply(function(x, y) "[<-"(x, "SegID", value = y),
                        multi.paths, IDs, SIMPLIFY = FALSE)
  
  #adjust the segment ID to group valley segments
  index.seq <- lapply(multi.paths, function(z) which(z[,"x"] == 1))
  index.seq <- index.seq[lapply(index.seq, function(x) length(x)) != 0]
  
  #start index sequence for valley segment merge
  starts <- lapply(index.seq, function(n) n[c(TRUE, diff(n) != 1)])
  #end index of sequence for valley segment merge
  ends <- lapply(index.seq, function (n) n[c(diff(n) != 1, TRUE)] + 1)
  
  #combine into dataframe
  start_stop <- mapply(function(x, y) data.frame(start = x, stop = y),
                       starts, ends, SIMPLIFY = FALSE)

  adj.id <- lapply(start_stop, function(a) 
    lapply(apply(a, 1, function(x) 
      list(seq(from = x[1], to = x[2]))), "[[", 1))
  
  #rename to index
  adj.id <- mapply(function(x, y) setNames(x, y), adj.id, starts, SIMPLIFY = F)
  
  
  adj.id <- lapply(adj.id, function(x) data.frame(ID = rep(names(x), sapply(x, length)),
                                                  index = unlist(x)))
  
  #Adjust names to valley group comids 
  for (i in names(multi.paths)){
    if(names(multi.paths[i]) %in% names(adj.id)){
      multi.paths[[i]][adj.id[[i]][,"index"], "SegID"] <- 
        paste0(as.character(multi.paths[[i]][adj.id[[i]][,"index"], "LevelPathI"]), 
               "_", as.character(adj.id[[i]][,"ID"]))
    }
  }
  
  #DISSOLVE LIST PATHS  
  multi.paths <- do.call(rbind.data.frame, multi.paths)
  rownames(multi.paths) <- NULL
  multi.paths <- multi.paths[,c("net_NHDPlusID", "NHDPlusID", "SegID")]
  
  #ADD IN THE SINGLE PATHS
  single.path <- data.frame(vaa.tbl[vaa.tbl[,"LevelPathI"] %in% single.path,
                                    c("net_NHDPlusID", "NHDPlusID")], 
                            SegID = paste0(single.path, "_1"))
  out <- rbind(multi.paths,single.path)
  out <- out[complete.cases(out),]
  
  wo_cat <- setdiff(out$net_NHDPlusID, netdelin$cat.SF$NHDPlusID)
  segid_wo_cat <- out[out$net_NHDPlusID %in% wo_cat, "SegID"]
  
  #these valley segs contain reaches that do not have catchemnt areas
  #remove them, often associated with ditches and interbasin transfers... 
  out <- out[!out$SegID %in% segid_wo_cat,]
  
  #update netdelin input
  netdelin$network <- merge(netdelin$network, out, by = c("net_NHDPlusID", "NHDPlusID"), all.x = T)
  netdelin$net.SF <- merge(netdelin$net.SF, out, by.x = "NHDPlusID", by.y = "net_NHDPlusID")
  #plot(st_geometry(st_zm(netdelin$net.SF, drop = T, what = "ZM")))
  netdelin$cat.SF <- merge(netdelin$cat.SF, out, by.x = "NHDPlusID", by.y = "net_NHDPlusID")
  
  #dissolve the line segments before running through the sample points- perhaps consider makeing this a functuion. 
  #the group_by summarize creates multilinestrings. st_cast expands the number of records making the dissolve useless 
  #st_line_merge works but only for mutliline geometries
  #note st_geometry_type and st_line_merge
  test_dissolve <- netdelin$net.SF %>% 
    st_zm(drop = T, what = "ZM") %>%
    group_by(SegID) %>% 
    summarise(v_length_KM = sum(LengthKM)) 
  
  lnstr <- test_dissolve[st_geometry_type(test_dissolve) != "MULTILINESTRING",]
  
  #some instances there are not multiline strings from test disolve. 
  if(any(st_geometry_type(test_dissolve) == "MULTILINESTRING")){
    mlnstr <- st_line_merge(test_dissolve[st_geometry_type(test_dissolve) == "MULTILINESTRING",])
    
    #gonna force the multilinestrings into line string
    #instance with LIPI where st_line_merge would not merge flowline
    #mlnstr[43,] as line string
    if(any(st_geometry_type(mlnstr) == "MULTILINESTRING")){
      forceln <- mlnstr[st_geometry_type(mlnstr) == "MULTILINESTRING", ]
      
      #function to force linestring 
      forceToLine <- function(x){ 
        crd <- st_coordinates(x)
        geometry <- st_sfc(st_linestring(matrix(crd[, c(1,2)], , 2)))
        sf <- st_as_sf(data.frame(st_set_geometry(x, NULL)), geometry, crs = st_crs(x))
        return(sf)
        }
  
      forceln <- lapply(split(forceln, forceln$SegID), function(x) forceToLine(x))
      forceln <- do.call(rbind, forceln)
      mutate(forceln, v_length_KM = set_units(st_length(forceln)/1000, NULL))
      mlnstr <- mlnstr[st_geometry_type(mlnstr) != "MULTILINESTRING", ]
      mlnstr <- rbind(mlnstr, forceln)
      
      warning(paste(forceln$SegID), " forced to LINESTRING", sep = " ")
    }
    
    test_dissolve <- rbind(lnstr, mlnstr)
    test_dissolve <- st_transform(test_dissolve, crs = 5070)
    } else {
      test_dissolve <- st_transform(lnstr, crs = 5070)  
    }
  
  netdelin$net.SF <- test_dissolve
  
  #dissolve catchments
  test_dissolvecat <- netdelin$cat.SF %>% 
    st_zm(drop = T, what = "ZM") %>%
    group_by(SegID) %>% 
    summarise(AreaSqKm = sum(AreaSqKm)) %>%
    st_transform(crs = 5070)
  netdelin$cat.SF<-test_dissolvecat
  
  return(netdelin)
}


#' @title down valley attributes 
#' @description calculate down valley slope, sinuosity 
#' @param v_segs dissolved reiver segments (confluence to confluence)
#' @param r is elevation raster cm 
#' @return Returns slope and sinuosity attributes appended to v_segs input
#' @importFrom add sf and rgdal libraries here
#' @export
#'
dwn_v_slope <- function(v_seg, r){

  v_coords.F <- data.frame(st_coordinates(v_seg)) %>%
  group_by(L1) %>%
  mutate(pt_seq = 1:n()) %>%
  group_by(L1) %>%
  top_n(-1, pt_seq) %>%
  ungroup()

# sequence (ppt sequence increases) such that the top 
# value is actually the bottom of the valley
v_coords.T <- data.frame(st_coordinates(v_seg)) %>%
  group_by(L1) %>%
  mutate(pt_seq = 1:n()) %>%
  group_by(L1) %>%
  top_n(c(1), pt_seq) %>%
  ungroup()

v_pts <- rbind(v_coords.F, v_coords.T)%>%
  arrange(L1)%>%
  sf::st_as_sf(coords = c("X","Y")) %>% 
  sf::st_set_crs(5070)%>%
  raster::extract(r, .)/100

rise <- sapply(seq(2, length(v_pts), by = 2), function(x) v_pts[x-1] - v_pts[x])
#change any negative values to zero - this is why we use hydroDEMs
#may cause issues down the road but seems alright for now
rise <- ifelse(rise < 0, 0, rise)

v_seg <- mutate(v_seg, 
                #v_length_m = round(st_length(v_seg), 3), 
                rise_m = rise, 
                v_slope_precent = round(rise/st_length(v_seg) * 100, 3))


#combine these into line - valley length
vl <- rbind(v_coords.F, v_coords.T)%>%
  arrange(L1)%>% 
  sf::st_as_sf(coords = c("X","Y")) %>% 
  sf::st_set_crs(5070)%>% 
  group_by(L1) %>% 
  summarize(pt_seq = max(pt_seq)) %>% 
  st_cast("LINESTRING")%>%
  mutate(v_length_m = round(st_length(.),3))

v_seg <- mutate(v_seg, 
                V_sinu = round(st_length(v_seg)/vl$v_length_m, 3))
v_seg <- data.frame(st_set_geometry(v_seg, NULL))
return(v_seg)  
}


#' @title Valley Floor Delineation 
#' @description derive floodplain using IDW 
#' @param ss stream skeleton raster
#' @param dem digital elevation model
#' @param fac flow accumulation grid
#' @param ws watershed boundary polygon
#' @param dtf depth to flood value (factor multiplied by bankfull height)
#' @return Returns sf object of floodplain
#' @importFrom add sf and rgdal libraries here
#' @export
#'
idw_fpln <- function(ss, dem, fac, ws, dtf){
  ss_ras <- clipraster(ss, ws) #clip stream skeleton by ws
  ss_pts <- st_as_sf(
    rasterToPoints(ss_ras, function(x) x == 1, T), crs = 5070) #create stream points
  ss_dem <- clipraster(dem, ws) #clip dem ws
  ss_ele <- extract(dem, ss_pts) #extract elevation values
  ss_fac <- extract(fac, ss_pts) %>% 
    sapply(function(x) (x * 10^2)/1000^2) #convert pixels to square km area
  
  #hydrologic geometry relationships form Bieger et al 2015
  #used parameter estimates for entire US rgion may 
  #consider adjusting for phisiographic region
  #added value to x to accound for FAC value starting at 0
  ss_bfh <- sapply(ss_fac, function(x) (0.30 *(x + 0.00001) ^ 0.213))  
  ss_dtf <- ss_bfh * dtf #multiply bankfull height by dft value 
  #multiplied dtf by 100 to convert from meters to centimeters
  #and match units of dem
  ss_idw <- (ss_dtf * 100) + ss_ele 
  ss_pts <- mutate(ss_pts, ss_fac = ss_fac, ss_ele = ss_ele, 
                   ss_bfh = ss_bfh, ss_dtf = ss_dtf, ss_idw = ss_idw)
  ss_pts <- data.frame(st_set_geometry(ss_pts, NULL), st_coordinates(ss_pts))
  
  #could do a sensitivity analysis on the nmax and idp values
  #arcgis has different default values for idw 
  idw.form <- gstat(id = "ss_idw", formula = ss_idw ~ 1, locations = ~X + Y, data = ss_pts, 
                    nmax = 7, set = list(idp = 0.5)) 
  idw <- interpolate(ss_ras, idw.form, xyNames = c("X","Y")) #idw flood elevation surface 
  
  #could maybe get some improvement if i do not convert to polygon 
  #use get Values instead 
  fp <- suppressWarnings(overlay(ss_dem, idw, fun = function(x,y) x <= y) %>% #select dem cells that are <= idw as flood plain
                           rasterToPolygons(., fun = function(x) x == 1, dissolve = T) %>% #converte to polygon
                           st_as_sf() %>% #convert to sf
                           st_cast("POLYGON")) #cast to polygon
  
  #remove flood plain occuring on ridge lines (i.e. not accosiated with a stream skeletion pixel)  
  a <- extract(ss_ras, fp)%>% #identify values associated with the polygon
    lapply(., function(x) sum(x, na.rm = T) >= 1)%>% #select those that have a stream pixel
    unlist() #unlist the lapply for T/F vector
  
  fp <- fp[a, ]
  fp <- st_as_sf(data.frame(GNIS_Name = as.character(ws$GNIS_Name), st_sfc(st_combine(fp))))
    
  return(fp)
}


#' @title Fill holes in flood plain resulting from small differences in topogrpy 
#' @description make flood plain look pretty with functions from the smoothr package
#' @param fpl flood plain object see idw_fpln
#' @param minarea minimum area filled in m^2 (1pixel  = 10m x 10m)
#' @return Returns sf object of filled floodplain
#' @importFrom add sf and rgdal libraries here
#' @export
#'
fillholes <- function(fpl, minarea = 1000){
  r_poly_dropped <- drop_crumbs(fpl, threshold = units::set_units(1000, m^2))%>%
    fill_holes(., units::set_units(maxarea, m^2))
  return(r_poly_dropped)
}
 
 
#' @title create sample reaches
#' @description splits valleys greater than max distance into equal lengths. appends smaller valleys
#' Modified form David Blodgett's "Split_lines" function
#' @param input_lines data.frame of class sf with LINESTRING sfc column.
#' @param max_length maximum reach length. all sample reaches will be less than this length
#' @param id name of ID column in input_lines
#' @return Sampling points along LINESTRING.
#' @importFrom dplyr group_by ungroup filter select mutate
#' @export
#'
sample_reach <- function(input_lines, max_length = 100, id = "ID") {
  
  if(grep("units", unlist(strsplit(st_crs(input_lines)$proj4string,"\\+")), value = T)!="units=m "){
    stop("use projected crs with units=m")
  }
  
  #identify geometry column in sf object 
  geom_column <- attr(input_lines, "sf_column")
  input_crs <- sf::st_crs(input_lines)
  
  #measure lengths of all input lines and set as new field 
  input_lines[["geom_len"]] <- sf::st_length(input_lines[[geom_column]])
  attr(input_lines[["geom_len"]], "units") <- NULL
  input_lines[["geom_len"]] <- as.numeric(input_lines[["geom_len"]])
  
  #drop all lines that are less that the min_length dist. 
  #sampling points will be constructed only for lines greater than this value 
  #those under the length should be bisected (one transect in the center)
  too_long <- filter(dplyr::select(input_lines, id, geom_column, geom_len), geom_len >= max_length)
  too_short <- filter(dplyr::select(input_lines, id, geom_column, geom_len), geom_len < max_length)
  #rm(input_lines) # just to control memory usage in case this is big.
  
  #adds new variables to "too_long"
  #pieces is the split by sampling frequency
  #fid is a unique id
  #then select all fields except -geom_len
  too_long <- mutate(too_long,
                     pieces = ceiling(geom_len/max_length),
                     fID = 1:nrow(too_long)) %>%
    dplyr::select(-geom_len)
  
  #drop the geometry, repeat row indicies based on pieces
  #essentially creates new id for each piece
  split_points <- sf::st_set_geometry(too_long, NULL)[rep(seq_len(nrow(too_long)), too_long[["pieces"]]),] %>%
    dplyr::select(-pieces)
  
  #add start and stop proportions to line
  split_points <- mutate(split_points, split_fID = row.names(split_points)) %>%
    group_by(fID) %>%
    mutate(piece = 1:n()) %>%
    mutate(start = (piece - 1) / n(),
           end = piece / n()) %>%
    ungroup()
  
  new_line <- function(i, f, t) {
    lwgeom::st_linesubstring(x = too_long[[geom_column]][i], from = f, to = t)[[1]]
  }
  
  split_lines <- apply(split_points[c("fID", "start", "end")], 1,
                       function(x) new_line(i = x[["fID"]], 
                                            f = x[["start"]], 
                                            t = x[["end"]]))
  
  split_lines <- st_sf(split_points[c(id, "split_fID")], 
                       geometry = st_sfc(split_lines, crs = input_crs))%>%
    mutate(geom_len = st_length(.))
  split_lines <- mutate(too_short, split_fID = seq((nrow(split_lines) + 1), 
                                       nrow(split_lines)+(nrow(too_short)), by = 1)) %>%
    rbind(.,split_lines)%>%
    arrange(split_fID)
  return(split_lines)
}


#' @title sample points and transects at midpoints of the sample reach
#' @description create sample points and transects at midpoints
#' @param reach lines to split at midpoint and to draw transects
#' @param t_dist length of lateral transect 
#' @return Sampling points and transects at midpoint of LINESTRING input.
#' @importFrom dplyr group_by ungroup filter select mutate
#' @export
#'
sample_pts<- function(reach, t_dist){   
  # take midpoints of line
  st_line_midpoints <- function(sf_lines = NULL) {
    
    g <- st_geometry(sf_lines)
    
    g_mids <- lapply(g, function(x) {
      
      coords <- as.matrix(x)
      
      # this is just a copypaste of View(maptools:::getMidpoints):
      get_mids <- function (coords) {
        dist <- sqrt((diff(coords[, 1])^2 + (diff(coords[, 2]))^2))
        dist_mid <- sum(dist)/2
        dist_cum <- c(0, cumsum(dist))
        end_index <- which(dist_cum > dist_mid)[1]
        start_index <- end_index - 1
        start <- coords[start_index, ]
        end <- coords[end_index, ]
        dist_remaining <- dist_mid - dist_cum[start_index]
        mid <- start + (end - start) * (dist_remaining/dist[start_index])
        return(mid)
      }
      
      mids <- st_point(get_mids(coords))
    })
    
    out <- st_sfc(g_mids, crs = st_crs(sf_lines))
    out <- st_sf(out)
    return(out)
  }
  
  midpts<-st_line_midpoints(reach)%>%
    mutate(SegID = reach$SegID, split_fID = reach$split_fID)
  
  #modification to extract points form end of split line
  #take the end points form the split lines as transect points
  coords <- data.frame(st_coordinates(midpts)) 
  
  #set transects perpendicular to line formed by the sampling point and one vertex above
  #start

  stpoint <- data.frame(st_coordinates(reach)) %>%
    group_by(L1) %>%
    mutate(pt_seq = 1:n()) %>%
    group_by(L1) %>%
    top_n(-1, pt_seq)
  #end
  endpoint <- data.frame(st_coordinates(reach)) %>%
    group_by(L1) %>%
    mutate(pt_seq = 1:n()) %>%
    group_by(L1) %>%
    top_n(1, pt_seq)
  
  tcoords <- rbind(stpoint,endpoint)%>%
    arrange(L1)
  
  #plot(st_as_sf(tcoords, coords = c("X","Y")), add = T, cex = 2)
  #calculate distance vector of the two sampling points, substract them
  vecs <- lapply(seq(2,nrow(tcoords), by = 2), function(x) 
    ((tcoords[(x), c("X","Y")] - tcoords[x-1, c("X","Y")]))) 
  
  #Opposite (perpendicular) distance vector
  vecs <- lapply(vecs, function(x) matrix(c(x[2], (-1 * x[1]))))
  
  #normalize distance vector by unit length (pythagoras theorm - each xy distance form orgin)
  vecs <- lapply(vecs, function(x) 
    matrix(c(x[[1]]/sqrt(x[[1]]^2 + x[[2]]^2), 
             x[[2]]/(sqrt(x[[1]]^2 + x[[2]]^2)))))
  
  #use point coordinate from coords and normalized distance vector to calculate point along line
  #follows form vector line equation
  #r (point vector) = a(starting point) + t_dist (distance along line)*normalized distance vector
  t_from <- data.frame(matrix(unlist(lapply(1:nrow(coords),            
             function (x) t(c(as.numeric(coords[x,"X"]), 
                              as.numeric(coords[x,"Y"])) + 
                              (t_dist * as.numeric(vecs[[x]]))))), ,2, byrow = T), 
             SegID = reach$SegID, split_fID = reach$split_fID)
  
  colnames(t_from)[c(1,2)] <- c("X", "Y")
  
  t_to <- data.frame(matrix(unlist(
    lapply(1:dim(coords)[1], 
           function (x) 
             t(c(as.numeric(coords[x,"X"]),
                 as.numeric(coords[x,"Y"])) - 
                 (t_dist * as.numeric(vecs[[x]]))))), ,2, byrow = T), 
     SegID = reach$SegID, split_fID = reach$split_fID)
  colnames(t_to)[c(1,2)] <- c("X", "Y")
  #make line form coordinates
  Lcoords <- rbind(t_from, t_to)%>%
    arrange(split_fID)
  
  #create transect line string
  t_lines<-st_as_sf(data.frame(Lcoords), coords = c("X","Y"), crs = 5070) %>%
    group_by(split_fID) %>% 
    summarize(SegID = unique(SegID)) %>% 
    st_cast("LINESTRING") 
    
  return(list(Points = midpts, Transects = t_lines))
  
}


#modified valley width function to include closest or furthese intersetion
#runs intersection between transect and catchment 
#which can create a multilinestring where it intersects feature 
#multiple times. fuuntion returns a continuous line (split at sampling point) 
#from the point where the transect 
#first crosses the catchment and the last time 
#this solves with issues with multipolygon features

#' @title Valley Widths
#' @description calculates widths for valley or valley floor widths based on intersection between line and polygon
#' @param x transects
#' @param y is polygon 
#' @param pt is point  
#' @param inner logical TRUE calculates distance to first intersection
#' @return linestrings representing valley widths for each sampling point 
#' @importFrom dplyr group_by ungroup filter select mutate
#' @export
#'
valley_w <- function(x, y, pt, inner = T) {
  #snap point to line incase any small deviations 
  #caused by computer rounding
  snappt <- st_snap(pt, x, tolerance = 1000)
  
  #insert snappt as vertex because st_split works to fins the closest vertex to list split
  x <- rbind(st_coordinates(x)[1, c(1,2)], st_coordinates(snappt),st_coordinates(x)[2, c(1,2)])%>%
    data.frame(split_fID = x$split_fID, SegID = x$SegID,.)%>%
    st_as_sf(coords = c("X","Y"), crs = 5070)%>%
    group_by(split_fID) %>% 
    summarize(SegID = unique(SegID)) %>% 
    st_cast("LINESTRING") 
  
  #intersect so each transect is separate in list
  int <- suppressWarnings(lwgeom::st_split(x, st_geometry(snappt)) %>%
                            st_collection_extract("LINESTRING") %>% 
                            st_intersection(., st_make_valid(y)) %>% 
                            split(.,rownames(.)) %>%
                            lapply(.,function(x) st_cast(x, "LINESTRING")))
  
  if (inner){ #determines which function to use to select coords form 
    #the minimum function assumes that if there 
    #are zeros then the ooccur at either first or last position
    #this should be true because point splits the transect
    #at either the start or end
    pt_fun <- function(x) {
      if(any(as.numeric(x)==0)){
        ifelse(which.min(x) == 1, (which.min(x)+1),(which.min(x)-1))
      } else {
        which.min(x)
      }
    }
  }else {
    pt_fun <- function(x) which.max(x)
  }
  
  if(length(int) > 0){
    ids <- names(int)
    #extract coordinates
    int <- lapply(int, function (x) st_coordinates(x)) %>%  
      lapply(., function (x) rbind(x[pt_fun(st_distance(st_as_sf(data.frame(x), 
                 coords = c("X","Y"), crs = 5070), snappt)), c(1,2)], st_coordinates(snappt)[,c(1,2)])) %>%
      lapply(., function(x) st_as_sf(data.frame(x), coords = c("X","Y"), crs = 5070)) %>%
      do.call(rbind,.) %>%
      data.frame(ID = rep(ids,each = 2), 
                 split_fID = rep(pt$split_fID,length(id)*2), 
                 SegID = rep(pt$SegID, length(id)*2), .) %>%
      st_as_sf() %>% 
      group_by(ID) %>% 
      summarize(split_fID = unique(split_fID), SegID = unique(SegID)) %>% 
      st_cast("LINESTRING") %>% 
      mutate(v_wid = round(st_length(.),3))
  } 
  else {
    int <- NULL
  }
  return(int)
}

#' @title Valley Side Slope
#' @description calculates valley sode from valley side slope transects 
#' @param vw valley width transects see valley_w
#' @param r elevation raster
#' @return linestrings valley widths with side slope for each transect  
#' @importFrom dplyr group_by ungroup filter select mutate
#' @export
#'
valley_s_slope <- function (vw, r){
  elev <- extract(r, st_coordinates(vw)[,-3])
  #verbose below: converte elevation to cm (divide by 100) 
  #then multply by 100 as percent   
  slope <- sapply(seq(2, length(elev), by = 2), 
                  function(k) max(elev[k], elev[k-1]) - 
                  min(elev[k], elev[k-1])) / 100 / vw$v_wid  
 
  vw <- mutate(vw, slope2 =  round(slope*100, 4))
  vw <- data.frame(st_set_geometry(vw, NULL)) %>% 
    dplyr::select(-"SegID")
  vw <- reshape2::dcast(split_fID ~ ID, data = vw, value.var = "slope2")
  names(vw)<-c("split_fID", "side_slope_1", "side_slope_2")
  vw$mean_slope <- (vw$side_slope_1+vw$side_slope_2)/2
  return(vw)
}
  
#' @title Remove spillover in basins - ensures valley floor does not exceed valley width
#' @description  constrains flood plain area by removeing fpl polygons that do not contain a stream pixel) also associates wach flp polygon with valley seg
#' @param flp_int result of arcgis intersection between flp and cat. In future will try to incorporate this step here
#' @param v_seg valley segments
#' @return floodplain polygons that are completely contained by catchments  
#' @importFrom dplyr group_by ungroup filter select mutate
#' @export
#'
flp_cat_trim <- function(fpl_int, v_seg){
  fpl_int <- suppressWarnings(st_cast(fpl_int,"POLYGON")) #cast multipart polygons
  fpl_int <- split(fpl_int, fpl_int$SegID) #split individual
  
  v_seg_pts <- suppressWarnings(st_cast(v_seg, "POINT")) #create points form valley lines 
  #determine whether point is contained by a polygon 
  s <- lapply(fpl_int, function (x) 
    st_contains(x, v_seg_pts[v_seg_pts$SegID == unique(x$SegID),], sparse = F))
  #change to logical to use as index
  s <- lapply(s, function(x) apply(x, 1, function(z) any(z==T)))

  #resolve split
  fpl_int <- do.call(rbind, fpl_int)
  s <- setNames(unlist(s), NULL)
  #select tru index 
  fpl_int <- fpl_int[s,]

  #calculate new area 
  fpl_int$AreaSqm <- round(st_area(fpl_int),3)

  return (fpl_int)
}

#' @title Create Channel belts
#' @description convex hull connects meander apex, with endpoints
#' @param reach_seg reach segments see sample reach 
#' @return channel belts that fully contain reach setments  
#' @importFrom dplyr group_by ungroup filter select mutate
#' @export
#'
channel_belt<-function(reach_seg){
  
  #order by split id for matching later on
  reach_seg <- reach_seg %>%
    arrange(split_fID)
  
  #select first and last coordinates to 
  #create centerline form line segment
  #input to function is a linestring
  centerline <- function(x) {
    coord <- st_coordinates(x)
    coord <- coord[c(1, nrow(coord)),]
    coord <- st_as_sf(data.frame(SegID= x$SegID, split_fID=x$split_fID, coord), coords = c("X","Y"))%>% 
      group_by(split_fID)%>%
      summarise(SegID = unique(SegID))%>%
      st_cast("LINESTRING")%>%
      mutate(CntLn_Len_m = st_length(.))%>%
      st_set_crs(5070)
  }
  
  #calculate centerlines for all reach segments
  #could be used later on to calculate reach sinuosity
  cntr_ln <- lapply(split(reach_seg,reach_seg$split_fID), function(x) centerline(x))%>%
    do.call(rbind,.)
  
  #off diag eleements are meander length (distance between successive points)
  #may want to count meanders later? 
  #input is a matrix
  offdiag <- function(x){
    x[row(x) == (col(x) - 1)]
  }
  
  #reach_seg that do not intersect are just center line lengths
  #every reach segment (that is not completely straight line)
  #had meander length
  meander_lengths <- suppressWarnings(
    lapply(split(cntr_ln, cntr_ln$split_fID), 
           function(x) st_intersection(x, reach_seg[reach_seg$split_fID == x$split_fID,])%>%
             st_cast("POINT")%>%st_distance()%>%
             offdiag()))
  
  #anything reach segment with >1 meander length intersects the centerline 
  ind <- lapply(meander_lengths,function(x) length(x) > 1) %>%
    unlist()
  
  #meander segments (i.e. portions of reach that intersect the center line) 
  #needed to find apex points 
  m_segs <- lapply(split(reach_seg[ind,], reach_seg$split_fID[ind]), 
                   function(x) 
                     st_split(st_geometry(x), st_geometry(cntr_ln[cntr_ln$split_fID == x$split_fID,])) %>%
                     st_collection_extract("LINESTRING")%>%
                     st_sfc(.)%>%
                     st_as_sf(data.frame(split_fID = x$split_fID, SegID = x$SegID),.))
  
  #coordinates of m_segs 
  apex <- lapply(m_segs, function(x)
    st_coordinates(x) %>% 
      data.frame() %>% 
      st_as_sf(.,coords = c("X","Y"), crs = 5070)%>%
      mutate(split_fID = unique(x$split_fID), SegID = unique(x$SegID)))
  
  #added reach segment points that do not intersect centerline
  #so that channel belts are created for reaches that only have one meander 
  #(i.e. do not intersect centerline)
  apex <- c(apex, suppressWarnings(lapply(split(reach_seg[ind==F,], reach_seg$split_fID[ind==F]), 
                                          function(x) st_cast(x, "POINT") %>% mutate(L1 = 1))))
  
  #identifies apex points that occur on the same side of a line
  #ln is the center line and points are points
  #output appends 1 or -1 to sites on the same side of line
  line_side <- function(ln, pts){
    #line (x1,y1) to (x2,y2)  and p is (x,y)
    ln_coords<-st_coordinates(ln)
    pt_coords<-st_coordinates(pts)
    x = pt_coords[,1];y=pt_coords[,2]
    x1 = ln_coords[1,1]; y1 = ln_coords[1,2]
    x2 = ln_coords[2,1]; y2 = ln_coords[2,2]
    #if values have same sign they are on same side of line
    #i really need to learn linear algebra and vectors
    side<-sign((x-x1)*(y2-y1)-(y-y1)*(x2-x1))
    
    return(mutate(pts, side = side))
  }
  
  #interested in finding the highest point associated with 
  #each meander on each side of the center line
  apex_pts <- lapply(apex, function(x)
    #split each line into a line select the max amplitude (i.e. furthest distance frome centerline)
    lapply(split(x, x$L1), function(y)  
      y[which.max(st_distance(cntr_ln[cntr_ln$split_fID == as.character(unique(y$split_fID)),], y)),]))%>%
    lapply(., function(x) do.call(rbind,x)) %>% #roll up sublists to split_fID with lapply do.call 
    #identify meander hights that are on same side of centerline
    lapply(., function(x) 
      line_side(ln = cntr_ln[cntr_ln$split_fID == as.character(unique(x$split_fID)),], pts = x))
  
  
  #draw line that is perpendicular (or parallel) 
  #to a straight line (ln) and passes through a point
  #function used to build tops and sides around reach segments 
  funline <- function(pt, ln, perpendicular = F){
    coords <- data.frame(st_coordinates(pt)) 
    
    #set transects perpendicular to line formed by the centerline
    tcoords <- rbind(st_coordinates(st_startpoint(ln)), 
                     st_coordinates(st_endpoint(ln)))
    
    #plot(st_as_sf(tcoords, coords = c("X","Y")), add = T, cex = 2)
    #calculate distance vector of the two sampling points, substract them
    vecs <- lapply(seq(2,nrow(tcoords), by = 2), function(x) 
      ((tcoords[(x), c("X","Y")] - tcoords[x-1, c("X","Y")])))%>% 
      #Opposite (perpendicular) distance vector #assiming parellel lines follow same distance vector
      ifelse(perpendicular,
             lapply(., function(x) matrix(c(x[2], (-1 * x[1])))),.) %>%
      #normalize distance vector by unit length (pythagoras theorm - each xy distance form orgin)
      lapply(., function(x) 
        matrix(c(x[[1]]/sqrt(x[[1]]^2 + x[[2]]^2), 
                 x[[2]]/(sqrt(x[[1]]^2 + x[[2]]^2)))))
    
    #use point coordinate from coords and normalized distance vector to calculate point along line
    #follows form vector line equation
    #r (point vector) = a(starting point) + t_dist (distance along line)*normalized distance vector
    t_from <- lapply(1:nrow(coords),function (x) 
      t(c(as.numeric(coords[x,"X"]), 
          as.numeric(coords[x,"Y"])) +
          (5000 * as.numeric(vecs[[1]]))))%>%
      unlist()%>%
      matrix(., ,2, byrow = T) %>% 
      data.frame()%>%
      mutate(Lid = seq_along(1:nrow(.))) 
    
    colnames(t_from)[c(1,2)] <- c("X", "Y")
    
    t_to <- lapply(1:nrow(coords),function (x) 
      t(c(as.numeric(coords[x,"X"]), 
          as.numeric(coords[x,"Y"])) -
          (1000000 * as.numeric(vecs[[1]]))))%>%
      unlist()%>%
      matrix(., ,2, byrow = T) %>% 
      data.frame()%>%
      mutate(Lid = seq_along(1:nrow(.)))
    
    colnames(t_to)[c(1,2)] <- c("X", "Y")
    #make line form coordinates
    Lcoords <- rbind(t_from, t_to) %>%
      arrange(Lid)
    
    #create transect line string
    t_lines<-st_as_sf(data.frame(Lcoords), coords = c("X","Y"), crs = 5070) %>%
      group_by(Lid) %>% 
      summarize(id = unique(Lid)) %>% 
      st_cast("LINESTRING")%>%
      mutate(split_fID = pt$split_fID, SegID = pt$SegID) 
    
    return(t_lines)
  }
  
  #returns endpoints of a line segment 
  endpts <- function(segs){
    seg_cords <- st_coordinates(segs)
    seg_cords <- st_as_sf(data.frame(seg_cords[c(1, nrow(seg_cords)),], 
                                     split_fID = unique(segs$split_fID), 
                                     SegID = unique(segs$SegID)), 
                          coords = c("X","Y"), crs = 5070)
  }
  
  #draw lines perpendicular to endpoints for all reach segments 
  tops <- lapply(split(reach_seg, reach_seg$split_fID), function(x) 
    funline(pt = endpts(x), 
            ln = cntr_ln[cntr_ln$split_fID == unique(as.character(x$split_fID)),],
            perpendicular = T))
  
  #list of list with (top and bottom of meander)
  #need ot select list elements it two steps
  #for matching names with list of list
  #ls is nested list. match names is character vestor
  select_list <- function(ls, matchname){
    temp<-ls[names(ls) == matchname]
    temp<-temp[[1]]
    return(temp)
  }
  
  #create channel belts
  cb <- suppressWarnings(
    lapply(apex_pts, function(x) 
      #select 1st and last row bc they will be amplitude closest to end of segment 
      lapply(split(x, x$side), 
             function(x) if(nrow(x)>1){x[c(1,nrow(x)),]} else{x[1,]})) %>%
      #Draw parallel lines at each apex point 
      lapply(., function(x) 
        lapply(x, function(x) #lapply transects over y - yeilds lines that are parellel to celter line and pass through apex
          funline(pt = x, ln = cntr_ln[cntr_ln$split_fID == as.character(unique(x$split_fID)),], perpendicular = F)))%>%
      #looking for points that intersect tops and sides
      lapply(., function(x) 
        lapply(x, function(x) 
          st_intersection(x, select_list(ls = tops, matchname = as.character(unique(x$split_fID))))))%>%
      #match intersect points with proper apex (id.1 is form tops and id is form apex_pts)
      #if there is only one point, apex needs to intersect both top and bottom
      lapply(., function(x) 
        lapply(x, function(x) 
          if(nrow(x)>2){ x[x$id == x$id.1,]} else {x[c(1:nrow(x)),]})) %>%
      #Dissolve nested list (i.e. roll up "side" sublists)
      lapply(., function(x) do.call(rbind, x)))
  
  #do.call separated here bc i added apex points and want to group them with the top and site points 
  cb <- lapply(c(cb, apex), function(x) x[,c("split_fID", "SegID")])%>%
    do.call(rbind, .) 
  
  cb <- lapply(split(cb,cb$split_fID), 
               function(x) st_as_sf(data.frame(split_fID = unique(x$split_fID), 
                                               SegID = unique(x$SegID), 
                                               st_sfc(st_convex_hull(st_combine(x)))))) %>% 
    do.call(rbind,.)
  
  #there were linestrings (probably form completely straight reach setments) which were causing error
  cb <- cb[(st_geometry_type(cb)=="POLYGON"),]
  
  return(cb)
}

#' @title Clip raster function. 
#' @description used extensively in floodplain delineaton
#' @param raster raster file to be clipped
#' @param polygon file to clip raster 
#' @return clip raster
#' @export
#' 
clipraster <- function(raster, polygon){
  cr <- crop(raster, extent(polygon), snap = "out")                    
  fr <- rasterize(polygon, cr)   
  out <- mask(x = cr, mask = fr)
  return(out)}

#' @title adjust catchment area for sampling points
#' @description adjust catchment area to accouunt for sampling location in middle of 
#' nhdPlusHR flowline
#' @param s_pts reach sampling points
#' @param fac flow accululation grid
#' @param cat catchment raster
#' @return catchment area adjustment - to be substracted from NHDPlus Total Area sqkm
#' @export
#' 
adj_cat_area <- function (s_pts, fac, cat, max_cat){ 
  
  #catseed <- raster(paste0(raster.dir,"/catseed.tif"))
  
  # check coordinate systems are the same
  #st_crs(cat) != st_crs(s_pts) &
  #st_crs(s_pts) != st_crs(fac)
  
  coords <- data.frame(st_coordinates(s_pts))
  #need to make sure sample point aligns with stream skeleton raster
  cells <- cellFromXY(cat, coords)
  gridcode <- cat[cells]
  cells <- data.frame(split_fID = s_pts$split_fID, cells, gridcode)
  
  ################
  # Identify heighest fac value in 3X3 
  #################
  
  rows <- rowFromY(fac, coords$Y)
  cols <- colFromX(fac, coords$X)
  
  ext <- data.frame(left = sapply(rows, function(x) x - 3), 
                    right = sapply(rows, function(x) x + 3),
                    up = sapply(cols, function(x) x + 3),
                    down = sapply(cols, function(x) x - 3))
  
  ext <- sapply(c(1:nrow(ext)), 
                function(x) 
                  extent(fac, 
                         ext$left[x], 
                         ext$right[x], 
                         ext$down[x], 
                         ext$up[x]))
  
  fac_ext <- lapply(ext, function(x) cellsFromExtent(fac, x))
  
  #extract on xy values is MUCH faster then extract by cells 
  fac_ext <- lapply(fac_ext, function(x) data.frame(xyFromCell(fac, x)))
  names(fac_ext) <- paste0("seq_",as.character(1:length(fac_ext)))
  fac_ext <- do.call(rbind, fac_ext)
  fac_ext <- data.frame(fac_ext, id = row.names(fac_ext))
  
  vals <- extract(fac, fac_ext[, c("x", "y")])
  fac_ext <- data.frame(fac_ext, fac_vals = vals)
  
  fac_ext$id <- unlist(lapply(strsplit(as.character(fac_ext$id), "\\."),"[[", 1))
  fac_ext$id <- as.numeric(unlist(lapply(strsplit(as.character(fac_ext$id), "_"),"[[", 2)))
  ss_fac <- lapply(split(fac_ext, f = fac_ext$id), function(x) max(x$fac_vals))
  ss_fac <- do.call(rbind, ss_fac)
  
  #############################  
  
  cells <- data.frame(cells, ss_fac)
  
  #identify max FAC in catchment zone. it is the catchment area at outlet
  #may not match the total area from NHD bc of divergent flow paths
  zone_max <- max_cat 
  
  cells <- merge(cells, 
                 zone_max, 
                 by ="gridcode", 
                 all.x = T)
  
  if(any(is.na(cells$ss_fac))){
    cells$ss_fac[is.na(cells$ss_fac)] <- 0
  }
  
  #head(cells)
  cells$area_adj <- cells$fac - cells$ss_fac
  
  #the FAC at a site can exceed FAC of catchment
  #because it was moved to the highest FAC within 3x3 area
  #these values are set to outlet area at outlet (no adjustment necessary)
  cells[cells$area_adj < 0, "area_adj"] <- 0
  cells$area_adj_sqkm <- (cells$area_adj * 10 ^ 2) / 1000 ^ 2
  
  return(cells)
}

#' @title attribute Reach sampling points usinv vaa catch area
#' @description temp fix to attr_dtf bc issues with estimates from fac
#' @param s_pts reach sampling points
#' @param fac flow accululation grid
#' @param dtf depth to flood constant multiplied by bankfull depth
#' @param ss_ras stream skeleton raster provided with NHDHR
#' @return attributed reach sampling points as SF object
#' @export
#' 
attr_dtf_vaa <- function(s_pts, dtf, dem, nhd_attr, cat_adjust){ 
  
  if(!any(names(nhd_attr) == "TotDASqKm")){
    stop("need TotDASqKm from nhd_attr")
  }
  
  nhd_attr <- Reduce(function(x,y) merge(x, y, by = "split_fID"), 
                     list(s_pts, nhd_attr, cat_adjust))
  
  ss_fac <- data.frame(nhd_attr$TotDASqKm - nhd_attr$area_adj_sqkm) 

  #cannot have negative drainage area; update to nhd area
  if(any(ss_fac[,1] < 0)){
    i <- which(ss_fac[,1] < 0)
    ss_fac[i, 1] <- st_set_geometry(nhd_attr[i, "TotDASqKm"], NULL)
  }
  
  ss_fac <- ss_fac[,1]
  ss_ele <- extract(dem, st_coordinates(nhd_attr)) #extract elevation values
  
  #hydrologic geometry relationships form Bieger et al 2015
  #used parameter estimates for entire US rgion may 
  #consider adjusting for phisiographic region
  
  #may need to adujust parameter estimates
  ss_bfh <- sapply(ss_fac, function(x) (0.30 *(x + 0.00001) ^ 0.213))  
  ss_dtf <- ss_bfh * dtf #multiply bankfull height by dft value 
  
  #multiplied dtf by 100 to convert from meters to centimeters
  #and match units of dem
  ss_idw <- (ss_dtf * 100) + ss_ele 

  s_pts <-nhd_attr %>% dplyr::select("split_fID")%>%
    mutate(fac_cat = ss_fac, ele = ss_ele, 
           bfh = ss_bfh, dtf = ss_dtf, fdep = ss_idw, 
           ss_flag = NA) %>% 
    dplyr::select(c("split_fID", "fac_cat", "ele", 
                    "bfh", "dtf", "fdep")) %>%
    st_set_geometry(., NULL)
  
  return(s_pts)
}


#' @title Attribute Valley Floor
#' @description attribute reach transect with valley floor
#' @param s_pts reach sampling points
#' @param dem digital elevation model 
#' @param v_tran reach transects
#' @return attributed reach transects points as SF object
#' @export
#' 
valleyfloor_w <- function(s_pts, dem, v_tran, attr_Points){
  
  if(any(names(v_tran) %in% c("fdep", "vf_wid"))){
    v_tran <- v_tran %>%
      dplyr::select(-(which(names(v_tran) %in% c("fdep", "vf_wid"))))
  }
  
  #pts_df <- st_set_geometry(s_pts, NULL)
  pts_df <- attr_Points
  #!names(pts_df)%in%c("split_ID", "fdep")
  v_tran <- merge(v_tran, pts_df[, c("split_fID", "fdep")], by = "split_fID")
  #vf_wid <- NULL 
 
  out <- data.frame()
  
  # add vertex to line to extract elevation values from dem
  for (i in 1:nrow(v_tran)){
    #i <- 1
    #insert points along transect line to speed up data extraction
    temp <- insert_vertex(v_tran[i,])
    temp <- data.frame(split_fID  = v_tran$split_fID[i], 
                       ID = v_tran$ID[i], 
                       fdep = setNames(v_tran$fdep[i], NULL), 
                       temp)
    
    out <- rbind(out,temp)
  }
  
  # extracting all elevation values at once saves a lot of time
  out <- data.frame(out, elevation = (extract(dem, out[,c("x","y")])))
    
  #names(out_2)[3] <- "fdep"
  #head(out_2)
  #flood = N is above water; the point that is closest to s_pts is fpwid
  #anything underwater is removed
  out$flood <- ifelse(out$elevation >= out$fdep, "N", "Y")
  out <- out[out$flood == "N", ]
  
  #create points for unflooded pixels 
  out <- st_as_sf(out, coords = c("x", "y"), crs = st_crs(s_pts))
  out <- split(out, factor(paste0(out$split_fID, "_", out$ID)))
  
  #identify shortest distance to unflooded point along transect
  out <- lapply(out, function(x) min(st_distance(x, s_pts[s_pts$split_fID == unique(x$split_fID),])))
  
  out <- do.call(rbind, out)
  colnames(out) <- "vf_wid"
  
  rename <- strsplit(row.names(out), split = "_")
  rename <- setNames(data.frame(do.call(rbind, rename)), c("split_fID", "ID"))
  
  out <- data.frame(rename, out)
    
  v_tran <- merge(v_tran, out, by = c("split_fID", "ID"), all.x = T)
  
  #if all pixels along transect are flooded vf_wid == v_wid (lenght of transect)
  v_tran[is.na(v_tran$vf_wid),"vf_wid"] <- st_length(v_tran[is.na(v_tran$vf_wid), ])
  v_tran <- st_set_geometry(v_tran,NULL) %>%
    dplyr::select(-"SegID")  
  
  v_tran <- lapply(split(v_tran, v_tran$split_fID), 
                   function(x) 
                     data.frame(split_fID = unique(x$split_fID),
                                v_wid = sum(x$v_wid), 
                                vf_wid = sum(x$vf_wid), 
                                v_vf_ratio = sum(x$vf_wid)/sum(x$v_wid)))
  
  v_tran <- do.call(rbind,v_tran)
  return(v_tran)
}


#' @title Attribute channel belt width
#' @description attribute reach transect with channel belt width
#' @param belts channel belts 
#' @param Reach_Transects reach transects
#' @return attributed reach transects points as SF object
#' @export
#' 
belt_w <- function(Reach_Transects, belts){
  belt_w <- suppressWarnings(lapply(split(Reach_Transects, Reach_Transects$split_fID),
                                    function(x) st_intersection(x, belts[as.character(belts$split_fID) == 
                                                                           unique(x$split_fID) ,]))%>%
                               lapply(., function(x) data.frame(split_fID = x$split_fID, ID = x$ID, belt_w = st_length(x))))
  belt_w <- do.call(rbind, belt_w) 
  Reach_Transects <- merge(Reach_Transects, belt_w, by = c("split_fID", "ID"), all.x = T)
  Reach_Transects$belt_w <- ifelse(is.na(Reach_Transects$belt_w), Reach_Transects$v_wid, Reach_Transects$belt_w)
  return(Reach_Transects)
}


#' @title Attribute channel shape
#' @description attribute reach segements with channel shapes
#' @param reach_seg reach segments
#' @return attributed reach segments as SF object
#' @export
#' 
channel_shp <- function(reach_seg = Reach_segments){
  #order by split id for matching later on
  reach_seg <- reach_seg %>% arrange(split_fID)
  
  #select first and last coordinates to 
  #create centerline form line segment
  #input to function is a linestring
  centerline <- function(x) {
    coord <- st_coordinates(x)
    coord <- coord[c(1, nrow(coord)),]
    coord <- st_as_sf(data.frame(SegID= x$SegID, split_fID=x$split_fID, coord), coords = c("X","Y"))%>% 
      group_by(split_fID)%>%
      summarise(SegID = unique(SegID))%>%
      st_cast("LINESTRING")%>%
      mutate(CntLn_Len_m = st_length(.))%>%
      st_set_crs(5070)
  }
  
  #calculate centerlines for all reach segments
  #could be used later on to calculate reach sinuosity
  cntr_ln <- lapply(split(reach_seg,reach_seg$split_fID), function(x) centerline(x))%>%
    do.call(rbind,.)
  
  cntr_ln$ch_sinu <- reach_seg$geom_len / cntr_ln$CntLn_Len_m
  #off diag eleements are meander length (distance between successive points)
  #may want to count meanders later? 
  #input is a matrix
  offdiag <- function(x){
    x[row(x) == (col(x) - 1)]
  }
  
  #reach_seg that do not intersect are just center line lengths
  #every reach segment (that is not completely straight line)
  #had meander length
  meander_lengths <- suppressWarnings(
    lapply(split(cntr_ln, cntr_ln$split_fID), 
           function(x) st_intersection(x, reach_seg[reach_seg$split_fID == x$split_fID,])%>%
             st_cast("POINT")%>%st_distance()%>%
             offdiag()))
  
  ml <- lapply(meander_lengths, 
               function(x) data.frame(meanders = length(x),
                                      wavelength = length(x)/2,
                                      mean_meander_length = mean(x)))
  ml <- do.call(rbind, ml)
  
  reach_seg <- merge(reach_seg, 
                     data.frame(split_fID = cntr_ln$split_fID, 
                                CntLn_Len_m = cntr_ln$CntLn_Len_m,
                                ch_sinu = cntr_ln$ch_sinu,
                                ml),  by = "split_fID")
  
  reach_seg <- st_set_geometry(reach_seg,NULL) %>% 
    dplyr::select( - "SegID")
  return(reach_seg)
}


#' @title Attribute sampling points with nhdplusHR Vaa
#' @description spatialy join sampling points to nhdplusHR network and query vaa attributes  
#' @param s_pts samling points 
#' @param nhdhr_net nhdHR shaplefile
#' @param vaa_tbl vaa table for 4-digit huc
#' @param vaa_attr column names from vaa table to add  
#' @return attributed reach segments as SF object
#' @export
#' 
attr_nhdplusHR <- function (s_pts, nhdhr_net, vaa_tbl, 
                            vaa_attr = c("StreamLeve", "StreamOrde", "StreamCalc")){

  #remove vaa names from 
  if(any(names(s_pts) %in% c(vaa_attr, "WBA_P_I"))){
    s_pts <- s_pts %>%
    dplyr::select( -which(names(s_pts)%in%c(vaa_attr, "WBA_P_I")))
  }
  
  #writing sf to shapefile shortens names to stupid shit
  if(any(names(s_pts) %in% c("StremLv", "StrmOrd", "StrmClc","TtDASqK", "AreaSqKm", "WBA_P_I"))){
    s_pts <- s_pts %>%
      dplyr::select(-which(names(s_pts) %in% c("StremLv", "StrmOrd", "StrmClc", "TtDASqK", "WBA_P_I")))
  }

  if(!is.character(vaa_tbl$NHDPlusID)){
    stop("NHDPlusID must be character")
  }
  
  #attribute points wiht NHDPlusID's
  # Spatial join to retreive NHDPlusID's
  nhdhr_net <- st_zm(nhdhr_net, drop = T, what = "ZM")
  nhdhr_net <- st_transform(nhdhr_net, crs = st_crs(Reach_Points))

  #small rounding errors require
  nhdhr_net_index <- st_is_within_distance(s_pts, nhdhr_net, dist = 0.001) 
  nhdhr_net_index <- do.call(rbind, nhdhr_net_index)

  if(ncol(nhdhr_net_index)>1){
    stop("try smaler distance")
  }
  
  NHDPlID <- st_set_geometry(nhdhr_net[nhdhr_net_index[,1], c("NHDPlID","WBA_P_I")], NULL)
  s_pts$NHDPlID <- NHDPlID[,1]
  s_pts$WBA_P_I <- NHDPlID[,2]

  s_pts <- merge(s_pts,
                 vaa_tbl[c("NHDPlusID", vaa_attr)], 
                 by.x ="NHDPlID", 
                 by.y = "NHDPlusID", all.x = T)

  s_pts <- st_set_geometry(s_pts, NULL)
  
  return(s_pts)
}


#' @title Range normalize function
#' @description normaize data form 0-1
#' @param x range of values to normalize
#' @return normalized data
#' @export
#' 
normal_range <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}


#' @title Functional Process Zone identification
#' @description agglomerative hierarchical clustering to dentify functional process zone 
#' @param x dataframe with FPZ values 
#' @parm FPZ_vars Variable string to cluster. Must be same names as column headings 
#' @parm normalize should normalize continuous variables
#' @parm max number of groups to calculate 
#' @return numeric vector of FPZ Class
#' @export
#'

FPZ <- function(x, FPZ_vars, normalize = T, nmax = 25, 
                saveplot = T, 
                savecluster = T, 
                write.dir = write.dir){
  
  if(any(is.na(x[,FPZ_vars]))){
    stop("na values in FPZ_vars not allowed")
  }
  
  if(any(!FPZ_vars%in%names(x))){
    stop("FPZ_vars missing from x")
  }
  
  if(normalize){
    norm_cols <- sapply(x[,FPZ_vars], function (x) class(x)) !="factor"
    x[,FPZ_vars[norm_cols]] <- apply(x[, FPZ_vars[norm_cols]], 2, function(x) normal_range(x))
  }
  
  if(nmax > nrow(x)){
    stop("too many groups; adjust nmax")
  }
  
  # generate distance matrix 
  d <- cluster::daisy(x[, FPZ_vars], metric = "gower")

  # agglomerative hierarchical clustering 
  hcluster <- cluster::agnes(d, method = "ward", 
                             keep.diss = F, 
                             keep.data = F)

  # average sillouette width optimal number of groups 
  SI <- numeric(nmax)
  
  for(k in 2:length(SI)){
    sil <- cluster::silhouette(cutree(hcluster, k), d)
    SI[k] <- summary(sil)$avg.width
  }
  
  # max silhouette is best grouping solution
  #select top three max silhouette, then choose one with most patches - its more intersting
  k.best <- order(SI, decreasing = T)[1:3]
  #select one with most clusters
  k.best.sil <- max(k.best) 
  
  #FPZ.sil <- sapply(1:k.best.sil, function(x) factor(cutree(hcluster, x)))
  #colnames(FPZ.sil) <- c(paste0("FPZ.sil", "_", 1:(k.best.sil - 1)), "FPZ.sil_best")
  
  FPZ.sil <- factor(cutree(hcluster, k.best.sil))
  
  #Comparisions betweeen the dissimilarity matrix and binary matrices representing partitions
  #kt <- data.frame(k = 1:nmax, r = 0)
  #for (i in 2:nmax){
   # gr <- cutree(hcluster, i)
    #veg <- as.data.frame(as.factor(gr))
    #distgr <- daisy(veg,"gower")
    #mt<-cor(d,distgr, method = "pearson")
    #kt[i, 2] <- mt
  #}
  
  #grouping with highest correlation coefficients
  #k.best.cor <- which.max(kt[, 2])
  #FPZ.cor <- factor(cutree(hcluster, k.best.cor))
  
  ########################
  if(saveplot){
  
    if(length(write.dir)==0){
      stop("provide write.dir")
    }
    
    n <- unlist(strsplit(write.dir, "/"))
    n <- n[length(n)]
    dend <- as.dendrogram(hcluster)
  
    # plot dendrogram with some cuts
    branch_heights <- get_branches_heights(dend,sort = T, decreasing = T)
  
    dend <- cut(dend, h = branch_heights[k.best.sil])$upper
    labels(dend) <- as.character(c(1:k.best.sil))
    
    windows()
    dend %>% 
      hang.dendrogram(hang = -1) %>% # un-hanging the leaves
      plot(main = n,  ylim = c(0, max(branch_heights)+1), 
           center = T, ylab = "height")
  
    savePlot(filename = paste0(write.dir,"/dendrogram.jpg"), type = "jpg")
    dev.off(which = dev.cur())
    
    
    dend_nodes <- group_dendrogram_nodes(dend = dend)
    write.csv(dend_nodes, paste0(write.dir, "/dendrogram_nodes.csv"), row.names = F)
}
  
  if(savecluster){
    save(hcluster, file = paste0(write.dir, "/hclust"))
  }
  ########################
  
  
  #return(data.frame(FPZ.sil, FPZ.cor))
  return(data.frame(FPZ.sil))
}


#' @title Functional Process Zone identification
#' @description agglomerative hierarchical clustering to dentify functional process zone 
#' @param Net_char dataframe with FPZ values 
#' @parm FPZ_vars Variable string to cluster. Must be same names as column headings 
#' @return FPZ Class attributed to network characteristic file
#' @export
#'  
  
cluster_fpz <- function(Net_Char, FPZ_vars){ 
  
  if(!all(complete.cases(Net_Char[, FPZ_vars]))){
    stop("NA values in FPZ_Vars")
  }
  
  # select complete cases, exclude WB - they are already classified
  all_net <- Net_Char[complete.cases(Net_Char[, FPZ_vars]) & 
                        is.na(Net_Char$WBA_P_I), ]
  
  # order to sample along catchment area gradient
  all_net <- all_net[order(all_net$fac_cat), ]
  
  # 30 equally spaced intervals along catchment gradient  
  sample_index <- all_net[all_net$fac_cat > 0, "split_fID"]
  sample_index <- sample_index[seq(1, length(sample_index), 
                                   by = ceiling(length(sample_index) / 30))]
  
  # indicate small scale sampling networks 
  all_net$sample <- NA
  all_net[all_net$splt_ID %in% sample_index, "sample"] <- 1
  
  #identify FPZ  
  nmax <- ifelse(nrow(all_net) > 25, 25, nrow(all_net) - 1)
  fpz <- FPZ(x = all_net, FPZ_vars = FPZ_vars, normalize = T, nmax = nmax)
  
  all_net$FPZ_all <- fpz
  all_net$FPZ_all_n <- nrow(all_net)
  all_net$FPZ_all_areasqkm <- max(all_net$TtDASqK) 
  all_net$FPZ_all_sub_netid <- "ALL"
  
  
  #add back in waterbodies
  WB <- Net_Char[complete.cases(Net_Char[, FPZ_vars]) & !is.na(Net_Char$WBA_P_I), ]
  
  if (nrow(WB)>0){
    levels(all_net$FPZ_all) <- c(levels(all_net$FPZ_all), "Waterbody")
    WB$sample <- NA
    WB$FPZ_all <- "Waterbody"
    WB$FPZ_all_n <- nrow(all_net)
    WB$FPZ_all_areasqkm <- max(all_net$TtDASqK) 
    WB$FPZ_all_sub_netid <- "ALL"
    all_net <- rbind(all_net, WB)
  }
  
  return(all_net)
}


#' @title Functional Process Zone Summary
#' @description calculate network scale summary metrics of functional Process Zone 
#' @param all_net output form cluster_fpz  
#' @param show_progress print count 
#' @return FPZ Class attributed to network characteristic file
#' @export
#'  

net_summary <- function (HPG_Sample = sample_index, 
                         HGP = FPZ, show_progress = T){
  
  NHDPlID <- sample_index$NHDPlID
  
  if(length(NHDPlID)==0){
    stop("need nhdplus ids")
  }
  
  FPZ_Summary <- data.frame()
  count<-1
  for (i in 1:nrow(sample_index)){
    # i <- 1
    # all flowlines above the root node
    z <- net_delin_HR(huc = x[site, "HUC4"], NHDPlusID = as.character(sample_index$NHDPlID[i]))
    
    # identify sub net  
    subnet <- HGP[HGP$NHDPlID %in% z$network$net_NHDPlusID, ]
    
    # calculate patch metrics here
    Root_ID <- as.character(sample_index$NHDPlID[i])
    areasqkm <- max(subnet[subnet$NHDPlID == as.character(sample_index$NHDPlID[i]), "fac_cat"])
    net_length <- sum(subnet[,c("geom_len")])
    
    FPZ_Richness <- length(unique(as.character(subnet$FPZ)))
    
    Shannon_index <- lapply(split(subnet, as.character(subnet$FPZ)), 
                            function (x) (sum(x$geom_len)/net_length) * 
                              (log(sum(x$geom_len)/net_length)))
    Shannon_index <- (-1) * sum(unlist(Shannon_index))
    
    temp <- data.frame(NHDPlID = unique(as.character(sample_index$NHDPlID[i])),
                       areasqkm, net_length, FPZ_Richness, Shannon_index)
    
    FPZ_Summary <- rbind(FPZ_Summary, temp)  
    
    if (show_progress){
      print(count)
    }
    
    count <- count + 1
  }
  
  return(FPZ_Summary)
}


#' @title insert vertex
#' @description place vertex to extracr vlaues 
#' @param ln sf lin to split  
#' @parm point_intervals intervals to insert
#' @return numeric vector of FPZ Class
#' @export
#'
insert_vertex <- function(ln, point_intervals = 10){
  
  split <- 1/ceiling(st_length(ln)/point_intervals)
  units(split) <- NULL
  ln <- data.frame(st_coordinates(ln))
  
  #direction
  direct <- c(ln$X[1] - ln$X[2], ln$Y[1] - ln$Y[2])
  
  #distance intervals for splitting line
  lambda <- seq(0, 1, by = split)

  # vector line equation: r = (point on line) + distance * direction 
  new.pt <- sapply(1:length(lambda), function(x) ln[1, -3] + -lambda[x] * direct)
  new.pt <- t(new.pt)
  
  # fixing list from sapply
  new.pt <- data.frame (x = do.call(rbind, new.pt[,1]), y  = do.call(rbind, new.pt[,2]))
  
return(new.pt)
}


#' @title identify node groupings of leaf in dendrogram
#' @description provides grouping table to investigate differenes between groups at each node 
#' @param ln sf lin to split  
#' @return dataframe with groupings
#' @export
#'
group_dendrogram_nodes <- function(dend){
  
  subtrees <- partition_leaves(dend)
  #init is focal group and will update
  init <- subtrees[1]
  compare.subtrees <- subtrees[-1]
  group.num <- 1
  group.out <- data.frame()
  
  while(length(init) > 0){
    group.table <- data.frame()
    for (j in 1:length(init)){
      #j <- 2 
      #print(j)
      #group.table <- data.frame()
      compare.subtrees <- compare.subtrees[unlist(lapply(compare.subtrees, function(x) !all(init[[j]] %in%x)))]
      #I need to iterate through all sub trees
      #focus <- init[[j]]
      for (i in c(1:length(compare.subtrees))){
        #i <- 1
        
        matched <- init[[j]][init[[j]] %in% compare.subtrees[[i]]]
        
        if(length(matched) > 0 & group.num == 1){
          gr1 <- matched
          init[[j]] <- init[[j]][!init[[j]] %in%  matched]
          group.num <- group.num + 1
        }
        
        if(length(matched) > 0 & group.num == 2){
          gr2 <- matched
          init[[j]] <- init[[j]][!init[[j]] %in%  matched]
        }
        
      }
      
      group.num <- 1
      group.table.temp <- data.frame(gr1 = paste0(gr1, collapse = " "), 
                                     gr2 = paste0(gr2, collapse = " "), 
                                     stringsAsFactors = F)  
      group.table <- rbind(group.table, group.table.temp)
    }
    
    #group table is going to grow so need to change by row...
    
    #Which values are greater than 1 - these ned to be split
    init <- unlist(group.table)
    init <- lapply(split(init, f = names(init)),function(x) unlist(strsplit(x," ")))
    init <- init[unlist(lapply(init, function(x) length(x) >= 2))]
    group.out <- rbind(group.out, group.table)
    
  }
  return(group.out)
}


#vertex <- insert_vertex (ln = v_tran[4367,], 10)
#a <- data.frame(x = 75, y = 100)
#b <- data.frame(x = 1, y = 1)
#ln <- rbind(a, b)
#ln <- st_linestring(as.matrix(ln))
#str(ln)

#plot(st_geometry(v_tran[4367,])) #[,c(1, 2)], col = c("red", "green"), pch = 19, cex = 2)
#points(vertex, col = "blue", pch = 19)

#st_length(st_linestring(as.matrix((vertex[c(5,6),]))))


#log sequence
lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

#Merging line sigements by geometries 
#function takes a series of flowlines and identifies indices with common geometries
#returns list with indecies of flowlines that are the same patch
#flowlines <- net_p1[[1]]

#' @title merge line segements 
#' @description used for merging line segments of the same type, used internally in class metrics 
#' @param flowline lines to merge  
#' @return list of merged lines
#' @export
#'
merge_lines <- function(flowlines){
  
  net_int <- st_intersection(flowlines)
  
  #merge the list elements that contain the same index
  inlist <- net_int$origins
  
  outlist <- list()
  count <- 1
  
  #inlist elements are iterariively removed 
  while(length(inlist) > 0){  
    j <- inlist[[1]] 
    #iteratively match across list elements 
    while(any(do.call(rbind, lapply(inlist, function (x) any(x %in% j))))){
      #index of matching elements
      ind <- do.call(rbind, lapply(inlist, function (x) any(x %in% j)))
      #new values to combine
      j <- unique(append(j, unique(unlist(inlist[ind])), length(j)))
      #update to remove elements that need to be combined
      inlist <- inlist[!ind]
    }
    
    #update new list with values
    outlist[[count]] <- j
    count <- count + 1
  }
  
  if (length(setdiff(1:nrow(flowlines), sort(unlist(outlist))))!=0|
      any(duplicated(sort(unlist(outlist))))){
    
    stop("something is wrong")
  }
  
  return(outlist)
}


#' @title sf to tidygraph 
#' @description adapted from https://www.r-spatial.org/r/2019/09/26/spatial-networks.html embedded within class metrics
#' @param x split lines object  
#' @return super useful graph
#' @export
#'
sf_to_tidygraph <- function(x, directed = TRUE) {
  #x<-splitln
  edges <- x %>%
    mutate(edgeID = c(1:n()))
  
  nodes <- edges %>%
    st_coordinates() %>%
    as_tibble() %>%
    rename(edgeID = L1) %>%
    group_by(edgeID) %>%
    slice(c(1, n())) %>%
    ungroup() %>%
    mutate(start_end = rep(c('start', 'end'), times = n()/2)) %>%
    mutate(xy = paste(.$X, .$Y)) %>% 
    mutate(nodeID = group_indices(., factor(xy, levels = unique(xy)))) %>%
    dplyr::select(-xy)
  
  source_nodes <- nodes %>%
    filter(start_end == 'start') %>%
    pull(nodeID)
  
  target_nodes <- nodes %>%
    filter(start_end == 'end') %>%
    pull(nodeID)
  
  edges = edges %>%
    mutate(from = source_nodes, to = target_nodes)
  
  nodes <- nodes %>%
    distinct(nodeID, .keep_all = TRUE) %>%
    dplyr::select(-c(edgeID, start_end)) %>%
    st_as_sf(coords = c('X', 'Y')) %>%
    st_set_crs(st_crs(edges))
  
  tbl_graph(nodes = nodes, edges = as_tibble(edges), directed = directed)
  
}


#' @title average distance between patches  
#' @description embedded within class metrics
#' @param plen a list with elements of flowlines of each patch type 
#' @param patch the patch type identified 
#' @param mu defines the rate of exponential distribution 
#' median parent offspring of 12km = 12,000m (Comte and Olden 2017)
#' mu = 1000; mean(rexp(100000, rate = 1/mu)) ~ 1,000
#' @return super useful graph
#' @export
#'
avgdist <- function(plen, patch = 2, mu = 12000){
  
  #identify matrix and patches  
  npatch <- do.call(rbind, plen[patch])
  npatch$patchid <- 1:nrow(npatch)
  names(npatch)[2] <- "geometry"
  st_geometry(npatch) <- "geometry"
  
  #plot(st_geometry(npatch))
  nmatrix <- do.call(rbind, plen[-patch])
  #nmatrix<-st_combine(nmatrix)
  w <- merge_lines(nmatrix)
  nmatrix <- lapply(w, function(x) st_sf(data.frame("matrix"), sf::st_combine(nmatrix[x,])))
  nmatrix <- do.call(rbind, nmatrix)
  nmatrix$matrixid <- 1:nrow(nmatrix)
  names(nmatrix)[2] <- "geometry"
  st_geometry(nmatrix) <- "geometry"
  
  intnodes <- sf::st_intersection(npatch, nmatrix)
  #dups <- intnodes$matrixid[duplicated(intnodes$matrixid)]
  #intnodes <- intnodes[intnodes$matrixid%in%dups, ]
  
  #create tidy graph
  #####
  full.net <- do.call(rbind, plen)
  splitln <- st_cast(full.net, "LINESTRING")
  names(splitln)[2] <- "geometry"
  st_geometry(splitln) <- "geometry"
  
  #function sf_to_tidygraph adapted from 
  #https://www.r-spatial.org/r/2019/09/26/spatial-networks.html
  graph <- sf_to_tidygraph(splitln, directed = F)
  
  #calculate length of edges
  graph <- graph %>% activate(edges) %>% mutate(length = st_length(geometry))
  
  #Match nodes from intersection
  nodeid <- graph %>% activate(nodes) %>%
    as_tibble() %>% st_as_sf() %>%  
    st_intersection(intnodes)
  ##################
  
  distances <- distances(graph = graph, 
                         weights = graph %>% activate(edges) %>% pull(length))
  
  #reduce distance matrix to notes connecting patches
  distances <- distances[nodeid$nodeID, nodeid$nodeID]
  
  #tosses error if not matrix... 
  if(!is.matrix(distances)){
    distances<-matrix(distances)
  }
  
  #format distance table
  ########
  rownames(distances) <- colnames(distances) <- nodeid$nodeID
  #melt into To/From table 
  distances <- reshape2::melt(distances)
  names(distances) <- c("To_nodeID", "From_nodeID", "dist_m")
  distances <- distances[!is.infinite(distances$dist_m), ]
  
  patchlength <- data.frame(patchid = npatch$patchid, len_m = st_length(npatch))
  #Add fields 
  distances$To_patchid <- nodeid$patchid[match(distances$To_nodeID, nodeid$nodeID)]
  distances$From_patchid <- nodeid$patchid[match(distances$From_nodeID, nodeid$nodeID)]
  distances$matrixid <- nodeid$matrixid[match(distances$To_nodeID, nodeid$nodeID)]
  distances$To_patch_len_m <- patchlength$len_m[match(distances$To_patchid, patchlength$patchid)]
  distances$From_patch_len_m <- patchlength$len_m[match(distances$From_patchid, patchlength$patchid)]
  
  #there is "flow thru" the patch needs to be removed 
  #such that the values flow to into patch
  distances <- distances[order(distances$From_patchid, distances$To_patchid, distances$dist_m), ]
  distances <- split(distances, list(distances$To_patchid, distances$From_patchid), drop = TRUE)
  distances <- do.call(rbind, lapply(distances, function(x) x[1,]))
  ##################  
  
  if(is.na(mean(unique(distances$dist_m)))){
    DCI <- NA
    } else {
    
      #DCI Calculation
      li <- (distances$From_patch_len_m / sum(unique(distances$From_patch_len_m)))
      lj <- (distances$To_patch_len_m / sum(unique(distances$From_patch_len_m)))    
  
      #this needs to add to 1 
      #if(!all.equal(sum(li*lj), 1)){
       # stop("remember me mother fucker")
      #}
      
      #probability is distance based on exponential distribution 
      #i.e. the probability of making it past a given distance
      p <- 1 - pexp(distances$dist_m, rate = 1/mu)
  
      #DCI calculation 
      DCI <- sum(p * li * lj) * 100
      
      #this needs to add to 1 
      #if(!all.equal(sum(li*lj), 1)){
        #stop("remember me mother fucker")
      #}
      
  #plotting examples - patch connetivity
  #######
  #pid <- st_set_geometry(z[z$matrixid == intnodes, "patchid"], NULL)
  #windows()
  #plot(st_geometry(nodeid[nodeid$nodeID %in% c(176, 40), ]), 
   #    col = "green", pch = 19, cex = 1.5, add = T)
  #plot(st_geometry(nodeid[nodeid$nodeID %in% c(175, 40), ]), 
   #   col = "red", pch = 19, cex = 1.5, add = T)
  #plot(st_geometry(npatch[npatch$patchid==4,]),lwd = 3, add = T)
  #plot(st_geometry(nmatrix), add = T, col = "purple")
  
  #plot(st_geometry(full.net))
  #plot(st_geometry(npatch), col = "purple", add = T)
  #plot(st_geometry(nmatrix), col = "red", add = T)
  #plot(st_geometry(npatch), add = T)
  
  #metric path
  #plot(st_geometry(nmatrix[nmatrix$matrixid == 3, ]), add = T)
  #plot(st_geometry(nodeid), add = T, col = "purple", pch = 19)
  #plot(st_geometry(nodeid[nodeid$nodeID == 42, ]), 
  #     add = T, col = "green", pch = 19, cex = 1.5)
  #head(distances[distances$From_nodeID == 40, ])
  #for (i in 1:nrow(distances[distances$From_nodeID == 200, ])){
    #i<-3
   # if(distances$From_nodeID[i]!=distances$To_nodeID[i]){
    #identify to/from nodes to navigate - used 1 and 2 for example
    #from_node <- graph %>% activate(nodes) %>% filter(nodeID == distances$From_nodeID[i]) %>% pull(nodeID)
    #to_node <- graph %>% activate(nodes) %>% filter(nodeID == distances$To_nodeID[i]) %>% pull(nodeID)
    
    #find shortest path between two nodes identified above
    #path <- shortest_paths(graph = graph, from = from_node, to = to_node,
     #                      output = 'both',weights = graph %>% activate(edges) %>% pull(length))
    
    #prepare to plot shortest path - summing length of all edges inbetween two nodes
    #path_graph <- graph %>% subgraph.edges(eids = path$epath %>% unlist()) %>% as_tbl_graph()
    #length
    #path_graph %>% activate(edges) %>% as_tibble() %>% summarise(length = sum(length))
    
    #plot path 
    #path_graph %>% activate(edges) %>% 
     # as_tibble() %>% st_as_sf() %>% st_geometry() %>%
      #plot(lwd = 2, col = 'black', add = T)
    #}
  #}
  #legend("left", c("patch", "matrix", "path"), 
   #      col = c("purple", "red", "green"), lty = 1, lwd = 2)
  
  #Plot nodes to from nodes
  #path_graph %>% activate(nodes) %>% filter(nodeID %in% c(from_node, to_node)) %>% 
  # as_tibble() %>% st_as_sf() %>% st_geometry() %>% plot(add = T, cex = 3)
  ###################
    }
  
  return(data.frame(DCI, meandist = mean(unique(distances$dist_m)), 
                    medianDist = median(unique(distances$dist_m))))
  #print(mean(unique(distances$dist_m)))
}

#' @title Class metrics 
#' @description calculate class level metrics
#' @param network network object
#' @param subnetID netids for subnetwork   
#' @return list of class scale metrics that can be aggregsted to landscape scale 
#' @export
#'

class_scale_metrics <- function(network, subnetID){
  #T1<-Sys.time()
  #select subnetwork from the full network 
  subnet <- network[network$NHDPlID %in% subnetID, ] 
  
  #split the subnetwork by the patches to create class scale metrics 
  net_p1 <- split(subnet, as.character(subnet$FPZ.sil))  
  
  #identify adjacent flowlines for eahc patch
  mergelines <- lapply(net_p1, function(x) merge_lines(x))
  
  #recondition network into patches
  #creates list of SF multilinestring representing patches
  
  #simlify = F; will not simplify outpur to array matcis or vector
  plen <- sapply(1:length(net_p1), function(q)
    mapply(function(x, y) st_combine(x[y, ]), net_p1[q], mergelines[[q]]), simplify = F)
  
  # a single patch patch
  if(length(plen) == 1) {
    
    plen <- st_sf(data.frame(patch = "patches"), do.call(c, plen))
    
    tlen.patch_m <- sum(st_length(plen))
    mean.patch_m <- mean(st_length(plen))
    No.Patch <- nrow(plen)
    mean.patch.dist_m <- data.frame(matrix(0, length(No.Patch), 3))
    names(mean.patch.dist_m) <- c("DCI",  "meandist", "medianDist")
    mean.patch.dist_m$DCI <- 100
  } else {
    plen <- lapply(plen, function(x) st_sf(data.frame(patch = "patches"), x))
    tlen.patch_m <- do.call(c, lapply(plen, function(x) sum(st_length(x))))
    mean.patch_m <- do.call(c, lapply(plen, function(x) mean(st_length(x))))
    No.Patch <- do.call(c, lapply(plen, function(x) nrow(x)))
    
    #only run patches with >1 ... if there is only single patch of a given type 
    #there is no distance between them... the vector catches assigns 
    #0 to any patch type no run 
    mean.patch.dist_m <- data.frame(matrix(0, length(No.Patch), 3))
    names(mean.patch.dist_m) <- c("DCI",  "meandist", "medianDist")
    mean.patch.dist_m$DCI <- 100
    if(any(No.Patch>1)){
      mean.patch <- suppressWarnings(sapply(which(No.Patch>1), 
                                            function(x) avgdist(plen = plen, patch = x), simplify = F))
      mean.patch.dist_m[which(No.Patch>1), ] <- do.call(rbind,mean.patch)
    }
  }
  
  #i can get mean patch length and number of patches for plen
  out <- data.frame(tlen.patch_m, mean.patch_m, No.Patch, mean.patch.dist_m)
  return(out)
}
