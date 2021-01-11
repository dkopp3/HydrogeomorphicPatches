# HydrogeomorphicPatches
Delineate Hydrogeomorphic Patches
Delineating Hydrogeomorphic Patches 
Hydrogeomorphic patches were identified using four modules written in R programming language. The modules divide processing into discrete steps (i.e. create, attribute, and cluster) that must be run in sequence. A final module, “landscape metrics” quantifies spatial pattern in each network at different spatial scales. The documentation for the functions is provided with the source code in “./Data/ NHDHR_Functions_example.r”. Do not change field names as many functions use names to call specific columns within larger datasets. We also provide an example river network, Ozarks Complex-2 site (HUC4 = 0316), to illustrate the output of the different modules.  For large networks (i.e. > 5,000km2) processing time for an individual module can exceed 18 hours.  Below are descriptions of the output files generated from each module.

Setup: 
•	Download NHDHR data and ancillary files. We reprojected the ancillary data to match the coordinate reference system of the NHDHR dataset. For this example, we only provide “reprojected” ancillary data to save processing time.  Since the ancillary data are consistent throughout the CONUS, reprojection only needs to be done once. The attribute module documents how the data were reprojected from their original source. Web links for download are provided in main text. See “.\Data\NHDPLUSHR_3_7_2019”. 
Contents: 
./NHDPLUS_H_0316_HU4_GDB – directory for NHDHR hydrography geodatabase (HUC 0316)
./NHDPLUS_H_0316_HU4_RASTER – directory for NHDHR rasters for (HUC 0316)
soilmu.shp – Reprojected soils data (Original source: SSURGO)
BDTICM_ras.tif – Reprojected depth to bedrock (Original Source: Soil Grids)  
ppt_ras.tif – Reprojected precipitation data (Original source: PRISM climate)  
temp_ras.tif – Reprojected temperature data (Original source: PRISM climate) 
 
•	Identify study network. See “./Data/NEON_outlet_example.csv”
Field Names:
site = Unique ID for the network
HUC4 = 4- digit hydrologic unit code. Identifies the directory containing the hydrography data.
NHDPLUSID = Unique identifier for root node of the network. 
DIVERGENCE = whether flowlines is divergent 
TOTDASQKM = upstream drainage areas in square kilometers. 
STREAMORDE = stream order. 

Create: 
•	The “create module” extracts NHDHR flowlines, identified reaches, sampling points, and transects. See “./Data/Create_Example.r”. See /

Output of create module:
"/NHDHR_Network.shp" – shapefile of stream network extracted form NHDHR
"/Valley_Segments.shp" – shapefile of valley segments 
"/Valley_Catchments.shp" – shapefile of catchments surrounding valley segments
"/Reach_segments.shp" – shapefile of reach segments 
"/Reach_Points.shp" – shapefile of sampling points 
"/Reach_Transects.shp" – shapefile of lateral transects for each reach

Attribute:
•	Reprojects ancillary data and uses shapefiles form the create module to assign hydrogeomorphic variables to each reach. See “./Data/Attribute_Example.r”

Output of Attribute module 
"/attr_Points.csv" – Hydrogrophrphic variables extracted at the sampling point of each reach
"/attr_valley.csv" – Hydrogeomorphic variables calculated at the valley scale (e.g. valley slope)
"/attr_transects.csv" – Hydrogeomorphic variables calculated along reach transects (e.g. valley width)
"/attr_reach.csv" – Hydrogeomorphic variables calculated at the reach scale (e.g. sinuosity)
"/Network_HGP.csv" – Merged all hydrogeomorphic variables for every site within the network

Cluster: 
•	A convenient wrapper for “daisy”, “agnes”, and “silhouette” functions from Maechler et al. (2019) used to identify reaches with similar hydrogreomorphic variables (ie. Hydrogeomorphic patches)

Output of cluster module
"/FPZ_2.csv" – classification of each reach from silhouette analysis 
"/FPZ_impute.csv" – assigns reaches with missing values to a HGP class using randomforest classification. Missing values were often a result of slope values in wide catchments. 
"/scaling_extent.csv" – randomly identify reaches associated with different catchment areas. Used to manipulate spatial extent
"/scaling_grain.csv" – Prunes network based on minimum path length (e.g. mainstem length). Used to manipulate spatial grain
"/FPZ_impute_all.csv" – assigns classes to reaches with missing data. Used to manipulate thematic resolution

Landscape Metrics:
•	The landscape metric module was used to calculate “landscape metrics” which quantify spatial pattern in river networks. 
Output of landscape metrics module
"/cL_metrics.csv" – class level metrics for each spatial extent within a network 
"/LS_metrics.csv" – landscape level metrics for each spatial extent within a network 
“/cL_metrics_grain.csv” - class level metrics for each grain size within a network
"/LS_metrics_grain.csv" - Landscape level metrics for each grain size within a network
"/cL_metrics_thematic.csv" - class level metrics for each thematic level within a network
"/LS_metrics_thematic.csv" – landscape level metrics for each thematic level within a network
