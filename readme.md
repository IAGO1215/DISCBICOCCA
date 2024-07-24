# Description
This github repository is built to calculate variogram and roman metrics, in order to evaluate the spatial representativeness of a group of field sites.   
Currently it only supports Sentinel-2 images. 
The spatial index utilized in this project is NIRv, of which the formula is NIRv = B8(radiance) * NDVI; while in the reference articles, surface albedo is the only spatial index being used instead. 

# Caveat
1. All the sites in this project are considered as point vectors. 
2. The term "distance" in this project refers to a squared area centered at the site whose side-length equal to the value of itself. For example, when it says "30m distance" of San Rossore, it refers to a squared area whose side-length is 30m, centered at the site San Rossore.   
3. In the article of Roman, the Roman metrics are only calculated based on a 1km area (internal area) and another 1.5km area (external area). But in this project, we don't adhere to the limitation that the external area should be 1.5 times of the internal area. 

# Required User Input
1. One L2A and one L1C Sentinel-2 image for each site, following the folder structure illustrated in the section below.  
2. The site coordinates in Lat-Long. 
3. The constant g = 2Htan(FOV°) for roman metrics.  
4. The desired distances of which the variograms will be calculated
5. The desired internal distance(s) and external distance(s) of which the roman metrics will be calculated

# Required Python Packages
1. Numpy
2. Pandas
3. Matplotlib
4. Shapely
5. GeoPandas
6. Rasterio
7. lxml (pip)
8. bs4
9. scipy
10. skgstat

# Folder Structure
DISC  
├── Python  
│   ├── My Functions  
│   │   ├── xxx.py  
│   │   ├── ...  
│   ├── NoteBooks  
│   │   ├── **Main 1**.ipynb  
│   │   ├── ...  
│   ├── NoteBooks Draft  
│   │   ├── xxx.ipynb  
│   │   ├── ...  
├── Sentinel-2 Images Raw  
│   ├── Site 1  
│   │   ├── L1C  
│   │   ├── L2A  
│   ├── Site 2  
│   ├── ...  
├── Sentinel-2 Images Processed  
│   ├── Site 1  
│   │   ├──   
│   ├── Site 2  
├── Results Merged
│   ├── Roman Metrics.csv  
│   ├── Variogram Parameters.csv  
│   ├── ...  

# References
1. Román, M.O., Schaaf, C.B., Woodcock, C.E., Strahler, A.H., Yang, X., Braswell, R.H., Curtis, P.S., Davis, K.J., Dragoni, D., Goulden, M.L., Gu, L., Hollinger, D.Y., Kolb, T. E., Meyers, T.P., Munger, J.W., Privette, J.L., Richardson, A.D., Wilson, T.B., Wofsy, S.C., 2009. The MODIS (collection v005) BRDF/albedo product: assessment of spatial representativeness over forested landscapes. Remote Sens. Environ. 113, 2476–2498. https://doi.org/10.1016/j.rse.2009.07.009.