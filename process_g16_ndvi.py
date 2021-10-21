#-----------------------------------------------------------------------------------------------------------
# Agrometeorology Course - Demonstration 2: Normalized Difference Vegetation Index (NDVI)
# Author: Diego Souza (INPE / Brazil)
#-----------------------------------------------------------------------------------------------------------
# Required modules
from netCDF4 import Dataset                    # Read / Write NetCDF4 files
import matplotlib.pyplot as plt                # Plotting library
from datetime import datetime                  # Basic Dates and time types
from datetime import timedelta, date, datetime # Basic Dates and time types
import cartopy, cartopy.crs as ccrs            # Plot maps
import cartopy.io.shapereader as shpreader     # Import shapefiles
import numpy as np                             # Scientific computing with Python
import os                                      # Miscellaneous operating system interfaces
from utilities import download_CMI             # Our own utilities
from utilities import geo2grid, latlon2xy, convertExtent2GOESProjection      # Our own utilities
#-----------------------------------------------------------------------------------------------------------
# Input and output directories
input = "/content/Samples"; os.makedirs(input, exist_ok=True)
output = "/content/Animation"; os.makedirs(output, exist_ok=True)

# Desired extent
extent = [-93.0, -60.00, -25.00, 15.00] # Min lon, Min lat, Max lon, Max lat

# Initial date and time to process
yyyymmddhhmn = '202107151800'

# Number of days to accumulate
ndays = 15

# Initial time and date
yyyy = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%Y')
mm = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%m')
dd = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%d')
hh = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%H')
mn = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%M')
date_ini = str(datetime(int(yyyy),int(mm),int(dd),int(hh),int(mn)))

# Choose number of days ahead
date_end = str(datetime(int(yyyy),int(mm),int(dd),int(hh),int(mn)) + timedelta(days=ndays))

#-----------------------------------------------------------------------------------------------------------
# NDVI colormap creation 
import matplotlib
colors = ["#8f2723", "#8f2723", "#8f2723", "#8f2723", "#af201b", "#af201b", "#af201b", "#af201b", "#ce4a2e", "#ce4a2e", "#ce4a2e", "#ce4a2e", 
          "#df744a", "#df744a", "#df744a", "#df744a", "#f0a875", "#f0a875", "#f0a875", "#f0a875", "#fad398", "#fad398", "#fad398", "#fad398",
          "#fff8ba",
          "#d8eda0", "#d8eda0", "#d8eda0", "#d8eda0", "#bddd8a", "#bddd8a", "#bddd8a", "#bddd8a", "#93c669", "#93c669", "#93c669", "#93c669", 
          "#5da73e", "#5da73e", "#5da73e", "#5da73e", "#3c9427", "#3c9427", "#3c9427", "#3c9427", "#235117", "#235117", "#235117", "#235117"]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
cmap.set_over('#235117')
cmap.set_under('#8f2723')
vmin = 0.1
vmax = 0.8
#-----------------------------------------------------------------------------------------------------------

# Variable to check if is the first iteration
first = True

while (date_ini <= date_end):

    # Date structure
    yyyymmddhhmn = datetime.strptime(date_ini, '%Y-%m-%d %H:%M:%S').strftime('%Y%m%d%H%M')
    print(yyyymmddhhmn)

    # Download the file for Band 02
    file_name_ch02 = download_CMI(yyyymmddhhmn, '2', input)

    # Download the file for Band 03
    file_name_ch03 = download_CMI(yyyymmddhhmn, '3', input)
    
    #-----------------------------------------------------------------------------------------------------------
    
    # If one the files hasn't been downloaded for some reason, skip the current iteration
    if not (os.path.exists(f'{input}/{file_name_ch02}.nc')) or not (os.path.exists(f'{input}/{file_name_ch03}.nc')):
      # Increment the date_ini
      date_ini = str(datetime.strptime(date_ini, '%Y-%m-%d %H:%M:%S') + timedelta(days=1))
      print("The file is not available, skipping this iteration.")
      continue

    # Open the GOES-R image (Band 02)
    file = Dataset(f'{input}/{file_name_ch02}.nc')
                      
    # Convert lat/lon to grid-coordinates
    lly, llx = geo2grid(extent[1], extent[0], file)
    ury, urx = geo2grid(extent[3], extent[2], file)
            
    # Get the pixel values
    data_ch02 = file.variables['CMI'][ury:lly, llx:urx][::2 ,::2]   
    
    #-----------------------------------------------------------------------------------------------------------
    
    # Open the GOES-R image (Band 03)
    file = Dataset(f'{input}/{file_name_ch03}.nc')
                      
    # Convert lat/lon to grid-coordinates
    lly, llx = geo2grid(extent[1], extent[0], file)
    ury, urx = geo2grid(extent[3], extent[2], file)
            
    # Get the pixel values
    data_ch03 = file.variables['CMI'][ury:lly, llx:urx]      
    
    #-----------------------------------------------------------------------------------------------------------
    
    # Make the arrays equal size
    cordX = np.shape(data_ch02)[0], np.shape(data_ch03)[0]
    cordY = np.shape(data_ch02)[1], np.shape(data_ch03)[1]

    minvalX = np.array(cordX).min()
    minvalY = np.array(cordY).min()

    data_ch02 = data_ch02[0:minvalX, 0:minvalY]
    data_ch03 = data_ch03[0:minvalX, 0:minvalY]
    
    #-----------------------------------------------------------------------------------------------------------
    
    # Calculate the NDVI
    data = (data_ch03 - data_ch02) / (data_ch03 + data_ch02)
    
    #-----------------------------------------------------------------------------------------------------------
    
    # Compute data-extent in GOES projection-coordinates
    img_extent = convertExtent2GOESProjection(extent)
    
    #-----------------------------------------------------------------------------------------------------------
    
    # If it's the first iteration, create the array that will store the max values
    if (first == True):
      first = False
      ndvi_max = np.full((data_ch02.shape[0],data_ch02.shape[1]),-9999)
    
    # Keep the maximuns, ignoring the nans
    ndvi_max = np.fmax(data,ndvi_max)
    # Remove low values from the accumulation
    #ndvi_max[ndvi_max < 0.1] = np.nan
    
    # Choose the plot size (width x height, in inches)
    plt.figure(figsize=(10,10))

    # Use the Geostationary projection in cartopy
    ax = plt.axes(projection=ccrs.Geostationary(central_longitude=-75.0, satellite_height=35786023.0))
   
    # Add a land mask
    ax.stock_img()       
    import cartopy.feature as cfeature
    land = ax.add_feature(cfeature.LAND, facecolor='gray', zorder=1)
  
    # Plot the image
    img = ax.imshow(ndvi_max, origin='upper', vmin=vmin, vmax=vmax, extent=img_extent, cmap=cmap, zorder=2)

    # Add an ocean mask
    ocean = ax.add_feature(cfeature.OCEAN, facecolor='lightsteelblue', zorder=3)

    # Add a shapefile
    shapefile = list(shpreader.Reader('ne_10m_admin_1_states_provinces.shp').geometries())
    ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='dimgray',facecolor='none', linewidth=0.3, zorder=4)

    # Add coastlines, borders and gridlines
    ax.coastlines(resolution='10m', color='black', linewidth=0.8, zorder=5)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5, zorder=6)
    ax.gridlines(color='gray', alpha=0.5, linestyle='--', linewidth=0.5, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), zorder=7)

    # Add a colorbar
    plt.colorbar(img, label='Normalized Difference Vegetation Index', extend='both', orientation='vertical', pad=0.05, fraction=0.05)

    # Extract date
    date = (datetime.strptime(file.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ'))

    # Add a title
    plt.title('GOES-16 NDVI ' + date.strftime('%Y-%m-%d %H:%M') + ' UTC', fontweight='bold', fontsize=10, loc='left')
    plt.title('Reg.: ' + str(extent) , fontsize=10, loc='right')
    #-----------------------------------------------------------------------------------------------------------
    # Save the image
    plt.savefig(f'{output}/G16_NDVI_{yyyymmddhhmn}.png', bbox_inches='tight', pad_inches=0, dpi=100)

    # Show the image
    plt.show()

    # Increment the date_ini
    date_ini = str(datetime.strptime(date_ini, '%Y-%m-%d %H:%M:%S') + timedelta(days=1))