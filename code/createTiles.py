## from https://github.com/DPIRD-DMA
# %% Import Packages
import os
import subprocess
from pathlib import Path

import math
from pathlib import Path
from tqdm.auto import tqdm
from multiprocess import Pool,cpu_count
%matplotlib inline

from osgeo import gdal,osr
import geopandas as gpd
from shapely import geometry
import pandas as pd
import numpy as np
import json

import rioxarray

import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
# %% Functions
def get_bounds(tif_path):
    # open file
    data = gdal.Open(tif_path)
    # grab bounds
    geoTransform = data.GetGeoTransform()
    left = geoTransform[0]
    top = geoTransform[3]
    right = left + geoTransform[1] * data.RasterXSize
    bottom = top + geoTransform[5] * data.RasterYSize
    # build dict to file bounds
    geo_tiff_bounds_dict = {'top':top,'left':left,'bottom':bottom,'right':right,'tif_path':tif_path}
    return geo_tiff_bounds_dict
def intersect_tile_with_geotiffs(tile_dict,geo_tiff_bounds):
    # setup set to collect results in, a set is used to avoid duplicates
    intersecting_geotiffs = set()
    # loop over each geotiff
    for geo_bounds in geo_tiff_bounds:
        # check is tile top or bottom is inside geotiff
        if (geo_bounds['top'] > tile_dict['top'] > geo_bounds['bottom'] or 
            geo_bounds['top'] > tile_dict['bottom'] > geo_bounds['bottom']):
            # check if left or right are inside a geotiff
            if geo_bounds['right'] > tile_dict['left'] > geo_bounds['left']:
                intersecting_geotiffs.add(geo_bounds['tif_path'])
            if geo_bounds['right'] > tile_dict['right'] > geo_bounds['left']:
                intersecting_geotiffs.add(geo_bounds['tif_path'])
    return intersecting_geotiffs
def make_polygons(row):
    tile_polygon_list = []
    tile_top = bound_y_max + y_move*row
    tile_bottom = tile_top + y_tile_size
    tile_left = bound_x_min

    for col in range(0,number_of_cols):
        tile_left = bound_x_min + col*x_move
        tile_right = tile_left + x_tile_size
        tile_dict = {'top':tile_top,'left':tile_left,'bottom':tile_bottom,'right':tile_right}
        tile_list = np.array([tile_top,tile_left,tile_bottom,tile_right])
        # check if valid tile
        intersect = intersect_tile_with_geotiffs(tile_dict,geo_tiff_bounds)
        raster_name = raster_name_append+str(row)+'_'+str(col)+'.tif'
        if len(intersect) > 0:
            polygon = {'geometry':geometry.Polygon([[tile_left, tile_top], [tile_right, tile_top], [tile_right, tile_bottom], [tile_left, tile_bottom]]),
                      'intersect':intersect, 'row':row, 'col':col, 'name':raster_name}
            tile_polygon_list.append(polygon)
    return tile_polygon_list
def intersector(geo_tiff):
    tiles_inside_geo_tiff = []
    # loop over each tile and check if the geotiff is the the intersect list
    for tile in tile_polygon_list:
        if geo_tiff in tile['intersect']:
            # count this so we know if the tile will be incomplete or not
            incomplete = len(tile['intersect'])>1
            # build dict with geom the current row and col for naming
            tiles_inside_geo_tiff.append({'geometry':tile['geometry'],'row':tile['row'],'col':tile['col'],'name':tile['name'],'incomplete':incomplete})
    return([geo_tiff,tiles_inside_geo_tiff])  
def cut_tiles(geotiff):
    # grab path to to file and open it
    geotiff_open = gdal.Open(geotiff[0])
    # grab the filename and strip the extension
    geo_tiff_filename = os.path.basename(geotiff[0]).replace(input_file_ext,'')
    incomplete_tile_list = []
    for tile in geotiff[1]:
        tile_geometry = tile['geometry']
        # shapely bounds returns "minx, miny, maxx, maxy" but we need minx, maxy, maxx, miny
        top = list(tile_geometry.bounds)[3]
        bottom = list(tile_geometry.bounds)[1]
        left = list(tile_geometry.bounds)[0]
        right =list(tile_geometry.bounds)[2]
        
        # make row folder path
        if output_to_row_folders:
            output_row_folder = os.path.join(output_folder,str(tile['row']))
            Path(output_row_folder).mkdir(parents=True, exist_ok=True)
        else:
            output_row_folder = output_folder
        # make row folder if necessary
        
        export_file_name = str(tile['name'])#str(tile['row'])+'_'+str(tile['col'])+'.tif'
        
        # check if tile is incomplete if so append the getiff name so that it is unique
        if tile['incomplete']:
            append_name = '-'+geo_tiff_filename+'_incomplete.tif'
            export_file_name = export_file_name.replace('.tif',append_name)
            # add tile to list so we dont need to refind them to compile incomplete tiles
            export_file_path = os.path.join(output_row_folder,export_file_name)
            incomplete_tile_list.append(export_file_path)
        else:
            export_file_path = os.path.join(output_row_folder,export_file_name)
        
        # check if already done
        if not os.path.isfile(export_file_path):

            # clip the data
            # make a string of tile dims to pass as a command line arg, this is kinda of hacky, would like a better option
            tile_clip_string = str(left) +' '+str(top) +' '+str(right) +' '+str(bottom)

            translate_options = gdal.TranslateOptions(gdal.ParseCommandLine("-projwin "+tile_clip_string)
                                                     ,creationOptions=['COMPRESS='+output_compression])

            tile_clip = gdal.Translate(export_file_path, geotiff_open, options = translate_options)
            # close the tile
            tile_clip = None
    return incomplete_tile_list
def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))
# %%
tile_size_px = [1024,1024]
tile_oxerlap_px = 0
input_file_ext = '.tif'
output_compression = 'LZW' 
output_to_row_folders = False
raster_name_append = ''
settings_file_path = os.path.join('/Users/joshuadimasaka/Desktop/PhD/GitHub/DeepTemporal/data/TILES/tiling_settings.json')
settings = {'tile_size_px':tile_size_px,
           'tile_oxerlap_px':tile_oxerlap_px,
           'raster_name_append':raster_name_append,
           'tile_format':'.tif'}
with open(settings_file_path, 'w') as fp:
    json.dump(settings, fp)

# geotiff_folder = os.path.join(parent_dir,'data/EMI/Geospatial References/Barangay Boundaries')
parent_dir = os.path.dirname(os.getcwd())
geo_tiff_list = [   os.path.join(parent_dir,'data/EMI/Geospatial References/Barangay Boundaries/brgyBoundary.tif'),
                    os.path.join(parent_dir,'data/openBLDGtemporal/2016_v3.tif'),
                    os.path.join(parent_dir,'data/openBLDGtemporal/2017_v3.tif'),
                    os.path.join(parent_dir,'data/openBLDGtemporal/2018_v3.tif'),
                    os.path.join(parent_dir,'data/openBLDGtemporal/2019_v3.tif'),
                    os.path.join(parent_dir,'data/openBLDGtemporal/2020_v3.tif'),
                    os.path.join(parent_dir,'data/openBLDGtemporal/2021_v3.tif'),
                    os.path.join(parent_dir,'data/openBLDGtemporal/2022_v3.tif'),
                    os.path.join(parent_dir,'data/openBLDGtemporal/2023_v3.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_W3.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_W2.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_W1.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_URM.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_URA.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_S4.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_S3.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_S2.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_S1.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_RM2.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_RM1.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_PC2.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_N.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_MWS.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_CWS.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_CHB.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_C4.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_C2.tif'),
                    os.path.join(parent_dir,'data/GMMA_Exposure_Compilation_August2013/prior_C1.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_worship.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_school.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_multihall.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_mrf.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_market.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_mall.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_hotel.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_hospital.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_healthcenter.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_gas.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_footbridges.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_factory.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_evacuation.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_daycare.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_coveredcourt.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_chemfactory.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_bridge.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_points_distance_from_brgyhall.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_lines_distance_from_uvroute.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_lines_distance_from_roadsecondary.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_lines_distance_from_roadprimary.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_lines_distance_from_roadothers.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_lines_distance_from_railway.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_lines_distance_from_jeepneyroute.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_lines_distance_from_citybusroute.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_lines_distance_from_bikelane.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_areas_distance_from_waterways.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_areas_distance_from_vacant.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_areas_distance_from_utility.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_areas_distance_from_socializedhousing.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_areas_distance_from_roads.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_areas_distance_from_residential.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_areas_distance_from_reservoir.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_areas_distance_from_recreational.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_areas_distance_from_openspace.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_areas_distance_from_military.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_areas_distance_from_institutional.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_areas_distance_from_informalsettlement.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_areas_distance_from_industrial.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_areas_distance_from_commercial.tif'),
                    os.path.join(parent_dir,'data/EMI/Built Environment/utility_areas_distance_from_cemetery.tif'),
                    os.path.join(parent_dir,'data/EMI/Hazards Vulnerability and Risk/Baseline Flood Simulation (QC DMP)/Baseline Flood Simulation - Flow Depth/100yr_flood_proximity_greaterthan_3m.tif'),
                    os.path.join(parent_dir,'data/EMI/Hazards Vulnerability and Risk/Baseline Flood Simulation (QC DMP)/Baseline Flood Simulation - Flow Depth/100yr_flood_proximity_greaterthan_1.5m.tif'),
                    os.path.join(parent_dir,'data/EMI/Hazards Vulnerability and Risk/Baseline Flood Simulation (QC DMP)/Baseline Flood Simulation - Flow Depth/100yr_flood_proximity_greaterthan_0.5m.tif'),
                    os.path.join(parent_dir,'data/EMI/Hazards Vulnerability and Risk/Baseline Flood Simulation (QC DMP)/Baseline Flood Simulation - Flow Depth/100yr_flood_proximity_greaterthan_0.2m.tif')]

# %% create folder
output_folder_list = []
for i in range(len(geo_tiff_list)):
    basename = os.path.basename(geo_tiff_list[i])
    output_folder = os.path.join(parent_dir,'data/TILES/',basename[:-4])
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    output_folder_list.append(output_folder)

# %% Valid for brgy boundary layer only   
with Pool() as pool:
    geo_tiff_bounds = list(tqdm(pool.imap(get_bounds, [geo_tiff_list[0]]), total=len([geo_tiff_list[0]])))
pure_bounds = []
for geo_tif_bounds in geo_tiff_bounds:
    pure_bounds.append([geo_tif_bounds['top'],geo_tif_bounds['left'],geo_tif_bounds['bottom'],geo_tif_bounds['right']])
pure_bounds_np = np.array(pure_bounds)
bound_y_max = float(pure_bounds_np[:,0].max()) #top
bound_x_min = float(pure_bounds_np[:,1].min()) #left
bound_y_min = float(pure_bounds_np[:,2].min()) #bottom
bound_x_max = float(pure_bounds_np[:,3].max()) #right

test_raster = gdal.Open(geo_tiff_list[0])
test_raster_gt =test_raster.GetGeoTransform()
pixel_size_x = test_raster_gt[1]
pixel_size_y = test_raster_gt[5]
proj = osr.SpatialReference(wkt=test_raster.GetProjection())
crs = 'EPSG:'+proj.GetAttrValue('AUTHORITY',1)

# calculate the geographical distance in each direction each tile must be from the last tile
x_move = pixel_size_x*(tile_size_px[0]-tile_oxerlap_px)
y_move = pixel_size_y*(tile_size_px[1]-tile_oxerlap_px)

# calculate the geographical size of each tile
x_tile_size = pixel_size_x*tile_size_px[0]
y_tile_size = pixel_size_y*tile_size_px[1]

# calculate the number of cols and rows so we can avoid using while loops
number_of_cols = math.ceil(abs((bound_x_max-bound_x_min)/x_move))
number_of_rows = math.ceil(abs((bound_y_max-bound_y_min)/y_move))

with Pool() as pool:
    tile_polygon_list = list(tqdm(pool.imap(make_polygons, range(0,number_of_rows)), total=len(range(0,number_of_rows))))

# %% this is returned as a list of list so it must be flattened
tile_polygon_list_old = list(np.concatenate(tile_polygon_list).ravel())

brgyboundarypath = geo_tiff_list[0]

tile_polygon_list = []
for i in range(len(tile_polygon_list_old)):
    polygon = tile_polygon_list_old[i]['geometry']
    xds = rioxarray.open_rasterio(brgyboundarypath, masked=True).rio.clip(geometries=[polygon], crs=4326, from_disk=True)
    if np.divide(np.count_nonzero(np.isnan(xds.to_numpy())), 1024*1024) <= 0.001:
        tile_polygon_list.append(tile_polygon_list_old[i])

%%time
#  convert into geodataframe
polygon_tiles_gpd = gpd.GeoDataFrame(tile_polygon_list,geometry='geometry',crs=crs)
del polygon_tiles_gpd['intersect']

polygon_export_loc = '/Users/joshuadimasaka/Desktop/PhD/GitHub/DeepTemporal/data/TILES/output_no_overlap.gpkg'
polygon_tiles_gpd.to_file(polygon_export_loc, driver="GPKG")   

# %%
for i in range(len(output_folder_list)):

    output_folder = output_folder_list[i]
    print(i, os.path.basename(output_folder))

    with Pool() as pool:
        geo_tiff_with_tiles = list(tqdm(pool.imap(intersector, [geo_tiff_list[0]]), total=len([geo_tiff_list[0]])))
    geo_tiff_with_tiles[0][0] = geo_tiff_list[i]

    concurrent_threads = cpu_count()
    tile_count = 0
    for raster,tiles in geo_tiff_with_tiles:
        tile_count+=len(tiles)
    average_tile_count = round(tile_count/len(geo_tiff_with_tiles))
    print(average_tile_count)

    # if you have more cores than rasters and more cores than average tile count then give every core a copy of the same raster to work on
    if concurrent_threads > len(geo_tiff_with_tiles) and concurrent_threads < average_tile_count:
        print('You have more threads than rasters so we each thread will get the same raster')
        incomplete_tile_list = []

        for raster,tiles in tqdm(geo_tiff_with_tiles):
            geo_tiff_with_tile_list = []

            for tile_split in split(tiles,concurrent_threads):
                geo_tiff_with_tile_list.append([raster,tile_split])

            pool = Pool(concurrent_threads)
            with pool:
                one_raster_incomplete_tile_list = list(pool.imap(cut_tiles,geo_tiff_with_tile_list))

            incomplete_tile_list.append(np.concatenate(one_raster_incomplete_tile_list))
        
    else:
        print('you have more rasters than threads so each thread will get its own raster')
        pool = Pool()
        with pool:
            incomplete_tile_list = list(tqdm(pool.imap(cut_tiles,geo_tiff_with_tiles), total=len(geo_tiff_with_tiles)))
