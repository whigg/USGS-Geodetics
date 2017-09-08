# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 15:59:48 2017

@author: cmcneil
"""

# Import arcpy module
import arcpy,os
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True
#feature_path="Q:\Project Data\GlacierData\GIS\Gulkana\Layers.gdb\coregistration_areas\coregistration_area_19570809"
#raster_dir="Q:\Project Data\GlacierData\GIS\Gulkana\Orthos\Raw\\2016.08.30"
#cell_size="5"
def clip(feature_path,raster_dir):
    '''function for clipping raster to shapefile extent'''
    temp_raster=raster_dir+"\\temp.tif"     
    if arcpy.Exists(temp_raster):
           arcpy.Delete_management(temp_raster)
    for file in os.listdir(raster_dir):
        if file.endswith(".tif"):
            raster=raster_dir+"\\"+file

    envelope = "in_memory\\envelope"
    buffer = "in_memory\\buffer"
    envelope=arcpy.FeatureEnvelopeToPolygon_management(feature_path,envelope, "SINGLEPART")
    buffer=arcpy.Buffer_analysis(envelope, buffer, 500, "FULL", "ROUND", "NONE", "", "PLANAR")
    rows = arcpy.SearchCursor(buffer)
    shapeName = arcpy.Describe(buffer).shapeFieldName
    for row in rows:
        feat = row.getValue(shapeName)
        extent = feat.extent
        bounds = str(extent.XMin) + ' ' + str(extent.YMin) + ' '+ str(extent.XMax) + ' ' + str(extent.YMax)
    arcpy.Clip_management(raster,bounds,temp_raster)
    arcpy.Delete_management("in_memory")
    return

def resample(cellsize,raster_dir):
    temp_raster=raster_dir+"\\temp.tif"     
    if arcpy.Exists(temp_raster):
           arcpy.Delete_management(temp_raster) 
    for file in os.listdir(raster_dir):
        if file.endswith(".tif"):
            raster=raster_dir+"\\"+file 
        
    resampled_cellsize=str(cellsize) +' '+ str(cellsize)
    arcpy.Resample_management(raster, temp_raster, resampled_cellsize, "BILINEAR")
    
def resample_and_clip(cellsize,feature_path,raster_dir):
    temp_raster=raster_dir+"\\temp.tif"
    if arcpy.Exists(temp_raster):
       arcpy.Delete_management(temp_raster)
    for file in os.listdir(raster_dir):
        if file.endswith(".tif"):
            raster=raster_dir+"\\"+file     
    resampled_raster="in_memory\\resampled"
    resampled_cellsize=str(cellsize) + ' ' + str(cellsize)
    arcpy.Resample_management(raster, resampled_raster, resampled_cellsize, "BILINEAR")  
    envelope = "in_memory\\envelope"
    buffer = "in_memory\\buffer"
    envelope=arcpy.FeatureEnvelopeToPolygon_management(feature_path,envelope, "SINGLEPART")
    buffer=arcpy.Buffer_analysis(envelope, buffer, 500, "FULL", "ROUND", "NONE", "", "PLANAR")
    rows = arcpy.SearchCursor(buffer)
    shapeName = arcpy.Describe(buffer).shapeFieldName
    for row in rows:
        feat = row.getValue(shapeName)
        extent = feat.extent
        envelope = str(extent.XMin) + ' ' + str(extent.YMin) + ' '+ str(extent.XMax) + ' ' + str(extent.YMax)
    arcpy.Clip_management(resampled_raster,envelope,temp_raster)
    arcpy.Delete_management("in_memory")
    return

def getmask(raster_path,mask_path,output_path):
        '''function to cinvert shapfeild to mask for dem'''
        raster=raster_path
        snap_tif = "in_memory\\snap.tif"
        snap_shp = "in_memory\\snap.shp"        
        nan_mask_shp = "in_memory\\nan_mask.shp"
        nan_mask_shp__3_ = nan_mask_shp
        coregistration_areas = mask_path        
        ca_tif = output_path+"mask.tif"        
        v_ca = coregistration_areas
        if arcpy.Exists(ca_tif):
                   arcpy.Delete_management(ca_tif)

        arcpy.gp.CreateConstantRaster_sa(snap_tif, "0", "INTEGER", DEM, DEM)
        arcpy.RasterToPolygon_conversion(snap_tif, snap_shp, "SIMPLIFY", "Value")
        arcpy.Erase_analysis(snap_shp, v_ca, nan_mask_shp)
        arcpy.DefineProjection_management(nan_mask_shp, "PROJCS['WGS_1984_UTM_Zone_8N',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-135.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]")
        arcpy.PolygonToRaster_conversion(nan_mask_shp__3_, "FID", ca_tif, "MAXIMUM_AREA", "NONE", raster)
        arcpy.Delete_management("in_memory")
        return
    