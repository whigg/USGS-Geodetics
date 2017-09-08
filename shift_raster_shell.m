clear all
close all
addpath functions


X=-2.6239;
Y=2.9789;
Z=0;
input_raster='Q:\Project Data\GlacierData\GIS\Wolverine\Orthos\Raw\2013.06.15\2013.06.15_ortho.tif';
ouput_raster='Q:\Project Data\GlacierData\GIS\Wolverine\Orthos\Processed\2013.06.15\2013.06.15_ortho.tif';
ouput_dir='Q:\Project Data\GlacierData\GIS\Wolverine\Orthos\Processed\2013.06.15\';
if ~exist(ouput_dir, 'dir')
    % Folder does not exist so create it.
    mkdir(ouput_dir);
end
    shift_rasterxy(X,Y,input_raster,ouput_raster)
