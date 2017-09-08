function [] = shift_rasterxy(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dbstop if error
if length(varargin)<4
    error('Need all input arguments')
else
%% Shift slave ortho image with DEM if it exists
shift_x=cell2mat(varargin(1));
shift_y=cell2mat(varargin(2));
input_raster=cell2mat(varargin(3));
ouput_raster=cell2mat(varargin(4));


Tinfo=imfinfo(input_raster);
copyfile(input_raster,ouput_raster)
t=Tiff(ouput_raster,'r+');
tagstructure.ModelTiepointTag=[Tinfo(1).ModelTiepointTag(1:3) Tinfo(1).ModelTiepointTag(4)+shift_x Tinfo(1).ModelTiepointTag(5)+shift_y Tinfo(1).ModelTiepointTag(6)];
t.setTag(tagstructure)
t.close
end
end

