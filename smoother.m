clear all
close all
addpath functions


%% DEFINE GLACIER and Select DEMs for coregistration
if ispc
    slash='\';
elseif isunix
    slash='/';
end
%search data directory for available DEMs
data_directory='C:\Users\cmcneil\Desktop\GIS\Geodetics\';
% data_directory=['Q:',slash,'Project Data',slash,'GlacierData',slash,'GIS',slash,];
d = dir(data_directory);
str = {d.name};
[s,v] = listdlg('PromptString','Select Master DEM for coregistration:',...
                'SelectionMode','single',...
                'ListString',str);
glacier=d(s).name;
raw_DEM_directory=[data_directory,glacier,slash,'DEMs',slash,'Raw'];%set path for raw DEMs
processed_DEM_directory=[data_directory,glacier,slash,'DEMs',slash,'Processed'];%set path for processed DEMs
%main_directory=['C:\Users\cmcneil\Desktop\Matlab\Taku and LemonCreeek glacier reanalysis code\2017.06.12\data\',glacier,'\DEMs\Raw'];%set path for raw DEMs
d = dir(raw_DEM_directory);
str = {d.name};
[s,v] = listdlg('PromptString','Select DEM for smoothing:',...
                'SelectionMode','single',...
                'ListString',str);
dem=d(s).name;
dem_path=[raw_DEM_directory,slash,'',dem,slash,dem,'_DEM.tif'];
coregistration_area_date=datestr(datenum(dem,'yyyy.mm.dd'),'yyyymmdd');
mask_path=[data_directory,glacier,'\Layers.gdb\coregistration_areas\coregistration_area_',coregistration_area_date];
file=dir(dem_path);

if file.bytes>1e+8
    
    feature_path=mask_path;
    dem_dir=[raw_DEM_directory,slash,'',dem];
    py.clip.clipRaster(feature_path,dem_dir)
    temp_dem_path=[raw_DEM_directory,slash,'',dem,slash,'temp.tif'];
    DEM = read_dem(temp_dem_path);
else
    DEM = read_dem(dem_path);
end
% DEM=read_dem(dem_path);
smoothing_radius=50;
n=round(smoothing_radius/DEM.cellsize);
slidingmean=conv2(DEM.Z,ones(n)/n^2,'same');
slidingstd=sqrt( conv2(DEM.Z.^2,ones(3)/9,'same') - slidingmean.^2 );

hillshade(DEM.Z,DEM.X(1,:),DEM.Y(:,1),'azimuth',45,'altitude',100,'plotit');

hillshade(slidingmean,DEM.X(1,:),DEM.Y(:,1),'azimuth',45,'altitude',100,'plotit');
% B=inpaint_nans(slidingmean,0)
% hillshade(B,DEM.X(1,:),DEM.Y(:,1),'azimuth',45,'altitude',100,'plotit');
