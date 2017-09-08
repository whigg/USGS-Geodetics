
close all
clear all

addpath functions

%% DEFINE GLACIER and Select DEMs for to stake and fill. Need to have at least two DEMs covering the full extend of the glacier
% because PC < all the rest
if ispc
    slash='\';
elseif isunix
    slash='/';
end
%silly computers

%define the root directory of DEMs
data_directory='C:\Users\cmcneil\Desktop\GIS\Geodetics\' 
%data_directory=['Q:',slash,'Project Data',slash,'GlacierData',slash,'GIS',slash,];
d = dir(data_directory);                                    %list all glaciers with data in directory
str = {d.name}; 
[s,v] = listdlg('PromptString','Select Master DEM for coregistration:',...
                'SelectionMode','single',...
                'ListString',str);                          %prompt dialog box for users to select which glacier to work on
glacier=d(s).name;
DEM_directory=[data_directory,glacier,slash,'DEMs',slash,'Processed'];%set path for processed DEMs
%main_directory=['C:\Users\cmcneil\Desktop\Matlab\Taku and LemonCreeek glacier reanalysis code\2017.06.12\data\',glacier,'\DEMs\Raw'];%set path for raw DEMs
d = dir(DEM_directory);                                     %list all processed (coregistered DEMs for directory/glacier
dems = {d.name};
[s,v] = listdlg('PromptString','Select DEMs to stack',...
                'SelectionMode','multiple',...
                'ListString',dems);                         %prompt dialog box for user to select DEMs to stack and interpolate
%%             
% for i=1:length(s)                                           %for all DEMs selected
%     raster_path=[DEM_directory,'\',cell2mat(dems(s(i))),'\',cell2mat(dems(s(i))),'_DEM.tif']; %set DEM path
%     dem_year=datestr(datenum(dems(s(i)),'yyyy.mm.dd'),'yyyy');%get year of DEM
%     mask_path=[data_directory,glacier,'\Layers.gdb\boundaries\C',dem_year,'_boundary'];%set path for glacier boundary in Layer.gdb
%     output_path=[data_directory,glacier,slash];             %define output path for python function that converts the boundary shapefile to a raster mask
%     py.Raster.getmask(raster_path,mask_path,output_path)         %python function that converts the boundary shapefile to a raster mask
%     mask_path=[output_path,'mask.tif'];                     %set path to raster mask of glacier boundary
%     boundary(i)=read_dem(mask_path); %import raster mask for area calculations
%     delete([output_path,'mask*'])
%     boundary(i).Z(isnan(boundary(i).Z))=1;                  %set raster mask to binary double for simplicity
%     cellsizes(i,1)=boundary(i).cellsize;
%     area(i,1)=(sum(boundary(i).Z(:)).*cellsizes(i,1)^2)./1000000;%determine glacier area for each bounadary                    
%     data_years(i,1)=str2num(dem_year);                      %year of dem data aquisition
% end

%%
area_table_path=[data_directory,glacier,'\Tables\',glacier,'Area.csv'];
data=importdata(area_table_path);
area_table=data.data;
[~,max_area_index]=max(area); %find max area of glacier during study interval
area=[];
%%
    for i=1:length(s)
        raster_path=[DEM_directory,'\',cell2mat(dems(s(i))),'\',cell2mat(dems(s(i))),'_DEM.tif'];%set DEM path
        raster_paths{i}=raster_path;                              %stash paths
        dem_year=datestr(datenum(dems(s(max_area_index)),'yyyy.mm.dd'),'yyyy');%get year of dem data aquistion
        data_years(i,1)=str2num(dem_year);  
        mask_path=[data_directory,glacier,'\Layers.gdb\boundaries\C',dem_year,'_boundary'];%set path for DEM
        output_path=[data_directory,glacier,slash];         %set output path for python function that converts shapefile to raster mask
        py.Raster.getmask(raster_path,mask_path,output_path)     %python function that converts shapefile to raster mask
        mask_path=[output_path,'mask.tif'];                 %set path for raster mask
        DEM(i) = read_dem(dem_path,mask_path);%import DEM
        delete([output_path,'mask*'])
        cellsizes(i,1)=DEM(i).cellsize;
        area(i,1)=area_table(area_table(:,1)==str2num(dem_year),2);%determine glacier area for each bounadary                    
        data_years(i,1)=str2num(dem_year);
    end
    

    [~,max_cellsize_index]=max(cellsizes); %find max area of glacier during study interval
    %define max glacier extent coordinates based on maximum glacier extent
    extent=[max(DEM(max_area_index).X((~isnan(DEM(max_area_index).Z))))+100 min(DEM(max_area_index).X((~isnan(DEM(max_area_index).Z))))-100 maxDEM(max_area_index).Y((~isnan(DEM(max_area_index).Z)))+100 min(DEM(max_area_index).Y((~isnan(DEM(max_area_index).Z))))-100];
    %Range of Eastings   
    Xs=DEM(max_cellsize_index).X(1,DEM(max_cellsize_index).X(1,:) <= extent(1,1) & DEM(max_cellsize_index).X(1,:) >= extent(1,2));
    %Range of Northings
    Ys=DEM(max_cellsize_index).Y(DEM(max_cellsize_index).Y(:,1) <= extent(1,3) & DEM(max_cellsize_index).Y(:,1) >= extent(1,4),1);
    %for each glacier
    
    for i=1:length(DEM)
        stack(i).X=repmat(Xs,length(Ys),1);                 %eastings for resampling
        stack(i).Y=repmat(Ys,1,length(Xs));                 %northings for resampling
        stack(i).Z=interp2(DEM(i).X,DEM(i).Y,DEM(i).Z_masked,stack(i).X,stack(i).Y, 'linear'); %resample DEM to the 
    end
%%

dz_1=stack(3).Z-stack(1).Z;
dz_2=stack(2).Z-stack(1).Z;
pixel_ratios=dz_2./dz_1;
dem_dz_ratio=nanmedian(abs(pixel_ratios(abs(pixel_ratios)<1)));

dz_4=stack(1).Z+(dz_1.*dem_dz_ratio);
DEM(2).Z(isnan(DEM(2).Z) & boundary(2).Z==1)=interp2(stack(1).X,stack(1).Y,dz_4,DEM(2).X(isnan(DEM(2).Z) & boundary(2).Z==1),DEM(2).Y(isnan(DEM(2).Z) & boundary(2).Z==1),'linear');

figure();hold on
imshow(DEM(2).Z,[])
colormap(jet)
colorbar

 original_dem=cell2mat(dem_paths(2));
 interpolated_dem=[original_dem(1:end-18),'interpolated',original_dem(end-18:end)];
 if ~exist([original_dem(1:end-18),'interpolated',slash], 'dir')
    % Folder does not exist so create it.
    mkdir([original_dem(1:end-18),'interpolated',slash]);
end
copyfile(original_dem,interpolated_dem)
image_info=DEM(2).Tinfo;
t=Tiff(interpolated_dem,'r+');
tagstructure.ModelTiepointTag=[image_info(1).ModelTiepointTag(1:3) image_info(1).ModelTiepointTag(4) image_info(1).ModelTiepointTag(5) image_info(1).ModelTiepointTag(6)];
tagstructure.BitsPerSample = 64;
tagstructure.SampleFormat = 3;
t.setTag(tagstructure)
t.write(DEM(2).Z);
t.close 


% [FileName,PathName,FilterIndex] = uigetfile(FilterSpec)