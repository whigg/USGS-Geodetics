clear all
close all

addpath functions
if ispc
    slash='\';
elseif isunix
    slash='/';
end
% data_directory='C:\Users\cmcneil\Desktop\GIS\Geodetics\'
data_directory=['Q:',slash,'Project Data',slash,'GlacierData',slash,'GIS',slash,];
d = dir(data_directory);
str = {d.name};
[s,v] = listdlg('PromptString','Select Master DEM for coregistration:',...
                'SelectionMode','single',...
                'ListString',str);
glacier=d(s).name;
Reference_DEM_directory=[data_directory,glacier,slash,'DEMs',slash,'Processed\dz_maps\'];%set path for processed DEMs
%main_directory=['C:\Users\cmcneil\Desktop\Matlab\Taku and LemonCreeek glacier reanalysis code\2017.06.12\data\',glacier,'\DEMs\Raw'];%set path for raw DEMs
d = dir(Reference_DEM_directory);
dems = {d.name};
[s,v] = listdlg('PromptString','Select DEMs to stack',...
                'SelectionMode','single',...
                'ListString',dems);
dz_map_directory= [Reference_DEM_directory,'\',cell2mat(dems(s)),'\']; 

d = dir(dz_map_directory);
dz_maps = {d.name};
[s,v] = listdlg('PromptString','Select DEMs to stack',...
                'SelectionMode','single',...
                'ListString',dz_maps);
dz_map_path=[dz_map_directory,'\',cell2mat(dz_maps(s))];

            %%
dzs=[];
map_name=cell2mat(dz_maps(s));
dem_years(1)=str2num(map_name(1:4));
dem_years(2)=str2num(datestr(datenum(map_name(1,end-16:end-7),'yyyy.mm.dd'),'yyyy'));
for i=1:2
    mask_path=[data_directory,glacier,'\Layers.gdb\boundaries\C',num2str(dem_years(i)),'_boundary'];
    output_path=[data_directory,glacier,slash];
    dem_path=dz_map_path;
    py.mask.getmask(dem_path,mask_path,output_path)
    mask_path=[output_path,'\mask.tif'];
    boundary(i)=read_dem(mask_path); 
    boundary(i).Z(isnan(boundary(i).Z))=1;
    cellsize(i)=boundary(i).cellsize;
    area(i,1)=(sum(boundary(i).Z(:)).*boundary(i).cellsize^2)./1000000;
end

[~,area_index]=max(area);
mask_path=[data_directory,glacier,'\Layers.gdb\boundaries\C',num2str(dem_years(area_index)),'_boundary'];
output_path=[data_directory,glacier,slash];
dem_path=dz_map_path;
py.mask.getmask(dem_path,mask_path,output_path)
mask_path=[output_path,'\mask.tif'];
dz_map = read_dem(dem_path,mask_path);
imshow(dz_map.Z_masked,[])
nanmean(dz_map.Z_masked(:))
