
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
% data_directory='C:\Users\cmcneil\Desktop\GIS\Geodetics\' 
data_directory=['Q:',slash,'Project Data',slash,'GlacierData',slash,'GIS',slash,];
d = dir(data_directory);                                    %list all glaciers with data in directory
str = {d.name}; 
[s,v] = listdlg('PromptString','Select glacier:',...
                'SelectionMode','single',...
                'ListString',str);                          %prompt dialog box for users to select which glacier to work on
glacier=d(s).name;
DEM_directory=[data_directory,glacier,slash,'DEMs',slash,'Processed'];%set path for processed DEMs
%main_directory=['C:\Users\cmcneil\Desktop\Matlab\Taku and LemonCreeek glacier reanalysis code\2017.06.12\data\',glacier,'\DEMs\Raw'];%set path for raw DEMs
d = dir(DEM_directory);                                     %list all processed (coregistered DEMs for directory/glacier
dems = {d.name};
[sa,v] = listdlg('PromptString','Select DEMs to difference',...
                'SelectionMode','multiple',...
                'ListString',dems);   

dz_directory=[DEM_directory,slash,'dz_maps'];
d = dir(dz_directory);                                     %list all processed (coregistered DEMs for directory/glacier
directory = {d.name};
[s,v] = listdlg('PromptString','Directory for dz_map (reference DEM for coregistration)',...
                'SelectionMode','multiple',...
                'ListString',directory);   
dz_map_path=[DEM_directory,slash,'dz_maps',slash,cell2mat(directory(s)),slash,cell2mat(dems(sa(2))),'_',cell2mat(dems(sa(1))),'_dz.tif'];
%%
% if ~exist([dz_map_directory,slash], 'dir')
%     % Folder does not exist so create it.
%     mkdir([dz_map_directory,slash]);
% end
if exist([DEM_directory,slash,cell2mat(dems(sa(1))),slash,'interpolated',slash], 'dir')
    dem1=[DEM_directory,slash,cell2mat(dems(sa(1))),slash,'interpolated',slash',cell2mat(dems(sa(1))),'_DEM.tif'];
else
     dem1=[DEM_directory,slash,cell2mat(dems(sa(1))),slash,cell2mat(dems(sa(1))),'_DEM.tif'];
end
if exist([DEM_directory,slash,cell2mat(dems(sa(2))),slash,'interpolated',slash], 'dir')
    dem2=[DEM_directory,slash,cell2mat(dems(sa(2))),slash,'interpolated',slash',cell2mat(dems(sa(2))),'_DEM.tif'];
else
     dem2=[DEM_directory,slash,cell2mat(dems(sa(2))),slash,cell2mat(dems(sa(2))),'_DEM.tif'];
end
[dz]=DEM_Difference(dem2,dem1,dz_map_path);

figure(); hold on
imshow(dz,[]);hold on
colormap(jet)
colorbar
box on