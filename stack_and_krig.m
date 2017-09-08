
close all
clear all

addpath functions

%% DEFINE GLACIER and Select DEMs for to stake and fill. Need to have at least two DEMs covering the full extend of the glacier
if ispc
    slash='\';
elseif isunix
    slash='/';
end
data_directory='C:\Users\cmcneil\Desktop\GIS\Geodetics\'
% data_directory=['Q:',slash,'Project Data',slash,'GlacierData',slash,'GIS',slash,];
d = dir(data_directory);
str = {d.name};
[s,v] = listdlg('PromptString','Select Master DEM for coregistration:',...
                'SelectionMode','single',...
                'ListString',str);
glacier=d(s).name;
DEM_directory=[data_directory,glacier,slash,'DEMs',slash,'Processed'];%set path for processed DEMs
%main_directory=['C:\Users\cmcneil\Desktop\Matlab\Taku and LemonCreeek glacier reanalysis code\2017.06.12\data\',glacier,'\DEMs\Raw'];%set path for raw DEMs
d = dir(DEM_directory);
dems = {d.name};
[s,v] = listdlg('PromptString','Select DEMs to stack',...
                'SelectionMode','multiple',...
                'ListString',dems);
%%             
for i=1:length(s)
    dem_path=[DEM_directory,'\',cell2mat(dems(s(i))),'\',cell2mat(dems(s(i))),'_DEM.tif'];
    dem_year=datestr(datenum(dems(s(i)),'yyyy.mm.dd'),'yyyy');
    mask_path=[data_directory,glacier,'\Layers.gdb\boundaries\C',dem_year,'_boundary'];
    output_path=[data_directory,glacier,slash];
    py.mask.getmask(dem_path,mask_path,output_path)
    mask_path=[output_path,'mask.tif'];
    boundary(i)=read_dem(mask_path); 
    boundary(i).Z(isnan(boundary(i).Z))=1;
    area(i,1)=(sum(boundary(i).Z(:)).*boundary(i).cellsize^2)./1000000;
    delete mask_path
    a=cell2mat(dems(s(i)));
    data_years(i,1)=str2num(a(1,1:4));
end

%%
    [~,area_index]=max(area);
    extent=[max(boundary(area_index).X((boundary(area_index).Z==1)))+100 min(boundary(area_index).X((boundary(area_index).Z==1)))-100 max(boundary(area_index).Y((boundary(area_index).Z==1)))+100 min(boundary(area_index).Y((boundary(area_index).Z==1)))-100];
    Xs=boundary(area_index).X(1,boundary(area_index).X(1,:) <= extent(1,1) & boundary(area_index).X(1,:) >= extent(1,2));
    Ys=boundary(area_index).Y(boundary(area_index).Y(:,1) <= extent(1,3) & boundary(area_index).Y(:,1) >= extent(1,4),1);
    for i=1:length(s)
    dem_path=[DEM_directory,'\',cell2mat(dems(s(i))),'\',cell2mat(dems(s(i))),'_DEM.tif'];
    dem_year=datestr(datenum(dems(s(area_index)),'yyyy.mm.dd'),'yyyy');
    mask_path=[data_directory,glacier,'\Layers.gdb\boundaries\C',dem_year,'_boundary'];
    output_path=[data_directory,glacier,slash];
    py.mask.getmask(dem_path,mask_path,output_path)
    mask_path=[output_path,'mask.tif'];
    DEM(i) = read_dem(dem_path,mask_path);
        stack(i).X=repmat(Xs,length(Ys),1);
        stack(i).Y=repmat(Ys,1,length(Xs));
        stack(i).Z=interp2(DEM(i).X,DEM(i).Y,DEM(i).Z,stack(i).X,stack(i).Y, 'linear');
        figure(i);hold on
        imshow(stack(i).Z,[])
    end

%%

figure(); hold on
imshow(mean_z,[]);hold on
colormap(jet)
colorbar

%%

dz_1=stack(3).Z-stack(1).Z;
dz_2=stack(2).Z-stack(1).Z;
dz_3=dz_2./dz_1;
dz_ratio=nanmedian(dz_3(:));
dz_4=stack(1).Z+(dz_1.*dz_ratio);
nanmean(dz_1(:).*dz_ratio)
figure();hold on
imshow(dz_4,[])
colormap(jet)
colorbar


% [FileName,PathName,FilterIndex] = uigetfile(FilterSpec)