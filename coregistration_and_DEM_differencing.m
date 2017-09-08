%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Coregistration and DEM differencing script 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INPUTS: 1) Master DEM - this is the DEM that will be held fixed in space
%        2) Slave DEM - this is the DEM that will be aligned to the master
%        DEM
%        3) A shapefile that defines stable, relatively smooth, and low
%        angle (<40*) terrain to be used for aligning DEMs

%general script outline
        




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
[s,v] = listdlg('PromptString','Select Master DEM for coregistration:',...
                'SelectionMode','single',...
                'ListString',str);
master_dem=d(s).name;
master_dem_path=[raw_DEM_directory,slash,'',master_dem,slash,master_dem,'_DEM.tif'];

[s,v] = listdlg('PromptString','Select Slave DEM for coregistration:',...
                'SelectionMode','single',...
                'ListString',str);
slave_dem=d(s).name; 
slave_dem_path=[raw_DEM_directory,slash,slave_dem,slash,slave_dem,'_DEM.tif'];
%% Import DEMs and find best fit alignment

tic
project_name = ['coregistration_of_',slave_dem,'_to_',master_dem];%Glacier Slaveyear Masteryear
addpath functions
name = strcat(project_name);%, timestamp);
max_iterations = 25; %max number of iterations
stop_iterations = 0.01; %the shift distance at which to stop iterations 
coregistration_area_date=datestr(datenum(slave_dem,'yyyy.mm.dd'),'yyyymmdd');
mask_path=[data_directory,glacier,'\Layers.gdb\coregistration_areas\coregistration_area_',coregistration_area_date];

w = waitbar(1/(max_iterations + 10), 'importing Master DEM');

file=dir(master_dem_path);
if file.bytes>1e+8
    
    feature_path=mask_path;
    dem_dir=[raw_DEM_directory,slash,'',master_dem];
    py.clip.clipRaster(feature_path,dem_dir)
    master_temp_dem_path=[raw_DEM_directory,slash,'',master_dem,slash,'temp.tif'];
    DEM1 = read_dem(master_temp_dem_path);
else
    DEM1 = read_dem(master_dem_path);
end

[DEM1.slope,DEM1.aspect] = CalcSlopeAspect(DEM1.Z,DEM1.cellsize);

file=dir(slave_dem_path);
if file.bytes>1e+8

    feature_path=mask_path;
    dem_dir=[raw_DEM_directory,slash,slave_dem];
    py.clip.clipRaster(feature_path,dem_dir)
    slave_temp_dem_path=[raw_DEM_directory,slash,'',slave_dem,slash,'temp.tif'];
    dem_path=slave_temp_dem_path;
    output_path=[data_directory,glacier,slash];
    py.mask.getmask(dem_path,mask_path,output_path)
    mask_path=[output_path,'\mask.tif'];
    DEM2 = read_dem(slave_temp_dem_path,mask_path);
else
    dem_path=slave_dem_path;
    output_path=[data_directory,glacier,slash];
    py.mask.getmask(dem_path,mask_path,output_path)
    mask_path=[output_path,'\mask.tif'];
    DEM2 = read_dem(slave_dem_path,mask_path);
end

cellsizes=[DEM1.cellsize  DEM2.cellsize];
[~,DEM_index]=max(cellsizes);


[DEM1.slope,DEM1.aspect] = CalcSlopeAspect(DEM1.Z,DEM1.cellsize);
[DEM2.slope,DEM2.aspect] = CalcSlopeAspect(DEM2.Z_masked,DEM2.cellsize);


% initialize
max_iterations = 25; %max number of iterations
stop_iterations = 0.01; %the shift distance at which to stop iterations  
shift.X = zeros(1,max_iterations); shift.Y = zeros(1,max_iterations); shift.Z = zeros(1,max_iterations);
mag = 1;
%%

for n = 1:max_iterations
    waitbar((n+3)/(max_iterations + 10),w, 'minimizing residuals');
    if mag >= stop_iterations; %r is the ratio of the last shift distance to the cellsize
        if DEM_index == 2
            good = isfinite(DEM2.Z_masked);
            dZ = DEM2.Z_masked - interp2(DEM1.X,DEM1.Y,DEM1.Z,DEM2.X + shift.X(n),DEM2.Y + shift.Y(n), 'linear') + shift.Z(n);
            slp = interp2(DEM1.X,DEM1.Y,DEM1.slope,DEM2.X + shift.X(n),DEM2.Y + shift.Y(n), 'linear');
            asp = interp2(DEM1.X,DEM1.Y,DEM1.aspect,DEM2.X + shift.X(n),DEM2.Y + shift.Y(n), 'linear');
            mytitle = ['XYshift #' num2str(n)];
            [myshift] = Fitting_routine( dZ(good), slp(good), asp(good), mytitle);
            shift.X(n+1) = shift.X(n) - myshift.x_adj; %note that iteration 1 is stored in column 2
            shift.Y(n+1) = shift.Y(n) - myshift.y_adj;
            shift.Z(n+1) = shift.Z(n) - myshift.z_adj;
            mag = sqrt(myshift.x_adj^2 + myshift.y_adj^2 + myshift.z_adj^2);
            iteration = n;
        elseif DEM_index == 1
            good = isfinite(DEM1.Z);
            dZ = interp2(DEM2.X + shift.X(n),DEM2.Y + shift.Y(n),DEM2.Z_masked + shift.Z(n),DEM1.X ,DEM1.Y , 'linear') - DEM1.Z;
            slp = interp2(DEM2.X + shift.X(n),DEM2.Y + shift.Y(n),DEM2.slope,DEM1.X ,DEM1.Y  , 'linear');
            asp = interp2(DEM2.X + shift.X(n),DEM2.Y + shift.Y(n),DEM2.aspect,DEM1.X ,DEM1.Y , 'linear');
            mytitle = ['XYshift #' num2str(n)];
            [myshift] = Fitting_routine( dZ(good), slp(good), asp(good), mytitle);
            shift.X(n+1) = shift.X(n) - myshift.x_adj; %note that iteration 1 is stored in column 2
            shift.Y(n+1) = shift.Y(n) - myshift.y_adj;
            shift.Z(n+1) = shift.Z(n) - myshift.z_adj;
            mag = sqrt(myshift.x_adj^2 + myshift.y_adj^2 + myshift.z_adj^2);
            iteration = n;
        end
    else  
    end
end



% the shifted DEM (off glacier);
waitbar((max_iterations + 5)/(max_iterations + 10),w, 'calculating final statistics');
DEM2.X_shifted = DEM2.X + shift.X(iteration + 1); 
DEM2.Y_shifted = DEM2.Y + shift.Y(iteration + 1); 
DEM2.Z_shifted = DEM2.Z_masked + shift.Z(iteration + 1);
if DEM_index ==1
    Error.dZ = interp2(DEM2.X_shifted,DEM2.Y_shifted,DEM2.Z_shifted,DEM1.X,DEM1.Y, 'linear') - DEM1.Z;
    Error.X=DEM1.X;
    Error.Y=DEM1.Y;
    Error.Z=DEM1.Z;
elseif DEM_index == 2
    Error.dZ = DEM2.Z_shifted - interp2(DEM1.X,DEM1.Y,DEM1.Z,DEM2.X_shifted,DEM2.Y_shifted, 'linear');
    Error.X=DEM2.X_shifted;
    Error.Y=DEM2.Y_shifted;
    Error.Z=DEM2.Z_shifted;
end
shift.Z(iteration + 2) = shift.Z(iteration + 1) - nanmean(Error.dZ(:)); %Final Z adjustment based on final XY pos
DEM2.Z_shifted = DEM2.Z_shifted - (shift.Z(iteration + 1) - shift.Z(iteration + 2)); 
Error.dZ = Error.dZ - (shift.Z(iteration + 1) - shift.Z(iteration + 2));
rmse = sqrt(nanmean(Error.dZ(:).^2)); %assuming no spatial covariance
se = rmse/(length(Error.dZ(:))); %estimate of vertical uncertainty between the final coregistered products assuming no spatial covariance

% plot the coregistration path
f(1) = figure;hold on
text(0,0,'start', 'color', [1 0 0], 'FontSize', 8)
title('Coregistration path')
ylabel('northing [m]') 
xlabel('easting [m]')
C = [flipud((1:iteration)'./iteration) ,zeros(iteration,1), (1:iteration)'./iteration]; C = [1,0,0;C];
scatter(shift.X(1:iteration+1)',shift.Y(1:iteration+1)', 200, C); 
for n = 2:iteration+1
    text(shift.X(n),shift.Y(n),num2str(n-1), 'color', C(n,:), 'FontSize', 8) 
end
text(shift.X(iteration)/2,shift.Y(iteration)/2,[num2str(iteration) ' iterations to convergence'], 'FontSize', 10)   
last_shift = mag*DEM1.cellsize;
text(shift.X(iteration)/2,shift.Y(iteration)/2-abs(shift.Y(iteration)/5),['Last shift ' num2str(last_shift, 2) ' m'], 'FontSize', 10)

% plot the dZ field to asses outliers, problem areas, and the spatial distribution

NMAD = nanmedian ( abs(dZ(:)));
%%
f(2) = figure; hold on
h=imagesc(Error.X(:,1),flipud(Error.Y(1,:)),Error.dZ);hold on
%h=imagesc(image(i).data);hold on
set(h,'alphadata',~isnan(Error.dZ))
colormap jet
colorbar
%axis image
% imshow(Error.dZ,[]);hold on
% colormap(jet)
% colorbar
title('Residuals')
ylabel('northing [m]') 
xlabel('easting [m]')
    %%
% calculate (semi) variance
waitbar((max_iterations + 6)/(max_iterations + 10),w, 'calculating the variogram');
f(3) = figure;  hold on
x=[reshape(Error.X(~isnan(Error.dZ)),[length(Error.X(~isnan(Error.dZ))),1]) reshape(Error.Y(~isnan(Error.dZ)),[length(Error.Y(~isnan(Error.dZ))),1]) reshape(Error.Z(~isnan(Error.dZ)),[length(Error.Z((~isnan(Error.dZ)))),1])];
y=reshape(Error.dZ(~isnan(Error.dZ)),[length(Error.dZ(~isnan(Error.dZ))),1]);
G = variogram(x,y, 'nrbins', 20, 'maxdist', 10000, 'type', 'gamma', 'plotit', true, 'subsample', 10000); %
% G = variogram([reshape(DEM2.X_shifted(~isnan(DEM2.error)),[length(DEM2.X_shifted(~isnan(DEM2.error))),1]) reshape(DEM2.Y_shifted(~isnan(DEM2.error)),[length(DEM2.Y_shifted(~isnan(DEM2.error))),1]) reshape(DEM2.Z_shifted(~isnan(DEM2.error)),[length(DEM2.Z_shifted((~isnan(DEM2.error)))),1])],reshape(DEM2.error(~isnan(DEM2.error)),[length(DEM2.error(~isnan(DEM2.error))),1]), 'nrbins', 20, 'maxdist', 10000, 'type', 'gamma', 'plotit', true, 'subsample', 10000); %


raw_dem=[raw_DEM_directory,slash,slave_dem,slash,slave_dem,'_DEM.tif'];
processed_dem=[processed_DEM_directory,slash,slave_dem,slash,slave_dem,'_DEM.TIF'];
%%
final_index(1)=find(shift.X~=0,1,'last');
final_index(2)=find(shift.Y~=0,1,'last');
final_index(3)=find(shift.Z~=0,1,'last');
Final_DEMs_Shift=[shift.X(final_index(1)) shift.Y(final_index(2)) shift.Z(final_index(3))];
shift_table=table({master_dem},Final_DEMs_Shift(1,1),Final_DEMs_Shift(1,2),Final_DEMs_Shift(1,3),NMAD,'VariableNames',{'master_dem' 'X_shift' 'Y_shift' 'Z_shift' 'NMAD'})

%% Shift slave DEM
file=dir(master_dem_path);
if file.bytes>1e+8
    delete_file=[raw_DEM_directory,slash,'',master_dem,slash,'temp*'];
    delete(delete_file);
end
if ~exist([processed_DEM_directory,slash,slave_dem,slash], 'dir')
        % Folder does not exist so create it.
        mkdir([processed_DEM_directory,slash,slave_dem,slash]);
end
clear DEM1 
file=dir(slave_dem_path);
if file.bytes>1e+8
    delete_file=[raw_DEM_directory,slash,'',slave_dem,slash,'temp*'];
    delete(delete_file);
    clear DEM2
    DEM2=read_dem(raw_dem);
else
end
copyfile(raw_dem,processed_dem)
t=Tiff(processed_dem,'r+');

tagstructure.ModelTiepointTag=[DEM2.Tinfo(1).ModelTiepointTag(1:3) DEM2.Tinfo(1).ModelTiepointTag(4)+Final_DEMs_Shift(1,1) DEM2.Tinfo(1).ModelTiepointTag(5)+Final_DEMs_Shift(1,2) DEM2.Tinfo(1).ModelTiepointTag(6)];
tagstructure.BitsPerSample = 64;
tagstructure.SampleFormat = 3;
t.setTag(tagstructure)
t.write(DEM2.Z+Final_DEMs_Shift(1,3));
t.close

writetable(shift_table,[processed_DEM_directory,slash,slave_dem,slash,slave_dem,'_shifts.csv']);
%% Shift slave ortho image with DEM if it exists

raw_Ortho_directory=['Q:',slash,'Project Data',slash,'GlacierData',slash,'GIS',slash,glacier,slash,'Orthos',slash,'Raw'];%set path for raw DEMs
processed_Ortho_directory=['Q:',slash,'Project Data',slash,'GlacierData',slash,'GIS',slash,glacier,slash,'Orthos',slash,'Processed'];
raw_ortho=[raw_Ortho_directory,slash,slave_dem,slash,slave_dem,'_Ortho.tif'];
processed_ortho=[processed_Ortho_directory,slash,slave_dem,slash,slave_dem,'_Ortho.TIF'];
if exist([raw_Ortho_directory,slash,slave_dem,slash], 'dir')
%      Ortho = read_dem(raw_ortho);
    if ~exist([processed_Ortho_directory,slash,slave_dem,slash], 'dir')
      
        % Folder does not exist so create it.
        mkdir([processed_Ortho_directory,slash,slave_dem,slash]);
    end
    Tinfo=imfinfo(raw_ortho);
    copyfile(raw_ortho,processed_ortho)
    t=Tiff(processed_ortho,'r+');
    tagstructure.ModelTiepointTag=[Tinfo(1).ModelTiepointTag(1:3) Tinfo(1).ModelTiepointTag(4)+Final_DEMs_Shift(1,1) Tinfo(1).ModelTiepointTag(5)+Final_DEMs_Shift(1,2) Tinfo(1).ModelTiepointTag(6)];
    t.setTag(tagstructure)
    %t.write();
    t.close
end

%% difference DEMs
clearvars -except master_dem slave_dem master_dem_path processed_dem processed_DEM_directory slash shift_table data_directory glacier

dz_map_directory=[processed_DEM_directory,slash,'dz_maps',slash,cell2mat(table2array(shift_table(1,1))),slash];
dz_map_path=[processed_DEM_directory,slash,'dz_maps',slash,cell2mat(table2array(shift_table(1,1))),slash,master_dem,'_',slave_dem,'_dz.tif'];
if ~exist([dz_map_directory,slash], 'dir')
    % Folder does not exist so create it.
    mkdir([dz_map_directory,slash]);
end
[map]=DEM_Difference(master_dem_path,processed_dem,dz_map_path);
%%
figure(); hold on
h=imagesc(map.X(1,:),map.Y(:,1),map.dz);hold on
set(h,'alphadata',~isnan(map.dz))
colormap jet
colorbar
axis image
box on
