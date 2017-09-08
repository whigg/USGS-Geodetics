% clear allspyd
% close all
% addpath functions
% tic

%% DEFINE GLACIER and Select DEMs for coregistration
if ispc
    slash='\';
elseif isunix
    slash='/';
end
%search data directory for available DEMs
% data_directory='C:\Users\cmcneil\Desktop\GIS\Geodetics\';
data_directory=['Q:',slash,'Project Data',slash,'GlacierData',slash,'GIS',slash,];
d = dir(data_directory);
str = {d.name};
[s,v] = listdlg('PromptString','Select Glacier',...
                'SelectionMode','single',...
                'ListString',str);
glacier=d(s).name;
ortho_directory=[data_directory,glacier,slash,'Orthos',slash,'Raw'];%set path for raw DEMs
d = dir(ortho_directory);
dems = {d.name};
[s,v] = listdlg('PromptString','Select Orthos for timelapse',...
                'SelectionMode','multiple',...
                'ListString',dems);  
            
d=d(s);  
boundary_date=d(1).name(1:4);
mask_path=[data_directory,glacier,'\Layers.gdb\\boundaries\C',boundary_date,'_boundary'];


py.area.UpdateAreaTable(glacier);
area_table_path=[data_directory,glacier,'\Tables\',glacier,'Area.csv'];
data=importdata(area_table_path);
years=data.data(:,1);
area=data.data(:,2);
area_table=[];
%%
for i=1:length(d)
    if i==1
        file=[ortho_directory,slash,slash,d(i).name];
        ortho_years(i,1)=str2num(d(i).name(1:4));
        year_index(i,1)=find(years(:,1)==ortho_years(i,1));
        area_table=[area_table;years(year_index(i,1),1) area(year_index(i,1),1)];
        feature_path=mask_path;
        raster_dir=file;
        cellsize=5;
        py.Raster.resample_and_clip(cellsize,feature_path,raster_dir);
        file=[ortho_directory,slash,d(i).name,slash,'temp.tif'];
        delete_file=[ortho_directory,slash,d(i).name,slash,'temp*'];
        a = read_dem(file);
        delete(delete_file)
        xs=a.X(1,1):5:a.X(1,end);
        ys=a.Y(end,1):5:a.Y(1,1);
        posts.X=repmat(xs,length(ys),1);
        posts.Y=flipud(repmat(ys,length(xs),1)');
        image(i).data=a.Z;
        image(i).data=interp2(a.X,a.Y,a.Z,posts.X,posts.Y,'linear');
        no_data=a.Tinfo.GDAL_NODATA;
        %make figure for snap shot
        figure(1);hold on
        % ortho image
        subplot(2,1,1)
        h=imagesc(posts.X(1,:),flipud(posts.Y(:,1)),image(i).data);hold on
        set(h,'alphadata',image(i).data~=0& ~isnan(image(i).data))
        colormap gray
        title([glacier,' Glacier: ',num2str(ortho_years(i))],'fontname','arial ','fontsize',14,'fontweight','bold');
        axis image
        axis off
        box on
        set(gcf, 'PaperPositionMode', 'auto');
        print -depsc2 gates_epoch2.eps
        % Area figure
        subplot(2,1,2)
        plot(area_table(1:i,1),area_table(1:i,2),'-s','linewidth',2,'color',[1 .5 0])
        ylabel('Area (km^2)','fontname','arial ','fontsize',14,'fontweight','bold')
        xlabel('Time (Years)','fontname','arial ','fontsize',14,'fontweight','bold');
        ylim([round(min(area)-1) round(max(area)+1)])
        ax=gca;
        ax.LineWidth=2;
        ax.FontName='Arial';
        ax.FontSize=14;
        ax.FontWeight='bold';
        xlim([1940 2020])
        
        axis square
        box on
        set(gcf, 'PaperPositionMode', 'auto');
        print -depsc2 gates_epoch2.eps 
        M(i) = getframe(gcf);
    else
        file=[ortho_directory,slash,slash,d(i).name];
        ortho_years(i,1)=str2num(d(i).name(1:4));
        year_index(i,1)=find(years==ortho_years(i,1));
        area_table=[area_table;years(year_index(i,1),1) area(year_index(i,1),1)];
        feature_path=mask_path;
        dem_dir=file;
         py.clip.clipRaster(feature_path,dem_dir);
        file=[ortho_directory,slash,d(i).name,slash,'temp.tif'];
        delete_file=[ortho_directory,slash,d(i).name,slash,'temp*'];
        a = read_dem(file);
        delete(delete_file)
        image(i).data=interp2(a.X,a.Y,a.Z,posts.X,posts.Y,'linear');
        no_data=a.Tinfo.GDAL_NODATA;
        figure(1);hold on
        subplot(2,1,1)
        h=imagesc(posts.X(1,:),flipud(posts.Y(:,1)),image(i).data);hold on
%         h=imagesc(image(i).data);hold on
        set(h,'alphadata',image(i).data~=0& ~isnan(image(i).data))
        colormap gray
        title([glacier,' Glacier: ',num2str(ortho_years(i))],'fontname','arial ','fontsize',14,'fontweight','bold');
        axis image
        axis off
        box on
        set(gcf, 'PaperPositionMode', 'auto');
        print -depsc2 gates_epoch2.eps
         % Area figure
        subplot(2,1,2)
        plot(area_table(1:i,1),area_table(1:i,2),'-s','linewidth',2,'color',[1 .5 0])
        ylabel('Area (km^2)','fontname','arial ','fontsize',14,'fontweight','bold')
        xlabel('Time (Years)','fontname','arial ','fontsize',14,'fontweight','bold');
        ylim([round(min(area)-1) round(max(area)+1)])
        ax=gca;
        ax.LineWidth=2;
        ax.FontName='Arial';
        ax.FontSize=14;
        ax.FontWeight='bold';
        xlim([1940 2020])
        axis square
        box on
        set(gcf, 'PaperPositionMode', 'auto');
        print -depsc2 gates_epoch2.eps 
       
        M(i) = getframe(gcf);
    end
end

%%
M=[M,M(end)];
v = VideoWriter(['Videos\',glacier,'.avi']); %creat blank movie file
v.FrameRate = 1; %set the frames per second rate. In this example we have 21 frames so at 1 frame per second we get a 21 second long movie
open(v) %open blank video file
writeVideo(v,M) %put in you movie frames
close(v) %close dat ish
% toc