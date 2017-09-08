function [map] = DEM_Difference(varargin)
dbstop if error
if length(varargin)<2
    error('Must have two DEMs to difference')
elseif length(varargin)<3
    error('Need path for output')
else
    %%
filename1=varargin{1};
filename2=varargin{2};
output=varargin{3};
DEM1=read_dem(filename1);
DEM2=read_dem(filename2);

cellsizes=[DEM1.cellsize  DEM2.cellsize];
dates=[datenum(filename1((length(filename1)-17):(length(filename1)-8)),'yyyy.mm.dd') datenum(filename2((length(filename2)-17):(length(filename2)-8)),'yyyy.mm.dd')];

[~,DEM_index]=max(cellsizes);
if DEM_index == 1 && dates(1)>dates(2)
    dz = DEM1.Z - interp2(DEM2.X,DEM2.Y,DEM2.Z,DEM1.X,DEM1.Y, 'linear');
elseif DEM_index == 1 && dates(1)<dates(2)
    dz =  interp2(DEM2.X,DEM2.Y,DEM2.Z,DEM1.X,DEM1.Y, 'linear') - DEM1.Z;
elseif DEM_index == 2 && dates(1)<dates(2)
    dz =  DEM2.Z- interp2(DEM1.X,DEM1.Y,DEM1.Z,DEM2.X,DEM2.Y, 'linear');
elseif DEM_index == 2 && dates(1)>dates(2)
    dz =  interp2(DEM1.X,DEM1.Y,DEM1.Z,DEM2.X,DEM2.Y, 'linear')-DEM2.Z;

end

if DEM_index==1
    copyfile(filename1,output)
    image_info=DEM1.Tinfo;
    t=Tiff(output,'r+');
    tagstructure.ModelTiepointTag=[image_info(1).ModelTiepointTag(1:3) image_info(1).ModelTiepointTag(4) image_info(1).ModelTiepointTag(5) image_info(1).ModelTiepointTag(6)];
    tagstructure.BitsPerSample = 64;
    tagstructure.SampleFormat = 3;
    t.setTag(tagstructure)
    
    t.write(dz);
    t.close
    map.dz=dz;
    map.X=DEM1.X;
    map.Y=DEM1.Y;
elseif DEM_index==2
   copyfile(filename2,output)
   image_info=DEM2.Tinfo;
    t=Tiff(output,'r+');
    tagstructure.ModelTiepointTag=[image_info(1).ModelTiepointTag(1:3) image_info(1).ModelTiepointTag(4) image_info(1).ModelTiepointTag(5) image_info(1).ModelTiepointTag(6)];
    tagstructure.BitsPerSample = 64;
    tagstructure.SampleFormat = 3;
    t.setTag(tagstructure)
    t.write(dz);
    t.close
    map.dz=dz;
    map.X=DEM2.X;
    map.Y=DEM2.Y;
end
end
end

