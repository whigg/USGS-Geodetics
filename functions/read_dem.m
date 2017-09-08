function [DEM]=read_dem(varargin)
% read_dem: read geotiff using imread and assign map info from infinfo.
%
% DEM =read_dem('filename');
% Reads whole images
% DEM = read_dem('filename','filename2');
% extract subset of the specified.

% output:
%DEM.X - matrix of cell easting coordinates 
%DEM.Y - matrix of cell Northing coordinates 
%DEM.Z - matrix of DEM elevations
dbstop if error
name = varargin{1};

Tinfo        = imfinfo(name);
info.samples = Tinfo(1).Width;
info.lines   = Tinfo(1).Height;
info.imsize  = Tinfo(1).Offset;
info.bands   = Tinfo(1).SamplesPerPixel;

sub = [1, info.samples, 1, info.lines];

info.map_info.dx = Tinfo(1).ModelPixelScaleTag(1);
info.map_info.dy = Tinfo(1).ModelPixelScaleTag(2);
info.map_info.mapx = Tinfo(1).ModelTiepointTag(4);
info.map_info.mapy = Tinfo(1).ModelTiepointTag(5);
%info.map_info.projection_name = Tinfo.GeoAsciiParamsTag;
%info.map_info.projection_info = Tinfo.GeoDoubleParamsTag;

minx = info.map_info.mapx;
maxy = info.map_info.mapy;
maxx = minx + (info.samples-1).*info.map_info.dx;
miny = maxy - (info.lines-1  ).*info.map_info.dy;

%info.CornerMap = [minx miny; maxx miny; maxx maxy; minx maxy; minx miny]; 

xm = info.map_info.mapx;
ym = info.map_info.mapy;
x_ = xm + ((0:info.samples-1).*info.map_info.dx);
y_ = ym - ((0:info.lines  -1).*info.map_info.dy);

I.x = x_(sub(1):sub(2));
I.y = y_(sub(3):sub(4));

[A] = imread(name);

I.z=A;
I.info = info;

I.cellsize=I.info.map_info.dx;
DEM.cellsize=I.cellsize;
% for i=1:length(I.y)
%     DEM.X(i,:)=I.x+(I.cellsize/2);
% end
x=I.x+(I.cellsize/2);
DEM.X=repmat(x,length(I.y),1);
% for i=1:length(I.x)
%     DEM.Y(:,i)=I.y-(I.cellsize/2);
% end
y=(I.y-(I.cellsize/2));
DEM.Y=repmat(y,length(I.x),1)';
DEM.Z=double(I.z);
no_data=str2num(Tinfo(1).GDAL_NODATA);
DEM.Z(DEM.Z==no_data)=nan;
DEM.Tinfo=Tinfo;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if mask was input import mask
if nargin == 2
    clearvars -except DEM varargin
    name = varargin{2};

    Tinfo        = imfinfo(name);
    info.samples = Tinfo(1).Width;
    info.lines   = Tinfo(1).Height;
    info.imsize  = Tinfo(1).Offset;
    info.bands   = Tinfo(1).SamplesPerPixel;

    sub = [1, info.samples, 1, info.lines];

    info.map_info.dx = Tinfo(1).ModelPixelScaleTag(1);
    info.map_info.dy = Tinfo(1).ModelPixelScaleTag(2);
    info.map_info.mapx = Tinfo(1).ModelTiepointTag(4);
    info.map_info.mapy = Tinfo(1).ModelTiepointTag(5);
    %info.map_info.projection_name = Tinfo.GeoAsciiParamsTag;
    %info.map_info.projection_info = Tinfo.GeoDoubleParamsTag;

    minx = info.map_info.mapx;
    maxy = info.map_info.mapy;
    maxx = minx + (info.samples-1).*info.map_info.dx;
    miny = maxy - (info.lines-1  ).*info.map_info.dy;

    %info.CornerMap = [minx miny; maxx miny; maxx maxy; minx maxy; minx miny]; 

    xm = info.map_info.mapx;
    ym = info.map_info.mapy;
    x_ = xm + ((0:info.samples-1).*info.map_info.dx);
    y_ = ym - ((0:info.lines  -1).*info.map_info.dy);

    I.x = x_(sub(1):sub(2));
    I.y = y_(sub(3):sub(4));

    [A] = double(imread(name));
    A(A==255)=1;
     B=DEM.Z.*A;
     B(B==0)=nan;
   DEM.Z_masked=B;
end

 end
       







