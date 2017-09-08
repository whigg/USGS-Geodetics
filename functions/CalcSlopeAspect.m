function [slope,aspect] = CalcSlopeAspect(DEM,cellsize)
%%  [slope] = Slope2(DEM,cellsize)
% This function calculates slopes of a DEM using the same methodology as
% ArcGIS. Border cells are set to NaN 
% Author: Christopher Nuth 

%% Initiate individual grids for matrix calculation 

[nrows ncols] = size(DEM);
    slope = nan(size(DEM));
     aspect = nan(size(DEM));
     

% DEM1 = DEM(1:end-2,1:end-2);
% DEM2 = DEM(1:end-2,2:end-1);
% DEM3 = DEM(1:end-2,3:end);
% DEM4 = DEM(2:end-1,1:end-2);
% DEM5 = DEM(2:end-1,2:end-1);
% DEM6 = DEM(2:end-1,3:end);
% DEM7 = DEM(3:end,1:end-2);
% DEM8 = DEM(3:end,2:end-1);
% DEM9 = DEM(3:end,3:end);


%% Calculate Slope

dx = ( (DEM(1:end-2,3:end) + 2.*DEM(2:end-1,3:end) + DEM(3:end,3:end)) - (DEM(1:end-2,1:end-2) + 2.*DEM(2:end-1,1:end-2) + DEM(3:end,1:end-2)) ) / (8 .* cellsize);
dy = ( (DEM(1:end-2,1:end-2) + 2.*DEM(1:end-2,2:end-1) + DEM(1:end-2,3:end)) - (DEM(3:end,1:end-2) + 2.*DEM(3:end,2:end-1) + DEM(3:end,3:end)) ) / (8 .* cellsize);

slope(2:nrows-1,2:ncols-1) = atand(sqrt( dx.^2 + dy.^2 ));
aspect(2:nrows-1,2:ncols-1) = 180 -  atand(dy./dx) + 90.*(dx./abs(dx));

% end
