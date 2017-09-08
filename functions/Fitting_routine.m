function [myshift] = Fitting_routine( dZ, slp, asp, mytitle)
%FITTING_ROUTINE uses the standard curve fitting toolbox to solve the
%co-registration equation. Uses a robust linear least squares approach.
%Also plots the solution of the cosinusoidal equation and the residuals.
% Based on a paper by Christopher Nuth (2011), and the original function written by the same  
%
%       Minor Revisions: Louis Sass (2014.01.27)
% 
%   [myshift ] = Fitting_routine( aspect, dhtan )
%       dZ = 1-D vector of the apparent change in surface elevation between
%           the base DEM and the comparison points 
%       asp = Aspect (1-D vector) of base DEM at the comparison points, in degrees from north
%       myshift = structured matrix of the returned fit parameters,
%           identical to the output from 'FIT' function. See help FIT for 
%           more details on the output. 

warning off

tanslp = tand(slp); %tanslp = 1-D vector, slope of base DEM at the comparison points
dZtan = dZ./tanslp; %dZtan = dZ / tan (slope) (1-D vector) in degrees
index = ~isnan(dZtan); dZtan = dZtan(index);asp = asp(index); %get rid of NaNs
index = ~isnan(asp); dZtan = dZtan(index);asp = asp(index); clear index

%% fitting using single curve fitting
f = fittype('a.*cosd(b-x) + c','coefficients',{'a','b','c'});
fopt = fitoptions('Method','LinearLeastSquares','Robust','LAR');
[cfun,gof,output] = fit(asp,dZtan,f,fopt);
param = coeffvalues(cfun);

%% prepare structured output

output.residuals = [];
output.Jacobian = [];

myshift.param = param;
myshift.param2 =cfun;
myshift.gof = gof;
myshift.output = output;

[myshift] = Add_shift_params(myshift,nanmean(tanslp));

%% Plotting Test

% create solution
myasp = 1:360;
mysol = param(1)*cosd(param(2)-myasp)+param(3);
figure (), hold on
        plot(asp,dZtan,'k.','MarkerSize',2);
        plot(myasp,mysol,'b-','LineWidth',2);
        plot(myasp,zeros(size(myasp)),'-','Color',[0.5 0.5 0.5],'LineWidth',2);
        axis([0 360 -300 300]);
        xlabel('Aspect [degrees]');
        ylabel('dZ / tan (slope)');
        title(mytitle,'Interpreter','none');
        text(45,250,['x shift = ' num2str(myshift.x_adj,2)] );
        text(45,220,['y shift = ' num2str(myshift.y_adj,2)] );
        text(45,190,['z shift = ' num2str(myshift.z_adj,2)] );
    hold off
snapnow
end

