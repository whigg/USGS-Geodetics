function [myshift] = Add_shift_params(myshift,slp_mean)
%ADJUST_SHIFT Uses the paramters given by 'Fitting_routine.m' to adjust the
%SLAVE dataset. 
%   -myshift = output from 'Fitting_routine.m'
%   -slp_mean = the mean slope of stable terrain used to determine the
%       co-registration parameters. Is an output from 'Preprocess.m'. 
%   -x, y = 1-D vectors describing the locations of the pixels in Z. 
%   -Z = 2-D DEM. 
%   -x_i,y_i,Z_i = co-registered DEM
%   -myshift = updated structured matrix with the x,y and z adjustments


%% get the x and y adjustment parameters from cosine fit parameters
% myshift.x_adj = myshift.param.a .* sind(myshift.param.b);
% myshift.y_adj = myshift.param.a .* cosd(myshift.param.b);
% myshift.z_adj = myshift.param.c .* (slp_mean);

myshift.x_adj = .5 .* myshift.param(1) .* sind(myshift.param(2));
myshift.y_adj = .5 .* myshift.param(1) .* cosd(myshift.param(2));
myshift.z_adj = .5 .* myshift.param(3) .* (slp_mean);


%% estimate the error if using 'fit'! TO BE UPDATED

mye = confint(myshift.param2,0.99); % 99% confidence interval
a = myshift.param(1);
b = myshift.param(2);
c = myshift.param(3);

e_a = diff(mye(:,1))/2;
    e_a_dx = abs(e_a.*sind(b));
    e_a_dy = abs(e_a.*cosd(b));
    
e_b = diff(mye(:,2))/2;
    e_b_dx = abs( (a.*sind(b+e_b))-(a.*sind(b-e_b)) );
    e_b_dy = abs( (a.*cosd(b+e_b))-(a.*cosd(b-e_b)) );

e_c = diff(mye(:,3))/2;

myshift.x_err = e_a_dx + e_b_dx;
myshift.y_err = e_a_dy + e_b_dy;
% myshift.z_err = e_c .* tand(slp_mean);
myshift.z_err = e_c .* tand(slp_mean);
myshift.slp_mean = slp_mean;
end

