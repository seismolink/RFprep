function [rad, tra] =  rad_tra (north, east, bazdeg)

% tranformation of north east to radial-transversal

% usage:
% north & east: time series, same length
% bazdeg: back azimuth in deg

% Copyright 2016 M.Reiss and G.Rümpker

baz = bazdeg*pi/180;

tra = cos(baz)*east-sin(baz)*north;
rad = sin(baz)*east+cos(baz)*north;

end