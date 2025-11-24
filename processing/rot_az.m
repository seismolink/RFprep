function [north1, east1] =  rot_az (north, east, rotdeg)

% rotation of north east components

% usage:
% north & east: time series, same length
% rotdeg: back azimuth in deg

% Copyright 2016 M.Reiss and G.Rümpker

rot = rotdeg*pi/180;

east1 = cos(rot)*east+sin(rot)*north;
north1 = -sin(rot)*east+cos(rot)*north;

end