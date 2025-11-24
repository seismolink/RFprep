function [phi,xlam1,xlam2] =  covar(north, east)

% covariace analysis to calculate axes of particle motion and
% BAZ

% usage: north: north component
%        east: east component

% Copyright 2016 M.Reiss and G.Rümpker

N  =  length(north);
co11 = 0;
co12 = 0;
co22 = 0;

% set up matrix
for k = 1:N
    co11 = co11+east(k)*east(k);
    co12 = co12+east(k)*north(k);
    co22 = co22+north(k)*north(k);
end

% calculate eigen values
xp = -co11-co22;
xq = co11*co22-co12^2;
rad = (xp/2)^2 - xq;

xlam1 = -xp/2 + sqrt(rad);
xlam2 = -xp/2 - sqrt(rad);

d22 = co22-xlam1;
d21 = co12;

yvec = -d21/d22;

% calculate BAZ
phi = 90-atan(yvec)*180/pi;

end