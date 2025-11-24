function [t_final]=interp_tt(phase,eq_dist,eq_depth)

% table look up for travel times
% usage: 
% phase: phase name
% eq_dist: distance between station and event
% eq_depth: event depth

% Copyright 2016 M.Reiss and G.Rümpker

% define parameters - same matrix size as stored phase tables (in defaults)
dist = 85:0.1:180;
depth = 0:10:700;

M = size(phase);
if M(1) > length(dist)
    dist = 0:0.1:180;
end

eq_dist = eq_dist + 0.0001;
eq_depth = eq_depth + 0.0001;

% find event parameters in vectors
[~, ind_dist] = sort(abs(dist-eq_dist));
[~, ind_depth] = sort(abs(depth-eq_depth));

dist_vec(1) = dist(ind_dist(1));
dist_vec(2) = dist(ind_dist(2));

depth_vec(1) = depth(ind_depth(1));
depth_vec(2) = depth(ind_depth(2));

% find travel times in phase matrix
t1_vec(1)=phase(ind_dist(1),ind_depth(1));
t1_vec(2)=phase(ind_dist(1),ind_depth(2));

t2_vec(1)=phase(ind_dist(2),ind_depth(1));
t2_vec(2)=phase(ind_dist(2),ind_depth(2));

t3_vec(1)=interp1(depth_vec,t1_vec,eq_depth);
t3_vec(2)=interp1(depth_vec,t2_vec,eq_depth);

t_final = interp1(dist_vec,t3_vec,eq_dist);

end