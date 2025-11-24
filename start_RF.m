function start_RF(path2opts,func,varargin)

% Initialization of RF prepare

% input
%   path2opts: path to file with options for running RFprep
%   func:
%       +download:  	   start download procedure for new data
%       +preprocessing:    start preprocessing and misorientation estimate
%       +makerf:           calculate multi taper receiver functions
%       +plotrf:           plot receiver functions
%       +isoHk:            do Hk stacking with receiver functions
%       +prepHD2D:         prepare input for HaRFE

% Copyright 2025 F.Link

close all
clc

if nargin < 2
    func = '';
end
if nargin < 1
    path2opts = '';
end

% clear all
clearvars -except path2opts func varargin
clear functions
clear global

warning('off','all')

txt = strvcat({'RFprep (auxilliary to HaRFE) - Copyright 2025 F.Link',...
    'This program is free software: you can redistribute it and/or ',...
    'modify it under the terms of the GNU General Public License as ',...
    'published by the Free Software Foundation, either version 3 of ',...
    'the License, or (at your option) any later version.', ' ',...
    'This program is distributed in the hope that it will be useful, ',...
    'but WITHOUT ANY WARRANTY; without even the implied warranty of ',...
    'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the ',...
    'GNU General Public License for more details. ', ' ',...
    'You should have received a copy of the GNU General Public License ',...
    'along with this program. ',...
    'If not, see <http://www.gnu.org/licenses/>.'});

disp(txt);

% get current folder name
cur_dir = pwd;
[~, df, ~] = fileparts(cur_dir); 

% add all sub directories to search path
addpath(genpath(['../',df]));

% Check which actions are requested
dwnld_flag = 0;
prpr_flag = 0;
mkrf_flag = 0;
pltrf_flag = 0;
stckrf_flag = 0;
phd2d_flag = 0;
pckhd2d_flag = 0;
if strcmp(func,'+download')
    dwnld_flag = 1;
end
if strcmp(func,'+preprocessing')
    prpr_flag = 1;
end
if strcmp(func,'+makerf')
    mkrf_flag = 1;
end
if strcmp(func,'+plotrf')
    pltrf_flag = 1;
end
if strcmp(func,'+isoHk')
    stckrf_flag = 1;
end
if strcmp(func,'+prepHD2D')
    phd2d_flag = 1;
end
if strcmp(func,'+pickHD2D')
    pckhd2d_flag = 1;
end

if ~dwnld_flag && ~prpr_flag && ~mkrf_flag && ~pltrf_flag && ~stckrf_flag && ~phd2d_flag && ~pckhd2d_flag
    disp(' ')
    disp(' ')
    disp('RFprep is a program for the calculation of Receiver Functions and Harmonic Decomposition.')
    disp('------------------------------------------------------------------------')
    disp('input')
    disp('  path2opts: path to file with options for running RFprep')
    disp('  func:')
    disp('      +download:  	   start download procedure for new data')
    disp('      +preprocessing:    start preprocessing and misorientation estimate')
    disp('      +makerf:           calculate multi taper receiver functions')
    disp('      +plotrf:           plot receiver functions')
    disp('      +isoHk:            do Hk stacking with receiver functions')
    disp('      +prepHD2D:         prepare input for HaRFE')
    disp('      +pickHD2D:         pick layer model for HaRFE')
end

if ~isempty(varargin)
    idx = 1:length(varargin);
    ind = idx(contains(varargin,'-stat'));
    if ~isempty(ind)
        sel_data.statstr = varargin{ind+1};
    else
        sel_data.statstr = '';
    end
else
    sel_data.statstr = '';
end


% perform requested actions
if dwnld_flag
    bdsr_download(path2opts,sel_data);
end
if prpr_flag
    bdsr_preprocess(path2opts,sel_data);
end
if mkrf_flag
    bdsr_makerf(path2opts,sel_data);
end
if pltrf_flag
    bdsr_plotrf(path2opts,sel_data);
end
if stckrf_flag
    bdsr_stackrf(path2opts,sel_data);
end
if phd2d_flag
    bdsr_prepHD2D(path2opts,sel_data);
end
if pckhd2d_flag
    bdsr_prepHD2Dpick(path2opts,sel_data);
end
end