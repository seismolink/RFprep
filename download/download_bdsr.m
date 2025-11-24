function download_bdsr(sel_data,nc,station,dci,svec0)

%% select events within specified range

sel_data.nc = nc;
sel_data.station = station;

% Copyright 2024 F.Link, J.Wolf and M.Reiss 

%read eventlist
fileID = fopen([sel_data.defaults_path,'/',sel_data.Eventlist]);
C = textscan(fileID,'%s %f %f %f %f %*[^\n]','HeaderLines',1);
fclose(fileID);

date = C{1,1};
events.lat = C{1,2};
events.lon = C{1,3};
events.depth = C{1,4};
events.mag = C{1,5};

% load travel time files
P = load('P.mat');
Pdiff = load('Pdiff.mat');
PKIKP = load('PKIKP.mat');
PKS = load('PKS.mat');
PKIKS = load('PKIKS.mat');
S = load('S.mat');
Sdiff = load('Sdiff.mat');
ScS = load('ScS.mat');
SKS = load('SKS.mat');
SKKS = load('SKKS.mat');
SKIKS = load('SKIKS.mat');
PcS = load('PcS.mat');
pPcS = load('pPcS.mat');
PS = load('PS.mat');
PPS = load('PPS.mat');

% Create station folder
stationfolder = strcat(sel_data.work_dir,'/mseed_files/',sel_data.nc,'.',sel_data.station);
if ~exist(stationfolder,'dir')
    mkdir(stationfolder)
% else
%     if exist([stationfolder '/channel_info.txt'],'file')
%         disp('Data already exists. skipping station.')
%         return
%     end
end
sel_data.save_dir = stationfolder;

if exist([stationfolder '/channel_info.txt'],'file')
    fileID2 = fopen([sel_data.defaults_path,'/data_centers.dat']);
    C2 = textscan(fileID2,'%s %s');
    fclose(fileID2);
    foundflag = 0;
    try
        [chan_info] = get_channel_info(strcat(stationfolder,'/channel_info.txt'));
        if ~foundflag
            stations.nc = chan_info(1).net;
            stations.name = chan_info(1).stat;
            stations.lat = chan_info(1).lat;
            stations.lon = chan_info(1).lon;
            stations.elevation = chan_info(1).ele;
            stations.done = 0;
            startvec = chan_info(1).start_vec;
            endvec = chan_info(end).end_vec;
            foundflag = 1;
        end
        sel_data.dci = dci;
    catch ME
        foundflag = 0;
    end
else
% read data center and check for applicable
fileID2 = fopen([sel_data.defaults_path,'/data_centers.dat']);
C2 = textscan(fileID2,'%s %s');
fclose(fileID2);
foundflag = 0;
if isempty(dci)
    n = 1;
    for ii = 1:length(C2{1,1})
        sel_data.dc_url = C2{1,2}{ii};
        getchannelinfo(sel_data);
        try
            [chan_info] = get_channel_info(strcat(stationfolder,'/tmp.txt'));
            if ~foundflag
                stations.nc = chan_info(1).net;
                stations.name = chan_info(1).stat;
                stations.lat = chan_info(1).lat;
                stations.lon = chan_info(1).lon;
                stations.elevation = chan_info(1).ele;
                stations.done = 0;
                startvec = chan_info(1).start_vec;
                endvec = chan_info(end).end_vec;
                foundflag = 1;
            end
            sel_data.dci(n) = ii;
            movefile(strcat(stationfolder,'/tmp.txt'),strcat(stationfolder,'/channel_info.txt'));
            n = n+1;
        catch ME
            delete (strcat(stationfolder,'/tmp.txt'));
        end
    end
else
    sel_data.dc_url = C2{1,2}{dci};
    getchannelinfo(sel_data);
    try
        [chan_info] = get_channel_info(strcat(stationfolder,'/tmp.txt'));
        if ~foundflag
            stations.nc = chan_info(1).net;
            stations.name = chan_info(1).stat;
            stations.lat = chan_info(1).lat;
            stations.lon = chan_info(1).lon;
            stations.elevation = chan_info(1).ele;
            stations.done = 0;
            startvec = chan_info(1).start_vec;
            endvec = chan_info(end).end_vec;
            foundflag = 1;
        end
        sel_data.dci = dci;
        movefile(strcat(stationfolder,'/tmp.txt'),strcat(stationfolder,'/channel_info.txt'));
    catch ME
        keyboard
        [chan_info] = get_channel_info(strcat(stationfolder,'/tmp.txt'));
        if ~foundflag
            stations.nc = chan_info(1).net;
            stations.name = chan_info(1).stat;
            stations.lat = chan_info(1).lat;
            stations.lon = chan_info(1).lon;
            stations.elevation = chan_info(1).ele;
            stations.done = 0;
            startvec = chan_info(1).start_vec;
            endvec = chan_info(end).end_vec;
            foundflag = 1;
        end
        sel_data.dci = dci;
        movefile(strcat(stationfolder,'/tmp.txt'),strcat(stationfolder,'/channel_info.txt'));
        % delete (strcat(stationfolder,'/tmp.txt'));
    end
end
end
if ~foundflag
    return
end
if isempty(startvec) || isempty(endvec)
    keyboard
end

% Create event folder
eventfolder = strcat(sel_data.work_dir,'/event_files/');
if ~exist(eventfolder,'dir')
    mkdir(eventfolder)
end
% check which events lie within specified range
if datetime([svec0(1:2) 1]) > datetime([startvec(1:2) 2])
    newstart_date = sprintf('%4d-%02d',svec0(1),svec0(2));
else
    newstart_date = sprintf('%4d-%02d',startvec(1),startvec(2));
end
newend_date = sprintf('%4d-%02d',endvec(1),endvec(2));

% find start and end rows according to specified dates
sd_find = strfind(date,char(newstart_date));
ed_find = strfind(date,char(newend_date));

start_row = find(~cellfun(@isempty,sd_find),1,'first');
end_row = find(~cellfun(@isempty,ed_find),1,'last');
if isempty(end_row)
    end_row = length(ed_find);
end

if ~isempty(start_row)

e_count = 0;
% loop for all events
for n = start_row:end_row
    e_count = e_count +1;
    events.year(n,:) = str2double(date{n,1}(1:4));
    events.month(n,:) = str2double(date{n,1}(6:7));
    events.day(n,:) = str2double(date{n,1}(9:10));
    events.hour(n,:) = str2double(date{n,1}(12:13));
    events.min(n,:) = str2double(date{n,1}(15:16));
    events.sec(n,:) = str2double(date{n,1}(18:23));
    events.ID(n,:) = strrep(strrep(sprintf('%4.0f%02.0f%02.0f%02.0f%02.0f%06.3fs%07.2fs%07.2fs%04.0fs%3.1f',...
        events.year(n,:),events.month(n,:),events.day(n,:),...
        events.hour(n,:),events.min(n,:),events.sec(n,:),events.lat(n),events.lon(n),events.depth(n),events.mag(n)),'.','d'),'-','n');
end
% end_row = start_row;
pp = gcp('nocreate');
if isempty(pp)
    if isempty(getenv('SLURM_CPUS_ON_NODE'))
        nWorkers = feature('NumCores');
    else
        nWorkers = str2double(getenv('SLURM_CPUS_ON_NODE'));
    end
    pp = parpool(nWorkers);
end
parfor n = start_row:end_row

    %check if event above specified magnitude
    if ge((events.mag(n)),sel_data.mag)

        %calculate distance between event and station
        [distance,~,baz]=delaz(events.lat(n),events.lon(n),stations.lat,stations.lon,0);

        if ge(distance,sel_data.min_dist) && le(distance,sel_data.max_dist)

            %get travel times
            [phases] = get_tp(distance,events.depth(n),P,Pdiff,PKIKP,PKS,PKIKS,S,Sdiff,ScS,SKS,SKKS,SKIKS,PcS,pPcS,PS,PPS);
            if ~isstruct(phases)
                continue;
            end
            
            timeref = datetime([events.year(n,:),events.month(n,:),events.day(n,:),events.hour(n,:),events.min(n,:),events.sec(n,:)]);
            flagN = 0;
            flagPKS = 0;
            flagPKKS = 0;
            flagPKIKS = 0;
            flagSKS = 0;
            flagSKKS = 0;
            flagSKIKS = 0;
            filename = strcat(stations.nc,...
                '_',stations.name,'_',...
                sprintf('%s', events.ID(n,:)));
            if exist([filename '_S.mseed'],'file')
                continue;
            end
            timestart2 = datetime('now');
            timeend2 = timeref;
            tt = 3600;
            for ii = 1:length(phases)
                if sel_data.prfflag
                    if distance < 115
                        if strcmp(phases(ii).name,'P') || strcmp(phases(ii).name,'Pdiff')
                            timestart2 = timeref+seconds(phases(ii).tt-200);
                            timeend2 = timeref+seconds(phases(ii).tt+700);
                            flagN = 1;
                            tt = phases(ii).tt-200;
                        end
                    else
                        if strcmp(phases(ii).name,'PKIKP')
                            timestart2 = timeref+seconds(phases(ii).tt-200);
                            timeend2 = timeref+seconds(phases(ii).tt+700);
                            flagN = 1;
                            tt = phases(ii).tt-200;
                        end
                    end
                elseif sel_data.srfflag
                    if distance > 55 && distance < 85
                        if strcmp(phases(ii).name,'S') || strcmp(phases(ii).name,'Sdiff')
                            timestart2 = min(timestart2,timeref+seconds(phases(ii).tt-700));
                            timeend2 = max(timeend2,timeref+seconds(phases(ii).tt+200));
                            flagN = 1;
                            tt = min(tt,phases(ii).tt-700);
                        end
                    end
                    if distance > 50 && distance < 80
                        if strcmp(phases(ii).name,'ScS')
                            timestart2 = min(timestart2,timeref+seconds(phases(ii).tt-700));
                            timeend2 = max(timeend2,timeref+seconds(phases(ii).tt+200));
                            flagN = 1;
                            tt = min(tt,phases(ii).tt-700);
                        end
                    end
                    if distance > 85 && distance < 120
                        if strcmp(phases(ii).name,'SKS')
                            timestart2 = min(timestart2,timeref+seconds(phases(ii).tt-700));
                            timeend2 = max(timeend2,timeref+seconds(phases(ii).tt+200));
                            flagN = 1;
                            tt = min(tt,phases(ii).tt-700);
                        end
                    end
                else
                    % if strcmp(phases(ii).name,'P') || strcmp(phases(ii).name,'Pdiff')
                    %     timestart = timeref+seconds(phases(ii).tt-200);
                    %     timeend = timeref+seconds(phases(ii).tt-50);
                    %     flagN = dwnlddata([filename '_N.mseed'],sel_data,C2,timestart,timeend);
                    %     if ~flagN
                    %         continue
                    %     end
                    % end
                    if strcmp(phases(ii).name,'PKS')
                        timestart2 = min(timestart2,timeref+seconds(phases(ii).tt-75));
                        timeend2 = max(timeend2,timeref+seconds(phases(ii).tt+75));
                        tt = min(tt,phases(ii).tt-75);
                        flagPKS = 1;
                        %flagPKS = dwnlddata([filename '_PKS.mseed'],sel_data,C2,timestart,timeend);
                    end
                    if strcmp(phases(ii).name,'PKKS')
                        timestart2 = min(timestart2,timeref+seconds(phases(ii).tt-75));
                        timeend2 = max(timeend2,timeref+seconds(phases(ii).tt+75));
                        tt = min(tt,phases(ii).tt-75);
                        flagPKKS = 1;
                        %flagPKKS = dwnlddata([filename '_PKKS.mseed'],sel_data,C2,timestart,timeend);
                    end
                    if strcmp(phases(ii).name,'PKIKS')
                        timestart2 = min(timestart2,timeref+seconds(phases(ii).tt-75));
                        timeend2 = max(timeend2,timeref+seconds(phases(ii).tt+75));
                        tt = min(tt,phases(ii).tt-75);
                        flagPKIKS = 1;
                        %flagPKIKS = dwnlddata([filename '_PKIKS.mseed'],sel_data,C2,timestart,timeend);
                    end
                    if strcmp(phases(ii).name,'SKS')
                        timestart2 = min(timestart2,timeref+seconds(phases(ii).tt-75));
                        timeend2 = max(timeend2,timeref+seconds(phases(ii).tt+75));
                        tt = min(tt,phases(ii).tt-75);
                        flagSKS = 1;
                        %flagSKS = dwnlddata([filename '_SKS.mseed'],sel_data,C2,timestart,timeend);
                    end
                    if strcmp(phases(ii).name,'SKKS')
                        timestart2 = min(timestart2,timeref+seconds(phases(ii).tt-75));
                        timeend2 = max(timeend2,timeref+seconds(phases(ii).tt+75));
                        tt = min(tt,phases(ii).tt-75);
                        flagSKKS = 1;
                        %flagSKKS = dwnlddata([filename '_SKKS.mseed'],sel_data,C2,timestart,timeend);
                    end
                    if strcmp(phases(ii).name,'SKIKS')
                        timestart2 = min(timestart2,timeref+seconds(phases(ii).tt-75));
                        timeend2 = max(timeend2,timeref+seconds(phases(ii).tt+75));
                        tt = min(tt,phases(ii).tt-75);
                        flagSKIKS = 1;
                        %flagSKIKS = dwnlddata([filename '_SKIKS.mseed'],sel_data,C2,timestart,timeend);
                    end
                end
            end
            
            if flagN+flagPKS+flagPKKS+flagPKIKS+flagSKS+flagSKKS+flagSKIKS>0
                dwnldflag = dwnlddata([filename '_S.mseed'],sel_data,C2,timestart2,timeend2);
                if ~dwnldflag
                    disp('No XKS-data downloaded')
                    continue
                end
                str2save = sprintf(['%4d %02d %02d %02d %02d %06.3f %8.2f %8.2f ',...
                    '%6.1f %4.1f %s %s %7.2f %7.2f %7.2f %7.2f %5.1f ',...
                    '%8.3f %1d %1d %1d %1d %1d %1d %1d'], ...
                    events.year(n) ,events.month(n),events.day(n),...
                    events.hour(n),events.min(n),events.sec(n),...
                    events.lat(n),events.lon(n),events.depth(n),events.mag(n), ...
                    stations.nc,stations.name,stations.lat,stations.lon,...
                    stations.elevation,distance,baz,tt,flagN,flagPKS,flagPKKS,...
                    flagPKIKS,flagSKS,flagSKKS,flagSKIKS);
                if ~exist([eventfolder '/' events.ID(n,:) '.txt'],'file')
                    fid = fopen([eventfolder '/' events.ID(n,:) '.txt'],'wt');
                else
                    comptxt = fileread([eventfolder '/' events.ID(n,:) '.txt']);
                    if contains(comptxt,str2save)
                        continue;
                    end
                    fid = fopen([eventfolder '/' events.ID(n,:) '.txt'],'at');
                end
                fprintf(fid,'%s\n',str2save);
                fclose(fid);
            end

        end
    end
end

end
% delete(pp);
if length(dir([stationfolder '/*.mseed'])) < 1
    rmdir(stationfolder,'s');
end

end