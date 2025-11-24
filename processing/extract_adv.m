function [trace]= extract_adv(file,chan_info)

% this function reads a miniseed file, checks for gaps, double entries
% a low pass of 1 s is administered to avoid aliasing effects

% usage:
% file: mseed file, chan_info (consist information of misalignment)

% Copyright 2021 F.Link, M.Reiss and G.Rümpker 

% read miniseed file
try
[X,I] = rdmseed2(file);
catch 
    disp('could not read mseed')
    trace = 0;
    return
end

% Check if correct station is in the header
statstr = [chan_info.net ':' chan_info.stat];
if ~strcmp(X(1).ChannelFullName(1:length(statstr)),statstr)
    disp('wrong file downloaded')
    chan_info.net = '';
    chan_info.stat = '';
    chckflag = 0;
else
    chckflag = 1;
end

for i = 1:length(chan_info.channel)
    locch{i} = strcat(chan_info.channel(i).location{1},chan_info.channel(i).channel_Z{1});
    locID(i) = str2double(chan_info.channel(i).location);
    locID(isnan(locID)) = -1;
    len(i) = chan_info.channel(i).chnum;
end
[locch_u,J] = unique(locch);
for i = 1:length(locch_u)
    loc_flag(i) = len(J(i))==3;
end
if isempty(locch_u(loc_flag))
    disp(['Less than three components measured in this period for station ' chan_info.stat])
    keyboard
    trace = 0;
    return
end
indxoi = J(loc_flag);


chanflag = 0;
for i_comp=1:length(I)
    
    clear comp_name
    clear comp
    clear trace0
    clear t
    
    comp_name = I(i_comp).ChannelFullName;
    if ~chckflag
        idx = strfind(comp_name,':');
        comp_name = ['::' comp_name(idx(2)+1:end)];
    end
    if ~chanflag
        for iloc = 1:length(indxoi)
            locstr = chan_info.channel(indxoi(iloc)).location;
            chanstr = [chan_info.net ':' chan_info.stat ':' char(locstr) ':'];
            if strcmp(comp_name(1:end-3),chanstr)
                chanflag = 1;
                break;
            end
        end
    end
    if ~strcmp(comp_name(1:end-3),chanstr)
        continue
    end
    comp = I(i_comp).XBlockIndex;
    
    trace0 = cat(1,X(comp).d);
    trace0 = trace0-mean(trace0);
    dt = 1/X(comp(1)).SampleRate;
    t = cat(1,X(comp).t);
    
    found_flag = 0;
    %% check if gaps exist
    if isempty(I(i_comp).GapBlockIndex)
        
         %% sort time and amp values, delete multiples
        
        [new_t, ind_sort] = sort(t);
        new_trace = trace0(ind_sort);
              
        [~,Index] = unique(new_t,'rows','first');
        t = datetime((new_t(Index))','ConvertFrom','datenum');
        trace0 = (new_trace(Index))';
               
        % low pass butterworth filter
        % trace0 = buttern_low(trace0,6,1,dt);
        % trace0 = buttern_filter(detrend(trace0).*tukeywin(length(trace0),0.1)',2,5,10,dt);
        
        % check if file contains appropriate amount of data points
        time_diff = datevec(t(end)-t(1));
        expec_no_samples = (time_diff(4)*3600 + time_diff(5)*60 + time_diff(6))/dt;
        
        if length(t) > expec_no_samples + 200
            disp('file corrupted')
            trace0 = 0;
            t = 0;
            dt = 0;
        end
        
    else
        disp('data points are missing')
        trace0 = 0;
        t = 0;
        dt = 0;
    end
    % use appropriate channel name
    comp_name = char(comp_name);
    % if ~found_flag
    for ii = 1:length(chan_info.channel)
        compref_Z = [chanstr char(chan_info.channel(ii).channel_Z)];
        compref_N = [chanstr char(chan_info.channel(ii).channel_N)];
        compref_E = [chanstr char(chan_info.channel(ii).channel_E)];
    if strcmp(comp_name,compref_Z)
        trace(3).comp = 'Z';
        trace(3).amp = trace0;
        trace(3).time = t;
        trace(3).dt = dt;
        trace(3).chID = ii;
    elseif strcmp(comp_name,compref_E)
        trace(1).comp = 'East';
        trace(1).amp = trace0;
        trace(1).time = t;
        trace(1).dt = dt;
        trace(1).chID = ii;
    elseif strcmp(comp_name,compref_N)
        trace(2).comp = 'North';
        trace(2).amp = trace0;
        trace(2).time = t;
        trace(2).dt = dt;
        trace(2).chID = ii;
    end
    end
    % else
    %     compref_Z = [chanstr(1:end-2) char(chan_info.channel(L).channel_Z)];
    %     compref_N = [chanstr(1:end-2) char(chan_info.channel(L).channel_N)];
    %     compref_E = [chanstr(1:end-2) char(chan_info.channel(L).channel_E)];
    %     if strcmp(comp_name,compref_Z)
    %         trace(3).comp = 'Z';
    %         trace(3).amp = trace0;
    %         trace(3).time = t;
    %         trace(3).dt = dt;
    %         trace(3).chID = L;
    %     elseif strcmp(comp_name,compref_E)
    %         trace(1).comp = 'East';
    %         trace(1).amp = trace0;
    %         trace(1).time = t;
    %         trace(1).dt = dt;
    %         trace(1).chID = L;
    %     elseif strcmp(comp_name,compref_N)
    %         trace(2).comp = 'North';
    %         trace(2).amp = trace0;
    %         trace(2).time = t;
    %         trace(2).dt = dt;
    %         trace(2).chID = L;
    %     end
    % end
end

if ~exist('trace','var')
    trace(1).comp = 'NotReadable';
    trace(1).amp = 0;
    trace(1).time = 0;
    trace(1).dt = 0;
    trace(1).chID = 0;
    close all
    return
end

close all
fclose('all');

end