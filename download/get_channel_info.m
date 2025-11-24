function [chan_info] = get_channel_info(txtfile)

% function to read infos about station, i.e. channel names & misalignment

% Copyright 2021 F.Link, M.Reiss and G.Rümpker 

% read channel file
fidCI = fopen(txtfile,'r');
formatSpec = '%s';
N = 17;
C_header = textscan(fidCI,formatSpec,N,'Delimiter','|');
C = textscan(fidCI,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','|');
fclose(fidCI);
net = C{1,1};
stat = C{1,2};
loc = C{1,3};
channel = C{1,4};
lat = str2double(C{1,5});
lon = str2double(C{1,6});
ele = str2double(C{1,7});
azimuth_txt = C{1,9};
dip_txt = C{1,10};
startdates = C{1,16};
enddates = C{1,17};

% find all different times spans in which the station settings (i.e. change
% in sensor / data logger / re alignment were changed
for j = 1:length(startdates)
    azimuth(j) = str2double(char(azimuth_txt{j}));
    dip(j) = str2double(char(dip_txt{j}));
    startvec{j} = datevec(char(startdates{j}(1:10)));
    startnum(j) = datenum(startvec{j});
    if isempty(enddates{j})
        endvec{j} = datevec(now);
        endnum(j) = datenum(endvec{j});
    else
        endvec{j} = datevec(char(enddates{j}(1:10)));
        endnum(j) = datenum(endvec{j});
    end
    if isempty(loc{j})
        loc{j} = 'NN';
    end
end

% sort those dates & assign channels
for k = 1:length(channel)
    chID{k} = strcat(loc{k},channel{k}(1:2));
end
check_ch = unique(chID);

jj = 1;
[fin_dates, I, ~] = unique(startdates);
fin_endnums = endnum(I);
% [fin_endnums,~,~] = unique(endnum);
for j = 1:length(fin_dates)
    n = 1;
    startv = datevec(fin_dates{j},'yyyy-mm-ddTHH:MM:SS');
    endv = datevec(fin_endnums(j));
    for i = 1:length(check_ch)
        chan_idx = find(contains(chID,check_ch{i}));
        clear dp az ch lo
        m = 1;
        for k = 1:length(chan_idx)
            if startnum(chan_idx(k))<=datenum((fin_dates{j}),'yyyy-mm-ddTHH:MM:SS') && endnum(chan_idx(k))>=((fin_endnums(j)))
                dp(m) = dip(chan_idx(k));
                az(m) = azimuth(chan_idx(k));
                ch(m) = channel(chan_idx(k));
                lo(m) = loc(chan_idx(k));
                st(m) = stat(chan_idx(k));
                nc(m) = net(chan_idx(k));
                lt(m) = lat(chan_idx(k));
                ln(m) = lon(chan_idx(k));
                m = m+1;
            end
        end
        if m == 1
            continue;
        end
        %    checks whether Z is vertical,
        indx = 1:length(dp);
        [~,I] = max(abs(dp));
        dipZ = dp(I);
        cordip = 0;
        if abs(abs(dipZ)-90) > 20 || std(dp) == 0
            foundflag = 0;
            for ic = 1:length(ch)
                if strcmp(ch{ic}(3),'Z')
                    I = ic;
                    dipZ = -90;
                    foundflag = 1;
                end
            end
            if ~foundflag
                continue
            end
        else
            if dipZ == 90
                cordip = 1;
            end
        end
        channel_Z = ch(I);
        if length(ch) == 3
            % check that N
            % &E components are aligned 90� to each other and check whether they ares
            % switched
            az_oi = az(indx~=I);
            ch_oi = ch(indx~=I);
            az_flag = round(mod(az_oi(1)-az_oi(2),360));
            foundflag = 0;
            switch az_flag
                case 90
                    Nidx = 2;
                    Eidx = 1;

                case 270
                    Nidx = 1;
                    Eidx = 2;
                otherwise
                    foundflag = 0;
                    for ic = 1:2
                        if strcmp(ch{ic}(3),'N')
                            Nidx = ic;
                            azN = 0;
                            foundflag = 1;
                        end
                        if strcmp(ch{ic}(3),'E')
                            Eidx = ic;
                            azE = 90;
                            foundflag = 1;
                        end
                    end
                    if ~foundflag
                        continue
                    end
            end
            if ~foundflag
                azN = az(Nidx);
                azE = az(Eidx);
            end
            channel_N = ch_oi(Nidx);
            channel_E = ch_oi(Eidx);

            % calculate whether components are misaligned
            cor_deg = azN;
            if cor_deg > 180
                cor_deg = cor_deg-360;
            elseif cor_deg < -180
                cor_deg = cor_deg+360;
            end
            if cor_deg == 0
                rot_flag = 0;
            else
                rot_flag = 1;
            end
            if length(channel_Z) > 1
                channel_Z = char(channel_Z{1});
            end
            if length(channel_N) > 1
                channel_N = char(channel_N{1});
            end
            if length(channel_E) > 1
                channel_E = char(channel_E{1});
            end
        end
    
    
        % save infos
        if strcmp(lo{1},'NN')
            lo{1} = '';
        end
        chan_info(jj).net = nc{1};
        chan_info(jj).stat = st{1};
        chan_info(jj).start_vec = startv; % take time vec from Z
        chan_info(jj).end_vec = endv;
        chan_info(jj).lat = lt(1);
        chan_info(jj).lon = ln(1);
        chan_info(jj).ele = ele(1);
        chan_info(jj).channel(n).location = lo(1);
        chan_info(jj).channel(n).chnum = length(ch);
        chan_info(jj).channel(n).channel_Z = channel_Z;
        if length(ch) == 3
            chan_info(jj).channel(n).channel_N = channel_N;
            chan_info(jj).channel(n).channel_E = channel_E;
            chan_info(jj).channel(n).rot_flag = rot_flag;
            chan_info(jj).channel(n).cor_deg = cor_deg;
            chan_info(jj).channel(n).cor_dip = cordip;
        end
        n = n+1;
    end
    if n > 1
        jj = jj+1;
    end
end
end