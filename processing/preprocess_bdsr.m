function preprocess_bdsr(sel_data,nc,station)

%% select events within specified range

sel_data.nc = nc;
sel_data.station = station;

% Copyright 2024 F.Link, J.Wolf and M.Reiss

if sel_data.prfflag
    dtn = 1/40;
elseif sel_data.srfflag
    dtn = 1/20;
else
    dtn = 100/256;
end

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

% check sanity of mseed file
% correct for orientation in metadata
% cut to 100 around phases
% downsample 
% calculate polarization angle from particle motion[
% estimate average misorientation

% get channel information

nomevents = 0;

preflag = 0;
if sel_data.manflag
    filename = [sel_data.save_dir,'/',sel_data.nc,'.',sel_data.station,'pre.mat'];
    if ~exist(filename,'file')
        preflag = 0;
    else
        preflag = 1;
    end
end

if ~preflag
[chan_info] = get_channel_info(strcat(sel_data.work_dir, ...
    '/mseed_files/', char(sel_data.nc),'.',char(sel_data.station), ...
    '/channel_info.txt'));

% prepare structure
data.stat = sel_data.station;
data.nc = sel_data.nc;
data.lat = chan_info(1).lat;
data.lon = chan_info(1).lon;
data.ele = chan_info(1).ele;

if exist([sel_data.work_dir '/' sel_data.misfile],'file')
    fid = fopen([sel_data.work_dir '/' sel_data.misfile],'rt');
    C = textscan(fid,'%s %s %s %f','Delimiter',',');
    fclose(fid);
    stat = C{1,1};
    t1 = C{1,2};
    t2 = C{1,3};
    phi = C{1,4};
    ix = 1:length(stat);
    ioi = ix(strcmp(stat,[nc '.' station]));
    if ~isempty(ioi)
        for j = 1:length(ioi)
            mis.start(j) = datetime(t1{ioi(j)},'InputFormat','yyyy-MM-dd''T''hh:mm:ss.SSSS');
            mis.end(j) = datetime(t2{ioi(j)},'InputFormat','yyyy-MM-dd''T''hh:mm:ss.SSSS');
            mis.phi(j) = phi(ioi(j));
        end
    end
    sel_data.mis = mis;
else
    if sel_data.srfflag
    disp('No misorientation file for this station (search for misorientation with SKS polarization or PRF first')
    return
    end
end

% find files for preprocessing
A = dir([sel_data.work_dir '/mseed_files/' sel_data.nc '.' sel_data.station '/*S.mseed']);
for i = 1:length(A)
    clear phases_to_analyze
    indx = strfind(A(i).name,'_');
    evID = A(i).name(indx(end-1)+1:indx(end)-1);

    % read event from text file & compare with channel info to find out
    % whether miniseed must be rotated because of misalignment
    eventdatevec = [str2double(evID(1:4)) str2double(evID(5:6)) str2double(evID(7:8)) 0 0 0];
    eventnum = datenum(eventdatevec);
    if length(chan_info) > 1
        found_ch = 0;
        for j = 2:length(chan_info)
            if eventnum<datenum(chan_info(j).start_vec)
                choi = j-1;
                found_ch = 1;
                break
            end
        end
        if ~found_ch
            choi = length(chan_info);
        end
    else
        choi = 1;
    end

    % check if noise and signal file exist
    file_readS = strcat(A(i).folder,'/',A(i).name);
    % file_readN = strrep(file_readS,'S.mseed','N.mseed');
    % if ~exist(file_readN,'file')
    %     continue
    % end

    % read miniseed
    % file_mseedS = extract_advv2(file_readS, chan_info(choi));
    file_mseedS = extract_adv(file_readS, chan_info(choi));
    % file_mseedN = extract_adv(file_readN, chan_info(choi));
    

    % check if more than 2 components exist
    if length(file_mseedS(1,:))>2% && length(file_mseedN(1,:))>2

    if isscalar(file_mseedS(1).amp) || isscalar (file_mseedS(2).amp) || isscalar(file_mseedS(3).amp)
        disp('file corrupted')
        continue
    end

    % first sanity check
    % sanchckN = max(sum(abs(file_mseedN(1).amp))./sum(abs(file_mseedN(2).amp)),sum(abs(file_mseedN(2).amp))./sum(abs(file_mseedN(1).amp)));
    sanchckS = max(sum(abs(file_mseedS(1).amp))./sum(abs(file_mseedS(2).amp)),sum(abs(file_mseedS(2).amp))./sum(abs(file_mseedS(1).amp)));
    nn = min(length(file_mseedS(1).amp),length(file_mseedS(2).amp));
    FSe = abs(fft(file_mseedS(2).amp(1:nn)));
    FSn = abs(fft(file_mseedS(1).amp(1:nn)));
    sanchck2S = max(mean(FSe(FSe>max(FSe)./100&FSn>max(FSn)./100)./FSn(FSe>max(FSe)./100&FSn>max(FSn)./100)),...
        mean(FSn(FSe>max(FSe)./100&FSn>max(FSn)./100)./FSe(FSe>max(FSe)./100&FSn>max(FSn)./100)));
    % nn = min(length(file_mseedN(1).amp),length(file_mseedN(2).amp));
    % FSe = abs(fft(file_mseedN(2).amp(1:nn)));
    % FSn = abs(fft(file_mseedN(1).amp(1:nn)));
    % sanchck2N = max(mean(FSe(FSe>max(FSe)./100&FSn>max(FSn)./100)./FSn(FSe>max(FSe)./100&FSn>max(FSn)./100)),...
    %     mean(FSn(FSe>max(FSe)./100&FSn>max(FSn)./100)./FSe(FSe>max(FSe)./100&FSn>max(FSn)./100)));

    

        if ~isempty(file_mseedS(1).comp) && ~isempty(file_mseedS(2).comp) && ...
            ~isempty(file_mseedS(3).comp) %&& ~isempty(file_mseedN(1).comp) && ...
            %~isempty(file_mseedN(2).comp) && ~isempty(file_mseedN(3).comp)

            %if sanchckS < 3 && sanchckN < 3 && sanchck2N < 2 && sanchck2S < 2
            % if sanchckS <3 && sanchck2S < 2 && nn > 100
            if nn > 100

                % prepare event structure
                [event] = read_event_file(sel_data,evID);
                if event.mag<sel_data.mag
                    continue
                end
    
                % resample data
                [sig] = re_sample(event,file_mseedS,dtn);
                % [noi] = re_sample(event,file_mseedN,dtn);
    
                % rotate if necessary
                if chan_info(choi).channel(file_mseedS(1).chID).rot_flag
                    [sig.n_amp,sig.e_amp] = rot_az(sig.n_amp,...
                        sig.e_amp, chan_info(choi).channel(file_mseedS(1).chID).cor_deg);
                    % [noi.n_amp,noi.e_amp] = rot_az(noi.n_amp,...
                    %     noi.e_amp, chan_info(choi).channel(file_mseedS(1).chID).cor_deg);
                end
                % flip Z component if necessary
                if chan_info(choi).channel(file_mseedS(1).chID).cor_dip == 1
                    sig.z_amp = -sig.z_amp;
                end
    
                %get travel times
                % try
                [phases,ttph,ppph,nomph] = get_tp(event.dist,event.depth,P,Pdiff,PKIKP,PKS,PKIKS,S,Sdiff,ScS,SKS,SKKS,SKIKS,PcS,pPcS,PS,PPS);
                % catch
                %     keyboard
                % end
    
                % calculate absolute travel times
                for iF = 1:length(phases)
                    phases(iF).tt_abs = ...
                        seconds(phases(iF).tt) + event.origin_time;
                end
                
                % detrend and filter traces
                % bhw = tukeywin(length(sig.z_amp),0.01);
                % sig.z_amp = detrend(sig.z_amp).*bhw';
                % sig.n_amp = detrend(sig.n_amp).*bhw';
                % sig.e_amp = detrend(sig.e_amp).*bhw';
                if sel_data.prfflag
                    wwin = tukeywin(length(sig.z_amp),0.01);
                    % sig.z_amp = buttern_filter(detrend(sig.z_amp).*wwin',4,0.02,20,dtn)';
                    % sig.n_amp = buttern_filter(detrend(sig.n_amp).*wwin',4,0.02,20,dtn)';
                    % sig.e_amp = buttern_filter(detrend(sig.e_amp).*wwin',4,0.02,20,dtn)';
                    sig.z_amp = buttern_filter(detrend(sig.z_amp).*wwin',4,0.05,20,dtn)';
                    sig.n_amp = buttern_filter(detrend(sig.n_amp).*wwin',4,0.05,20,dtn)';
                    sig.e_amp = buttern_filter(detrend(sig.e_amp).*wwin',4,0.05,20,dtn)';
                    sig.time = sig.time';
                elseif sel_data.srfflag
                    wwin = tukeywin(length(sig.z_amp),0.01);
                    sig.z_amp = buttern_filter(detrend(sig.z_amp).*wwin',4,0.01,20,dtn)';
                    sig.n_amp = buttern_filter(detrend(sig.n_amp).*wwin',4,0.01,20,dtn)';
                    sig.e_amp = buttern_filter(detrend(sig.e_amp).*wwin',4,0.01,20,dtn)';
                    sig.time = sig.time';
                else
                    wwin = tukeywin(length(sig.z_amp),0.01);
                    sig.z_amp = buttern_filter(detrend(sig.z_amp).*wwin',2,1/sel_data.p2,1./sel_data.p1,sig.sr);
                    sig.n_amp = buttern_filter(detrend(sig.n_amp).*wwin',2,1/sel_data.p2,1./sel_data.p1,sig.sr);
                    sig.e_amp = buttern_filter(detrend(sig.e_amp).*wwin',2,1/sel_data.p2,1./sel_data.p1,sig.sr);
                end
                % sig.z_amp = buttern_filter(sig.z_amp,2,...
                %     1/sel_data.p2,1/sel_data.p1,sig.sr);
                % sig.n_amp = buttern_filter(sig.n_amp,2,...
                %     1/sel_data.p2,1/sel_data.p1,sig.sr);
                % sig.e_amp = buttern_filter(sig.e_amp,2,...
                %     1/sel_data.p2,1/sel_data.p1,sig.sr);
                % bhw = tukeywin(length(noi.z_amp),0.01);
                % noi.z_amp = detrend(noi.z_amp).*bhw';
                % noi.n_amp = detrend(noi.n_amp).*bhw';
                % noi.e_amp = detrend(noi.e_amp).*bhw';
                % noi.z_amp = buttern_filter(noi.z_amp,2,...
                %     1/sel_data.p2,1/sel_data.p1,noi.sr);
                % noi.n_amp = buttern_filter(noi.n_amp,2,...
                %     1/sel_data.p2,1/sel_data.p1,noi.sr);
                % noi.e_amp = buttern_filter(noi.e_amp,2,...
                %     1/sel_data.p2,1/sel_data.p1,noi.sr);
                % first quality check
                if sel_data.prfflag
                    % disp(num2str(event.dist))
                    try
                        phases = rc_phasesRF(phases,sig,event,sel_data);
                    catch
                        disp('problem during assessing quality')
                        keyboard
                    end
                elseif sel_data.srfflag
                    phases = rc_phasesSRF(phases,sig,event,sel_data);
                else
                    phases = rc_phases(phases,sig,sel_data.p1,sel_data.p2);
                end
                to_save = 0;
    
                % check if phases fulfill set criteria
                for i_phase = 1:length(phases)
                    if isfield(phases(i_phase),'snr')
                        if ~isempty(phases(i_phase).snr)
                            to_save = to_save +1;
                            phases_to_analyze(to_save) = i_phase;
                        end
                    end
                end
    
                % save phase
                if to_save > 0
                    if sel_data.prfflag
                        data.(['ev' evID]) = phases;
                    elseif sel_data.srfflag
                        for iph = 1:length(phases)
                            data.(['ev' evID '_' phases(iph).name]) = phases(iph);
                        end
                    else
                        event.ttph = ttph;
                        event.pp = ppph;
                        event.nomph = nomph;
                        event.trace = sig;
                        % event.noise = noi;
                        event.phases = phases(phases_to_analyze);
                        data.(['ev' evID]) = event;
                    end
                    nomevents = nomevents+1;
                    disp([num2str(nomevents) '/' num2str(i) ' of ' num2str(length(A))])
                end
            else
                disp('traces too short')
            end

        end
    else
        disp('Less than three components')
    end

end
end
if nomevents>0
    % [clean_mean_miss_baz,clean_std] = baz_stat(data);
    % data.cor_deg = clean_mean_miss_baz;
    % data.cor_deg_err = clean_std;
    if sel_data.manflag
        save([sel_data.save_dir,'/',sel_data.nc,'.',sel_data.station,'pre.mat'],'-struct','data','-v7.3');
    else
        save([sel_data.save_dir,'/',sel_data.nc,'.',sel_data.station,'.mat'],'-struct','data','-v7.3');
    end
else
    if ~preflag
        disp('there are no events for this station')
    end
end

if sel_data.manflag
    clear data
    datain = load([sel_data.save_dir,'/',sel_data.nc,'.',sel_data.station,'pre.mat']);
    data = qc_rfmanual(datain,sel_data);
    save([sel_data.save_dir,'/',sel_data.nc,'.',sel_data.station,'.mat'],'-struct','data','-v7.3');
end

end