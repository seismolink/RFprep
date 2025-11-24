function makesrf_bdsr(sel_data,nc,station)

% calculate receiver functions from processed data

sel_data.nc = nc;
sel_data.station = station;
sel_data.drate=20;
sel_data.tpre=20;
sel_data.tdur=160;
sel_data.tdur2=160;
sel_data.fmax=6.0;
% sel_data.fmax = 24.0;

[sel_data] = calcconst(sel_data);

% Copyright 2024 F.Link, J.Wolf and M.Reiss

% load preprocessed data
if exist([sel_data.save_dir,'/',sel_data.nc,'.',sel_data.station,'.mat'],'file')
    disp('RF already exists skipping station.')
    return
end
data = load([sel_data.load_dir,'/',sel_data.nc,'.',sel_data.station,'.mat']);

% define reference parameters for station
sel_data.lat = data.lat;
sel_data.lon = data.lon;

% read event ID
fn = fieldnames(data);
fn(~contains(fn,'ev')) = [];
% misflag = 0;
vpsflag = 0;
% if exist([sel_data.save_dir '/' sel_data.misfile],'file')
%     fid = fopen([sel_data.save_dir '/' sel_data.misfile],'rt');
%     C = textscan(fid,'%s %s %s %f','Delimiter',',');
%     fclose(fid);
%     stat = C{1,1};
%     t1 = C{1,2};
%     t2 = C{1,3};
%     phi = C{1,4};
%     ix = 1:length(stat);
%     ioi = ix(strcmp(stat,[nc '.' station]));
%     if ~isempty(ioi)
%         misflag = 1;
%         for i = 1:length(ioi)
%             mis.start(i) = datetime(t1{ioi(i)},'InputFormat','yyyy-MM-dd''T''hh:mm:ss.SSSS');
%             mis.end(i) = datetime(t2{ioi(i)},'InputFormat','yyyy-MM-dd''T''hh:mm:ss.SSSS');
%             mis.phi(i) = phi(ioi(i));
%         end
%     end
% end
if exist([sel_data.save_dir '/' sel_data.vpsfile],'file')
    fid = fopen([sel_data.save_dir '/' sel_data.vpsfile],'rt');
    C = textscan(fid,'%s %f','Delimiter',',');
    fclose(fid);
    stat = C{1,1};
    temp = C{1,2};
    ix = 1:length(stat);
    ioi = ix(strcmp(stat,[nc '.' station]));
    if ~isempty(ioi)
        vpsflag = 1;
        for i = 1:length(ioi)
            vpf = temp(ioi(i));
        end
    end
end
% if ~misflag
%     mis = findmisorientation(sel_data,data,fn);
% end
if ~vpsflag
    vsf = findsurfacevs2(sel_data,data,fn);
end

disp("Calculate receiver functions")
% z_targets = [-50,0,50,100,150,200]; % Targeting Depth
z_targets = linspace(-100,250,10);
if sum(z_targets==0) == 0
    z_targets(end+1) = 0;
    z_targets = sort(z_targets);
end
% 1-D velocity model for mapping to depth
ldepth = [20, 35, 300]; % depths of interfaces
velp = [5.8, 6.5, 8]; % P wave velocity
vels = [3.4, 3.75, 4.5]; % S wave velocity

for iev = 1:length(fn)
    phase = data.(fn{iev});
    slow = phase.pp./111.11;
    if length(mis.start)>1
        phi = mis.phi(mis.start<=phase.tt_abs&mis.end>=phase.tt_abs);
    else
        phi = mis.phi;
    end
    if isempty(phi)
        phi = mis.phi(end);
    end
    baz = phase.event.baz-phi;
    rayp = phase.pp./111.11;
    if vpf*rayp > 1
        incl = 89;
    else
        incl = asind(vpf*rayp);
    end
    % phase0 = phase;
    % for i = 1:3
    % phase.tr(i).data = detrend(phase.tr(i).data).*tukeywin(length(phase.tr(i).data),0.005)';
    % phase.tr(i).data = buttern_high(phase.tr(i).data,2,0.1,phase.dt)';
    % end
    % keyboard
    [phase.tr(1).data,phase.tr(2).data,phase.tr(3).data]=rotVNE2VRT(phase.tr(1).data,phase.tr(2).data,phase.tr(3).data,baz);
    [phase.tr(1).data,phase.tr(2).data,phase.tr(3).data]=rotVRT2PSvSh(phase.tr(1).data,phase.tr(2).data,phase.tr(3).data,incl);
    phase.comp = 'lqt';
    for j = 1:length(z_targets)
        z_target = z_targets(j);
        tshift = z2t(z_target,ldepth,velp,vels,slow);
        % [hr,ht,stdevr,stdevt] = MTC(phase,tshift,sel_data);
        [hr,ht,stdevr,stdevt] = MTC2(phase,tshift,sel_data);
        rf.(fn{iev}).rrf(j).data = hr;
        rf.(fn{iev}).trf(j).data = ht;
        rf.(fn{iev}).rstd(j).data = stdevr;
        rf.(fn{iev}).tstd(j).data = stdevt;
        rf.(fn{iev}).tshift(j) = tshift;
        rf.(fn{iev}).zt(j) = z_target;
    end
    rf.(fn{iev}).event = data.(fn{iev}).event;
    rf.(fn{iev}).tt_abs = data.(fn{iev}).tt_abs;
    rf.(fn{iev}).pp = data.(fn{iev}).pp;
    % rf.(fn{iev}).phase = phase;
end
rf.lat = sel_data.lat;
rf.lon = sel_data.lon;
rf.stat = sel_data.station;
rf.net = sel_data.nc;

% data = datao;
save([sel_data.save_dir,'/',sel_data.nc,'.',sel_data.station,'.mat'],'-struct','rf','-v7.3');
disp('Receiver function calculation finished.')

end