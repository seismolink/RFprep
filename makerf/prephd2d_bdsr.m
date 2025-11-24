function prephd2d_bdsr(sel_data)

% calculate receiver functions from processed data

% Copyright 2025 F.Link

% constants for harmonics
ldepth = [20, 35, 300]; % depths of interfaces
velp = [5.8, 6.5, 8]; % P wave velocity
vels = [3.4, 3.75, 4.5]; % S wave velocity
z_out = -25:0.25:150;
zpts = length(z_out);
Nb = 100;
weightflag = sel_data.wflag;

A = dir([sel_data.load_dir '/*.*_rf.mat']);
skipflag = 0;
if isempty(A)
    A = dir([sel_data.save_dir '/*.*orig.mat']);
    statlist = [sel_data.save_dir 'stationlst.txt'];
    fid = fopen(statlist,'rt');
    D = textscan(fid,'%f %f %s %s','Delimiter',',','Headerlines',0);
    fclose(fid);
    lata = D{1,1};
    lona = D{1,2};
    stata = D{1,4};
    neta = D{1,3};
    skipflag = 1;
end
for i = 1:length(A)
    if ~skipflag
        temp = load([A(i).folder '/' A(i).name]);
        todel = sel_data.min_dist>temp.epi|sel_data.max_dist<temp.epi;
        temp.epi(todel) = [];
        temp.baz(todel) = [];
        temp.N(todel) = [];
        temp.slow(todel) = [];
        temp.tr(todel) = [];
    else
        temp.lat = lata(i);
        temp.lon = lona(i);
        temp.stat = stata{i};
        temp.net = neta{i};
    end
    lat(i) = temp.lat;
    lon(i) = temp.lon;
    statall{i} = temp.stat;
    net{i} = temp.net;
    if exist([sel_data.save_dir '/' net{i} '.' statall{i} 'orig.mat'],'file')
        data(i) = load([sel_data.save_dir '/' net{i} '.' statall{i} 'orig.mat']);
    else
        % prepare harmonics
        [hm,errhm] = hm_decon_depth3(temp,z_out,ldepth,velp,vels,Nb,weightflag);
        % [hm2,errhm2] = hm_decon_depth3(temp,z_out,ldepth,velp,vels,Nb,0);
        data(i).hm = hm;
        data(i).errhm = errhm;
        % prepare rf
        tb = 0:15:360;
        iall = 1:length(temp.N);
        n = 1;
        for j = 1:length(tb)-1
            ioi = iall(temp.baz>=tb(j)&temp.baz<=tb(j+1));
            if isempty(ioi)
                continue
            end
            [nflag,ii] = sort(temp.N(ioi),'descend');
            data(i).RF(n).slow = temp.slow(ioi(ii(1)));
            data(i).RF(n).baz = temp.baz(ioi(ii(1)));
            data(i).RF(n).N = nflag(1);
            data(i).RF(n).time = temp.time;
            data(i).RF(n).rad = temp.tr(ioi(ii(1))).rrf;
            data(i).RF(n).tra = temp.tr(ioi(ii(1))).trf;
            n = n+1;
        end
        % [nflag,ii] = sort(temp.N,'descend');
        % keyboard
        % for j = 1:min(18,length(temp.N))
        %     data(i).RF(j).slow = temp.slow(ii(j));
        %     data(i).RF(j).baz = temp.baz(ii(j));
        %     data(i).RF(j).N = nflag(j);
        %     data(i).RF(j).time = temp.time;
        %     data(i).RF(j).rad = temp.tr(ii(j)).rrf;
        %     data(i).RF(j).tra = temp.tr(ii(j)).trf;
        % end
        temp = data(i);
        save([sel_data.save_dir '/' net{i} '.' statall{i} 'orig.mat'],'-struct','temp','-v7.3');
    end
end

% Define starting point of the profile
fig = figure; 
subplot(2,1,1)
plot(lon,lat,'.','MarkerSize',20)
for i = 1:length(statall)
    text(lon(i),min(lat)+0.05.*(mean(lat)-min(lat)),[net{i} '.' statall{i}],'rotation',90)
end
[lonref,latref] = ginput(1);




% disp('The following stations are available as reference: ')
% for i = 1:length(lat)
%     disp([num2str(i) ' ' statall{i} ' ' net{i} ' ' num2str(lat(i)) ' ' num2str(lon(i))])
% end
% ii = input('Station identifier for reference point: ');
% % ii = 1;
% disp(['Station ' statall{ii} ' selected as reference.']);
% latref = lat(ii);
% lonref = lon(ii);

stats = statall;
dd = sqrt((latref-lat).^2+(lonref-lon).^2);
% dd = sqrt((latref-lat).^2);
[dd,I] = sort(dd);
stats = stats(I);
lat = lat(I);
lon = lon(I);
net = net(I);
data = data(I);

dt = 0.025;

zmax1 = 90;
zmax2 = 30;
% profname = 'NESTn';
profname = '';

% Define profile and depth coordinates
zz = -25:0.25:150;
ddeg = 0.001;

n = 1;
i1 = 1;
for i = 1:length(stats)
    HD = data(i).hm;
    % HD = csvread( [filedir,'/',filenameHD] );
    
    % hdt = HD(:,i1:end);
    hdt = HD(1,:);
    if isnan(mean(abs(hdt(:))))
        disp([stats{i} ' bad data'])
        continue
    end
    statr{n} = stats{i};
    statn{n} = stats{i};%sprintf('M%02d',i);

    amp1(n,:) = HD(1,:)./max(abs(hdt(:)));
    amp2(n,:) = HD(2,:)./max(abs(hdt(:)));
    amp3(n,:) = HD(3,:)./max(abs(hdt(:)));
    amp4(n,:) = HD(4,:)./max(abs(hdt(:)));
    amp5(n,:) = HD(5,:)./max(abs(hdt(:)));
    amp1b(n,:) = HD(6,:)./max(abs(hdt(:)));
    amp2b(n,:) = HD(7,:)./max(abs(hdt(:)));
    amp3b(n,:) = HD(8,:)./max(abs(hdt(:)));
    amp4b(n,:) = HD(9,:)./max(abs(hdt(:)));
    amp5b(n,:) = HD(10,:)./max(abs(hdt(:)));
    n = n+1;    
end
subplot(2,1,2)
imagesc(amp1')
for i = 1:length(stats)
    text(i,length(HD(1,:))-20,[net{i} '.' stats{i}],'rotation',90)
end
hold on
button = 0;
uflag = [];
n = 1;
% while button ~= 3
%     [x,~,button] = ginput(1);
%     if button == 1
%         if sum(uflag==round(x)) == 0
%             uflag(end+1) = round(x);
%         else
%             uflag(uflag==round(x)) = [];
%         end
%         if exist('pl','var')
%             delete(pl)
%         end
%         pl = plot(uflag,100,'r.','MarkerSize',20);
%     end
% end
while button ~= 3
    [x,~,button] = ginput(1);
    if button == 1
        if sum(uflag==round(x)) == 0
            uflag(end+1) = round(x);
        else
            uflag(uflag==round(x)) = [];
        end
        if exist('pl','var')
            delete(pl)
        end
        pl = plot(uflag,100,'r.','MarkerSize',20);
    end
end

% stats = stats(sort(uflag));
% stat = stats;
% net = net(sort(uflag));
% lat = lat(sort(uflag));
% lon = lon(sort(uflag));
% data = data(sort(uflag));
% dd = dd(sort(uflag));
% amp1 = amp1(sort(uflag),:);
% amp2 = amp2(sort(uflag),:);
% amp3 = amp3(sort(uflag),:);
% amp4 = amp4(sort(uflag),:);
% amp5 = amp5(sort(uflag),:);
% amp1b = amp1b(sort(uflag),:);
% amp2b = amp2b(sort(uflag),:);
% amp3b = amp3b(sort(uflag),:);
% amp4b = amp4b(sort(uflag),:);
% amp5b = amp5b(sort(uflag),:);

stats = stats((uflag));
stat = stats;
net = net((uflag));
lat = lat((uflag));
lon = lon((uflag));
data = data((uflag));
dd = dd((uflag));
amp1 = amp1((uflag),:);
amp2 = amp2((uflag),:);
amp3 = amp3((uflag),:);
amp4 = amp4((uflag),:);
amp5 = amp5((uflag),:);
amp1b = amp1b((uflag),:);
amp2b = amp2b((uflag),:);
amp3b = amp3b((uflag),:);
amp4b = amp4b((uflag),:);
amp5b = amp5b((uflag),:);

M = size(amp1);
indx = 1:M(1);
ix = round(linspace(1,M(2),20));
amp1o = amp1;
amp2o = amp2;
amp3o = amp3;
amp4o = amp4;
amp5o = amp5;
amp1 = zeros(size(amp1));
amp2 = amp1;
amp3 = amp1;
amp4 = amp1;
amp5 = amp1;

wtot = zeros(size(amp1(1,:)));
for jjj = 1:length(ix)-2
    amp1oo = amp1o(:,ix(jjj):ix(jjj+2));
    amp2oo = amp2o(:,ix(jjj):ix(jjj+2));
    amp3oo = amp3o(:,ix(jjj):ix(jjj+2));
    amp4oo = amp4o(:,ix(jjj):ix(jjj+2));
    amp5oo = amp5o(:,ix(jjj):ix(jjj+2));
    % ampp = [amp1(:,ii:end)./mean(abs(amp1(:,ii:end)),2) amp2(:,ii:end)./mean(abs(amp2(:,ii:end)),2) amp3(:,ii:end)./mean(abs(amp3(:,ii:end)),2) amp4(:,ii:end)./mean(abs(amp4(:,ii:end)),2) amp5(:,ii:end)./mean(abs(amp5(:,ii:end)),2)];
    ampp = [amp1o(:,ix(jjj):ix(jjj+2))./mean(abs(amp1o(:,ix(jjj):ix(jjj+2))),2) amp2o(:,ix(jjj):ix(jjj+2))./mean(abs(amp2o(:,ix(jjj):ix(jjj+2))),2) ...
        amp3o(:,ix(jjj):ix(jjj+2))./mean(abs(amp3o(:,ix(jjj):ix(jjj+2))),2) amp4o(:,ix(jjj):ix(jjj+2))./mean(abs(amp4o(:,ix(jjj):ix(jjj+2))),2) ...
        amp5o(:,ix(jjj):ix(jjj+2))./mean(abs(amp5o(:,ix(jjj):ix(jjj+2))),2)];
    %ampp = [amp2./mean(abs(amp2),2) amp3./mean(abs(amp3),2) amp4./mean(abs(amp4),2) amp5(:,120:end)./mean(abs(amp5),2)];
    %ampp = [amp2./mean(abs(amp2),2) amp3./mean(abs(amp3),2)];% amp4./mean(abs(amp4),2) amp5(:,120:end)./mean(abs(amp5),2)];
    if jjj == 1
        temp = tukeywin(length(amp1oo(1,:))*2,0.5)';
        wwin = temp(length(temp)/2+1:end);
    elseif jjj == length(ix)-2
        temp = tukeywin(length(amp1oo(1,:))*2,0.5)';
        wwin = temp(1:length(temp)/2);
    else
        wwin = tukeywin(length(amp1oo(1,:)),0.5)';
    end
    clear todel
    n = 1;
    for i = 1:M(1)
        % ii = indx(abs(dd-dd(i))<0.135);
        ii = indx(abs(dd-dd(i))<sel_data.wdist);
        temp = (mean(ampp(ii,:),1)+ampp(i,:))./2;
        clear temp3
        for j = 1:length(ii)
            temp2 = corrcoef(temp,ampp(ii(j),:));
            temp3(j) = temp2(1,2);
        end
        w = temp3;
        iw = ii;
        if w(iw==i) < 0
            w = -w;
        end
        w(w<0) = 0;
        w2 = ones(size(amp1oo(1,:)));
        amp1(i,ix(jjj):ix(jjj+2)) = amp1(i,ix(jjj):ix(jjj+2))+(w*(w2.*amp1oo(iw,:))).*wwin;
        amp2(i,ix(jjj):ix(jjj+2)) = amp2(i,ix(jjj):ix(jjj+2))+(w*amp2oo(iw,:)).*wwin;
        amp3(i,ix(jjj):ix(jjj+2)) = amp3(i,ix(jjj):ix(jjj+2))+(w*amp3oo(iw,:)).*wwin;
        amp4(i,ix(jjj):ix(jjj+2)) = amp4(i,ix(jjj):ix(jjj+2))+(w*amp4oo(iw,:)).*wwin;
        amp5(i,ix(jjj):ix(jjj+2)) = amp5(i,ix(jjj):ix(jjj+2))+(w*amp5oo(iw,:)).*wwin;
        if i == 1
            wtot(ix(jjj):ix(jjj+2)) = wtot(ix(jjj):ix(jjj+2))+wwin;
        end
    end
end
for i = 1:M(1)
amp1(i,:) = amp1(i,:)./wtot;
amp2(i,:) = amp2(i,:)./wtot;
amp3(i,:) = amp3(i,:)./wtot;
amp4(i,:) = amp4(i,:)./wtot;
amp5(i,:) = amp5(i,:)./wtot;
end


fnorm = mean(abs([amp1(:,1:end) amp2(:,1:end) amp3(:,1:end) amp4(:,1:end) amp5(:,1:end)]),2);
amp1nn = amp1./(fnorm*ones(1,M(2)));
amp2nn = amp2./(fnorm*ones(1,M(2)));
amp3nn = amp3./(fnorm*ones(1,M(2)));
amp4nn = amp4./(fnorm*ones(1,M(2)));
amp5nn = amp5./(fnorm*ones(1,M(2)));

% Plot exemplary damped harmonics
nn = round(M(1)/2);
nn = 6;
figure; 
hold on
plot(amp1nn(nn,:)./max(abs(amp1nn(nn,:)))*0.7,'Color','blue','LineWidth',1.1); 
plot(amp2nn(nn,:)./max(abs(amp1nn(nn,:)))*0.7+1,'Color','blue','LineWidth',1.1); 
plot(amp3nn(nn,:)./max(abs(amp1nn(nn,:)))*0.7+2,'Color','blue','LineWidth',1.1); 
plot(amp4nn(nn,:)./max(abs(amp1nn(nn,:)))*0.7+3,'Color','blue','LineWidth',1.1); 
plot(amp5nn(nn,:)./max(abs(amp1nn(nn,:)))*0.7+4,'Color','blue','LineWidth',1.1)

plot(amp1o(nn,:)./max(abs(amp1o(nn,:)))*0.7,'Color','red'); 
plot(amp2o(nn,:)./max(abs(amp1o(nn,:)))*0.7+1,'Color','red'); 
plot(amp3o(nn,:)./max(abs(amp1o(nn,:)))*0.7+2,'Color','red'); 
plot(amp4o(nn,:)./max(abs(amp1o(nn,:)))*0.7+3,'Color','red'); 
plot(amp5o(nn,:)./max(abs(amp1o(nn,:)))*0.7+4,'Color','red')

% create red waveform
nk = round(256/6);
nw = round(nk/10);
r = zeros(256, 1);
r(1:nk) = linspace(0.3,1,length(1:nk));
r(nk+1:128-nw) = 1;
r(128-nw+1:128+nw) = 1;
r(129+nw:256-nk) = linspace(1,0,length(129+nw:256-nk));
r(256-(nk-1):256) = 0;
% Create green waveform:
g = zeros(256, 1);
g(1:nk) = 0;
g(nk+1:128-nw) = linspace(0,1,length(nk+1:128-nw));
g(128-nw+1:128+nw) = 1;
g(129+nw:256-nk) = linspace(1,0,length(129+nw:256-nk));
g(256-(nk-1):256) = 0;
% Create blue waveform:
b = zeros(256, 1);
b(1:nk) = 0;
b(nk+1:128-nw) = linspace(0,1,length(nk+1:128-nw));
b(129-nw:128+nw) = 1;
b(129-nw:256-nk) = 1;
b(256-(nk-1):256) = linspace(1,0.3,length(256-(nk-1):256));
% % Now that the individual color channel mappings have been created,
% % stitch the red, green, and blue vectors together to create the 256-by-3 colormap matrix.
customColorMap = [r, g, b];


fnorm = max(abs(amp1nn(:)));
fnorm0 = max(abs(amp1o(:)));
amp1o = amp1o./fnorm0;
amp2o = amp2o./fnorm0;
amp3o = amp3o./fnorm0;
amp4o = amp4o./fnorm0;
amp5o = amp5o./fnorm0;
amp1nn = amp1nn./fnorm;
amp2nn = amp2nn./fnorm;
amp3nn = amp3nn./fnorm;
amp4nn = amp4nn./fnorm;
amp5nn = amp5nn./fnorm;
amp1b = amp1b./fnorm0;
amp2b = amp2b./fnorm0;
amp3b = amp3b./fnorm0;
amp4b = amp4b./fnorm0;
amp5b = amp5b./fnorm0;

% store results
if ~skipflag
fid = fopen([sel_data.save_dir '/stationlst.txt'],'wt');
for i = 1:length(lat)
    fprintf(fid,'%f,%f,%s,%s\n',lat(i),lon(i),net{i},stat{i});
end
fclose(fid);
end
for i = 1:length(data)
    data(i).hm(1,:) = amp1nn(i,:);
    data(i).hm(2,:) = amp2nn(i,:);
    data(i).hm(3,:) = amp3nn(i,:);
    data(i).hm(4,:) = amp4nn(i,:);
    data(i).hm(5,:) = amp5nn(i,:);
    data(i).hm(6,:) = amp1b(i,:);
    data(i).hm(7,:) = amp2b(i,:);
    data(i).hm(8,:) = amp3b(i,:);
    data(i).hm(9,:) = amp4b(i,:);
    data(i).hm(10,:) = amp5b(i,:);
    temp = data(i);
    save([sel_data.save_dir '/' net{i} '.' stat{i} 'weighted.mat'],'-struct','temp','-v7.3');
end

dn = linspace(min(dd)-0.05,max(dd)+0.05,1000);
% for i = 1:M(2)
%     amp1nnn(:,i) = interp1(dd,amp1(:,i),dn);
% end
[x,y] = meshgrid(dd,1:M(2));
[X,Y] = meshgrid(dn,1:M(2));
fnorm = max(abs(amp1nn(:)));
amp1nnn = interp2(x,y,amp1nn'./fnorm,X,Y,'makima')';
fnorm2 = max(abs(amp1nnn),[],2);
amp1nnn = amp1nnn./fnorm2;
amp2nnn = interp2(x,y,amp2nn'./fnorm,X,Y,'makima')'./fnorm2;
amp3nnn = interp2(x,y,amp3nn'./fnorm,X,Y,'makima')'./fnorm2;
amp4nnn = interp2(x,y,amp4nn'./fnorm,X,Y,'makima')'./fnorm2;
amp5nnn = interp2(x,y,amp5nn'./fnorm,X,Y,'makima')'./fnorm2;

fnormo = max(abs(amp1o(:)));
oamp1nnn = interp2(x,y,amp1o'./fnormo,X,Y,'makima')';
fnormo2 = max(abs(oamp1nnn),[],2);
oamp1nnn = oamp1nnn./fnormo2;
oamp2nnn = interp2(x,y,amp2o'./fnormo,X,Y,'makima')'./fnormo2;
oamp3nnn = interp2(x,y,amp3o'./fnormo,X,Y,'makima')'./fnormo2;
oamp4nnn = interp2(x,y,amp4o'./fnormo,X,Y,'makima')'./fnormo2;
oamp5nnn = interp2(x,y,amp5o'./fnormo,X,Y,'makima')'./fnormo2;


FZ = 11;
% Plot all harmonics (sorted in profile order)
clim1 = max(abs(amp1nnn(:))).*1;
clim2 = clim1; clim3 = clim1; clim4 = clim1; clim5 = clim1;
pos = [50 50 850 350];

figcbar = figure; colorbar; colormap(customColorMap); clim([-clim1 clim1]);

% First, print station sorted results
fig1o = figure('Position',pos); 
contourf(dn*111.11,zz,oamp1nnn',40,'LineColor','none'); hold on; plot([min(dn) max(dn)]*111.11,[0 0],'k','LineWidth',1.2); for i = 1:4; plot([min(dn) max(dn)]*111.11,[i*20 i*20],'--','Color',[0 0 0 0.5],'LineWidth',1.2); end
colormap(customColorMap); set(gca,'FontSize',FZ,'YTickLabel',[],'YDir','reverse'); clim([-clim1 clim1]); axis([min(dn)*111.11 max(dn)*111.11 -15 zmax1]); ylabel('Depth in [km]'); xlabel('Profile Coordinate in [km]')
plot(dd*111.11,zeros(size(dd))-2.5,'kv','MarkerFaceColor','k','MarkerSize',10)
text(dd*111.11,zeros(size(dd))-7,stat,'Rotation',30);
fig1b = figure('Position',pos); 
contourf(dn*111.11,zz,amp1nnn',40,'LineColor','none'); hold on; plot([min(dn) max(dn)]*111.11,[0 0],'k','LineWidth',1.2); for i = 1:4; plot([min(dn) max(dn)]*111.11,[i*20 i*20],'--','Color',[0 0 0 0.5],'LineWidth',1.2); end
colormap(customColorMap); set(gca,'FontSize',FZ,'YTickLabel',[],'YDir','reverse'); clim([-clim1 clim1]); axis([min(dn)*111.11 max(dn)*111.11 -15 zmax1]); ylabel('Depth in [km]'); xlabel('Profile Coordinate in [km]')
plot(dd*111.11,zeros(size(dd))-2.5,'kv','MarkerFaceColor','k','MarkerSize',10)
text(dd*111.11,zeros(size(dd))-7,stat,'Rotation',30);

fig2o = figure('Position',pos); 
contourf(dn*111.11,zz,oamp2nnn',40,'LineColor','none'); hold on; plot([min(dn) max(dn)]*111.11,[0 0],'k','LineWidth',1.2); for i = 1:4; plot([min(dn) max(dn)]*111.11,[i*20 i*20],'--','Color',[0 0 0 0.5],'LineWidth',1.2); end
colormap(customColorMap); set(gca,'FontSize',FZ,'YTickLabel',[],'YDir','reverse'); clim([-clim2 clim2]); axis([min(dn)*111.11 max(dn)*111.11 -15 zmax1]); ylabel('Depth in [km]'); xlabel('Profile Coordinate in [km]')
plot(dd*111.11,zeros(size(dd))-2.5,'kv','MarkerFaceColor','k','MarkerSize',10)
text(dd*111.11,zeros(size(dd))-7,stat,'Rotation',30);
fig2b = figure('Position',pos);
contourf(dn*111.11,zz,amp2nnn',40,'LineColor','none'); hold on; plot([min(dn) max(dn)]*111.11,[0 0],'k','LineWidth',1.2); for i = 1:4; plot([min(dn) max(dn)]*111.11,[i*20 i*20],'--','Color',[0 0 0 0.5],'LineWidth',1.2); end
colormap(customColorMap); set(gca,'FontSize',FZ,'YTickLabel',[],'YDir','reverse'); clim([-clim2 clim2]); axis([min(dn)*111.11 max(dn)*111.11 -15 zmax1]); ylabel('Depth in [km]'); xlabel('Profile Coordinate in [km]')
plot(dd*111.11,zeros(size(dd))-2.5,'kv','MarkerFaceColor','k','MarkerSize',10)
text(dd*111.11,zeros(size(dd))-7,stat,'Rotation',30);

fig3o = figure('Position',pos); 
contourf(dn*111.11,zz,oamp3nnn',40,'LineColor','none'); hold on; plot([min(dn) max(dn)]*111.11,[0 0],'k','LineWidth',1.2); for i = 1:4; plot([min(dn) max(dn)]*111.11,[i*20 i*20],'--','Color',[0 0 0 0.5],'LineWidth',1.2); end
colormap(customColorMap); set(gca,'FontSize',FZ,'YTickLabel',[],'YDir','reverse'); clim([-clim3 clim3]); axis([min(dn)*111.11 max(dn)*111.11 -15 zmax1]); ylabel('Depth in [km]'); xlabel('Profile Coordinate in [km]')
plot(dd*111.11,zeros(size(dd))-2.5,'kv','MarkerFaceColor','k','MarkerSize',10)
text(dd*111.11,zeros(size(dd))-7,stat,'Rotation',30);
fig3b = figure('Position',pos);
contourf(dn*111.11,zz,amp3nnn',40,'LineColor','none'); hold on; plot([min(dn) max(dn)]*111.11,[0 0],'k','LineWidth',1.2); for i = 1:4; plot([min(dn) max(dn)]*111.11,[i*20 i*20],'--','Color',[0 0 0 0.5],'LineWidth',1.2); end
colormap(customColorMap); set(gca,'FontSize',FZ,'YTickLabel',[],'YDir','reverse'); clim([-clim3 clim3]); axis([min(dn)*111.11 max(dn)*111.11 -15 zmax1]); ylabel('Depth in [km]'); xlabel('Profile Coordinate in [km]')
plot(dd*111.11,zeros(size(dd))-2.5,'kv','MarkerFaceColor','k','MarkerSize',10)
text(dd*111.11,zeros(size(dd))-7,stat,'Rotation',30);

fig4o = figure('Position',pos); 
contourf(dn*111.11,zz,oamp4nnn',40,'LineColor','none'); hold on; plot([min(dn) max(dn)]*111.11,[0 0],'k','LineWidth',1.2); for i = 1:4; plot([min(dn) max(dn)]*111.11,[i*20 i*20],'--','Color',[0 0 0 0.5],'LineWidth',1.2); end
colormap(customColorMap); set(gca,'FontSize',FZ,'YTickLabel',[],'YDir','reverse'); clim([-clim4 clim4]); axis([min(dn)*111.11 max(dn)*111.11 -15 zmax1]); ylabel('Depth in [km]'); xlabel('Profile Coordinate in [km]')
plot(dd*111.11,zeros(size(dd))-2.5,'kv','MarkerFaceColor','k','MarkerSize',10)
text(dd*111.11,zeros(size(dd))-7,stat,'Rotation',30);
fig4b = figure('Position',pos); 
contourf(dn*111.11,zz,amp4nnn',40,'LineColor','none'); hold on; plot([min(dn) max(dn)]*111.11,[0 0],'k','LineWidth',1.2); for i = 1:4; plot([min(dn) max(dn)]*111.11,[i*20 i*20],'--','Color',[0 0 0 0.5],'LineWidth',1.2); end
colormap(customColorMap); set(gca,'FontSize',FZ,'YTickLabel',[],'YDir','reverse'); clim([-clim4 clim4]); axis([min(dn)*111.11 max(dn)*111.11 -15 zmax1]); ylabel('Depth in [km]'); xlabel('Profile Coordinate in [km]')
plot(dd*111.11,zeros(size(dd))-2.5,'kv','MarkerFaceColor','k','MarkerSize',10)
text(dd*111.11,zeros(size(dd))-7,stat,'Rotation',30);

fig5o = figure('Position',pos); 
contourf(dn*111.11,zz,oamp5nnn',40,'LineColor','none'); hold on; plot([min(dn) max(dn)]*111.11,[0 0],'k','LineWidth',1.2); for i = 1:4; plot([min(dn) max(dn)]*111.11,[i*20 i*20],'--','Color',[0 0 0 0.5],'LineWidth',1.2); end
colormap(customColorMap); set(gca,'FontSize',FZ,'YTickLabel',[],'YDir','reverse'); clim([-clim5 clim5]); axis([min(dn)*111.11 max(dn)*111.11 -15 zmax1]); ylabel('Depth in [km]'); xlabel('Profile Coordinate in [km]')
plot(dd*111.11,zeros(size(dd))-2.5,'kv','MarkerFaceColor','k','MarkerSize',10)
text(dd*111.11,zeros(size(dd))-7,stat,'Rotation',30);
fig5b = figure('Position',pos); 
contourf(dn*111.11,zz,amp5nnn',40,'LineColor','none'); hold on; plot([min(dn) max(dn)]*111.11,[0 0],'k','LineWidth',1.2); for i = 1:4; plot([min(dn) max(dn)]*111.11,[i*20 i*20],'--','Color',[0 0 0 0.5],'LineWidth',1.2); end
colormap(customColorMap); set(gca,'FontSize',FZ,'YTickLabel',[],'YDir','reverse'); clim([-clim5 clim5]); axis([min(dn)*111.11 max(dn)*111.11 -15 zmax1]); ylabel('Depth in [km]'); xlabel('Profile Coordinate in [km]')
plot(dd*111.11,zeros(size(dd))-2.5,'kv','MarkerFaceColor','k','MarkerSize',10)
text(dd*111.11,zeros(size(dd))-7,stat,'Rotation',30);

print(fig1o,[sel_data.save_dir '/orig_ConstantTerm' profname '.jpg'],'-r600','-djpeg');
print(fig2o,[sel_data.save_dir '/orig_FirstOrderTerm_sin' profname '.jpg'],'-r600','-djpeg');
print(fig3o,[sel_data.save_dir '/orig_FirstOrderTerm_cos' profname '.jpg'],'-r600','-djpeg');
print(fig4o,[sel_data.save_dir '/orig_SecondOrderTerm_sin' profname '.jpg'],'-r600','-djpeg');
print(fig5o,[sel_data.save_dir '/orig_SecondOrderTerm_cos' profname '.jpg'],'-r600','-djpeg');

print(fig1b,[sel_data.save_dir '/weight_ConstantTerm' profname '.jpg'],'-r600','-djpeg');
print(fig2b,[sel_data.save_dir '/weight_FirstOrderTerm_sin' profname '.jpg'],'-r600','-djpeg');
print(fig3b,[sel_data.save_dir '/weight_FirstOrderTerm_cos' profname '.jpg'],'-r600','-djpeg');
print(fig4b,[sel_data.save_dir '/weight_SecondOrderTerm_sin' profname '.jpg'],'-r600','-djpeg');
print(fig5b,[sel_data.save_dir '/weight_SecondOrderTerm_cos' profname '.jpg'],'-r600','-djpeg');

print(figcbar,[sel_data.save_dir '/colorbar' profname '.jpg'],'-r600','-djpeg');

% print(fig1o,[sel_data.save_dir '/orig_ConstantTerm' profname '.pdf'],'-r600','-dpdf','-vector','-bestfit');
% print(fig2o,[sel_data.save_dir '/orig_FirstOrderTerm_sin' profname '.pdf'],'-r600','-dpdf','-vector','-bestfit');
% print(fig3o,[sel_data.save_dir '/orig_FirstOrderTerm_cos' profname '.pdf'],'-r600','-dpdf','-vector','-bestfit');
% print(fig4o,[sel_data.save_dir '/orig_SecondOrderTerm_sin' profname '.pdf'],'-r600','-dpdf','-vector','-bestfit');
% print(fig5o,[sel_data.save_dir '/orig_SecondOrderTerm_cos' profname '.pdf'],'-r600','-dpdf','-vector','-bestfit');
% 
% print(fig1b,[sel_data.save_dir '/weight_ConstantTerm' profname '.pdf'],'-r600','-dpdf','-vector','-bestfit');
% print(fig2b,[sel_data.save_dir '/weight_FirstOrderTerm_sin' profname '.pdf'],'-r600','-dpdf','-vector','-bestfit');
% print(fig3b,[sel_data.save_dir '/weight_FirstOrderTerm_cos' profname '.pdf'],'-r600','-dpdf','-vector','-bestfit');
% print(fig4b,[sel_data.save_dir '/weight_SecondOrderTerm_sin' profname '.pdf'],'-r600','-dpdf','-vector','-bestfit');
% print(fig5b,[sel_data.save_dir '/weight_SecondOrderTerm_cos' profname '.pdf'],'-r600','-dpdf','-vector','-bestfit');
% 
% print(figcbar,[sel_data.save_dir '/colorbar' profname '.pdf'],'-r600','-dpdf','-vector');

close all


end