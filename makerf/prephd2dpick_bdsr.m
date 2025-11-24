function prephd2dpick_bdsr(sel_data)

% calculate receiver functions from processed data

% Copyright 2025 F.Link

if sel_data.origflag==0
    addname = 'orig';
elseif sel_data.origflag == 1
    addname = 'weighted';
else
    addname = '';
end

% load station information
fid = fopen([sel_data.load_dir '/stationlst.txt'],'rt');
D = textscan(fid,'%f %f %s %s','Delimiter',',','Headerlines',0);
fclose(fid);
% station = C{1,1};
lat = D{1,1};
lon = D{1,2};
net = D{1,3};
stat = D{1,4};
% indx = 1:length(stat);

% load harmonics from path and station list
for ii = 1:length(stat)
    A = dir([sel_data.load_dir '/' net{ii} '.' stat{ii} addname '.mat']);
    
    load([sel_data.load_dir '/' A(1).name]);
    
    hm(6:end,:) = [];
    bootm0 = hm'./max(abs(hm(1,:)));
    hdt = bootm0(:);
    
    amp1(ii,:) = bootm0(:,1)./max(abs(hdt(:)));
    amp2(ii,:) = bootm0(:,2)./max(abs(hdt(:)));
    amp3(ii,:) = bootm0(:,3)./max(abs(hdt(:)));
    amp4(ii,:) = bootm0(:,4)./max(abs(hdt(:)));
    amp5(ii,:) = bootm0(:,5)./max(abs(hdt(:)));
end

% preparation of parameters
ampa = (abs(amp1)./max(abs(amp1(:)))).^0.25+(abs(amp2)./max(abs(amp2(:)))).^0.25+(abs(amp3)./max(abs(amp3(:)))).^0.25+(abs(amp4)./max(abs(amp4(:)))).^0.25+(abs(amp5)./max(abs(amp5(:)))).^0.25;
ampa = smoothdata(ampa,2);
ampa = ampa./max(abs(ampa),[],2);

xvec = 1:length(lat);
zvec = -25:0.25:150;
%zvec = -2:0.025:12;
xlim = [xvec(1)-0.5 xvec(end)+0.5];
zlim = [zvec(1) zvec(end)];
selmat = zeros(size(amp1));
M = size(selmat);
mohovec = zeros(M(1),1);

preflag = 0;
if exist([sel_data.save_dir '/LayerModel.txt'],'file')
    preflag = 1;
    fid = fopen([sel_data.save_dir '/LayerModel.txt'],'rt');
    line = fgetl(fid);
    while ischar(line)
        ix = strfind(line,',');
        format = ['%s ' repmat('%f ',1,length(ix)-1), '%f'];
        C = textscan(line,format,'Delimiter',',');

        temps = char(C{1,1});
        tempm = C{1,2};
        templ = [];
        for i = 2:length(ix)-2
            if C{1,i+1} == 0
                break;
            end
            templ(i-1) = C{1,i+1};
        end
        mohovec(contains(stat,temps)) = tempm;
        for i = 1:length(templ)
            [~,I] = min((templ(i)-zvec).^2);
            selmat(contains(stat,temps),I) = 1;
        end
        line = fgetl(fid);
    end
    fclose(fid);
end


sm = "hexagram";
cm = [1 1 1];
ss = ["o","square","diamond","^","v",">","<","pentagram","+","*",".","x","_","|"];
cs = [   0         0    1.0000;
         0    0.7008         0;
    1.0000         0    0.4444;
    0.5873         0         0;
    0.5873    1.0000         0;
    0.1587    0.5512    1.0000;
         0         0    0.4762;
    1.0000    0.6772    0.0476;
    0.8413         0    1.0000;
    0.5873    0.4882    0.4286;
    1.0000    0.6142    0.9683;
    0.1111    0.8583    0.8730;
    0.1746    0.0866         0;
    0.4603    0.0709    0.6190;
         0    0.3701         0];
   
 % create red waveform
nk = round(256/12);
nw = round(nk/3);
r = zeros(256, 1);
r(1:nk) = linspace(0.6,1,length(1:nk));
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
b(256-(nk-1):256) = linspace(1,0.6,length(256-(nk-1):256));
% % Now that the individual color channel mappings have been created,
% % stitch the red, green, and blue vectors together to create the 256-by-3 colormap matrix.
customColorMap = [r, g, b];

% produce initial guess
if ~preflag
minzi = find(zvec>10,1);
for i = 1:(M(1))
    [~,iz] = max((amp1(i,minzi:end)));
    mohovec(i) = zvec(iz+minzi);
end
end

% Interactive Plot and Selection
finished_flag = 0;
while ~finished_flag
    anyinput = input('PI>>','s');
    if strcmp(anyinput,'quit')
        finished_flag = 1;
        continue
    end

    if strcmp(anyinput,'plotall') || strcmp(anyinput,'pa')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,ampa')
        colormap(jet);
        caxis([min(ampa(:)) max(ampa(:))])
        set(gca,'Ydir','reverse')
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
        % plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','r','LineWidth',1.1)
        % plot(xvec(indxn),zeros(length(indxn),1),'w*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:M(1)
            temp = selmat(i,:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','w','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec);
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(mpl),sm,'Color','w','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = 0;
        continue
    end

    if strcmp(anyinput,'plotamp1') || strcmp(anyinput,'p1')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,amp1')
        set(gca,'Ydir','reverse')
        colormap(customColorMap);
        caxis([-max(abs(amp1(:))) max(abs(amp1(:)))])
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
        % plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','g','LineWidth',1.1)
        % plot(xvec(indxn),zeros(length(indxn),1),'g*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:M(1)
            temp = selmat(i,:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec);
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(mpl),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = 1;
        continue
    end

    if strcmp(anyinput,'plotamp2') || strcmp(anyinput,'p2')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,amp2')
        set(gca,'Ydir','reverse')
        colormap(customColorMap);
        caxis([-max(abs(amp2(:))) max(abs(amp2(:)))])
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
        % plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','g','LineWidth',1.1)
        % plot(xvec(indxn),zeros(length(indxn),1),'g*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:M(1)
            temp = selmat(i,:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec);
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(mpl),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = 2;
        continue
    end

    if strcmp(anyinput,'plotamp3') || strcmp(anyinput,'p3')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,amp3')
        set(gca,'Ydir','reverse')
        colormap(customColorMap);
        caxis([-max(abs(amp3(:))) max(abs(amp3(:)))])
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
        % plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','g','LineWidth',1.1)
        % plot(xvec(indxn),zeros(length(indxn),1),'g*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:M(1)
            temp = selmat(i,:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec);
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(mpl),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = 3;
        continue
    end

    if strcmp(anyinput,'plotamp4') || strcmp(anyinput,'p4')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,amp4')
        set(gca,'Ydir','reverse')
        colormap(customColorMap);
        caxis([-max(abs(amp4(:))) max(abs(amp4(:)))])
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
        % plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','g','LineWidth',1.1)
        % plot(xvec(indxn),zeros(length(indxn),1),'g*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:M(1)
            temp = selmat(i,:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec);
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(mpl),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = 4;
        continue
    end

    if strcmp(anyinput,'plotamp5') || strcmp(anyinput,'p5')
        if ~exist('figmain','var')
            figmain = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
        else
            figure(figmain);
        end
        hold on
        imagesc(xvec,zvec,amp5')
        set(gca,'Ydir','reverse')
        colormap(customColorMap);
        caxis([-max(abs(amp5(:))) max(abs(amp5(:)))])
        plot([xvec(1)-0.5 xvec(end)+0.5],[0 0],'k','LineWidth',1.1)
        axis([xlim zlim])
        % plot([xvec(ig)' xvec(ig)']'+0.5,[zvec(1).*ones(size(xvec(ig)')) zvec(end).*ones(size(xvec(ig)'))]','g','LineWidth',1.1)
        % plot(xvec(indxn),zeros(length(indxn),1),'g*','MarkerSize',8,'LineWidth',1.1)
        n = 1;
        for i = 1:M(1)
            temp = selmat(i,:);
            it = find(temp);
            if ~isempty(it)
                for j = 1:length(it)
                    pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                    n = n+1;
                end
            end
        end
        mpl = find(mohovec);
        if ~isempty(mpl)
            plm = plot(mpl,mohovec(mpl),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
        end
        hold off
        opt = 5;
        continue
    end
    
    if strcmp(anyinput,'zoom') || strcmp(anyinput,'z')
        xlimt = input('Provie new limits for x (two entries): ','s');
        zlimt = input('Provide new limits for z (two entries): ','s');
        I = strfind(xlimt,' ');
        J = strfind(zlimt,' ');
        if isempty(I) || isempty(J)
            disp('more than two entries in the format: lim1 lim2')
            continue
        end
        xlim(1) = str2double(xlimt(1:I-1));
        xlim(2) = str2double(xlimt(I+1:end));
        zlim(1) = str2double(zlimt(1:J-1));
        zlim(2) = str2double(zlimt(J+1:end));
        if round(xlim(1)-0.1) == xlim(1)
            xlim(1) = xlim(1)-0.5;
        end
        if round(xlim(2)+0.1) == xlim(2)
            xlim(2) = xlim(2)+0.5;
        end
        figure(figmain);
        axis([xlim zlim])
        continue
    end

    if strcmp(anyinput,'depthzoom') || strcmp(anyinput,'dz')
        zlimt = input('Provide new limits for z (two entries): ','s');
        J = strfind(zlimt,' ');
        if isempty(J)
            disp('more than two entries in the format: lim1 lim2')
            continue
        end
        zlim(1) = str2double(zlimt(1:J-1));
        zlim(2) = str2double(zlimt(J+1:end));
        figure(figmain);
        axis([xlim zlim])
        continue
    end

    if strcmp(anyinput,'reset') || strcmp(anyinput,'r')
        zlim = [min(zvec) max(zvec)];
        xlim = [min(xvec)-0.5 max(xvec)+0.5];
        figure(figmain);
        axis([xlim zlim])
        continue
    end

    if strcmp(anyinput,'guizoom') || strcmp(anyinput,'gz')
        figure(figmain);
        hold on
        [xlim(1),zlim(1)] = ginput(1);
        pl2 = plot(xlim(1),zlim(1),'g+','LineWidth',1.1);
        [xlim(2),zlim(2)] = ginput(1);
        xlim = sort(xlim);
        zlim = sort(zlim);
        axis([xlim zlim])
        delete(pl2);
        hold off
        continue
    end

    if strcmp(anyinput,'setmoho') || strcmp(anyinput,'m')
        exitflag = 0;
        nn = 1;
        mohovecnew = mohovec;
        figure(figmain);
        hold on
        while ~exitflag
            [xtemp,ztemp,button] = ginput(1);
            if button == 1
                xn(nn) = round(xtemp);
                [~,I] = min((zvec-ztemp).^2);
                zn(nn) = zvec(I);
                temp = mohovec(xn(nn));
                if temp==0
                    pl2(nn) = plot(xn(nn),zn(nn),'g+');
                    mohovecnew(xn(nn)) = zn(nn);
                else
                    pl2(nn) = plot(xn(nn),mohovecnew(xn(nn)),'kx');
                    mohovecnew(xn(nn)) = zn(nn);
                end
                nn = nn+1;
            elseif button == 3
                exitflag = 1;
            end
        end
        if exist('pl2','var')
            acceptflag = input('Accept new layer assignments? 1=yes,0=no');
            if ~acceptflag
                delete(pl2)
            else
                delete(pl2)
                if exist('plm','var')
                delete(plm)
                end
                mohovec = mohovecnew;
                mpl = find(mohovec);
                if opt == 0
                    if ~isempty(mpl)
                        plm = plot(mpl,mohovec(mpl),sm,'Color','w','MarkerFaceColor',cm,'MarkerSize',6);
                    end
                else
                    if ~isempty(mpl)
                        plm = plot(mpl,mohovec(mpl),sm,'Color','k','MarkerFaceColor',cm,'MarkerSize',6);
                    end
                end
            end
        end
        hold off
        continue
    end

    if strcmp(anyinput,'setlayers') || strcmp(anyinput,'s')
        exitflag = 0;
        nn = 1;
        selmatnew = selmat;
        figure(figmain);
        hold on
        while ~exitflag
            [xtemp,ztemp,button] = ginput(1);
            if button == 1
                xn(nn) = round(xtemp);
                [~,I] = min((zvec-ztemp).^2);
                zn(nn) = zvec(I);
                iz = find(zvec==zn(nn),1);
                temp = selmat(xn(nn),:);
                izt = find(temp(iz-5:min(length(zvec),iz+5)));
                if isempty(izt)
                    pl2(nn) = plot(xn(nn),zn(nn),'g+');
                    selmatnew(xn(nn),iz) = 1;
                else
                    izt = izt(1);
                    pl2(nn) = plot(xn(nn),zvec(iz-5+izt-1),'kx');
                    selmatnew(xn(nn),iz-5+izt-1) = 0;
                end
                nn = nn+1;
            elseif button == 3
                exitflag = 1;
            end
        end
        if exist('pl2','var')
        acceptflag = input('Accept new layer assignments? 1=yes,0=no');
        if ~acceptflag
            delete(pl2)
        else
            delete(pl2)
            if exist('pl','var')
            delete(pl)
            end
            selmat = selmatnew;
            if opt == 0
                n = 1;
                for i = 1:M(1)
                    temp = selmat(i,:);
                    it = find(temp);
                    if ~isempty(it)
                        for j = 1:length(it)
                            pl(n) = plot(i,zvec(it(j)),ss(j),'Color','w','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                            n = n+1;
                        end
                    end
                end
            else
                n = 1;
                for i = 1:M(1)
                    temp = selmat(i,:);
                    it = find(temp);
                    if ~isempty(it)
                        for j = 1:length(it)
                            pl(n) = plot(i,zvec(it(j)),ss(j),'Color','k','MarkerFaceColor',cs(j,:),'MarkerSize',6);
                            n = n+1;
                        end
                    end
                end
            end
        end
        end
        hold off
        continue
    end
    
    if strcmp(anyinput,'help')
        disp('List of commands:');
        disp('    - quit          : close program')
        disp('    - plotall/pa       : plot figure with squared and summed harmonics')
        disp('    - plotamp1/p1      : plot figure with constant term')
        disp('    - plotamp2/p2      : plot figure with cos(baz) term')
        disp('    - plotamp3/p3      : plot figure with sin(baz) term')
        disp('    - plotamp4/p4      : plot figure with cos(2*baz) term')
        disp('    - plotamp5/p5      : plot figure with sin(2*baz) term')
        disp('    - zoom/z           : zoom into the figure by entering x and z limits')
        disp('    - guizoom/gz       : zoom into the figure by graphically setting boundary box')
        disp('    - setmoho/m        : assign or correct moho picks')
        disp('    - setlayer/s       : assign layers or delete assignments')
        disp('    - reset/r          : reset zoom')
        continue
    end

    try 
        eval(anyinput)
    catch
        disp('No matlab or PI command')
    end
end

nn = 1;
nlay = [];
for i = 1:M(1)
    temp = selmat(i,:);
    if mohovec(i) == 0
        continue
    end
    if sum(temp) == 0 
        continue
    else
        stats{nn} = stat{(i)};
        mohoz(nn) = mohovec((i));
        ix = find(temp);
        layt{nn} = zvec(ix);
        if isempty(nlay)
            nlay = length(ix);
        else
            nlay = max(nlay,length(ix));
        end
        nn = nn+1;
    end
end
lays = zeros(length(stats),nlay);
for i = 1:length(stats)
    lays(i,1:length(layt{i})) = layt{i};
end
filename = 'LayerModel.txt';
fid = fopen([sel_data.save_dir '\' filename],'wt');
format = ['%s,' repmat('%f,',1,nlay) '%f,%f,%f\n'];
for i = 1:length(stats)
    fprintf(fid,format,stats{i},mohoz(i),lays(i,:),1,1);
end
fclose(fid);

close all
end