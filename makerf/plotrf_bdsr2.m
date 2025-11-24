function plotrf_bdsr2(sel_data,nc,station)


plotbazflag = sel_data.plbaz;
plotepiflag = sel_data.plepi;
plothdflag = sel_data.plhd;
sel_data.nc = nc;
sel_data.station = station;
% Copyright 2024 F.Link, J.Wolf and M.Reiss

% load preprocessed data
rf = load([sel_data.load_dir,'/',sel_data.nc,'.',sel_data.station,'_rf.mat']);

% read event ID
% fn = fieldnames(rf);
% fn(~contains(fn,'ev')) = [];

if plotbazflag
% back azimuthal sweep
% compute 36 composite RFs spaced:  baz = 5,355,10
ddeg = 10;
% ddeg = 25;

tt = rf.time;
mint = -5;
maxt = 10;
ntot = 0;
fig = figure;
subplot(1,2,1)
for ii = 1:round(360/ddeg)%36
    fnorm(ii) = 0;
    % baz = 5+10*(ii-1);
    baz = ddeg/2+ddeg*(ii-1);
    jrec = 0;
    rff = zeros(size(rf.tr(1).rrf));
    for jj = 1:length(rf.baz)
        if rf.N(jj) < 2
            continue
        end
        ebaz = rf.baz(jj);
        test = abs(baz-ebaz);
        if test<ddeg || test >360-ddeg
            if isnan(max(abs(rf.tr(jj).rrf)))
                continue
            end
            ntot = ntot+1;
            jrec = jrec+1;
            rff = rff + rf.tr(jj).rrf;
        end
    end
    if jrec > 0
        fnorm(ii) = max(abs(rff));
        bbb = baz.*ones(size(rff))+10*rff./fnorm(ii);
        plot(tt,bbb,'Color','black','LineWidth',0.5);
        hold on
        localzero = baz;
        fill_trace(tt',bbb',localzero,0,'b');
        fill_trace(tt',bbb',localzero,1,'r');
    end
end
plot([0 0],[-10.,370],'k','LineWidth',0.1);
title('Radial RF')
xlabel("Delay Time (sec)")
ylabel("Back Azimuth")
xlim([mint,maxt])
ylim([-10.,370])
disp(ntot);

subplot(1,2,2)
for ii = 1:round(360/ddeg)
    % baz = 5+10*(ii-1);
    baz = ddeg/2+ddeg*(ii-1);
    jrec = 0;
    rff = zeros(size(rf.tr(1).rrf));
    for jj = 1:length(rf.baz)
        if rf.N(jj) < 2
            continue
        end
        ebaz = rf.baz(jj);
        test = abs(baz-ebaz);
        if test<ddeg || test >360-ddeg
            if isnan(max(abs(rf.tr(jj).trf)))
                continue
            end
            ntot = ntot+1;
            jrec = jrec+1;
            rff = rff + rf.tr(jj).trf;
        end
    end
    if jrec > 0
        bbb = baz.*ones(size(rff))+10*rff./fnorm(ii);
        plot(tt,bbb,'Color','black','LineWidth',0.5);
        hold on
        localzero = baz;
        fill_trace(tt',bbb',localzero,0,'b');
        fill_trace(tt',bbb',localzero,1,'r');
    end
end
plot([0 0],[-10.,370],'k','LineWidth',0.1);
title('Transverse RF')
xlabel("Delay Time (sec)")
% ylabel("Back Azimuth")
xlim([mint,maxt])
ylim([-10.,370])
fname=[sel_data.save_dir '/' sel_data.nc '.' sel_data.station '-' '-rfbaz2.png'];
print(fig,fname,'-dpng','-r600');
disp("Azimuthal moveout printed")
close(fig);
end

if plotepiflag
% epicentral sweep
% compute composite RFs spaced:  epi = 25,155,10
depi = 10;
epis = 25:depi:155;

tt = rf.time;
mint = -5;
maxt = 10;
ntot = 0;
fig = figure;
subplot(1,2,1)
for ii = 1:length(epis)
    fnorm(ii) = 0;
    epi = 25+depi*(ii-1);
    jrec = 0;
    nev = 0;
    rff = zeros(size(rf.tr(1).rrf));
    rft = rff;
    for jj = 1:length(rf.epi)
        if rf.N(jj) < 2
            continue
        end
        eepi = rf.epi(jj);
        test = abs(epi-eepi);%abs(baz-ebaz);
        if test<depi
            if isnan(max(abs(rf.tr(jj).rrf)))
                continue
            end
            ntot = ntot+1;
            jrec = jrec+1;
            rff = rff + rf.tr(jj).rrf;
        end
        if ii == 1
            if isnan(max(abs(rf.tr(jj).rrf)))
                continue
            end
            rft = rft + rf.tr(jj).rrf;
        end
    end
    if jrec > 0
        fnorm(ii) = max(abs(rff));
        bbb = epi.*ones(size(rff))+10*rff./fnorm(ii);
        plot(tt,bbb,'Color','black','LineWidth',0.5);
        hold on
        localzero = epi;
        fill_trace(tt',bbb',localzero,0,'b');
        fill_trace(tt',bbb',localzero,1,'r');
    end
    if ii == 1
        fnormsum = max(abs(rft));
        bbb = 170.*ones(size(rft))+10*rft./fnormsum;
        plot(tt,bbb,'Color','black','LineWidth',0.5);
        hold on
        localzero = 170;
        fill_trace(tt',bbb',localzero,0,'b');
        fill_trace(tt',bbb',localzero,1,'r');
    end
end
plot([0 0],[15.,180],'k','LineWidth',0.1);
title('Radial RF')
xlabel("Delay Time (sec)")
ylabel("Epicentral Distance")
xlim([mint,maxt])
ylim([15.,180])
disp(ntot);

subplot(1,2,2)
for ii = 1:length(epis)
    epi = 25+depi*(ii-1);
    jrec = 0;
    nev = 0;
    rff = zeros(size(rf.tr(1).rrf));
    rft = rff;
    for jj = 1:length(rf.epi)
        if rf.N(jj) < 2
            continue
        end
        eepi = rf.epi(jj);
        test = abs(epi-eepi);%abs(baz-ebaz);
        if test<depi
            if isnan(max(abs(rf.tr(jj).trf)))
                continue
            end
            ntot = ntot+1;
            jrec = jrec+1;
            rff = rff + rf.tr(jj).trf;
        end
        if ii == 1
            if isnan(max(abs(rf.tr(jj).trf)))
                continue
            end
            rft = rft + rf.tr(jj).trf;
        end
    end
    if jrec > 0
        bbb = epi.*ones(size(rff))+10*rff./fnorm(ii);
        plot(tt,bbb,'Color','black','LineWidth',0.5);
        hold on
        localzero = epi;
        fill_trace(tt',bbb',localzero,0,'b');
        fill_trace(tt',bbb',localzero,1,'r');
    end
    if ii == 1
        bbb = 170.*ones(size(rft))+10*rft./fnormsum;
        plot(tt,bbb,'Color','black','LineWidth',0.5);
        hold on
        localzero = 170;
        fill_trace(tt',bbb',localzero,0,'b');
        fill_trace(tt',bbb',localzero,1,'r');
    end
end
plot([0 0],[15.,180],'k','LineWidth',0.1);
title('Transverse RF')
xlabel("Delay Time (sec)")
% ylabel("Back Azimuth")
xlim([mint,maxt])
ylim([15,180])

fname=[sel_data.save_dir '/' sel_data.nc '.'  sel_data.station '-' '-rfepi2.png'];
print(fig,fname,'-dpng','-r600');
disp("Epicentral moveout printed")
close(fig);
end

if plothdflag
% 1-D velocity model for mapping to depth
ldepth = [20, 35, 300]; % depths of interfaces
velp = [5.8, 6.5, 8]; % P wave velocity
vels = [3.4, 3.75, 4.5]; % S wave velocity
z_out = -25:0.25:150;
zpts = length(z_out);
Nb = 100;
% weightflag = 1;
weightflag = sel_data.wflag;

% hm = hm_decon_depth2(rf,z_out,ldepth,velp,vels);
[hm,errhm] = hm_decon_depth3(rf,z_out,ldepth,velp,vels,Nb,weightflag);
disp(['Harmonics calculated'])
hm = hm./max(abs(hm(:)));

fac = 1;%0.33;

fig = figure;
subplot(1,2,1)
for ii = 1:5
    localzero = (4-ii+1)*10.0;
    rft = hm(ii,:);
    erft = errhm(ii,:);
    bbb = localzero*ones(1,zpts)+10*rft.*0.5;

    bbbmin = localzero*ones(1,zpts)+10*(rft-fac*erft).*0.5;
    bbbmax = localzero*ones(1,zpts)+10*(rft+fac*erft).*0.5;

    plot(z_out,bbb,'Color','black','LineWidth',0.01);
    hold on
    fill_trace(z_out,bbb,localzero,0,'b');
    fill_trace(z_out,bbb,localzero,1,'r');
    patch([z_out flip(z_out)],[bbbmin flip(bbbmax)],[0.5 0.5 0.5],'EdgeColor','none')
    plot(z_out,bbb,'Color','black','LineWidth',0.1);
end
plot([0.,0.],[0.,45.0],'k')
zz = 10:10:160;
for i = 1:length(zz)
    plot([zz(i),zz(i)],[0.,45.],'k--','LineWidth',0.5)
end
title(['Anisotropy/Dip RFs'])
xlabel('Depth (km)')
xlim([-20,100])
xticks(-20:20:120)
ylim([-5.0,50.0])
yticks([0,10,20,30,40])
yticklabels(["sin(2*baz)","cos(2*baz)","sin(baz)","cos(baz)","constant"])


subplot(1,2,2)
for iii = 1:5
    localzero = (4-iii+1)*10.0;
    ii=iii+5;
    rft = hm(ii,:);
    erft = errhm(ii,:);
    bbb = localzero*ones(1,zpts)+10.*rft*0.5;

    bbbmin = localzero*ones(1,zpts)+10*(rft-fac*erft).*0.5;
    bbbmax = localzero*ones(1,zpts)+10*(rft+fac*erft).*0.5;

    plot(z_out,bbb,'Color','black','LineWidth',0.01);
    hold on
    fill_trace(z_out,bbb,localzero,0,'b');
    fill_trace(z_out,bbb,localzero,1,'r');
    patch([z_out flip(z_out)],[bbbmin flip(bbbmax)],[0.5 0.5 0.5],'EdgeColor','none')
    plot(z_out,bbb,'Color','black','LineWidth',0.5);
    hold on
end
plot([0.,0.],[0.,45.0],'k')
zz = 10:10:160;
for i = 1:length(zz)
    plot([zz(i),zz(i)],[0.,45.],'k--','LineWidth',0.5)
end
title(['Unmodelled RFs'])
xlabel('Depth (km)')
xlim([-20,100])
xticks(-20:20:120)
ylim([-5.0,50.0])
yticks([0,10,20,30,40])
yticklabels([])
% keyboard

fname=[sel_data.save_dir '/' sel_data.nc '.' sel_data.station '-' '-Harmonics2.png'];
print(fig,fname,'-dpng','-r600');
disp("Azimuthal moveout printed")
close(fig);
end


end