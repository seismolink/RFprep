function Do_isoHk_main2(sel_data,nc,station)

sel_data.nc = nc;
sel_data.station = station;
% h=sel_data.h_min:sel_data.d_h:sel_data.h_max;
% rat=sel_data.k_min:sel_data.d_k:sel_data.k_max;
h=sel_data.hmin:0.1:sel_data.hmax;
rat = sel_data.kmin:0.0025:sel_data.kmax;
% Copyright 2024 F.Link, J.Wolf and M.Reiss

Nb = 100;

sel_data.vp = 6.5;
sel_data.weights = [1 1 1];

vp = sel_data.vp;

% load preprocessed data
rf = load([sel_data.load_dir,'/',sel_data.nc,'.',sel_data.station,'_rf.mat']);
statpos = load([sel_data.load_dir,'/',sel_data.nc,'.',sel_data.station,'.mat'],'lat','lon');

indxoi = 1:length(rf.N);

weight = zeros(size(indxoi));
pray_cmp=linspace(min(rf.slow),max(rf.slow),15);
for i = 2:length(pray_cmp)
    clear indexsub
    indexsub = indxoi(rf.slow(indxoi)<pray_cmp(i)&rf.slow(indxoi)>=pray_cmp(i-1));
    if numel(indexsub)>=2
        weight(indexsub) = 1./length(indexsub);
    end
end

vs = vp./rat;

result.FSF = zeros(length(rat),length(h),1,'single');
FSFb = zeros([size(result.FSF),Nb]);
for ib = 1:Nb
SFps = result.FSF;
SFppps = result.FSF;
SFppss = result.FSF;
Wdps = result.FSF;
Wdppps = result.FSF;
Wdppss = result.FSF;
idx = randi(length(indxoi),length(indxoi),1);
for i = 1:length(idx)
    if rf.N(idx(i)) < 2
        continue
    end
    ind = idx(i);

    tps = (sqrt((1./vs.^2-rf.slow(ind).^2))-sqrt((1./vp.^2-rf.slow(ind).^2)))'*(h);
    FSFps = interp1(rf.time,rf.tr(ind).rrf,tps);
    tppps = (sqrt((1./vs.^2-rf.slow(ind).^2))+sqrt((1./vp.^2-rf.slow(ind).^2)))'*(h);
    FSFppps = interp1(rf.time,rf.tr(ind).rrf,tppps);
    tppss = 2*sqrt((1./vs.^2-rf.slow(ind).^2))'*(h);
    FSFppss = interp1(rf.time,rf.tr(ind).rrf,tppss);

    ISFps = sign(FSFps).*FSFps.^2;
    ISFppps = sign(FSFppps).*FSFppps.^2;
    ISFppss = -sign(FSFppss).*FSFppss.^2;
    
    ISFps(isnan(ISFps)) = 0;
    ISFppps(isnan(ISFppps)) = 0;
    ISFppss(isnan(ISFppss)) = 0;

    SFps = SFps+(weight(indxoi(i)).*ISFps);
    SFppps = SFppps+(weight(indxoi(i)).*ISFppps);
    SFppss = SFppss+(weight(indxoi(i)).*ISFppss);

    % Wnps = SFps.^2;
    % Wnppps = SFppps.^2;
    % Wnppss = SFppss.^2;
    % 
    % Wdps = Wdps+((weight(indxoi(i)).*ISFps).^2);
    % Wdppps = Wdppps+((weight(indxoi(i)).*ISFppps).^2);
    % Wdppss = Wdppss+((weight(indxoi(i)).*ISFppss).^2);

    Wnps = abs(SFps);
    Wnppps = abs(SFppps);
    Wnppss = abs(SFppss);

    Wdps = Wdps+abs((weight(indxoi(i)).*ISFps));
    Wdppps = Wdppps+abs((weight(indxoi(i)).*ISFppps));
    Wdppss = Wdppss+abs((weight(indxoi(i)).*ISFppss));

    temp = Wnps./Wdps.*sel_data.weights(1).*SFps+Wnppps./Wdppps.*sel_data.weights(2).*SFppps+Wnppss./Wdppss.*sel_data.weights(3).*SFppss;
    [~,position]=max(temp(:));
    [II,JJ] = ind2sub([length(rat) length(h)],position);
    mrat=rat(II);mh=h(JJ);
    % fprintf('H = %.1f km    k = %.2f    for step %i/%i\n',mh,mrat,i,length(indxoi));
end
FSFb(:,:,ib) = temp;
mratb(ib) = mrat;
mhb(ib) = mh;
fprintf('Step %i of %i\n',ib,Nb)
end
result.FSF = mean(FSFb,3);
result.FSF = result.FSF./max(result.FSF(:));
mrat = mean(mratb);
mh = mean(mhb);

P = 0.9545;
hbstrp = h;
ratbstrp = rat;
h_b = mh;
rat_b = mrat;
dvlim(1) = sqrt(2/pi)*(h(2)-h(1));
dvlim(2) = sqrt(2/pi)*(rat(2)-rat(1));
SS(1).var = mhb(mhb>mean(mhb)-5&mhb<mean(mhb)+5);
SS(2).var = mratb(mratb>mean(mratb)-0.1&mratb<mean(mratb)+0.1);
for i = 1:length(SS)
    SS(i).mvar = mean(SS(i).var);
    SS(i).dmv = 2*sqrt(var(SS(i).var));
    SS(i).nomv = length(SS(i).var);
    if SS(i).dmv < dvlim(i)
        SS(i).dmv = dvlim(i);
    end
end
Sres.mh = mhb;
Sres.stat.mh = SS(1).mvar;
Sres.stat.nomh = SS(1).nomv;
Sres.stat.dmh = SS(1).dmv;
Sres.mmh = SS(1).mvar;
Sres.dmh = SS(1).dmv;
Sres.mrat = mratb;
Sres.stat.mrat = SS(2).mvar;
Sres.stat.nomrat = SS(2).nomv;
Sres.stat.dmrat = SS(2).dmv;
Sres.mmrat = SS(2).mvar;
Sres.dmrat = SS(2).dmv;

toplot.stat.mx = round(Sres.stat.mh*10)/10;
toplot.stat.dx = round(Sres.stat.dmh*10)/10;
toplot.stat.my = round(Sres.stat.mrat*100)/100;
toplot.stat.dy = round(Sres.stat.dmrat*100)/100;
toplot.stat.xcount = Sres.mh;
toplot.stat.ycount = Sres.mrat;
toplot.stat.nomx = Sres.stat.nomh;
toplot.stat.nomy = Sres.stat.nomrat;
toplot.stat.x = linspace(min(h),max(h),1000);
toplot.stat.y = linspace(min(rat),max(rat),1000);
toplot.stat.periodx = 0;
toplot.stat.periody = 0;
mm = 1;
for kk = 1:length(toplot.stat.mx)
    toplot.stat.Px(kk) = round(toplot.stat.nomx(kk)./length(toplot.stat.xcount)*P*100*100)/100;
    toplot.resultstr{mm} = ['H: ' num2str(toplot.stat.mx(kk)) ' + ' num2str(toplot.stat.dx(kk)) ' km; [' num2str(toplot.stat.Px(kk)) ' %]'];
    mm = mm+1;
end
for kk = 1:length(toplot.stat.my)
    toplot.stat.Py(kk) = round(toplot.stat.nomy(kk)./length(toplot.stat.ycount)*P*100*100)/100;
    toplot.resultstr{mm} = ['k: ' num2str(toplot.stat.my(kk)) ' + ' num2str(toplot.stat.dy(kk)) '; [' num2str(toplot.stat.Py(kk)) ' %]'];
    mm = mm+1;
end
[~,I] = max(toplot.stat.Px);
toplot.href = toplot.stat.mx(I);
[~,I] = max(toplot.stat.Py);
toplot.kref = toplot.stat.my(I);


title_fontsize2 = 11;
title_fontsizec = num2str(title_fontsize2);
axes_fontsizec = 9;

% epicentral sweep
% compute composite RFs spaced:  epi = 25,155,10
depi = 10;
epis = 25:depi:155;
ds = (pray_cmp(2)-pray_cmp(1)).*1.5;

tps = (sqrt((1./vs(II).^2-pray_cmp.^2))-sqrt((1./vp.^2-pray_cmp.^2)))'*(h(JJ));
tppps = (sqrt((1./vs(II).^2-pray_cmp.^2))+sqrt((1./vp.^2-pray_cmp.^2)))'*(h(JJ));
tppss = 2*sqrt((1./vs(II).^2-pray_cmp.^2))'*(h(JJ));

tt = rf.time;
mint = -5;
% maxt = 25;
ntot = 0;
fig = figure;
set(fig,'Position',[100 50 900 700])
% hp3_1_hk_rf=subplot(9,2,[3 5 7 9 11 13 15 17]);hold on
% set(hp3_1_hk_rf,'TickDir','Out')
object = mysubplot(fig,0.05,0,0.5,1);
delete(object.c)
axes(object.a)
for ii = 1:length(pray_cmp)
    slow = pray_cmp(ii);
    jrec = 0;
    nev = 0;
    rff = zeros(size(rf.tr(1).rrf));
    rft = rff;
    for jj = 1:length(rf.epi)
        if rf.N(jj) < 2
            continue
        end
        pray = rf.slow(jj);
        test = abs(slow-pray);%abs(baz-ebaz);
        if test<ds
            ntot = ntot+1;
            jrec = jrec+1;
            rff = rff + rf.tr(jj).rrf;
        end
        if ii == 1
            rft = rft + rf.tr(jj).rrf;
        end
    end
    if jrec > 0
        fnorm = max(abs(rff));
        bbb = slow.*ones(size(rff))+ds/1.5*rff./fnorm;
        plot(tt,bbb,'Color','black','LineWidth',0.5);
        hold on
        localzero = slow;
        fill_trace(tt',bbb',localzero,0,'b');
        fill_trace(tt',bbb',localzero,1,'r');
    end
    if ii == 1
        fnormsum = max(abs(rft));
        % localzero = max(pray_cmp)+ds*1.5;
        localzerosum = 0;
        bbbsum = localzerosum.*ones(size(rft))+ds/1.5*rft./fnormsum;
        % plot(tt,bbb,'Color','black','LineWidth',0.5);
        % hold on
        % fill_trace(tt',bbb',localzerosum,0,'b');
        % fill_trace(tt',bbb',localzerosum,1,'r');
    end
end
plot(tps,pray_cmp,'k','LineWidth',1.2)
plot(tppps,pray_cmp,'k','LineWidth',1.2)
plot(tppss,pray_cmp,'k','LineWidth',1.2)
% title('Radial RF')
xlabel("Delay Time (sec)")
ylabel("Slowness (sec/km)")
xlim([mint,ceil(max(tppss)/5)*5])
ylim([min(pray_cmp)-ds,max(pray_cmp)+ds])

% hp1_1_hk_rf=subplot(9,2,1);hold on
% set(hp1_1_hk_rf,'TickDir','Out')
axes(object.b)
plot(tt,bbbsum,'Color','black','LineWidth',0.5);
hold on
fill_trace(tt',bbbsum',localzerosum,0,'b');
fill_trace(tt',bbbsum',localzerosum,1,'r');
ylabel('\Sigma RF_i')
Rtitle=title({['\fontsize{' title_fontsizec '} Radial Q-Component']});
set(Rtitle,'FontWeight','Bold')
min_y = min(bbbsum);
if min_y > 0
    min_y = -0.0001;
end
max_y = max(bbbsum);
if min_y*1.1 > max_y*1.1
    min_y = -max([abs(min_y) abs(max_y)]);
    max_y = max([abs(min_y) abs(max_y)]);
end
xlim([mint,ceil(max(tppss)/5)*5])
ylim([min_y*1.1 max_y*1.1])
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
set(gca,'FontSize',axes_fontsizec);

% hp2_1_hk_rf=subplot(9,2,[2 4 6 8 10 12 14 16 18]);
% pos = get(hp2_1_hk_rf,'Position');
% object = mysubplot(fig,pos(1),pos(2),pos(3),pos(4));
object = mysubplot(fig,0.45,0,0.5,1);
axes(object.a)

% object = mysubplot(fig,xpos,1-jj*hh,1/n_tot,hh);
axes(object.b);
title('H-k Stack','FontSize',title_fontsize2)
axes(object.a)
axpos = get(object.a,'Position');
axpos2 = get(object.c,'Position');
clr_bar = colorbar('EastOutSide','FontSize',axes_fontsizec);
clr_bar.Position= ([axpos2(1)+axpos2(3)+0.01 axpos2(2) 0.025 axpos2(4)]);
clr_bar.AxisLocation = 'out';
zlab = get(clr_bar,'ylabel');
set(zlab,'String','Normalized Stack','FontSize',axes_fontsizec); 
object.a.Position = axpos;

axes(object.a);
[X,Y] = meshgrid(h,rat);
contourf(X,Y,result.FSF,50,'Linestyle','none');
CC=load('pal_hk');
CC2=interp_pal(CC,10);
colormap(CC2);
clim([0 max(result.FSF(:))]);
% plot([min(h) max(h)],[mrat mrat],'k-.','LineWidth',1.0);
% plot([mh mh],[min(rat) max(rat)],'k-.','LineWidth',1.0);
% plot(mh,mrat,'k*','MarkerSize',30,'LineWidth',1.0)
ylabel('vp/vs','FontSize',axes_fontsizec)
xlabel('Thickness in [km]','FontSize',axes_fontsizec)
xt = get(object.a,'XTick');
yt = get(object.a,'YTick');
for i = 2:length(xt)
    plot([xt(i) xt(i)],[min(rat) max(rat)],'Color',[1,1,1]*0.25,'Linestyle',':');
end
for j = 2:length(yt)
    plot([min(h) max(h)],[yt(j) yt(j)], 'Color',[1,1,1]*0.25,'Linestyle',':');
end
set(object.a,'FontSize',axes_fontsizec)

axes(object.b);
xcount = toplot.stat.xcount;
dx = (toplot.stat.x(2)-toplot.stat.x(1))*10;
mx = toplot.stat.mx;
ex = toplot.stat.dx;
x = toplot.stat.x;
nomx = toplot.stat.nomx;
edges = (min(xcount)-dx/2):dx:(max(xcount)+dx/2);
histogram(xcount,edges)
ymax = max(hist(xcount,ceil((max(xcount)-min(xcount))/dx)+1));
hold all
for i = 1:length(mx)
    sigma = ex(i)/2;
    bt = x;
    aa = exp(-(bt-mx(i)).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*nomx(i).*dx;
    plot(bt,aa,'LineWidth',1.2)
    ymax = max([ymax max(aa)]);
%     text(S.stat.mh(i)+sel_data.d_h,max(aa)*1.05,[num2str(round(S.stat.mh(i)*100)/100) ' ' char(177) ' ' num2str(round(S.stat.dmh(i)*100)/100) ' km ' char(10) 'P = ' num2str(round(S.stat.nomh(i)./length(S.mh)*P*100*100)/100) ' %'],'VerticalAlignment','bottom')
end
for i = 1:length(mx)
    plot([mx(i) mx(i)],[0 1.15*ymax],'r-.','LineWidth',1.2)
end
axis([min(x) max(x) 0 1.15*ymax])
set(object.b,'XTickLabel',[],'YAxisLocation','right');%,'XTickLabel',[]);
set(object.b,'YGrid','on','XGrid','on');
clear yticks
yticks_temp = get(object.b,'YTick');
for i = 1:(length(yticks_temp)-1)
    yticks{i} = '';
end
yticks{end+1} = num2str(yticks_temp(end));
set(object.b,'YTickLabel',yticks,'FontSize',axes_fontsizec)
set(object.b,'FontSize',axes_fontsizec)
%         ylabel('Nb','FontSize',AxesFontSize)
box on
%         box off
axes(object.a);
for i = 1:length(mx)
    plot([mx(i) mx(i)],[min(toplot.stat.y) max(toplot.stat.y)],'r-.','LineWidth',1.2)
end

axes(object.c);
xcount = toplot.stat.ycount;
dx = (toplot.stat.y(2)-toplot.stat.y(1))*10;
mx = toplot.stat.my;
ex = toplot.stat.dy;
x = toplot.stat.y;
nomx = toplot.stat.nomy;
edges = (min(xcount)-dx/2):dx:(max(xcount)+dx/2);
histogram(xcount,edges,'Orientation','horizontal');
ymax = max(hist(xcount,ceil((max(xcount)-min(xcount))/dx)+1));
hold all
for i = 1:length(mx)
    sigma = ex(i)/2;
    bt = x;
    aa = exp(-(bt-mx(i)).^2./(2*sigma*sigma))./(sqrt(2*pi)*sigma).*nomx(i).*dx;
    plot(aa,bt,'LineWidth',1.2)
    ymax = max([ymax max(aa)]);
%     text(S.stat.mrat(i)+sel_data.d_k,max(aa)*1.05,[num2str(round(S.stat.mrat(i)*100)/100) ' ' char(177) ' ' num2str(round(S.stat.dmrat(i)*100)/100) ' ' char(10) 'P = ' num2str(round(S.stat.nomrat(i)./length(S.mrat)*P*100*100)/100) ' %'],'VerticalAlignment','bottom')
end
for i = 1:length(mx)
    plot([0 1.15*ymax],[mx(i) mx(i)],'r-.','LineWidth',1.2)
end
axis([0 1.15*ymax min(x) max(x)])
%         set(object.c,'YTickLabel',[],'XTickLabel',[]);
set(object.c,'YTickLabel',[],'XAxisLocation','top');
xticks = get(object.c,'XTickLabel');
for i = 1:(length(xticks)-1)
xticks{i} = [];
end
set(object.c,'XTickLabel',xticks,'FontSize',axes_fontsizec)
set(object.c,'YGrid','on','XGrid','on');
box on
%         box off
%         xlabel('N','FontSize',AxesFontSize)
axes(object.a);
for i = 1:length(mx)
    plot([min(toplot.stat.x) max(toplot.stat.x)],[mx(i) mx(i)],'r-.','LineWidth',1.2);
end
for i = 1:length(toplot.stat.my)
    for j = 1:length(toplot.stat.mx)
        plot(toplot.stat.mx(j),toplot.stat.my(i),'r*','MarkerSize',20,'LineWidth',1.2)
    end
end

fname=[sel_data.save_dir '/' sel_data.nc '.' sel_data.station '-' '-isoHk.png'];
print(fig,fname,'-dpng','-r600');

fid = fopen([sel_data.save_dir '/' sel_data.Hkfile],'at');
fprintf(fid,'%s,%f,%f,%f,%f,%f,%f\n',[sel_data.nc '.' sel_data.station],statpos.lat,statpos.lon,Sres.mmh,Sres.dmh,Sres.mmrat,Sres.dmrat);
fclose(fid);

close(fig);