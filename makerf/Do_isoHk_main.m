function Do_isoHk_main(sel_data,nc,station)

sel_data.nc = nc;
sel_data.station = station;
% h=sel_data.h_min:sel_data.d_h:sel_data.h_max;
% rat=sel_data.k_min:sel_data.d_k:sel_data.k_max;
h=20:0.1:60;
rat = 1.5:0.01:2.1;
% Copyright 2024 F.Link, J.Wolf and M.Reiss

sel_data.vp = 6.5;
sel_data.weights = [0.2 1 1];

vp = sel_data.vp;

% load preprocessed data
rf = load([sel_data.load_dir,'/',sel_data.nc,'.',sel_data.station,'_rf.mat']);

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
SFps = result.FSF;
SFppps = result.FSF;
SFppss = result.FSF;
Wdps = result.FSF;
Wdppps = result.FSF;
Wdppss = result.FSF;
for i = 1:length(indxoi)
    if rf.N(indxoi(i)) < 2
        continue
    end
    ind = indxoi(i);

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

    Wnps = SFps.^2;
    Wnppps = SFppps.^2;
    Wnppss = SFppss.^2;

    Wdps = Wdps+((weight(indxoi(i)).*ISFps).^2);
    Wdppps = Wdppps+((weight(indxoi(i)).*ISFppps).^2);
    Wdppss = Wdppss+((weight(indxoi(i)).*ISFppss).^2);

    result.FSF = Wnps./Wdps.*sel_data.weights(1).*SFps+Wnppps./Wdppps.*sel_data.weights(2).*SFppps+Wnppss./Wdppss.*sel_data.weights(3).*SFppss;
    [~,position]=max(result.FSF(:));
    [II,JJ] = ind2sub([length(rat) length(h)],position);
    mrat=rat(II);mh=h(JJ);
    fprintf('H = %.1f km    k = %.2f    for step %i/%i\n',mh,mrat,i,length(indxoi));
end


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
maxt = 25;
ntot = 0;
fig = figure;
subplot(1,2,1)
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
        localzero = max(pray_cmp)+ds*1.5;
        bbb = localzero.*ones(size(rft))+ds/1.5*rft./fnormsum;
        plot(tt,bbb,'Color','black','LineWidth',0.5);
        hold on
        fill_trace(tt',bbb',localzero,0,'b');
        fill_trace(tt',bbb',localzero,1,'r');
    end
end
plot(tps,pray_cmp,'k','LineWidth',1.2)
plot(tppps,pray_cmp,'k','LineWidth',1.2)
plot(tppss,pray_cmp,'k','LineWidth',1.2)
title('Radial RF')
xlabel("Delay Time (sec)")
ylabel("Slowness (sec/km)")
xlim([mint,maxt])
ylim([min(pray_cmp)-ds,max(pray_cmp)+ds*1.5*2])

subplot(1,2,2)
imagesc(h,rat,result.FSF);
hold on
plot(mh,mrat,'r.','MarkerSize',20)
plot([mh mh],[rat(1) rat(end)],'r-.','LineWidth',0.5)
plot([h(1) h(end)],[mrat mrat],'r-.','LineWidth',0.5)
set(gca,'YDir','normal'); 
colorbar; 
caxis([0 max(result.FSF(:))])
axis([h(1) h(end) rat(1) rat(end)])
xlabel('Thickness [km]');
ylabel('vp/vs ratio')

keyboard


% Searching for the maximum
% % -------------------------
% disp('Searching for maximum');
% [~,position]=max(result.FSF(:));
% [II,JJ] = ind2sub(size(result.FSF),position);
% result.mrat=rat(II);result.mh=h(JJ);
% disp(sprintf('H = %.1f km    k = %.2f  ',result.mh,result.mrat));
% result.baz = rf.baz(indxoi);
% result.pray = rf.slow(indxoi);
% result.filename = filename(indxoi);
% % save(strcat(sel_data.res_dir,'/results_',sel_data.modmethod,'_',num2str(round(sel_data.cutoff.*100)),'.mat'),'sel_data','result');
% 
% 
% 
% % Create graphical representation
% % -------------------------------
% disp('Create graphics for Stacking function.')
% graphics_name = strcat(sel_data.res_dir,'/HkcombRF_',sel_data.modmethod,'_',num2str(round(sel_data.cutoff.*100)));
% create_graphics_iso(result,sel_data,graphics_name);
% 
% close all