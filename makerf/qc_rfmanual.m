function dataout = qc_rfmanual(datain,sel_data)


dataout.stat = datain.stat;
dataout.nc = datain.nc;
dataout.lat = datain.lat;
dataout.lon = datain.lon;
dataout.ele = datain.ele;

nomevents = 0;

fig = figure('Units','normalized','Position',[0.25 0.1 0.5 0.75]);
fn = fieldnames(datain);
fn(~contains(fn,'ev')) = [];
for i = 1:length(fn)
    phase = datain.(fn{i});
    if sel_data.srfflag
        [s1,L1,snr1] = snrcheck(phase,sel_data.snrlim,0.02,1,'r');
        [s2,L2,snr2] = snrcheck(phase,sel_data.snrlim,0.02,0.1,'r');
        [s3,L3,snr3] = snrcheck(phase,sel_data.snrlim,0.075,0.5,'r');
        [~,snr4a,snr4b,rf,tt,test] = snrsrfcheck(phase,0.5,2);
    else
        [~,L1,snr1] = snrcheck(phase,5,0.02,5,'z');
        [~,L2,snr2] = snrcheck(phase,5,0.02,0.15,'z');
        [~,L3,snr3] = snrcheck(phase,5,0.5,5,'z');
        [~,snr4a,snr4b,rf,tt,trf,~,test] = snrrfcheck(phase,1,2);
    end
    % fig = figure('Units','normalized','Position',[0.25 0.1 0.5 0.75]);
    subplot(5,1,1)
    pl1 = plot(seconds(phase.time-phase.time(1)),L1);
    hold on;
    pl2 = plot([seconds(phase.tt_abs-phase.time(1)) seconds(phase.tt_abs-phase.time(1))],[min(L1) max(L1)],'r');
    % text(10,0,num2str(snr1),'FontSize',16)
    legend({'0.02-5Hz'},'Location','northwest')
    title([dataout.nc '.' dataout.stat ' ' fn{i}])
    hold off
    subplot(5,1,2)
    pl3 = plot(seconds(phase.time-phase.time(1)),L2);
    hold on;
    pl4 = plot([seconds(phase.tt_abs-phase.time(1)) seconds(phase.tt_abs-phase.time(1))],[min(L2) max(L2)],'r');
    % text(10,0,num2str(snr2),'FontSize',16)
    legend({'0.02-0.15Hz'},'Location','northwest')
    hold off
    subplot(5,1,3)
    pl5 = plot(seconds(phase.time-phase.time(1)),L3);
    hold on;
    pl6 = plot([seconds(phase.tt_abs-phase.time(1)) seconds(phase.tt_abs-phase.time(1))],[min(L3) max(L3)],'r');
    % text(10,0,num2str(snr3),'FontSize',16)
    legend({'0.5-5Hz'},'Location','northwest')
    hold off
    subplot(5,1,4)
    pl7 = plot(tt(1:length(rf)),rf,'LineWidth',1.2);
    hold on
    if sel_data.prfflag
        pl8 = plot(tt(1:length(trf)),trf,'LineWidth',1.2);
    else
        pl8 = 0;
    end
    axis([-20 20 -max(abs(rf)) max(abs(rf))])
    % text(-19,0,[num2str(snr4a) ' ' num2str(snr4b)],'FontSize',16)
    legend({'rad RF','tra RF'},'Location','northwest')
    hold off
    subplot(5,1,5)
    pl9 = plot(tt(1:length(rf)),test(1:length(rf)),'LineWidth',1.2);
    hold on
    pl10 = plot([-20 20],[2 2],'r-.','LineWidth',0.1);
    axis([-20 20 0 1.05*max(abs(test))])
    % text(-19,1,num2str(max(test)),'FontSize',16)
    legend({'STA/LTA','limit'})
    hold off
    figure(fig);
    sm = input('keep trace for analysis? [y]es or [enter]','s');
    % close(fig);
    if ~strcmp(sm,'y')
        disp(['trace discarded - ' num2str(nomevents) '/' num2str(i) ' of ' num2str(length(fn))])
    else
        dataout.(fn{i}) = phase;
        nomevents = nomevents+1;
        disp([num2str(nomevents) '/' num2str(i) ' of ' num2str(length(fn))])
    end
    delete([pl1 pl2 pl3 pl4 pl5 pl6 pl7 pl8 pl9 pl10])
end
close(fig)

end