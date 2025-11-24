function phasesn = rc_phasesSRF(phases,trace,event,sel_data)

% first quality control / initial preparation of data
% usage: event: struct with event parameters, channels, phases etc.

mis = sel_data.mis;

% check all phases;
ii = [];
for i_phase = 1:length(phases)
    s_str = char(phases(i_phase).name);
    if event.dist > 55 && event.dist < 85
        find_S = strcmp(s_str,'S');
        if ~find_S
            find_S = strcmp(s_str,'Sdiff');
        end
        if find_S
            ii = [ii,i_phase];
        end
    end
    if event.dist > 50 && event.dist < 80
        find_S = strcmp(s_str,'ScS');
        if find_S
            ii = [ii,i_phase];
        end
    end
    if event.dist > 85 && event.dist < 120
        find_S = strcmp(s_str,'SKS');
        if find_S
            ii = [ii,i_phase];
        end
    end
end
if isempty(ii)
    phasesn = 0;
    disp('phase not found')
    return
end
n = 1;
for jj = 1:length(ii)
    phase = phases(ii(jj));
    tp = phase.tt_abs;
    
    %cutting an appropriate time window if possible
    cut = (tp-seconds(300));
    [~, index] = min(abs(trace.time-cut));
    if seconds(cut-trace.time(1)) < 0
        phasesn = 0;
        disp('time too short')
        return
    end
    
    index_end = index+(450/trace.sr);
    if index_end > length(trace.time)
        phasesn = 0;
        disp('time too short')
        return
    end
    find_cut = trace.time(index:index_end);
    bhw = tukeywin(length(find_cut),0.005);

    if length(mis.start)>1
        phi = mis.phi(mis.start<=phase.tt_abs&mis.end>=phase.tt_abs);
    else
        phi = mis.phi;
    end
    if isempty(phi)
        phi = mis.phi(end);
    end
    baz = event.baz-phi;

    phase.tr(1).data = detrend(trace.z_amp(index:index_end)).*bhw';
    phase.tr(2).data = detrend(trace.n_amp(index:index_end)).*bhw';
    phase.tr(3).data = detrend(trace.e_amp(index:index_end)).*bhw';
    [phase.tr(1).data,phase.tr(2).data,phase.tr(3).data]=rotVNE2VRT(phase.tr(1).data,phase.tr(2).data,phase.tr(3).data,baz);
    % [phase.tr(2).data, phase.tr(3).data] =  rad_tra (phase.tr(2).data, phase.tr(3).data, baz);
    phase.comp = 'zrt';
    phase.time = find_cut;
    phase.dt = trace.sr;
    phase.event = event;
    [s1,L1,snr1] = snrcheck(phase,sel_data.snrlim,0.02,1,'r');
    [s2,L2,snr2] = snrcheck(phase,sel_data.snrlim,0.02,0.1,'r');
    [s3,L3,snr3] = snrcheck(phase,sel_data.snrlim,0.075,0.5,'r');
        % if s1+s2+s3 > 0
        %     [s4,snr4a,snr4b,rf,tt,test] = snrsrfcheck(phase,0.3,2);
        %     % [s4b,rf2,tt2,trf2] = snrrfcheck2(phase,1,3);
        % else
        %     rf = L1;
        %     tt = linspace(-20,20,length(rf));
        %     s4 = 0;
        % end

        % if s1+s2+s3 > 0
        % fig = figure;
        % subplot(4,1,1)
        % plot(phase.time,L1)
        % hold on;
        % if s1
        %     plot([tp tp],[min(L1) max(L1)],'g')
        % else
        %     plot([tp tp],[min(L1) max(L1)],'r')
        % end
        % % text(100,0,num2str(s1),'FontSize',16)
        % subplot(4,1,2)
        % plot(phase.time,L2)
        % hold on;
        % if s2
        %     plot([tp tp],[min(L1) max(L1)],'g')
        % else
        %     plot([tp tp],[min(L1) max(L1)],'r')
        % end
        % % text(100,0,num2str(s2),'FontSize',16)
        % subplot(4,1,3)
        % plot(phase.time,L3)
        % hold on;
        % if s3
        %     plot([tp tp],[min(L1) max(L1)],'g')
        % else
        %     plot([tp tp],[min(L1) max(L1)],'r')
        % end
        % % text(100,0,num2str(s3),'FontSize',16)
        % subplot(4,1,4)
        % plot(tt(1:length(rf)),rf)
        % hold on
        % % plot(tt(1:length(trf)),trf)
        % axis([-20 20 -max(abs(rf)) max(abs(rf))])
        % % text(-19,0,num2str(s4),'FontSize',16)
        % text(-19,0,num2str(max(snr4a,snr4b)),'FontSize',16)
        % pause
        % close(fig);
        % end
        
        % if s1+s2+s3 == 3
        %     keyboard
        % end
    
        if ~(s1+s2+s3>1)% && s4>0)
            phase = 0;
            disp('low snr')
        else
            % phase.snr = s1+s2+s3+s4;
            phase.snr = s1+s2+s3;
        end
    % else
    %     % if traces are too short
    %     % disp('cut traces are less than 240s')
    %     % event_faulty_comp = event_faulty_comp+1;
    % end
    if isstruct(phase)
    phasesn(n) = phase;
    n = n+1;
    end
end
if n == 1
    phasesn = 0;
end
    



end