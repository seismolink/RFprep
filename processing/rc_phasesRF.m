function phase = rc_phasesRF(phases,trace,event,sel_data)

% first quality control / initial preparation of data
% usage: event: struct with event parameters, channels, phases etc.

% check all phases;
ii = 0;
for i_phase = 1:length(phases)
    p_str = char(phases(i_phase).name);
    if event.dist < 115
        % select ..KS  phases only
        find_P = strfind(p_str(end), 'P');
        if isempty(find_P)
            find_P = strfind(p_str, 'Pdiff');
        end
    else
        find_P = strfind(p_str, 'PKIKP');
    end
    if find_P >= 1
        ii = i_phase;
        break
    end
end
if ii == 0
    phase = 0;
    disp('phase not found')
    return
end
phase = phases(ii);
tp = phase.tt_abs;

%cutting an appropriate time window if possible
cut = (tp-seconds(150));
[~, index] = min(abs(trace.time-cut));
if seconds(cut-trace.time(1)) < 0
    phase = 0;
    disp('time too short')
    return
end

index_end = index+(450/trace.sr);
if index_end > length(trace.time)
    phase = 0;
    disp('time too short')
    return
end

find_cut = trace.time(index:index_end);
bhw = tukeywin(length(find_cut),0.005);
% if length(find_cut)>239/trace.sr
    phase.tr(1).data = detrend(trace.z_amp(index:index_end)).*bhw';
    phase.comp = 'zne';
    phase.tr(2).data = detrend(trace.n_amp(index:index_end)).*bhw';
    phase.tr(3).data = detrend(trace.e_amp(index:index_end)).*bhw';
    phase.time = find_cut;
    phase.dt = trace.sr;
    phase.event = event;
    [s1,L1,snr1] = snrcheck(phase,sel_data.snrlim,0.02,5,'z');
    [s2,L2,snr2] = snrcheck(phase,sel_data.snrlim,0.02,0.15,'z');
    [s3,L3,snr3] = snrcheck(phase,sel_data.snrlim,0.5,5,'z');
    % if s1+s3+s3 > 1
        [s4,snr4a,snr4b,rf,tt,trf,~,test] = snrrfcheck(phase,1,2);
        % [s4b,rf2,tt2,trf2] = snrrfcheck2(phase,1,3);
    % else
    %     s4 = 0;
    % end
    % if event.mag >= 5.8
    % fig = figure;
    % subplot(4,1,1)
    % plot(L1)
    % text(100,0,num2str(s1),'FontSize',16)
    % subplot(4,1,2)
    % plot(L2)
    % text(100,0,num2str(s2),'FontSize',16)
    % subplot(4,1,3)
    % plot(L3)
    % text(100,0,num2str(s3),'FontSize',16)
    % subplot(4,1,4)
    % plot(tt(1:length(rf)),rf)
    % hold on
    % plot(tt(1:length(trf)),trf)
    % axis([-20 20 -max(abs(rf)) max(abs(rf))])
    % text(-19,0,num2str(s4),'FontSize',16)
    % keyboard
    % close(fig);
    % end
    
    % if s1+s2+s3 == 3
    %     keyboard
    % end

    if ~(s1+s2+s3>1 && s4>0)
        phase = 0;
        disp('low snr')
    else
        phase.snr = s1+s2+s3+s4;
        % disp(num2str(snr4a))
    end
% else
%     % if traces are too short
%     % disp('cut traces are less than 240s')
%     % event_faulty_comp = event_faulty_comp+1;
% end
    



end