function [snrflag,trL_f,SNR_L] = snrcheck(phase,thresL,fmin,fmax,comp)

trL_f = buttern_filter(phase.tr(strfind(phase.comp,comp)).data,4,fmin,fmax,phase.dt);
ref = phase.tt_abs;

[~, Nstart] = min(abs(phase.time-((ref-seconds(100)))));
[~, Nend] = min(abs(phase.time-((ref-seconds(20)))));
[~, Pstart] = min(abs(phase.time-((ref-seconds(5)))));
[~, Pend] = min(abs(phase.time-((ref+seconds(10)))));

% [~, Nstart] = min(abs(phase.time-((ref-seconds(100)))));
% [~, Nend] = min(abs(phase.time-((ref-seconds(20)))));
% [~, Pstart] = min(abs(phase.time-((ref-seconds(5)))));
% [~, Pend] = min(abs(phase.time-((ref+seconds(100)))));

% [~, Nstart] = min(abs(phase.time-((ref-seconds(20)))));
% [~, Nend] = min(abs(phase.time-((ref-seconds(10)))));
% [~, Pstart] = min(abs(phase.time-((ref-seconds(5)))));
% [~, Pend] = min(abs(phase.time-((ref+seconds(10)))));

trL_noise = trL_f(Nstart:Nend);
trL_signal = trL_f(Pstart:Pend);

nL = mean(abs(trL_noise));
pL = mean(abs(trL_signal));
SNR_L = 10*log10((pL/nL).^2);

if round(SNR_L)>=thresL
    snrflag = 1;
else
    snrflag = 0;
end

end
