function [newtrace]=re_sample(event,trace,new_sr)

% this function resamples event traces
% the interp1 function is used to have the data points starting exactly at
% origin time

% resample
dt=datenum(seconds(new_sr));

% define new time vector
start_t = max([trace(1).time(1) trace(2).time(1) trace(3).time(1)]); 
end_t = min([trace(1).time(end) trace(2).time(end) trace(3).time(end)]);
new_dt =start_t:dt:end_t;

for n=1:length(trace)
    trace(n).new_amp = interp1(trace(n).time,trace(n).amp,...
        new_dt);
end

newtrace.sr = new_sr;
newtrace.z_amp = trace(3).new_amp;
newtrace.n_amp = trace(2).new_amp;
newtrace.e_amp = trace(1).new_amp;
newtrace.time = new_dt;

end