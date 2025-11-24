function [clean_mean_miss_baz,clean_std] = baz_stat(data)

% calculates mean misorientation value from previously selected events

% input
% data: struct with events

events = fieldnames(data);
indx = 1:length(events);
len = cellfun(@length,events);
its = indx(len==45);
no_up = 0;

% read only information of phases fulfilling previous quality criteria

for iEvents = 1: length(its)
    fn = events{its(iEvents)};
    for i_phase = 1:length(data.(fn).phases)
        if data.(fn).phases(i_phase).snr>2.5
            no_up = no_up+1;
            cordeg(no_up) = data.(fn).phases(i_phase).cordeg;
            baz(no_up) = data.(fn).baz;
            time_vec(no_up) = data.(fn).origin_time;
        end
    end
end

if ~exist('cordeg','var')
    cordeg = 0;
    clean_mean_miss_baz = 0;
    clean_std = 0;
end

% calculate mean miss aligment
if length(cordeg)>1
    
    % calc mean
    mean_miss_baz = mean(cordeg);
    % calc standard deviation
    std_miss_baz = std(cordeg);
    % find values within one standard deviation
    abs_val = abs(cordeg-mean_miss_baz);
    miss_baz_find=find(abs_val<std_miss_baz);
    clean_miss_baz=cordeg(miss_baz_find);
    %calc new standard deviation without outliers
    clean_mean_miss_baz = mean(clean_miss_baz);
    clean_std = std(clean_miss_baz);
       
elseif isscalar(cordeg)
    
    clean_mean_miss_baz = cordeg;
    clean_std = 0;

end


end