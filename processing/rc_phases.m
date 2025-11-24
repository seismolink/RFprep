function phases = rc_phases(phases,trace,p1,p2)

% first quality control / initial preparation of data
% usage: event: struct with event parameters, channels, phases etc.

% no phase selected yet
i_xks = 0;

% check all phases;
for i_phase = 1:length(phases)
    
    % select ..KS  phases only
    p_str = char(phases(i_phase).name);
    find_xks = strfind(p_str, 'KS');
    if isempty(find_xks)
        find_xks = strfind(p_str, 'PcS');
    end
    
    if find_xks >= 1
        i_xks = i_xks+1;
        
        % get arrival time of phase
        t_xks = phases(i_phase).tt_abs;
        if (i_xks == 1)
            t_xks_previous = t_xks;
        end
        
        % check if window contains phase & ...
        
        %cutting an appropriate time window is possible
        cut = datenum((t_xks-seconds(50)));
        [~, index] = min(abs(trace.time-cut));
        
        index_end = index+(100/trace.sr);
        if index_end > length(trace.time)
            continue
        end
        find_cut = trace.time(index:index_end);
        
        if length(find_cut)>99/trace.sr
            
            % check for overlapping phases
            i_overlap=0;
            if ((i_xks > 1) && (seconds(t_xks-t_xks_previous) < 10.0))
                i_overlap=1;
                % for next check
                t_xks_previous = t_xks;
            end
            
            % avoid analysis of overlapping phases
            if (i_overlap == 0)
                
                % cut selected time windows 100 s length
                % time_cut = trace.time(index:index_end);

                north_filt = buttern_filter(trace.n_amp,2,...
                    1/p2,1/p1,trace.sr);
                east_filt = buttern_filter(trace.e_amp,2,...
                    1/p2,1/p1,trace.sr);
                
                north_cut = north_filt(index:index_end);
                east_cut = east_filt(index:index_end);
                
                % calculate signal to noise ratio of XKS phase on 
                % effective horizontal component
                h_trace = sqrt((north_cut).^2 + (east_cut).^2);
                snr_h = snr(h_trace,trace.sr,20,25);
                
                phases(i_phase).snr = snr_h;
                
                % % cut analysis window
                % cut2 = datenum(t_xks - seconds(5));
                % [~, index2] = min(abs(time_cut-cut2));
                % 
                % index_end2 = index2+round(30/trace.sr);
                % east_cut2 = east_cut(index2:index_end2);
                % north_cut2 = north_cut(index2:index_end2);
                % 
                % % calculate long to short axis
                % [~,xlam1,xlam2] = covar(north_cut2,east_cut2);
                % 
                % % particle motion at long periods
                % 
                % north_cut_pm = buttern_filter(north_cut,2,...
                %     1/50,1/15,trace.sr);
                % 
                % east_cut_pm = buttern_filter(east_cut,2,...
                %     1/50,1/15,trace.sr);
                % 
                % east_cut2_filt = east_cut_pm(index2:index_end2);
                % north_cut2_filt = north_cut_pm(index2:index_end2);
                % 
                % maxn = max(abs(north_cut2_filt));
                % maxe = max(abs(east_cut2_filt));
                % maxne = max(maxe,maxn);
                % north_cut2_filt = north_cut2_filt/maxne;
                % east_cut2_filt = east_cut2_filt/maxne;
                % [baz_long,~,~] = covar(north_cut2_filt,east_cut2_filt);
                % 
                % % calculate miss alignment
                % if baz > baz_long
                %     baz_miss = mod(baz-baz_long,360);
                % elseif baz < baz_long
                %     baz_miss = mod(baz-baz_long,-360);
                % elseif baz == baz_long
                %     baz_miss = 0;
                % end
                % 
                % if ~exist('baz_miss')
                %     return
                % end
                % % change baz_miss to correct value
                % if baz_miss > 90
                %     baz_miss = baz_miss-180;
                % end
                % 
                % if baz_miss < -90
                %     baz_miss = baz_miss +180;
                % end
                % 
                % % save misalignment
                % phases(i_phase).cordeg = baz_miss;
                % 
                % % short/long axis ratio
                % phases(i_phase).xlam_ratio = xlam2/xlam1;
                % 
            end
            
            % take next '..KS' phase
            
        else
            % if traces are too short
            disp('cut traces are less than 100s')
            % event_faulty_comp = event_faulty_comp+1;
        end
    end
    
    % take next phase
    
end


end