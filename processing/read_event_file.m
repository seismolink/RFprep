function [event,phflag] = read_event_file(sel_data,evID)

% extract event data from textfile
fid = fopen(strcat(sel_data.work_dir,'/event_files/', evID, '.txt'),'rt');
C2 = textscan(fid,['%d %d %d %d %d %f %f %f ',...
                    '%f %f %s %s %f %f %f %f %f ',...
                    '%f %d %d %d %d %d %d %d']);
fclose(fid);
net = C2{1,11};
station = C2{1,12};
comb = strcat(net,'.',station);
ind = strcmp(comb,strcat(sel_data.nc,'.',sel_data.station));
nn = 1:length(C2{1,1});
n = nn(ind);

% read origin time
date_val = [C2{1,1}(n),C2{1,2}(n),C2{1,3}(n),C2{1,4}(n),C2{1,5}(n),C2{1,6}(n)];
event.origin_time = datetime(date_val);
event.ot_num=datenum(event.origin_time);

event.lat = C2{1,7}(n);
event.lon = C2{1,8}(n);
event.depth = C2{1,9}(n);
event.mag = C2{1,10}(n);

event.dist = C2{1,16}(n);
event.baz = C2{1,17}(n);

end