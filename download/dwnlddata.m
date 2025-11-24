function [foundflag] = dwnlddata(filename,sel_data,C2,timestart,timeend)

% Copyright 2024 F.Link, J.Wolf and M.Reiss 

time1 = datevec(timestart);
time2 = datevec(timeend);

% built string for fdsn request

%start time
st_time = strcat(num2str(time1(1)),'-',...
    sprintf('%02d',time1(2)),'-',...
    sprintf('%02d',time1(3)),'T',...
    sprintf('%02d',time1(4)),':',...
    sprintf('%02d',time1(5)),':',...
    sprintf('%06.3f',time1(6)),'Z');

%end time
e_time = strcat(num2str(time2(1)),'-',...
    sprintf('%02d',time2(2)),'-',...
    sprintf('%02d',time2(3)),'T',...
    sprintf('%02d',time2(4)),':',...
    sprintf('%02d',time2(5)),':',...
    sprintf('%06.3f',time2(6)),'Z');

name_mseed = strcat(sel_data.save_dir,'/',filename);
            
if exist(char(name_mseed),'file')
    disp(['mseed ' char(filename) ' already exists. Skipping request.'])
    foundflag = 1;
    return
end
            
foundflag = 0;
streams = {'HH','BH'};
opts = 2;
ndc = length(sel_data.dci);
nt = 1;
nt2 = 1;
while ~foundflag && nt <= opts && nt2 <= ndc
%check whether location code was chosen
    disp(char(strcat({'Request data for '},char(filename),{' and channel '},char(streams{nt}))))

    dc_url = C2{1,2}{sel_data.dci(nt2)};

    if ~isempty(sel_data.auth)
        req_line1 = ['curl -L --digest --user ' sel_data.auth];
        req_line1 = strcat(req_line1,{' "'},dc_url,...
            'queryauth?net=',...
            sel_data.nc, '&sta=',sel_data.station, ...
            '&cha=',streams{nt}, ...
            '*&starttime=',st_time,...
            '&endtime=',e_time,'" -o ');
    else
        req_line1 = 'curl -L --digest';
        req_line1 = strcat(req_line1,{' "'},dc_url,...
            'query?net=',...
            sel_data.nc, '&sta=',sel_data.station, ...
            '&cha=',streams{nt}, ...
            '*&starttime=',st_time,...
            '&endtime=',e_time,'" -o ');
    end
    

    req_line = ...
        sprintf('%s %s',char(req_line1),char(name_mseed));

    % sent fdsn request
    [~, cmdout] = system(req_line);
    disp(cmdout);

    %check if data is downloaded
    MyFileInfo = dir(char(name_mseed));

    if isempty(MyFileInfo) || (isstruct(MyFileInfo) && MyFileInfo.bytes < 2000)
        if nt2 < ndc
            nt2 = nt2+1;
        else
            nt = nt+1;
        end
    else
        disp('Downloaded data for requested file.')
        foundflag = 1;
    end

end
if ~foundflag
    disp('request returned no data')
    delete(char(name_mseed));
end

end