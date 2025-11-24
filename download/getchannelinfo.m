function getchannelinfo(sel_data)

% Copyright 2024 F.Link, J.Wolf and M.Reiss 

% get information on station
    find_ind = strfind(sel_data.dc_url, '/dataselect/1/');

% check whether location code was used
    request_url = strcat('curl -L --digest "',char(sel_data.dc_url(1:find_ind)),...
        'station/1/query?network=',char(sel_data.nc), '&station=',...
        char(sel_data.station), ...
        '&channel=');
    request_url = strcat(request_url,'HH*,BH*');
    request_url = strcat(char(request_url),'&format=text&level=channel" -o ');

    output = strcat(sel_data.save_dir,'/tmp.txt');
    req_line =  sprintf('%s %s', char(request_url),char(output));
    system(req_line);
end