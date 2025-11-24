function bdsr_download(path2opts,sel_data)

% control function to gather nescessary parameters and invoke the download
% of new data
%
if isfolder(path2opts)
    sel_data.work_dir = path2opts;
    sel_data.optsfile = [path2opts '/SR4BigData_options.txt'];
else
    sel_data.work_dir = fileparts(path2opts);
    sel_data.optsfile = path2opts;
end

cur_dir = pwd;
sel_data.defaults_path = strcat(cur_dir,'/defaults/');

% read relevant information of options file
sel_data = readoptsfile(sel_data);

% make download folder
if ~exist([sel_data.work_dir '/mseed_files/'],'dir')
    mkdir([sel_data.work_dir '/mseed_files/']);
end

dci = [];
if ~sel_data.minlat && ~sel_data.maxlat && ~sel_data.minlon && ~sel_data.maxlon
    disp('no limits read from file')
    dwnldlist = dir([sel_data.work_dir '/mseed_files/' sel_data.statlist]);
    disp(['Read stations from ' sel_data.statlist])
% read info 
    n = 1;
    m = inf;
    fileID = fopen([sel_data.work_dir '/mseed_files/' dwnldlist.name]);
    line = fgetl(fileID);
    while ischar(line)
        if length(line)<10
            dcitest = str2double(line);
            if ~isnan(dcitest)
                m = 0;
            end
        end
        if length(line)>10 && ~contains(line,'DATACENTER') && ~contains(line,'#NETWORK','IgnoreCase',true)
            C = textscan(line,'%s %s %f %f %f %s %s %s','Delimiter','|');
            try
                starttime = datetime(C{1,7}{1}(1:19),'InputFormat','yyyy-MM-dd''T''HH:mm:ss');
            catch ME
                starttime = datetime(C{1,7}{1});
            end
            if isempty(C{1,8})
                endtime = datetime('now');
            else
                try
                    endtime = datetime(C{1,8}{1}(1:19),'InputFormat','yyyy-MM-dd''T''HH:mm:ss');
                catch ME
                    endtime = datetime(C{1,8}{1});
                end
            end
            disp(num2str(n))
            if days(endtime-starttime) < 60
                line = fgetl(fileID);
                continue
            end
            nc{n} = C{1,1}{1};
            station{n} = C{1,2}{1};
            if m < 5
                dci(n) = dcitest;
            end
            svec{n} = datevec(starttime);
            n = n+1;
        else
            m = m+1;
        end
        line = fgetl(fileID);
    end
    fclose(fileID);
else
% Check for existing stations in lat/lon box
    fileID2 = fopen([sel_data.defaults_path,'/data_centers.dat']);
    C2 = textscan(fileID2,'%s %s');
    fclose(fileID2);
    station = {};
    nc = {};
    svec = {};
    n = 1;
    for ii = 1:length(C2{1,1})
        dc_url = C2{1,2}{ii};
        find_ind = strfind(dc_url, '/dataselect/1/');
        if ~isempty(sel_data.auth)
            request_url = strcat('curl -L -s --digest "',char(dc_url(1:find_ind)),...
                'station/1/query?minlatitude=',num2str(sel_data.minlat), '&maxlatitude=',...
                num2str(sel_data.maxlat),'&minlongitude=',num2str(sel_data.minlon), ...
                '&maxlongitude=',num2str(sel_data.maxlon),'&includerestricted=true&channel=');
        else
            request_url = strcat('curl -L -s --digest "',char(dc_url(1:find_ind)),...
                'station/1/query?minlatitude=',num2str(sel_data.minlat), '&maxlatitude=',...
                num2str(sel_data.maxlat),'&minlongitude=',num2str(sel_data.minlon), ...
                '&maxlongitude=',num2str(sel_data.maxlon),'&includerestricted=false&channel=');
        end
        request_url = strcat(request_url,'HH*,BH*');
        request_url = strcat(char(request_url),'&format=text&level=station"');
        [a,b] = system(request_url);
        disp(a)
        disp(b)
        if ~a
            headerstr = '#Network';
            indhead = strfind(b,headerstr);
            if isempty(indhead)
                continue
            end
            fid = fopen([sel_data.work_dir '/mseed_files/my-stations.txt'],'at');
            fprintf(fid,'%2d\n%s',ii,b(indhead:end));
            fclose(fid);
            C = textscan(b(indhead:end),'%s %s %f %f %f %s %s %s','Delimiter','|','HeaderLines',1);
            nc = [nc C{1,1}'];
            station = [station C{1,2}'];
            try
                starttime = datetime(C{1,7}{1}(1:19),'InputFormat','yyyy-MM-dd''T''HH:mm:ss');
            catch ME
                starttime = datetime(C{1,7}{1});
            end
            svec{end+1} = datevec(starttime);
            dci = [dci ii*ones(size(C{1,1}))'];
        end
    end
end

startstat = 1;
endstat = length(station);
for i = startstat:endstat
    dwnld = tic; 
    download_bdsr(sel_data,nc{i},station{i},dci(i),svec{i});
    disp([nc{i},'.',station{i} ' in ' num2str(toc(dwnld)) 's'])
end

end