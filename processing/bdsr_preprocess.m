function bdsr_preprocess(path2opts,sel_data)

% control function to gather nescessary parameters and invoke the
% processing of data
%
if isfolder(path2opts)
    sel_data.work_dir = path2opts;
    sel_data.optsfile = [path2opts '/RFprep_config.txt'];
else
    sel_data.work_dir = fileparts(path2opts);
    sel_data.optsfile = path2opts;
end

cur_dir = pwd;
sel_data.defaults_path = strcat(cur_dir,'/defaults/');

% read relevant information of options file
sel_data = readoptsfile(sel_data);

% make folder for processed data
sel_data.save_dir = [sel_data.work_dir '/processed/'];
if ~exist(sel_data.save_dir,'dir')
    mkdir(sel_data.save_dir);
end

% produce list of stations to be processed
if isempty(sel_data.statstr)
    A = dir([sel_data.work_dir '/mseed_files/*.*']);
    if isempty(A)
        A = dir([sel_data.work_dir '/processed/*.*pre.mat']);
        preflag = 1;
    else
        preflag = 0;
    end
    n = 1;
    for i = 1:length(A)
        if length(A(i).name) < 4
            continue
        end
        ii = strfind(A(i).name,'.');
        if ii(1) ~= 3
            continue
        end
        if ~preflag
        if ~exist([sel_data.save_dir A(i).name '.mat'],'file')
            station{n} = A(i).name(4:end);
            nc{n} = A(i).name(1:2);
            n = n+1;
        end
        else
            stattemp = A(i).name(4:end-7);
            nettemp = A(i).name(1:2);
            if ~exist([sel_data.save_dir nettemp '.' stattemp '.mat'],'file')
                station{n} = stattemp;
                nc{n} = nettemp;
                n = n+1;
            end
        end
    end
else
    station{1} = sel_data.statstr(4:end);
    nc{1} = sel_data.statstr(1:2);
end

% pp = gcp('nocreate');
% if isempty(pp)
%     if isempty(getenv('SLURM_CPUS_ON_NODE'))
%         nWorkers = feature('NumCores');
%     else
%         nWorkers = str2double(getenv('SLURM_CPUS_ON_NODE'));
%     end
%     pp = parpool(nWorkers);
% end
% parfor n = 1:length(station)
for n = 1:length(station)
    if exist([sel_data.save_dir nc{n} '.' station{n} '.mat'],'file')
        disp('Station already processed')
        continue
    end
    prep = tic;
    preprocess_bdsr(sel_data,nc{n},station{n});
    disp([nc{n},'.',station{n} ' in ' num2str(toc(prep)) 's'])
end

end