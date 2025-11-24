function bdsr_stackrf(path2opts,sel_data)

% control function to gather nescessary parameters and invoke the
% quality check of processed data
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

% make folder for processed data
sel_data.load_dir = [sel_data.work_dir '/rf/'];
sel_data.save_dir = [sel_data.work_dir '/rf/'];
if ~exist(sel_data.save_dir)
    mkdir(sel_data.save_dir);
end

% produce list of stations to be processed
if isempty(sel_data.statstr)
    A = dir([sel_data.load_dir '/*.*.mat']);
else
    A = dir([sel_data.load_dir '/' sel_data.statstr '.mat']);
end
n = 1;
for i = 1:length(A)
    if length(A(i).name) < 4
        continue
    end
    if contains(A(i).name,'_rf')
        continue
    end
    ii = strfind(A(i).name,'.');
    if ii(1) ~= 3
        continue
    end
    if exist([sel_data.load_dir A(i).name],'file')
        station{n} = strrep(A(i).name(4:end),'.mat','');
        nc{n} = A(i).name(1:2);
        n = n+1;
    end
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
    % if ~strcmp(station{n},'CS07')
    %     continue
    % end
% for n = 12:12%length(station)
    rfchk = tic;
    Do_isoHk_main2(sel_data,nc{n},station{n})
    disp([nc{n},'.',station{n} ' in ' num2str(toc(rfchk)) 's'])
end
% if sel_data.stackflag
%     for n = 1:length(station)
%         stack_bdsr(sel_data,nc{n},station{n});
%         disp([nc{n},'.',station{n}])
%     end
% end