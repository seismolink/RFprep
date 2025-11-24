function [sel_data] = readoptsfile(sel_data)

sel_data.statlist = 'gmap-stations.txt';
sel_data.minlat = 0;
sel_data.maxlat = 0;
sel_data.minlon = 0;
sel_data.maxlon = 0;
sel_data.Eventlist = '';
sel_data.mag = 0;
sel_data.min_dist = 0;
sel_data.max_dist = 0;
sel_data.auth = '';
sel_data.p1 = 0;
sel_data.p2 = 0;
sel_data.snrlim = 3;
sel_data.prfflag = 1;
sel_data.srfflag = 0;
sel_data.misfile = 0;
sel_data.vpsfile = 0;
sel_data.fcut = 1;
sel_data.Hkfile = 'results_isoHk.txt';
sel_data.hmin = 20;
sel_data.hmax = 80;
sel_data.kmin = 1.5;
sel_data.kmax = 2.1;
sel_data.wflag = 1;
sel_data.plepi = 1;
sel_data.plbaz = 1;
sel_data.plhd = 1;
sel_data.wdist = 0.135;

fid = fopen(sel_data.optsfile,'rt');
line = fgetl(fid);
while ischar(line)
    % line = fgetl(fid);
    ix = strfind(line,'keyword');
    if ~isempty(ix)
        sel_data.keyword = str2double(line(ix+lenkw+1:end));
    end
    ix = strfind(line,'statlist');
    if ~isempty(ix)
        sel_data.statlist = line(ix+8+1:end);
    end
    ix = strfind(line,'minlat');
    if ~isempty(ix)
        sel_data.minlat = str2double(line(ix+6+1:end));
    end
    ix = strfind(line,'maxlat');
    if ~isempty(ix)
        sel_data.maxlat = str2double(line(ix+6+1:end));
    end
    ix = strfind(line,'minlon');
    if ~isempty(ix)
        sel_data.minlon = str2double(line(ix+6+1:end));
    end
    ix = strfind(line,'maxlon');
    if ~isempty(ix)
        sel_data.maxlon = str2double(line(ix+6+1:end));
    end
    ix = strfind(line,'Eventlist');
    if ~isempty(ix)
        sel_data.Eventlist = line(ix+9+1:end);
    end
    ix = strfind(line,'mag');
    if ~isempty(ix)
        sel_data.mag = str2double(line(ix+3+1:end));
    end
    ix = strfind(line,'min_dist');
    if ~isempty(ix)
        sel_data.min_dist = str2double(line(ix+8+1:end));
    end
    ix = strfind(line,'max_dist');
    if ~isempty(ix)
        sel_data.max_dist = str2double(line(ix+8+1:end));
    end
    ix = strfind(line,'auth');
    if ~isempty(ix)
        sel_data.auth = line(ix+4+1:end);
    end
    ix = strfind(line,'min_per');
    if ~isempty(ix)
        sel_data.p1 = str2double(line(ix+7+1:end));
    end
    ix = strfind(line,'max_per');
    if ~isempty(ix)
        sel_data.p2 = str2double(line(ix+7+1:end));
    end
    ix = strfind(line,'snrlim');
    if ~isempty(ix)
        sel_data.snrlim = str2double(line(ix+6+1:end));
    end
    ix = strfind(line,'prfflag');
    if ~isempty(ix)
        sel_data.prfflag = str2double(line(ix+7+1:end));
    end
    ix = strfind(line,'srfflag');
    if ~isempty(ix)
        sel_data.srfflag = str2double(line(ix+7+1:end));
    end
    ix = strfind(line,'misfile');
    if ~isempty(ix)
        sel_data.misfile = (line(ix+7+1:end));
    end
    ix = strfind(line,'vpsfile');
    if ~isempty(ix)
        sel_data.vpsfile = (line(ix+7+1:end));
    end
    ix = strfind(line,'manflag');
    if ~isempty(ix)
        sel_data.manflag = str2double(line(ix+7+1:end));
    end
    ix = strfind(line,'fcut');
    if ~isempty(ix)
        sel_data.fcut = str2double(line(ix+4+1:end));
    end
    ix = strfind(line,'Hkfile');
    if ~isempty(ix)
        sel_data.Hkfile = (line(ix+6+1:end));
    end
    ix = strfind(line,'hmin');
    if ~isempty(ix)
        sel_data.hmin = str2double(line(ix+4+1:end));
    end
    ix = strfind(line,'hmax');
    if ~isempty(ix)
        sel_data.hmax = str2double(line(ix+4+1:end));
    end
    ix = strfind(line,'kmin');
    if ~isempty(ix)
        sel_data.kmin = str2double(line(ix+4+1:end));
    end
    ix = strfind(line,'kmax');
    if ~isempty(ix)
        sel_data.kmax = str2double(line(ix+4+1:end));
    end
    ix = strfind(line,'wflag');
    if ~isempty(ix)
        sel_data.wflag = str2double(line(ix+5+1:end));
    end
    ix = strfind(line,'plepi');
    if ~isempty(ix)
        sel_data.plepi = str2double(line(ix+5+1:end));
    end
    ix = strfind(line,'plbaz');
    if ~isempty(ix)
        sel_data.plbaz = str2double(line(ix+5+1:end));
    end
    ix = strfind(line,'plhd');
    if ~isempty(ix)
        sel_data.plhd = str2double(line(ix+4+1:end));
    end
    ix = strfind(line,'wdist');
    if ~isempty(ix)
        sel_data.wdist = str2double(line(ix+5+1:end));
    end
    line = fgetl(fid);
end
fclose(fid);
end