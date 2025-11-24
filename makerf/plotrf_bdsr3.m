function plotrf_bdsr3(sel_data,nc,station)

if exist([sel_data.save_dir '/' nc '.' station '_rf.mat'],'file')
    disp('File exists. Only preparing plots.')
    return
end

% calculate receiver functions from processed data

sel_data.nc = nc;
sel_data.station = station;
sel_data.drate=40;
sel_data.tpre=20;
sel_data.tdur=80;
sel_data.tdur2=80;
% sel_data.fmax=6.0;
sel_data.fmax=24.0;

[sel_data] = calcconst(sel_data);

%nfz = [int(nf/6), int(nf/3), int(nf/2), int(nf/1.5), int(nf/1)]
%fcutz = ["0.5","1.0","1.5","2.0","3.0"]
fcutz = sel_data.fcut; % ["0.5","1.0","1.5","2.0","3.0"]
%nfz = round(sel_data.nf/(3/(fcutz))); % [int(nf/6), int(nf/3), int(nf/2), int(nf/1.5), int(nf/1)]
fcut = fcutz;
% f2nf = round(sel_data.nf/(3/fcut));
f2nf = round(sel_data.nf/((sel_data.fmax/2)/fcut));

% Copyright 2024 F.Link, J.Wolf and M.Reiss

% load preprocessed data
rf = load([sel_data.load_dir,'/',sel_data.nc,'.',sel_data.station,'.mat']);

rfn.lat = rf.lat;
rfn.lon = rf.lon;
rfn.stat = rf.stat;
rfn.net = rf.net;

% read event ID
fn = fieldnames(rf);
fn(~contains(fn,'ev')) = [];

% back azimuthal sweep
% compute 36 composite RFs spaced:  baz = 5,355,10
ddeg = 10;
depi = 20;
epis = 25:depi/2:155;
% epis = 25:depi:155;

rramp = pi*((1:f2nf)-1)./f2nf;
ctap = zeros(sel_data.nf,1);
ctaper = cos(rramp)+1;
ctap(1:f2nf) = ctaper./2;
nnmt = 1;
wwin = tukeywin((1+nnmt)*sel_data.nmt,0.5);
rfn.time = sel_data.dt*((1:sel_data.dpts2)-sel_data.npre)';
n = 1;
for ii = 1:72%36
    baz = 5+ddeg/2*(ii-1);
    % baz = 5+ddeg*(ii-1);
    for iii = 1:length(epis)
        epi = epis(iii);
        jrec = 0;
        rff = zeros(length(rf.(fn{1}).zt),sel_data.nf);
        drff = rff;
        rfft = rff;
        drfft = rff;
        tsh = zeros(length(rf.(fn{1}).zt),1);
        slow = 0;
        bazf = 0;
        % bazf2 = 0;
        epif = 0;
        for jj = 1:length(fn)
            phase = rf.(fn{jj});
            ebaz = mod(phase.event.baz,360);
            test = abs(baz-ebaz);
            eepi = phase.event.dist;
            % if eepi > 115
            %     continue
            % end
            test2 = abs(epi-eepi);%abs(baz-ebaz);
            if (test<ddeg || test >360-ddeg) && test2<depi
                jrec = jrec+1;
                for iz = 1:length(phase.zt)
                    sigsq = abs(phase.rstd(iz).data).*abs(phase.rstd(iz).data);
                    rff(iz,:) = rff(iz,:) + (phase.rrf(iz).data./sigsq)';
                    drff(iz,:) = drff(iz,:) + (ones(sel_data.nf,1)./sigsq)';
                    sigsqt = abs(phase.tstd(iz).data).*abs(phase.tstd(iz).data);
                    rfft(iz,:) = rfft(iz,:) + (phase.trf(iz).data./sigsqt)';
                    drfft(iz,:) = drfft(iz,:) + (ones(sel_data.nf,1)./sigsqt)';
                    tsh(iz) = tsh(iz) + phase.tshift(iz);
                end
                slow = slow+rf.(fn{(jj)}).pp/111.11;
                if abs(baz-ebaz)>180
                    bazf = bazf+ebaz-360;
                    % bazf2 = bazf2+ebaz;
                else
                    bazf = bazf+ebaz;
                    % bazf2 = bazf2+ebaz+360;
                end
                epif = epif+eepi;
            end
        end
        if jrec > 0
            rft = zeros(sel_data.dpts2,1);
            rftt = rft;
            ww = rft;
            slow = slow./jrec;
            bazf = mod(bazf./jrec,360);
            % bazf2 = bazf2./jrec;
            % disp([num2str(bazf) ' ' num2str(bazf2)])
            epif = epif./jrec;
            for iz = 1:length(tsh)
                tsn = tsh(iz)/jrec;
                rff(iz,:) = rff(iz,:)./drff(iz,:);
                ftfft = ctap.*rff(iz,:)';
                bb = real(ifft(ftfft,sel_data.dpad));
                rfft(iz,:) = rfft(iz,:)./drfft(iz,:);
                ftfftt = ctap.*rfft(iz,:)';
                bbt = real(ifft(ftfftt,sel_data.dpad));
                m1 = sel_data.npre+round(tsn/sel_data.dt)-sel_data.nmt+1;
                m2 = sel_data.npre+round(tsn/sel_data.dt)+nnmt*sel_data.nmt;
                if m1 <= 0
                    dm = -m1+2;
                    m1 = 1;
                else
                    dm = 1;
                end
                if m2>length(rft)
                    dm2 = m2-length(rft);
                else
                    dm2 = 0;
                end
                if (-sel_data.nmt+dm)>=0
                    rft(m1:m2-dm2) = rft(m1:m2-dm2)+bb((end-sel_data.nmt+dm):(nnmt*sel_data.nmt-dm2)).*wwin(dm:(end-dm2));
                    rftt(m1:m2-dm2) = rftt(m1:m2-dm2)+bbt((end-sel_data.nmt+dm):(nnmt*sel_data.nmt-dm2)).*wwin(dm:(end-dm2));
                else
                    temp = [bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(nnmt*sel_data.nmt))']';
                    rft(m1:m2-dm2) = rft(m1:m2-dm2)+temp(1:end-dm2).*wwin(dm:(end-dm2));
                    tempt = [bbt((end-sel_data.nmt+dm):sel_data.dpad)',bbt(1:(nnmt*sel_data.nmt))']';
                    rftt(m1:m2-dm2) = rftt(m1:m2-dm2)+tempt(1:end-dm2).*wwin(dm:(end-dm2));
                end
                ww(m1:m2-dm2) = ww(m1:m2-dm2)+wwin(dm:(end-dm2));
            end
            ww(ww<1) = 1;
            rft = rft./ww;
            rftt = rftt./ww;
            fnorm = max(abs(rft));
            rfn.tr(n).rrf = rft/fnorm;
            rfn.tr(n).trf = rftt/fnorm;
            rfn.N(n) = jrec;
            rfn.baz(n) = bazf;
            rfn.slow(n) = slow;
            rfn.epi(n) = epif;
            n = n+1;
        end
    end
end
if ispc
    save([sel_data.save_dir '\' nc '.' station '_rf.mat'],'-struct','rfn','-v7.3')
else
    save([sel_data.save_dir '/' nc '.' station '_rf.mat'],'-struct','rfn','-v7.3')
end
end