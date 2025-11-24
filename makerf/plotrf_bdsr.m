function plotrf_bdsr(sel_data,nc,station)

% calculate receiver functions from processed data

sel_data.nc = nc;
sel_data.station = station;
sel_data.drate=40;
sel_data.tpre=20;
sel_data.tdur=80;
sel_data.tdur2=80;
sel_data.fmax=6.0;

[sel_data] = calcconst(sel_data);

plotbazflag = 0;
plotepiflag = 0;

%nfz = [int(nf/6), int(nf/3), int(nf/2), int(nf/1.5), int(nf/1)]
%fcutz = ["0.5","1.0","1.5","2.0","3.0"]
fcutz = 1; % ["0.5","1.0","1.5","2.0","3.0"]
%nfz = round(sel_data.nf/(3/(fcutz))); % [int(nf/6), int(nf/3), int(nf/2), int(nf/1.5), int(nf/1)]
fcut = fcutz;
f2nf = round(sel_data.nf/(3/fcut));

% Copyright 2024 F.Link, J.Wolf and M.Reiss

% load preprocessed data
rf = load([sel_data.load_dir,'/',sel_data.nc,'.',sel_data.station,'.mat']);

% read event ID
fn = fieldnames(rf);
fn(~contains(fn,'ev')) = [];

if plotbazflag
% back azimuthal sweep
% compute 36 composite RFs spaced:  baz = 5,355,10
ddeg = 10;

rramp = pi*((1:f2nf)-1)./f2nf;
ctap = zeros(sel_data.nf,1);
ctaper = cos(rramp)+1;
ctap(1:f2nf) = ctaper./2;
nnmt = 1;
wwin = tukeywin((1+nnmt)*sel_data.nmt,0.5);
tt = sel_data.dt*((1:sel_data.dpts2)-sel_data.npre)';
mint = -5;
maxt = 10;
ntot = 0;
fig = figure;
subplot(1,2,1)
for ii = 1:36
    fnorm(ii) = 0;
    baz = 5+10*(ii-1);
    jrec = 0;
    rff = zeros(length(rf.(fn{1}).zt),sel_data.nf);
    drff = rff;
    tsh = zeros(length(rf.(fn{1}).zt),1);
    for jj = 1:length(fn)
        phase = rf.(fn{jj});
        ebaz = phase.event.baz;
        test = abs(baz-ebaz);
        if test<ddeg || test >360-ddeg
            temp = real(ifft(ctap.*phase.rrf(phase.zt==0).data,sel_data.dpad));
            temp2 = real(ifft(ctap.*phase.trf(phase.zt==0).data,sel_data.dpad));
            rfflag = rfcrit(temp,temp2);
            if ~rfflag
                disp('problematic trace')
                % keyboard
                % continue
            end
            % if sum(abs(temp(1:20)))/sum(abs(temp(end-29:end-10))) < 1 || sum(abs(temp(1:100)))/sum(abs(temp(end-99:end))) < 1 ...
            %         || sum(abs(temp2(1:100)))/sum(abs(temp2(end-99:end))) < 1 || max(abs(temp2(1:300)))/max(abs(temp2(end-29:end-10))) < 1
            %     continue
            % end
            ntot = ntot+1;
            jrec = jrec+1;
            for iz = 1:length(phase.zt)
                sigsq = abs(phase.rstd(iz).data).*abs(phase.rstd(iz).data);
                % sigsq = ones(size(sigsq));
                rff(iz,:) = rff(iz,:) + (phase.rrf(iz).data./sigsq)';
                drff(iz,:) = drff(iz,:) + (ones(sel_data.nf,1)./sigsq)';
                tsh(iz) = tsh(iz) + phase.tshift(iz);
            end
        end
    end
    if jrec > 0
        rft = zeros(sel_data.dpts2,1);
        ww = rft;
        for iz = 1:length(tsh)
            tsn = tsh(iz)/jrec;
            rff(iz,:) = rff(iz,:)./drff(iz,:);
            ftfft = ctap.*rff(iz,:)';
            bb = real(ifft(ftfft,sel_data.dpad));
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
            else
                temp = [bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(nnmt*sel_data.nmt))']';
                % rft(m1:m2-dm2) = rft(m1:m2-dm2)+[bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(3*sel_data.nmt-dm2))']'.*wwin(dm:(end-dm2));
                rft(m1:m2-dm2) = rft(m1:m2-dm2)+temp(1:end-dm2).*wwin(dm:(end-dm2));
            end
            ww(m1:m2-dm2) = ww(m1:m2-dm2)+wwin(dm:(end-dm2));
        end
        ww(ww<1) = 1;
        rft = rft./ww;
        fnorm(ii) = max(abs(rft));
        bbb = baz.*ones(sel_data.dpts2,1)+10*rft./fnorm(ii);
        plot(tt,bbb,'Color','black','LineWidth',0.5);
        hold on
        localzero = baz;
        fill_trace(tt',bbb',localzero,0,'b');
        fill_trace(tt',bbb',localzero,1,'r');
    end
end
title(['Radial RF @ ' num2str(fcut) ' Hz'])
xlabel("Delay Time (sec)")
ylabel("Back Azimuth")
xlim([mint,maxt])
ylim([-10.,370])
disp(ntot);

subplot(1,2,2)
for ii = 1:36
    baz = 5+10*(ii-1);
    jrec = 0;
    rff = zeros(length(rf.(fn{1}).zt),sel_data.nf);
    drff = rff;
    tsh = zeros(length(rf.(fn{1}).zt),1);
    for jj = 1:length(fn)
        phase = rf.(fn{jj});
        ebaz = phase.event.baz;
        test = abs(baz-ebaz);
        if test<ddeg || test >360-ddeg
            % temp = real(ifft(ctap.*phase.rrf(phase.zt==0).data,sel_data.dpad));
            % temp2 = real(ifft(ctap.*phase.trf(phase.zt==0).data,sel_data.dpad));
            % rfflag = rfcrit(temp,temp2);
            % if ~rfflag
            %     continue
            % end
            % if sum(abs(temp(1:20)))/sum(abs(temp(end-29:end-10))) < 1 || sum(abs(temp(1:100)))/sum(abs(temp(end-99:end))) < 1 ...
            %         || sum(abs(temp2(1:100)))/sum(abs(temp2(end-99:end))) < 1 || max(abs(temp2(1:300)))/max(abs(temp2(end-29:end-10))) < 1
            %     continue
            % end
            jrec = jrec+1;
            for iz = 1:length(phase.zt)
                sigsq = abs(phase.tstd(iz).data).*abs(phase.tstd(iz).data);
                % sigsq = ones(size(sigsq));
                rff(iz,:) = rff(iz,:) + (phase.trf(iz).data./sigsq)';
                drff(iz,:) = drff(iz,:) + (ones(sel_data.nf,1)./sigsq)';
                tsh(iz) = tsh(iz) + phase.tshift(iz);
            end
        end
    end
    if jrec > 0
        rft = zeros(sel_data.dpts2,1);
        ww = rft;
        for iz = 1:length(tsh)
            tsn = tsh(iz)/jrec;
            rff(iz,:) = rff(iz,:)./drff(iz,:);
            ftfft = ctap.*rff(iz,:)';
            bb = real(ifft(ftfft,sel_data.dpad));
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
            else
                temp = [bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(nnmt*sel_data.nmt))']';
                % rft(m1:m2-dm2) = rft(m1:m2-dm2)+[bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(3*sel_data.nmt-dm2))']'.*wwin(dm:(end-dm2));
                rft(m1:m2-dm2) = rft(m1:m2-dm2)+temp(1:end-dm2).*wwin(dm:(end-dm2));
            end
            ww(m1:m2-dm2) = ww(m1:m2-dm2)+wwin(dm:(end-dm2));
        end
        ww(ww<1) = 1;
        rft = rft./ww;
        bbb = baz.*ones(sel_data.dpts2,1)+10*rft./fnorm(ii);
        plot(tt,bbb,'Color','black','LineWidth',0.5);
        hold on
        localzero = baz;
        fill_trace(tt',bbb',localzero,0,'b');
        fill_trace(tt',bbb',localzero,1,'r');
    end
end
title(['Transverse RF @ ' num2str(fcut) ' Hz'])
xlabel("Delay Time (sec)")
% ylabel("Back Azimuth")
xlim([mint,maxt])
ylim([-10.,370])

fname=[sel_data.save_dir '/' sel_data.station '-' '-rfbaz.png'];
print(fig,fname,'-dpng','-r600');
disp("Azimuthal moveout printed")
close(fig);
end

if plotepiflag
% epicentral sweep
% compute composite RFs spaced:  epi = 25,155,10
depi = 10;
epis = 25:depi:155;

rramp = pi*((1:f2nf)-1)./f2nf;
ctap = zeros(sel_data.nf,1);
ctaper = cos(rramp)+1;
ctap(1:f2nf) = ctaper./2;
nnmt = 1;
wwin = tukeywin((1+nnmt)*sel_data.nmt,0.5);
tt = sel_data.dt*((1:sel_data.dpts2)-sel_data.npre)';
mint = -5;
maxt = 10;
ntot = 0;
fig = figure;
subplot(1,2,1)
for ii = 1:length(epis)
    fnorm(ii) = 0;
    epi = 25+depi*(ii-1);
    jrec = 0;
    nev = 0;
    rff = zeros(length(rf.(fn{1}).zt),sel_data.nf);
    drff = rff;
    rffsum = zeros(length(rf.(fn{1}).zt),sel_data.nf);
    drffsum = rffsum;
    tsh = zeros(length(rf.(fn{1}).zt),1);
    tshsum = tsh;
    for jj = 1:length(fn)
        phase = rf.(fn{jj});
        % ebaz = phase.event.baz;
        eepi = phase.event.dist;
        test = abs(epi-eepi);%abs(baz-ebaz);
        if test<depi
            % temp = real(ifft(ctap.*phase.rrf(phase.zt==0).data,sel_data.dpad));
            % temp2 = real(ifft(ctap.*phase.trf(phase.zt==0).data,sel_data.dpad));
            % rfflag = rfcrit(temp,temp2);
            % if ~rfflag
            %     continue
            % end
            % if sum(abs(temp(1:20)))/sum(abs(temp(end-29:end-10))) < 1 || sum(abs(temp(1:100)))/sum(abs(temp(end-99:end))) < 1 ...
            %         || sum(abs(temp2(1:100)))/sum(abs(temp2(end-99:end))) < 1 || max(abs(temp2(1:300)))/max(abs(temp2(end-29:end-10))) < 1
            %     continue
            % end
            ntot = ntot+1;
            jrec = jrec+1;
            for iz = 1:length(phase.zt)
                sigsq = abs(phase.rstd(iz).data).*abs(phase.rstd(iz).data);
                % sigsq = ones(size(sigsq));
                rff(iz,:) = rff(iz,:) + (phase.rrf(iz).data./sigsq)';
                drff(iz,:) = drff(iz,:) + (ones(sel_data.nf,1)./sigsq)';
                tsh(iz) = tsh(iz) + phase.tshift(iz);
            end
        end
        if ii == 1
            temp = real(ifft(ctap.*phase.rrf(phase.zt==0).data,sel_data.dpad));
            temp2 = real(ifft(ctap.*phase.trf(phase.zt==0).data,sel_data.dpad));
            rfflag = rfcrit(temp,temp2);
            if ~rfflag
                continue
            end
            nev = nev+1;
            for iz = 1:length(phase.zt)
                sigsq = abs(phase.rstd(iz).data).*abs(phase.rstd(iz).data);
                rffsum(iz,:) = rffsum(iz,:) + (phase.rrf(iz).data./sigsq)';
                drffsum(iz,:) = drffsum(iz,:) + (ones(sel_data.nf,1)./sigsq)';
                tshsum(iz) = tshsum(iz) + phase.tshift(iz);
            end
        end
    end
    if jrec > 0
        rft = zeros(sel_data.dpts2,1);
        ww = rft;
        for iz = 1:length(tsh)
            tsn = tsh(iz)/jrec;
            rff(iz,:) = rff(iz,:)./drff(iz,:);
            ftfft = ctap.*rff(iz,:)';
            bb = real(ifft(ftfft,sel_data.dpad));
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
            else
                temp = [bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(nnmt*sel_data.nmt))']';
                % rft(m1:m2-dm2) = rft(m1:m2-dm2)+[bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(3*sel_data.nmt-dm2))']'.*wwin(dm:(end-dm2));
                rft(m1:m2-dm2) = rft(m1:m2-dm2)+temp(1:end-dm2).*wwin(dm:(end-dm2));
            end
            ww(m1:m2-dm2) = ww(m1:m2-dm2)+wwin(dm:(end-dm2));
        end
        ww(ww<1) = 1;
        rft = rft./ww;
        fnorm(ii) = max(abs(rft));
        bbb = epi.*ones(sel_data.dpts2,1)+10*rft./fnorm(ii);
        plot(tt,bbb,'Color','black','LineWidth',0.5);
        hold on
        localzero = epi;
        fill_trace(tt',bbb',localzero,0,'b');
        fill_trace(tt',bbb',localzero,1,'r');
    end
    if ii == 1
        rft = zeros(sel_data.dpts2,1);
        ww = rft;
        for iz = 1:length(tshsum)
            tsn = tshsum(iz)/nev;
            rffsum(iz,:) = rffsum(iz,:)./drffsum(iz,:);
            ftfft = ctap.*rffsum(iz,:)';
            bb = real(ifft(ftfft,sel_data.dpad));
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
            else
                temp = [bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(nnmt*sel_data.nmt))']';
                % rft(m1:m2-dm2) = rft(m1:m2-dm2)+[bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(3*sel_data.nmt-dm2))']'.*wwin(dm:(end-dm2));
                rft(m1:m2-dm2) = rft(m1:m2-dm2)+temp(1:end-dm2).*wwin(dm:(end-dm2));
            end
            ww(m1:m2-dm2) = ww(m1:m2-dm2)+wwin(dm:(end-dm2));
        end
        ww(ww<1) = 1;
        rft = rft./ww;
        fnormsum = max(abs(rft));
        bbb = 170.*ones(sel_data.dpts2,1)+10*rft./fnormsum;
        plot(tt,bbb,'Color','black','LineWidth',0.5);
        hold on
        localzero = 170;
        fill_trace(tt',bbb',localzero,0,'b');
        fill_trace(tt',bbb',localzero,1,'r');
    end
end
title(['Radial RF @ ' num2str(fcut) ' Hz'])
xlabel("Delay Time (sec)")
ylabel("Epicentral Distance")
xlim([mint,maxt])
ylim([15.,180])
disp(ntot);

subplot(1,2,2)
for ii = 1:length(epis)
    epi = 25+depi*(ii-1);
    jrec = 0;
    nev = 0;
    rff = zeros(length(rf.(fn{1}).zt),sel_data.nf);
    drff = rff;
    rffsum = zeros(length(rf.(fn{1}).zt),sel_data.nf);
    drffsum = rffsum;
    tsh = zeros(length(rf.(fn{1}).zt),1);
    tshsum = tsh;
    for jj = 1:length(fn)
        phase = rf.(fn{jj});
        % ebaz = phase.event.baz;
        eepi = phase.event.dist;
        test = abs(epi-eepi);%abs(baz-ebaz);
        if test<depi
            % temp = real(ifft(ctap.*phase.rrf(phase.zt==0).data,sel_data.dpad));
            % temp2 = real(ifft(ctap.*phase.trf(phase.zt==0).data,sel_data.dpad));
            % rfflag = rfcrit(temp,temp2);
            % if ~rfflag
            %     continue
            % end
            % if sum(abs(temp(1:20)))/sum(abs(temp(end-29:end-10))) < 1 || sum(abs(temp(1:100)))/sum(abs(temp(end-99:end))) < 1 ...
            %         || sum(abs(temp2(1:100)))/sum(abs(temp2(end-99:end))) < 1 || max(abs(temp2(1:300)))/max(abs(temp2(end-29:end-10))) < 1
            %     continue
            % end
            jrec = jrec+1;
            for iz = 1:length(phase.zt)
                sigsq = abs(phase.tstd(iz).data).*abs(phase.tstd(iz).data);
                % sigsq = ones(size(sigsq));
                rff(iz,:) = rff(iz,:) + (phase.trf(iz).data./sigsq)';
                drff(iz,:) = drff(iz,:) + (ones(sel_data.nf,1)./sigsq)';
                tsh(iz) = tsh(iz) + phase.tshift(iz);
            end
        end
        if ii == 1
            temp = real(ifft(ctap.*phase.rrf(phase.zt==0).data,sel_data.dpad));
            temp2 = real(ifft(ctap.*phase.trf(phase.zt==0).data,sel_data.dpad));
            rfflag = rfcrit(temp,temp2);
            if ~rfflag
                continue
            end
            nev = nev+1;
            for iz = 1:length(phase.zt)
                sigsq = abs(phase.tstd(iz).data).*abs(phase.tstd(iz).data);
                rffsum(iz,:) = rffsum(iz,:) + (phase.trf(iz).data./sigsq)';
                drffsum(iz,:) = drffsum(iz,:) + (ones(sel_data.nf,1)./sigsq)';
                tshsum(iz) = tshsum(iz) + phase.tshift(iz);
            end
        end
    end
    if jrec > 0
        rft = zeros(sel_data.dpts2,1);
        ww = rft;
        for iz = 1:length(tsh)
            tsn = tsh(iz)/jrec;
            rff(iz,:) = rff(iz,:)./drff(iz,:);
            ftfft = ctap.*rff(iz,:)';
            bb = real(ifft(ftfft,sel_data.dpad));
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
            else
                temp = [bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(nnmt*sel_data.nmt))']';
                % rft(m1:m2-dm2) = rft(m1:m2-dm2)+[bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(3*sel_data.nmt-dm2))']'.*wwin(dm:(end-dm2));
                rft(m1:m2-dm2) = rft(m1:m2-dm2)+temp(1:end-dm2).*wwin(dm:(end-dm2));
            end
            ww(m1:m2-dm2) = ww(m1:m2-dm2)+wwin(dm:(end-dm2));
        end
        ww(ww<1) = 1;
        rft = rft./ww;
        bbb = epi.*ones(sel_data.dpts2,1)+10*rft./fnorm(ii);
        plot(tt,bbb,'Color','black','LineWidth',0.5);
        hold on
        localzero = epi;
        fill_trace(tt',bbb',localzero,0,'b');
        fill_trace(tt',bbb',localzero,1,'r');
    end
    if ii == 1
        rft = zeros(sel_data.dpts2,1);
        ww = rft;
        for iz = 1:length(tshsum)
            tsn = tshsum(iz)/nev;
            rffsum(iz,:) = rffsum(iz,:)./drffsum(iz,:);
            ftfft = ctap.*rffsum(iz,:)';
            bb = real(ifft(ftfft,sel_data.dpad));
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
            else
                temp = [bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(nnmt*sel_data.nmt))']';
                % rft(m1:m2-dm2) = rft(m1:m2-dm2)+[bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(3*sel_data.nmt-dm2))']'.*wwin(dm:(end-dm2));
                rft(m1:m2-dm2) = rft(m1:m2-dm2)+temp(1:end-dm2).*wwin(dm:(end-dm2));
            end
            ww(m1:m2-dm2) = ww(m1:m2-dm2)+wwin(dm:(end-dm2));
        end
        ww(ww<1) = 1;
        rft = rft./ww;
        bbb = 170.*ones(sel_data.dpts2,1)+10*rft./fnormsum;
        plot(tt,bbb,'Color','black','LineWidth',0.5);
        hold on
        localzero = 170;
        fill_trace(tt',bbb',localzero,0,'b');
        fill_trace(tt',bbb',localzero,1,'r');
    end
end
title(['Transverse RF @ ' num2str(fcut) ' Hz'])
xlabel("Delay Time (sec)")
% ylabel("Back Azimuth")
xlim([mint,maxt])
ylim([15,180])

fname=[sel_data.save_dir '/' sel_data.station '-' '-rfepi.png'];
print(fig,fname,'-dpng','-r600');
disp("Azimuthal moveout printed")
close(fig);
end


% 1-D velocity model for mapping to depth
ldepth = [20, 35, 300]; % depths of interfaces
velp = [5.8, 6.5, 8]; % P wave velocity
vels = [3.4, 3.75, 4.5]; % S wave velocity
z_out = -25:0.25:150;
zpts = length(z_out);

[hm,rfn] = hm_decon_depth(rf,fn,sel_data,z_out,ldepth,velp,vels,fcut);
disp(['Harmonics calculated for cut-off frequency ',num2str(fcut)])
hm = hm./max(abs(hm(:)));
save([sel_data.save_dir '\' nc '.' station '_rf.mat'],'-struct','rfn','-v7.3')

fig = figure;
subplot(1,2,1)
for ii = 1:5
    localzero = (4-ii+1)*10.0;
    rft = hm(ii,:);
    bbb = localzero*ones(1,zpts)+10*rft.*0.5;
    plot(z_out,bbb,'Color','black','LineWidth',0.5);
    hold on
    fill_trace(z_out,bbb,localzero,0,'b');
    fill_trace(z_out,bbb,localzero,1,'r');
end
plot([0.,0.],[0.,45.0],'k')
zz = 10:10:160;
for i = 1:length(zz)
    plot([zz(i),zz(i)],[0.,45.],'k--','LineWidth',0.5)
end
title(['Anisotropy/Dip RFs @ ' num2str(fcut) ' Hz corner'])
xlabel('Depth (km)')
xlim([-20,100])
xticks(-20:20:120)
ylim([-5.0,50.0])
yticks([0,10,20,30,40])
yticklabels(["sin(2*baz)","cos(2*baz)","sin(baz)","cos(baz)","constant"])


subplot(1,2,2)
for iii = 1:5
    localzero = (4-iii+1)*10.0;
    ii=iii+5;
    rft = hm(ii,:);
    bbb = localzero*ones(1,zpts)+10.*rft*0.5;
    plot(z_out,bbb,'Color','black','LineWidth',0.5);
    hold on
    fill_trace(z_out,bbb,localzero,0,'b');
    fill_trace(z_out,bbb,localzero,1,'r');
end
plot([0.,0.],[0.,45.0],'k')
zz = 10:10:160;
for i = 1:length(zz)
    plot([zz(i),zz(i)],[0.,45.],'k--','LineWidth',0.5)
end
title(['Unmodelled RFs @ ' num2str(fcut) ' Hz corner'])
xlabel('Depth (km)')
xlim([-20,100])
xticks(-20:20:120)
ylim([-5.0,50.0])
yticks([0,10,20,30,40])
yticklabels([])
% keyboard

fname=[sel_data.save_dir '/' sel_data.station '-' '-Harmonics.png'];
print(fig,fname,'-dpng','-r600');
disp("Azimuthal moveout printed")
close(fig);


end