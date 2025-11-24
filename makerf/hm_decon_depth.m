function [hm,rfn] = hm_decon_depth(rf,fn,sel_data,z_out,ldepth,velp,vels,fcut)
% Harmonic decomposition after mapping to depth
rfn.lat = rf.lat;
rfn.lon = rf.lon;
rfn.stat = rf.stat;
rfn.net = rf.net;
z_targets = rf.(fn{1}).zt;
zpts = length(z_out);
f2nf = round(sel_data.nf/(3/fcut));

% nf2Hz = round(nf/3*fmax);
rramp = pi*((1:f2nf)-1)./f2nf;
ctap = zeros(sel_data.nf,1);
ctaper = cos(rramp)+1;
ctap(1:f2nf) = ctaper;
tt = sel_data.dt*((1:sel_data.dpts2)-sel_data.npre)';
nnmt = 1;
wwin = tukeywin((1+nnmt)*sel_data.nmt,0.5);

rfn.time = tt;

nev = 0;
drrf_sum = zeros(length(z_targets),sel_data.nf);
dtrf_sum = zeros(length(z_targets),sel_data.nf);
for i = 1:length(fn)
    % if rf.(fn{i}).event.dist>115
    %     continue
    % end
    % temp = real(ifft(ctap.*rf.(fn{i}).rrf(rf.(fn{i}).zt==0).data,sel_data.dpad));
    % temp2 = real(ifft(ctap.*rf.(fn{i}).trf(rf.(fn{i}).zt==0).data,sel_data.dpad));
    % % rfflag = rfcrit(temp,temp2);
    % [rfflag,r1(i),r2(i),r3(i),r4(i),r5(i)] = rfcrit(temp,temp2);
    % if ~rfflag
    %     continue
    % end
    for iz = 1:length(z_targets)
        drrf_sum(iz,:) = drrf_sum(iz,:)+(ones(sel_data.nf,1)./(abs(rf.(fn{i}).rstd(iz).data).^2))';
        dtrf_sum(iz,:) = dtrf_sum(iz,:)+(ones(sel_data.nf,1)./(abs(rf.(fn{i}).tstd(iz).data).^2))';
    end
    nev = nev+1;
end
disp(num2str(nev));
drrf_sum = drrf_sum./(nev);
dtrf_sum = dtrf_sum./(nev);
bootg = zeros(10,10,zpts);
bootd = zeros(10,zpts);

for i = 1:length(fn)
    % if rf.(fn{i}).event.dist>115
    %     continue
    % end
    % temp = real(ifft(ctap.*rf.(fn{i}).rrf(rf.(fn{i}).zt==0).data./abs(rf.(fn{i}).rstd(rf.(fn{i}).zt==0).data),sel_data.dpad));
    % temp2 = real(ifft(ctap.*rf.(fn{i}).trf(rf.(fn{i}).zt==0).data./abs(rf.(fn{i}).tstd(rf.(fn{i}).zt==0).data),sel_data.dpad));
    % temp = real(ifft(ctap.*rf.(fn{i}).rrf(rf.(fn{i}).zt==0).data,sel_data.dpad));
    % temp2 = real(ifft(ctap.*rf.(fn{i}).trf(rf.(fn{i}).zt==0).data,sel_data.dpad));
    % rfflag = rfcrit(temp,temp2);
    % if ~rfflag
    %     continue
    % end
    
    ebaz = rf.(fn{(i)}).event.baz;
    slow = rf.(fn{(i)}).pp/111.11;
    
    cbaz = cosd(ebaz);
    sbaz = sind(ebaz);
    c2baz = cbaz*cbaz-sbaz*sbaz;
    s2baz = 2*cbaz*sbaz;

    bfr = [1 cbaz sbaz c2baz s2baz 0 cbaz sbaz c2baz s2baz];
    bft = [0 sbaz -cbaz s2baz -c2baz 1 -sbaz cbaz -s2baz c2baz];
    
    bfdyadr = bfr'*bfr;
    bfdyadt = bft'*bft;

    rft = zeros(sel_data.dpts2,1);
    ww = zeros(sel_data.dpts2,1);
    
    for iz = 1:length(z_targets)
        rff = rf.(fn{(i)}).rrf(iz).data./(abs(rf.(fn{(i)}).rstd(iz).data).^2);
        rff = rff./drrf_sum(iz,:)';
        ftfft = ctap.*rff;
        tsn = rf.(fn{i}).tshift(iz);
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
            rft(m1:m2-dm2) = rft(m1:m2-dm2)+bb((-sel_data.nmt+dm):(nnmt*sel_data.nmt-dm2)).*wwin(dm:(end-dm2));
        else
            temp = [bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(nnmt*sel_data.nmt))']';
            % rft(m1:m2-dm2) = rft(m1:m2-dm2)+[bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(3*sel_data.nmt-dm2))']'.*wwin(dm:(end-dm2));
            rft(m1:m2-dm2) = rft(m1:m2-dm2)+temp(1:end-dm2).*wwin(dm:(end-dm2));
        end
        ww(m1:m2-dm2) = ww(m1:m2-dm2)+wwin(dm:(end-dm2));
    end
    ww(ww<1) = 1;
    rft = rft./ww;

    % Depth mapping
    out = depthmapping(tt,rft,z_out,ldepth,velp,vels,slow);
    fnorm = max(abs(out));
    bootg = bootg+bsxfun(@times,bfdyadr,reshape(ones(size(out)),1,1,numel(out)));
    bootd = bootd+bfr'*out./fnorm;
    rfn.tr(i).rrf = rft/fnorm;
    rfn.baz(i) = ebaz;
    rfn.slow(i) = slow;
    rfn.event(i) = rf.(fn{i}).event;

    out2 = out;
    ww2 = ww;

    rft = zeros(sel_data.dpts2,1);
    ww = zeros(sel_data.dpts2,1);

    for iz = 1:length(z_targets)
        rff = rf.(fn{(i)}).trf(iz).data./(abs(rf.(fn{(i)}).tstd(iz).data).^2);
        rff = rff./dtrf_sum(iz,:)';
        ftfft = ctap.*rff;
        tsn = rf.(fn{i}).tshift(iz);
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
            rft(m1:m2-dm2) = rft(m1:m2-dm2)+bb((-sel_data.nmt+dm):(nnmt*sel_data.nmt-dm2)).*wwin(dm:(end-dm2));
        else
            temp = [bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(nnmt*sel_data.nmt))']';
            % rft(m1:m2-dm2) = rft(m1:m2-dm2)+[bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(3*sel_data.nmt-dm2))']'.*wwin(dm:(end-dm2));
            rft(m1:m2-dm2) = rft(m1:m2-dm2)+temp(1:end-dm2).*wwin(dm:(end-dm2));
        end
        ww(m1:m2-dm2) = ww(m1:m2-dm2)+wwin(dm:(end-dm2));
    end
    ww(ww<1) = 1;
    rft = rft./ww;
    rfn.tr(i).trf = rft./fnorm;

    % Depth mapping
    out = depthmapping(tt,rft,z_out,ldepth,velp,vels,slow);
    bootg = bootg+bsxfun(@times,bfdyadt,reshape(ones(size(out)),1,1,numel(out)));
    bootd = bootd+bft'*out./fnorm;
    % keyboard
end

hm = zeros(10,zpts);
for n = 1:zpts
    bootgg = bootg(:,:,n);
    bootdd = bootd(:,n);
    bootmm = bootgg\bootdd;
    hm(:,n) = hm(:,n)+bootmm;
end

end