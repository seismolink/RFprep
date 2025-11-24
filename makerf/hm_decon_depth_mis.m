function [bma_list,stt_list,edt_list] = hm_decon_depth_mis(rf,fn,evtime,sel_data,z_out,ldepth,velp,vels,fcut,z_targets)
% Harmonic decomposition after mapping to depth
%def hm_decon_depth_mis(rrf,trf,drrf,dtrf,ttshft,dataset_save,z_out,ldepth, velp, vels,f_corner,z_targets):

% [dt,dpts,dpts2,tprior,nmt,npre,npost,dln2,dpad,tpad,df,nf] = calcconst(drate,tpre,tdur,tdur2,fmax);

zpts = length(z_out);
f2nf = round(sel_data.nf/((sel_data.fmax/2)/fcut));

% nf2Hz = round(nf/3*fmax);
rramp = pi*((1:f2nf)-1)./f2nf;
ctap = zeros(sel_data.nf,1);
ctaper = cos(rramp)+1;
ctap(1:f2nf) = ctaper;
tt = sel_data.dt*((1:sel_data.dpts2)-sel_data.npre)';
nnmt = 1;
wwin = tukeywin((1+nnmt)*sel_data.nmt,0.5);

drrf_sum = zeros(sel_data.nf,1);
dtrf_sum = zeros(sel_data.nf,1);
for i = 1:length(fn)
    for iz = 1:length(z_targets)
        drrf_sum = drrf_sum+ones(sel_data.nf,1)./(abs(rf.(fn{i}).rstd(iz).data).^2);
        dtrf_sum = dtrf_sum+ones(sel_data.nf,1)./(abs(rf.(fn{i}).tstd(iz).data).^2);
    end
end
drrf_sum = drrf_sum./(length(fn)+length(z_targets));
dtrf_sum = dtrf_sum./(length(fn)+length(z_targets));

[~,ii] = sort(evtime,'ascend');
for i = 1:length(evtime)
    bootg = zeros(10,10,zpts);
    bootd = zeros(10,zpts);
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
        rff = rf.(fn{(i)}).rrf(iz).data;%./(abs(rf.(fn{(i)}).rstd(iz).data).^2);
        %rff = rff./drrf_sum;
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

    % Depth mapping
    out = depthmapping(tt,rft,z_out,ldepth,velp,vels,slow);
    fnorm = max(abs(out));
    bootg = bootg+bsxfun(@times,bfdyadr,reshape(ones(size(out)),1,1,numel(out)));
    bootd = bootd+bfr'*out./fnorm;

    rft = zeros(sel_data.dpts2,1);
    ww = zeros(sel_data.dpts2,1);

    for iz = 1:length(z_targets)
        rff = rf.(fn{(i)}).trf(iz).data;%./(abs(rf.(fn{(i)}).tstd(iz).data).^2);
        %rff = rff./dtrf_sum;
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

    % Depth mapping
    out = depthmapping(tt,rft,z_out,ldepth,velp,vels,slow);
    bootg = bootg+bsxfun(@times,bfdyadt,reshape(ones(size(out)),1,1,numel(out)));
    bootd = bootd+bft'*out./fnorm;

    bd_list{i} = bootd;
    bg_list{i} = bootg;
    et_list(i) = evtime(ii(i));
end

idx = 1:length(bd_list);
nflag = 0;
m = 1;
for i = 1:length(et_list)
    starttime = et_list(i);
    endtime = starttime+days(12*30);
    if endtime>max(evtime) && nflag == 0
        nflag = 1;
    elseif endtime > max(evtime)
        break
    end
    bootg = zeros(10,10,zpts);
    bootd = zeros(10,zpts);
    ix = idx(et_list>=starttime&et_list<=endtime);
    if length(ix)<6
        continue
    end
    for j = 1:length(ix)
        bootg = bootg+bg_list{ix(j)};
        bootd = bootd+bd_list{ix(j)};
    end
    bootma = zeros(10,zpts);
    for n = 1:zpts
        bootgg = bootg(:,:,n);
        bootdd = bootd(:,n);
        bootmm = bootgg\bootdd;
        bootma(:,n) = bootma(:,n)+bootmm;
    end
    bma_list{m} = real(bootma);
    stt_list(m) = starttime;
    edt_list(m) = et_list(ix(end));
    m = m+1;
end

end