function [hm] = hm_decon_depth(rf,fn,sel_data,z_out,ldepth,velp,vels,fcut)
% Harmonic decomposition after mapping to depth

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

drrf_sum = zeros(length(z_targets),sel_data.nf);
dtrf_sum = zeros(length(z_targets),sel_data.nf);
for i = 1:length(fn)
    temp = real(ifft(ctap.*rf.(fn{i}).rrf(rf.(fn{i}).zt==0).data,sel_data.dpad));
    temp2 = real(ifft(ctap.*rf.(fn{i}).trf(rf.(fn{i}).zt==0).data,sel_data.dpad));
    rfflag = rfcrit(temp,temp2);
    if ~rfflag
        continue
    end
    for iz = 1:length(z_targets)
        drrf_sum(iz,:) = drrf_sum(iz,:)+(ones(sel_data.nf,1)./(abs(rf.(fn{i}).rstd(iz).data).^2))';
        dtrf_sum(iz,:) = dtrf_sum(iz,:)+(ones(sel_data.nf,1)./(abs(rf.(fn{i}).tstd(iz).data).^2))';
    end
end
drrf_sum = drrf_sum./(length(fn));
dtrf_sum = dtrf_sum./(length(fn));

bootg = zeros(length(z_targets),10,10,sel_data.nf);
bootd = zeros(length(z_targets),10,sel_data.nf);
tsn = zeros(length(z_targets),1);
nev = 0;
for i = 1:length(fn)
    temp = real(ifft(ctap.*rf.(fn{i}).rrf(rf.(fn{i}).zt==0).data,sel_data.dpad));
    temp2 = real(ifft(ctap.*rf.(fn{i}).trf(rf.(fn{i}).zt==0).data,sel_data.dpad));
    rfflag = rfcrit(temp,temp2);
    if ~rfflag
        continue
    end
    
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
    
    for iz = 1:length(z_targets)
        rff = (rf.(fn{(i)}).rrf(iz).data./(abs(rf.(fn{(i)}).rstd(iz).data).^2))';
        rff = (rff./drrf_sum(iz,:));
        tsn(iz) = tsn(iz)+rf.(fn{i}).tshift(iz);
        % bootg(iz,:,:,:) = squeeze(bootg(iz,:,:,:))+bsxfun(@times,bfdyadr,reshape(ones(size(rff')),1,1,numel(rff')));
        % bootd(iz,:,:) = squeeze(bootd(iz,:,:))+bfr'*(rff);

        rfft = (rf.(fn{(i)}).trf(iz).data./(abs(rf.(fn{(i)}).tstd(iz).data).^2))';
        rfft = (rfft./dtrf_sum(iz,:));
        tsn(iz) = tsn(iz)+rf.(fn{i}).tshift(iz);
        % bootg(iz,:,:,:) = squeeze(bootg(iz,:,:,:))+bsxfun(@times,bfdyadt,reshape(ones(size(rff')),1,1,numel(rff')));
        % bootd(iz,:,:) = squeeze(bootd(iz,:,:))+bft'*(rff);
        bootg(iz,:,:,:) = squeeze(bootg(iz,:,:,:)) + repmat(bfdyadr,1,1,sel_data.nf).*permute(repmat(ones(sel_data.nf,1),1,10,10),[3 2 1]) + repmat(bfdyadt,1,1,sel_data.nf).*permute(repmat(ones(sel_data.nf,1),1,10,10),[3 2 1]);
        bootd(iz,:,:) = squeeze(bootd(iz,:,:)) + bfr'*rff + bft'*rfft;
    end
    nev = nev+1;


    % for iz = 1:length(z_targets)
    %     rff = rf.(fn{(i)}).trf(iz).data./(abs(rf.(fn{(i)}).tstd(iz).data).^2);
    %     rff = rff./dtrf_sum;
    %     ftfft = ctap.*rff;
    %     tsn = rf.(fn{i}).tshift(iz);
    %     bb = real(ifft(ftfft,sel_data.dpad));
    %         m1 = sel_data.npre+round(tsn/sel_data.dt)-sel_data.nmt+1;
    %         m2 = sel_data.npre+round(tsn/sel_data.dt)+nnmt*sel_data.nmt;
    %         if m1 <= 0
    %             dm = -m1+2;
    %             m1 = 1;
    %         else
    %             dm = 1;
    %         end
    %         if m2>length(rft)
    %             dm2 = m2-length(rft);
    %         else
    %             dm2 = 0;
    %         end
    %         if (-sel_data.nmt+dm)>=0
    %             rft(m1:m2-dm2) = rft(m1:m2-dm2)+bb((end-sel_data.nmt+dm):(nnmt*sel_data.nmt-dm2)).*wwin(dm:(end-dm2));
    %         else
    %             temp = [bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(nnmt*sel_data.nmt))']';
    %             % rft(m1:m2-dm2) = rft(m1:m2-dm2)+[bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(3*sel_data.nmt-dm2))']'.*wwin(dm:(end-dm2));
    %             rft(m1:m2-dm2) = rft(m1:m2-dm2)+temp(1:end-dm2).*wwin(dm:(end-dm2));
    %         end
    %         ww(m1:m2-dm2) = ww(m1:m2-dm2)+wwin(dm:(end-dm2));
    % end
    % ww(ww<1) = 1;
    % rft = rft./ww;
    % 
    % % Depth mapping
    % out = depthmapping(tt,rft,z_out,ldepth,velp,vels,slow);
    % bootg = bootg+bsxfun(@times,bfdyadt,reshape(ones(size(out)),1,1,numel(out)));
    % bootd = bootd+bft'*out;
end
tsna = tsn./nev;

bootma = zeros(iz,10,sel_data.nf);
for i = 1:length(z_targets)
    for j = 1:sel_data.nf
        bootgg = squeeze(bootg(i,:,:,j));
        bootdd = squeeze(bootd(i,:,j));
        bootmm = bootgg\bootdd';
        bootma(i,:,j) = bootma(i,:,j)+bootmm';
    end
end

hmt = zeros(sel_data.dpts2,10);
ww = zeros(sel_data.dpts2,1);
for iz = 1:length(z_targets)
    for id = 1:10
        ftfft = ctap.*(squeeze(bootma(iz,id,:)));
        tsn = tsna(iz);
        bb = real(ifft(ftfft,sel_data.dpad));
        m1 = sel_data.npre+round(tsn/sel_data.dt)-sel_data.nmt+1;
        m2 = sel_data.npre+round(tsn/sel_data.dt)+nnmt*sel_data.nmt;
        if m1 <= 0
            dm = -m1+2;
            m1 = 1;
        else
            dm = 1;
        end
        if m2>sel_data.dpts2
            dm2 = m2-sel_data.dpts2;
        else
            dm2 = 0;
        end
        if (-sel_data.nmt+dm)>=0
            hmt(m1:m2-dm2,id) = hmt(m1:m2-dm2,id)+bb((-sel_data.nmt+dm):(nnmt*sel_data.nmt-dm2)).*wwin(dm:(end-dm2));
        else
            temp = [bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(nnmt*sel_data.nmt))']';
            % rft(m1:m2-dm2) = rft(m1:m2-dm2)+[bb((end-sel_data.nmt+dm):sel_data.dpad)',bb(1:(3*sel_data.nmt-dm2))']'.*wwin(dm:(end-dm2));
            hmt(m1:m2-dm2,id) = hmt(m1:m2-dm2,id)+temp(1:end-dm2).*wwin(dm:(end-dm2));
        end
        if id == 1
            ww(m1:m2-dm2) = ww(m1:m2-dm2)+wwin(dm:(end-dm2));
        end
    end
end
ww(ww<1) = 1;
for i = 1:10
hmt(:,i) = hmt(:,i)./ww;
hm(:,i) = depthmapping(tt,hmt(:,i),z_out,ldepth,velp,vels,slow);
end
hm = hm';
keyboard
% hm = zeros(10,zpts);
% for n = 1:zpts
%     bootgg = bootg(:,:,n);
%     bootdd = bootd(:,n);
%     bootmm = bootgg\bootdd;
%     hm(:,n) = hm(:,n)+bootmm;
% end


end