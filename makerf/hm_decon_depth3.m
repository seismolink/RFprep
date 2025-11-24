function [res,err] = hm_decon_depth3(rf,z_out,ldepth,velp,vels,Nb,wflag)
% Harmonic decomposition after mapping to depth

zpts = length(z_out);
tt = rf.time;

% prepare data and calculate harmonic decomposition
trhm = zeros(length(rf.tr),10,zpts);
rhs = zeros(length(rf.tr),10,10);
for i = 1:length(rf.tr)
    if isnan(max(abs(rf.tr(i).rrf))) || isnan(max(abs(rf.tr(i).trf)))
        continue
    end
    ebaz = rf.baz(i);
    slow = rf.slow(i);
    
    cbaz = cosd(ebaz);
    sbaz = sind(ebaz);
    c2baz = cbaz*cbaz-sbaz*sbaz;
    s2baz = 2*cbaz*sbaz;

    bfr = [1 cbaz sbaz c2baz s2baz 0 cbaz sbaz c2baz s2baz];
    bft = [0 sbaz -cbaz s2baz -c2baz 1 -sbaz cbaz -s2baz c2baz];

    bfdyadr = bfr'*bfr;
    bfdyadt = bft'*bft;

    rhs(i,:,:) = bfdyadr+bfdyadt;

    outr = depthmapping(tt,rf.tr(i).rrf,z_out,ldepth,velp,vels,slow);
    fnorm = max(abs(outr));
    outt = depthmapping(tt,rf.tr(i).trf,z_out,ldepth,velp,vels,slow);

    trhm(i,:,:) = bfr'*outr./fnorm+bft'*outt./fnorm;

    Tr(i,:) = bfr;
    Tt(i,:) = bft;
    Dr(i,:) = outr;
    Dt(i,:) = outt;
end
T = [Tr' Tt']';
D = [Dr' Dt']';
hm = zeros(Nb,10,zpts);
% bootstrap loop
for i = 1:Nb
    tic
    if Nb == 1
        idx = 1:length(rf.tr);
    else
        idx = randi(length(rf.tr),length(rf.tr),1);
    end
    % Ndis = round(length(rf)/3);
    % idx = 1:length(rf.tr);
    % idx(randperm(length(idx),Ndis)) = [];
    % idx(end+1:length(rf)) = idx(randperm(length(idx),Ndis));

    % Wn = squeeze(sum(trhm(idx,:,:),1)).^2;
    % Wd = squeeze(sum(trhm(idx,:,:).^2,1));
    if wflag
        Wn = abs(squeeze(sum(trhm(idx,:,:),1)));
        Wd = squeeze(sum(abs(trhm(idx,:,:)),1));
    %     % Wn2 = abs(squeeze(sum(rhs(idx,:,:),1)));
    %     % Wd2 = squeeze(sum(abs(rhs(idx,:,:)),1));
    % 
    %     bootd = squeeze(sum(trhm(idx,:,:),1)).*(Wn./Wd);
    %     % bootg = bsxfun(@times,squeeze(sum(rhs(idx,:,:),1)),reshape(ones(zpts,1),1,1,zpts));
    %     % rhs2 = (squeeze(sum(rhs(idx,:,:),1)).*(Wn2./Wd2));
    else
        Wn = 1;
        Wd = 1;
    %     bootd = squeeze(sum(trhm(idx,:,:),1));
    %     % rhs2 = squeeze(sum(rhs(idx,:,:),1));
    end
    bootd = squeeze(sum(trhm(idx,:,:),1));

    % A_reshaped = reshape(bootg, 10, 10, 1, zpts);
    % B_reshaped = reshape(bootd, 10, 1, zpts);

    % keyboard

    % hm(i,:,:) = pagefun(@mldivide, squeeze(sum(rhs(idx,:,:),1)), B_reshaped);
    if wflag
        temp = ((bootd'/squeeze(sum(rhs(idx,:,:),1)).')');
        test2 = 1-abs(trhm(idx,:,:)-permute(pagemtimes(permute(rhs(idx,:,:),[3,2,1]),squeeze(temp)),[3,1,2]))./abs(trhm(idx,:,:));
        test2(test2<0) = 0;
        w2 = squeeze(mean(test2,1));
        w2 = (w2./max(abs(w2(:))));
        % w2(w2>0.8) = w2(w2>0.8).^0.5;
        % w2(w2<0.8) = w2(w2<0.8).^2;
        w2(w2>0.7) = 1;
        w2 = movmean(w2,round(zpts/50),2);
        hm(i,:,:) = temp;%.*w2;
        % hm2(i,:,:) = temp.*w2;
        temp3 = (((bootd.*w2)'/squeeze(sum(rhs(idx,:,:),1)).')');
        hm2(i,:,:) = temp3./max(abs(temp3(:)));
        % hm2(i,:,:) = temp.*Wn./Wd;
        TTinv_T = (T([idx' idx'+length(rf.tr)]',:)'*T([idx' idx'+length(rf.tr)]',:))\T([idx' idx'+length(rf.tr)]',:)';
        mest = TTinv_T*D([idx' idx'+length(rf.tr)]',:);
        r = D([idx' idx'+length(rf.tr)]',:)'-(T([idx' idx'+length(rf.tr)]',:)*mest)';
        % rr = movmean(sum(r'.^2,1),round(zpts/30));
        rr = sum(r'.^2,1);
        r2 = D'-(T*mest)';
        TTinv_T2 = (T'*T)\T';
        dm = movmean(TTinv_T2*r2',round(zpts/30),2);
        % dm = TTinv_T2*r2';
        res_norm = (rr-min(rr(:)))./(max(rr(:))-min(rr(:)));
        dm_norm = (dm.^2-min(dm(:).^2))./(max(dm(:).^2)-min(dm(:).^2));
        w3 = 1./(1+ones(10,1)*res_norm./(dm_norm+eps));
        % w3 = 1./(1+ones(10,1)*res_norm);
        % alpha = 2;  % sensitivity to residuals
        % beta = 1;   % sensitivity to perturbation energy
        % w4 = exp(-alpha * res_norm) .* (dm_norm).^beta;
        % w3 = ones(10,1)*(1-sqrt(sum(r'.^2)));
        % w3 = w3./max(abs(w3(:)));
        hm3(i,:,:) = temp.*w3;
    else
        hm(i,:,:) = ((bootd'/squeeze(sum(rhs(idx,:,:),1)).')').*Wn./Wd;
    end
    % hm(i,:,:) = (bootd'/rhs2.')';

    % for n = 1:zpts
    %     bootgg = bootg(:,:,n);
    %     bootdd = bootd(:,n);
    %     bootmm = bootgg\bootdd;
    %     hm(i,:,n) = bootmm;
    % end

    % hm(i,:,:) = bsxfun(@ldivide,bootg,reshape(bootd,1,10,zpts));
    disp([num2str(i) ' ' num2str(toc)])
end
% if wflag
%     % test = pagemtimes(rhs(idx,:,:),squeeze(hm(1,:,:)));%
%     % test = 1-sqrt((trhm-permute(pagemtimes(permute(rhs,[3,2,1]),squeeze(hm(1,:,:))),[3,1,2])).^2);
%     test = 1-abs(trhm-permute(pagemtimes(permute(rhs,[3,2,1]),squeeze(hm(1,:,:))),[3,1,2]))./abs(permute(pagemtimes(permute(rhs,[3,2,1]),squeeze(hm(1,:,:))),[3,1,2]));
%     test(test<0) = 0;
%     w = squeeze(mean(test,1));
%     test2 = 1-abs(trhm-permute(pagemtimes(permute(rhs,[3,2,1]),squeeze(hm(1,:,:))),[3,1,2]))./abs(trhm);
%     test2(test2<0) = 0;
%     w2 = squeeze(mean(test2,1));
% end
% if wflag
%     Wn = abs(squeeze(sum(hm,1)));
%     Wd = squeeze(sum(abs(hm),1));
%     res =squeeze(mean(hm,1)).*Wn./Wd;
%     err = squeeze(std(hm,1)).*Wn./Wd;
% else
res = squeeze(mean(hm,1));
err = squeeze(std(hm,[],1));
if wflag
    res = squeeze(mean(hm2,1));
    % res2 = squeeze(mean(hm2,1));
    err = squeeze(std(hm2,[],1));
end
% end

% bootg = zeros(10,10,zpts);
% bootd = zeros(10,zpts);
% 
% for i = 1:length(rf.tr)
%     % if rf.N(i) < 2
%     %     continue
%     % end
% 
%     ebaz = rf.baz(i);
%     slow = rf.slow(i);
% 
%     cbaz = cosd(ebaz);
%     sbaz = sind(ebaz);
%     c2baz = cbaz*cbaz-sbaz*sbaz;
%     s2baz = 2*cbaz*sbaz;
% 
%     bfr = [1 cbaz sbaz c2baz s2baz 0 cbaz sbaz c2baz s2baz];
%     bft = [0 sbaz -cbaz s2baz -c2baz 1 -sbaz cbaz -s2baz c2baz];
% 
%     bfdyadr = bfr'*bfr;
%     bfdyadt = bft'*bft;
% 
%     % Depth mapping
%     out = depthmapping(tt,rf.tr(i).rrf,z_out,ldepth,velp,vels,slow);
%     fnorm = max(abs(out));
%     bootg = bootg+bsxfun(@times,bfdyadr,reshape(ones(size(out)),1,1,numel(out)));
%     bootd = bootd+bfr'*out./fnorm;
% 
%     % Depth mapping
%     out = depthmapping(tt,rf.tr(i).trf,z_out,ldepth,velp,vels,slow);
%     bootg = bootg+bsxfun(@times,bfdyadt,reshape(ones(size(out)),1,1,numel(out)));
%     bootd = bootd+bft'*out./fnorm;
%     % keyboard
% end
% 
% hm = zeros(10,zpts);
% for n = 1:zpts
%     bootgg = bootg(:,:,n);
%     bootdd = bootd(:,n);
%     bootmm = bootgg\bootdd;
%     hm(:,n) = hm(:,n)+bootmm;
% end

end