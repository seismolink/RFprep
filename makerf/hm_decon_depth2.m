function hm = hm_decon_depth2(rf,z_out,ldepth,velp,vels)
% Harmonic decomposition after mapping to depth

zpts = length(z_out);
tt = rf.time;

bootg = zeros(10,10,zpts);
bootd = zeros(10,zpts);

for i = 1:length(rf.tr)
    if rf.N(i) < 2
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

    % Depth mapping
    out = depthmapping(tt,rf.tr(i).rrf,z_out,ldepth,velp,vels,slow);
    fnorm = max(abs(out));
    bootg = bootg+bsxfun(@times,bfdyadr,reshape(ones(size(out)),1,1,numel(out)));
    bootd = bootd+bfr'*out./fnorm;

    % Depth mapping
    out = depthmapping(tt,rf.tr(i).trf,z_out,ldepth,velp,vels,slow);
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