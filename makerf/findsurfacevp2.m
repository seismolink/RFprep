function vsout = findsurfacevp2(sel_data,data,fn,mis)

disp("Search for surface velocity")
vpc = 6.5;
H = 1:2:150;
fac = 1;%4;% 0.33;
fcutz = sort(unique([vpc./(fac.*H),1]));
f2nfz = round(sel_data.nf./((sel_data.fmax/2)./fcutz));
% f2nfz = round(sel_data.nf./((sel_data.fmax/fac)./fcutz));
[f2nfz,ii] = unique(f2nfz);
fcutz = fcutz(ii);
[fcutz,ii] = sort(fcutz,'descend');
f2nfz = f2nfz(ii);
N = 500; % 100
% fcutz = 3; % ["0.5","1.0","1.5","2.0","3.0"]
% fcutz = [0.3,0.4,0.5,0.75,1.0];
% fcutz = linspace(0.05,9,12);
% for jj = 1:length(fcutz)
% % nfz = round(sel_data.nf/(3/(fcutz))); % [int(nf/6), int(nf/3), int(nf/2), int(nf/1.5), int(nf/1)]
% fcut = fcutz(jj);
% f2nf = round(sel_data.nf/((sel_data.fmax/2)/fcut));
% 
% % nf2Hz = round(nf/3*fmax);
% rramp = pi*((1:f2nf)-1)./f2nf;
% ctap = zeros(sel_data.nf,1);
% ctaper = cos(rramp)+1;
% ctap(1:f2nf) = ctaper;

% vpt = 1.5:0.5:8.5;
% vpt = 3:0.1:8;
vst = 1:0.1:5;
% Find correct upper layer velocity
for ii = 1:length(vst)%vpt)
    disp(['Testing velocity ' num2str(vst(ii))]) % vpt(ii))])
    drrf_sum = zeros(sel_data.nf,1);
    ftfft = zeros(sel_data.nf,1);
    for ij = 1:min([100,length(fn)])
        phase = data.(fn{ij});
        if phase.event.dist<25 || phase.event.dist>95
            continue
        end
        if length(mis.start)>1
            phi = mis.phi(mis.start<=phase.tt_abs&mis.end>=phase.tt_abs);
        else
            phi = mis.phi;
        end
        if isempty(phi)
            phi = mis.phi(end);
        end
        baz = phase.event.baz-phi;
        % Vpi = vpt(ii); % Near surface velocity used for LQT rotation
        rayp = phase.pp./111.11;
        % if Vpi*rayp > 1
        %     incl = 89;
        % else
        %     incl = asind(Vpi*rayp);
        % end
        Vsi = vst(ii);
        rho = 3;
        [incl] = get_theo_inc(baz,rayp,Vsi.*1.73,Vsi,rho);
        for i = 1:3
        phase.tr(i).data = detrend(phase.tr(i).data).*tukeywin(length(phase.tr(i).data),0.005)';
        % phase.tr(i).data = buttern_high(phase.tr(i).data,2,0.1,phase.dt)';
        end
        [phase.tr(1).data,phase.tr(2).data,phase.tr(3).data]=rotVNE2VRT(phase.tr(1).data,phase.tr(2).data,phase.tr(3).data,baz);
        [phase.tr(1).data,phase.tr(2).data,phase.tr(3).data]=rotVRT2PSvSh(phase.tr(1).data,phase.tr(2).data,phase.tr(3).data,incl);
        phase.comp = 'lqt';
        [hr,~,rstd,~] = MTC2(phase,0,sel_data);
        
        % rstd = ones(size(rstd));
        sigsq = (abs(rstd).*abs(rstd));
        ftfft = ftfft+(hr./sigsq);
        drrf_sum = drrf_sum+(ones(sel_data.nf,1)./(sigsq));
    end
    
    for jj = 1:length(fcutz)
        % fcut = fcutz(jj);
        % f2nf = round(sel_data.nf/((sel_data.fmax/2)/fcut));
        f2nf = f2nfz(jj);
        rramp = pi*((1:f2nf)-1)./f2nf;
        ctap = zeros(sel_data.nf,1);
        ctaper = cos(rramp)+1;
        ctap(1:f2nf) = ctaper;
        % temp = tukeywin(2*f2nf,0.5);
        % ctap2 = zeros(sel_data.nf,1);
        % ctap2(1:f2nf) = temp(f2nf+1:end);
        % ctest{jj} = ctap;
        % ctest2{jj} = ctap2;
        % ctap = ctap2;
        bb = real(ifft(ctap.*(ftfft./drrf_sum),sel_data.dpad));
        nnorm = sel_data.dpts2./sel_data.nf;
        bb = nnorm*bb;

        amp(ii,jj) = sum(abs(bb((end-20):sel_data.dpad).^2));
        % amp(ii) = sum(abs([bb((end-100):sel_data.dpad)',bb(1:50)']'.^2));
        amp4pl{ii,jj} = [bb((end-N):sel_data.dpad)',bb(1:50)']';
    end
end
for jj = 1:length(fcutz)
amp0(jj) = min(amp(:,jj));
% vpf(jj) = vpt(amp0(jj)==amp(:,jj));
% disp(vpf(jj));
vsf(jj) = vst(amp0(jj)==amp(:,jj));
disp(vsf(jj));
end

fname=[sel_data.save_dir '/' sel_data.nc '.' sel_data.station '-' sel_data.vpsfile];
fid = fopen(fname,'wt');
for i = 1:length(fcutz)
fprintf(fid,'%f %f %f\n',fcutz(i),vsf(i),(vsf(i).*1.73)./(fac.*fcutz(i)));
end
fclose(fid);

fig = figure; 
plot((vsf.*1.73)./(fac.*fcutz),vsf,'.')
hold on
plot((vsf.*1.73)./(fac.*fcutz),vsf.*1.73,'r.')
xlabel('Pseudodepth in [km]')
ylabel('Velocity in [km/s]')
legend({'vs','vp_e'})

fname=[sel_data.save_dir '/' sel_data.nc '.' sel_data.station '-VelocityPseudoProfile.png'];
print(fig,fname,'-dpng','-r600');
close(fig);

foi = sel_data.fcut;
[~,ioi] = min(abs(fcutz-foi));

vsout = vsf(ioi);

fig = figure;
plot(vst,amp(:,ioi),'color','k','LineWidth',1.2)
hold on
plot([vsf(ioi) vsf(ioi)],[0.9.*min(amp(:,ioi)) 1.1.*max(amp(:,ioi))],'Color','red','LineWidth',1.2,'LineStyle','--');
title(['Surface Velocity ' num2str(vsf(ioi)) ' km/s at ' num2str(foi) ' Hz'])
xlabel('vs [km/s]');
xlim([min(vst) max(vst)])
ylim([0.9*min(amp(:,ioi)) 1.1*max(amp(:,ioi))])

fname=[sel_data.save_dir '/' sel_data.nc '.' sel_data.station '-SurfaceVelocity_' num2str(foi) 'Hz.png'];
print(fig,fname,'-dpng','-r600');
close(fig);

cols = jet(length(vst));
fig = figure;
hold on
for i = 1:length(vst)
    plot((-N:50).*phase.dt,amp4pl{i,ioi},'LineWidth',1.2,'Color',cols(i,:))
end
colormap(cols);
if datetime(version('-date')) >= datetime('March 9, 2022') % Matlab 2022a release
    clim([min(vst) max(vst)]);
else
    caxis([min(vst) max(vst)]);
end
    
colorbar
title(['Surface Velocity ' num2str(vsf(ioi)) ' km/s at ' num2str(foi) ' Hz'])
xlabel('time in [s]');
xlim([min((-100:50).*phase.dt) max((-100:50).*phase.dt)])

fname=[sel_data.save_dir '/' sel_data.nc '.' sel_data.station '-SurfaceVelocity_' num2str(foi) 'Hz_amp.png'];
print(fig,fname,'-dpng','-r600');
close(fig);

disp("Surface velocity identified") 
% 
% fid = fopen([sel_data.save_dir '/' sel_data.vpsfile],'at');
% fprintf(fid,'%s,%f\n',[sel_data.nc '.' sel_data.station],vsf);
% fclose(fid);

end