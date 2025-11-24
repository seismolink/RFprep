function vpf = findsurfacevp(sel_data,data,fn,mis)

disp("Search for surface velocity")

fcutz = 3; % ["0.5","1.0","1.5","2.0","3.0"]
% nfz = round(sel_data.nf/(3/(fcutz))); % [int(nf/6), int(nf/3), int(nf/2), int(nf/1.5), int(nf/1)]
fcut = fcutz;
f2nf = round(sel_data.nf/(3/fcut));

% nf2Hz = round(nf/3*fmax);
rramp = pi*((1:f2nf)-1)./f2nf;
ctap = zeros(sel_data.nf,1);
ctaper = cos(rramp)+1;
ctap(1:f2nf) = ctaper;

vpt = 1.5:0.5:8.5;
% Find correct upper layer velocity
for ii = 1:length(vpt)
    disp(['Testing velocity ' num2str(vpt(ii))])
    drrf_sum = zeros(sel_data.nf,1);
    ftfft = zeros(sel_data.nf,1);
    for ij = 1:min([100,length(fn)])
        phase = data.(fn{ij});
        if length(mis.start)>1
            phi = mis.phi(mis.start<=phase.tt_abs&mis.end>=phase.tt_abs);
        else
            phi = mis.phi;
        end
        if isempty(phi)
            phi = mis.phi(end);
        end
        baz = phase.event.baz-phi;
        Vpi = vpt(ii); % Near surface velocity used for LQT rotation
        rayp = phase.pp./111.11;
        if Vpi*rayp > 1
            incl = 89;
        else
            incl = asind(Vpi*rayp);
        end
        % for i = 1:3
        % phase.tr(i).data = detrend(phase.tr(i).data).*tukeywin(length(phase.tr(i).data),0.005)';
        % phase.tr(i).data = buttern_high(phase.tr(i).data,2,0.1,phase.dt)';
        % end
        [phase.tr(1).data,phase.tr(2).data,phase.tr(3).data]=rotVNE2VRT(phase.tr(1).data,phase.tr(2).data,phase.tr(3).data,baz);
        [phase.tr(1).data,phase.tr(2).data,phase.tr(3).data]=rotVRT2PSvSh(phase.tr(1).data,phase.tr(2).data,phase.tr(3).data,incl);
        phase.comp = 'lqt';
        [hr,~,rstd,~] = MTC(phase,0,sel_data);
        
        % rstd = ones(size(rstd));
        sigsq = (abs(rstd).*abs(rstd));
        ftfft = ftfft+(hr./sigsq);
        drrf_sum = drrf_sum+(ones(sel_data.nf,1)./(sigsq));
        
    end
    bb = real(ifft(ctap.*(ftfft./drrf_sum),sel_data.dpad));
    nnorm = sel_data.dpts2./sel_data.nf;
    bb = nnorm*bb;

    amp(ii) = sum(abs(bb((end-20):sel_data.dpad).^2));
    % amp(ii) = sum(abs([bb((end-100):sel_data.dpad)',bb(1:50)']'.^2));
    amp4pl{ii} = [bb((end-100):sel_data.dpad)',bb(1:50)']';
end
amp0 = min(amp);
vpf = vpt(amp0==amp);
disp(vpf);

fig = figure;
plot(vpt,amp,'color','k','LineWidth',1.2)
hold on
plot([vpf vpf],[0.9.*min(amp) 1.1.*max(amp)],'Color','red','LineWidth',1.2,'LineStyle','--');
title(['Surface Velocity ' num2str(vpf) ' km/s'])
xlabel('vp [km/s]');
xlim([min(vpt) max(vpt)])
ylim([0.9*min(amp) 1.1*max(amp)])

fname=[sel_data.save_dir '/' sel_data.nc '.' sel_data.station '-SurfaceVelocity.png'];
print(fig,fname,'-dpng','-r600');
close(fig);

cols = jet(length(amp4pl));
fig = figure;
hold on
for i = 1:length(amp4pl)
    plot((-100:50).*phase.dt,amp4pl{i},'LineWidth',1.2,'Color',cols(i,:))
end
colormap(cols);
if datetime(version('-date')) >= datetime('March 9, 2022') % Matlab 2022a release
    clim([min(vpt) max(vpt)]);
else
    caxis([min(vpt) max(vpt)]);
end
    
colorbar
title(['Surface Velocity ' num2str(vpf) ' km/s'])
xlabel('time in [s]');
xlim([min((-100:50).*phase.dt) max((-100:50).*phase.dt)])

fname=[sel_data.save_dir '/' sel_data.nc '.' sel_data.station '-SurfaceVelocity_amp.png'];
print(fig,fname,'-dpng','-r600');
close(fig);

disp("Surface velocity identified") 

fid = fopen([sel_data.save_dir '/' sel_data.vpsfile],'at');
fprintf(fid,'%s,%f\n',[sel_data.nc '.' sel_data.station],vpf);
fclose(fid);

end