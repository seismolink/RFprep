function misorient = findmisorientation(sel_data,data,fn)

% identify misorientation from P wave polarization in harmonics

z_targets = [0]; % Targeting Depth
% 1-D velocity model for mapping to depth
ldepth = [20, 35, 300]; % depths of interfaces
velp = [5.8, 6.5, 8]; % P wave velocity
vels = [3.4, 3.75, 4.5]; % S wave velocity

% Define the grid of depth output
z_out = -25:0.25:150;
% fcutz = 1; % ["0.5","1.0","1.5","2.0","3.0"]
% % nfz = round(sel_data.nf/(3/(fcutz))); % [int(nf/6), int(nf/3), int(nf/2), int(nf/1.5), int(nf/1)]
% fcut = fcutz;
fcut = sel_data.fcut;

% f2nf = round(sel_data.nf/((sel_data.fmax/2)/fcut));
f2nf = round(sel_data.nf/((sel_data.fmax/4)/fcut));

% nf2Hz = round(nf/3*fmax);
rramp = pi*((1:f2nf)-1)./f2nf;
ctap = zeros(sel_data.nf,1);
ctaper = cos(rramp)+1;
ctap(1:f2nf) = ctaper;

nrf = 1;
for iev = 1:length(fn)
    phase = data.(fn{iev});
    slow = phase.pp./111.11;
    if phase.event.dist<25 || phase.event.dist>95
        continue
    end
    % for i = 1:3
    % phase.tr(i).data = detrend(phase.tr(i).data).*tukeywin(length(phase.tr(i).data),0.005)';
    % phase.tr(i).data = buttern_high(phase.tr(i).data,4,0.2,phase.dt)';
    % end
    [phase.tr(1).data,phase.tr(2).data,phase.tr(3).data]=rotVNE2VRT(phase.tr(1).data,phase.tr(2).data,phase.tr(3).data,phase.event.baz);
    phase.comp = 'zrt';

    for j = 1:length(z_targets)
        z_target = z_targets(j);
        tshift = z2t(z_target,ldepth,velp,vels,slow);
        [hr,ht,stdevr,stdevt] = MTC2(phase,tshift,sel_data);
        % temp = real(ifft(ctap.*hr,sel_data.dpad));
        % temp2 = real(ifft(ctap.*ht,sel_data.dpad));
        % rfflag = rfcrit(temp,temp2);
        % if ~rfflag
        %     disp('problematic trace')
        %     keyboard
        %     continue
        % end
        rf.(fn{iev}).rrf(j).data = hr;
        rf.(fn{iev}).trf(j).data = ht;
        rf.(fn{iev}).rstd(j).data = stdevr;
        rf.(fn{iev}).tstd(j).data = stdevt;
        rf.(fn{iev}).tshift(j) = tshift;
        rf.(fn{iev}).zt(j) = z_target;
    end
    % if ~rfflag
    %     continue
    % end
    rf.(fn{iev}).event = data.(fn{iev}).event;
    rf.(fn{iev}).tt_abs = data.(fn{iev}).tt_abs;
    rf.(fn{iev}).pp = data.(fn{iev}).pp;
    fn2{nrf} = fn{iev};
    evtime(nrf) = data.(fn{iev}).tt_abs;
    nrf = nrf+1;
end
disp("Receiver Functions calculated for preprocessing")

% Space domain harmonic regression
[bma_list,stt_list,edt_list] = hm_decon_depth_mis(rf,fn2,evtime,sel_data,z_out,ldepth,velp,vels,fcut,z_targets);

disp('Search misorientation')
% Search for misorientation
misx = 0:0.1:180;
PHI = zeros(length(misx),length(bma_list));
phi1 = zeros(length(bma_list),1);
phit = zeros(length(bma_list),1);

for ii = 1:length(bma_list)
    bootma = bma_list{ii};
    % fr = sqrt(sum((-sind(misx)'*bootma(1,z_out>-2&z_out<2)+cosd(misx)'*bootma(6,z_out>-2&z_out<2)).^2,2));
    fr = sqrt(sum((-sind(misx)'*bootma(1,z_out>-5&z_out<5)+cosd(misx)'*bootma(6,z_out>-5&z_out<5)).^2,2));
    phi0 = misx(fr==min(fr));
    phi0 = phi0(1);
    PHI(:,ii) = fr;
    % aa = min(fr);
    % hmr = cosd(phi0)*bootma(1,:) + sind(phi0)*bootma(6,:); 
    % hmr2 = cosd(phi0+180)*bootma(1,:) + sind(phi0+180)*bootma(6,:); 
    if phi0 > 90
        det1 = -1;
        det2 = 1;
    else
        det1 = 1;
        det2 = -1;
    end
    if det1 > 0
        phi = round(phi0,1);
        Tflag(ii) = false;
    end
    if det2 > 0
        phi = round(phi0-180,1);
        Tflag(ii) = true;
    end
    phi1(ii) = abs(phi);
    phit(ii) = phi;
end

startdate = stt_list(1);
enddate = edt_list(end);
% store the dates between two dates in a list
dates = startdate:days(1):enddate;
for ii = 1:length(dates)
    a0 = inf;
    for j = 1:length(stt_list)
        if dates(ii) >= stt_list(j) && dates(ii) <= edt_list(j)
            if a0 > min(PHI(:,j))
                a0 = min(PHI(:,j));
                column = PHI(:,j);
                phioi = phi1(j);
                phitoi = phit(j);
                Tflagoi = Tflag(j);
            end
        end
    end
    if isinf(a0)
        column = PHI(:,end);
        phioi = phi1(end);
        phitoi = phit(end);
        Tflagoi = Tflag(end);
    end
    PHI2(:,ii) = column;
    phi2(ii) = phioi;
    phit2(ii) = phitoi;
    Tflag2(ii) = Tflagoi;
end

dphlim = 10;
dph = abs(diff(phit2));
istep = find(dph>dphlim);
if isempty(istep)
    phires = mean(phit2);
    stdt = dates(1);
    eddt = dates(end);
else
    phires(1) = mean(phit2(1:istep(1)));
    stdt(1) = dates(1);
    eddt(1) = dates(istep(1));
    for ii = 1:length(istep)-1
        phires(ii+1) = mean(phit2(istep(ii):istep(ii+1)));
        stdt(ii+1) = dates(istep(ii));
        eddt(ii+1) = dates(istep(ii+1));
    end
    phires(end+1) = mean(phit2(istep(end):end));
    stdt(end+1) = dates(istep(end));
    eddt(end+1) = dates(end);
end

% % Calculate ticks for plots later
% z = linspace(10,160,10);
% t_ticks = z2t(z,ldepth,velp,vels,0);

misorient.start = stdt;
misorient.end = eddt;
misorient.phi = phires;

fid = fopen([sel_data.save_dir '/' sel_data.misfile],'at');
for i = 1:length(misorient.start)
    fprintf(fid,'%s,%s,%s,%f\n',[sel_data.nc '.' sel_data.station],char(misorient.start(i),'yyyy-MM-dd''T''hh:mm:ss.SSSS'),char(misorient.end(i),'yyyy-MM-dd''T''hh:mm:ss.SSSS'),misorient.phi(i));
end
fclose(fid);

fig = figure;
imagesc(dates,misx,PHI2);
set(gca,'YDir','normal');
hold on
scatter(dates(Tflag2),phi2(Tflag2),0.5,'r')
scatter(dates(~Tflag2),phi2(~Tflag2),0.5,'b')
for i = 1:length(misorient.start)
    plot([misorient.start(i) misorient.end(i)],[misorient.phi(i) misorient.phi(i)],'w--','LineWidth',1.2)
    text(misorient.start(i)+(misorient.end(i)-misorient.start(i))/2,abs(misorient.phi(i))+5,[num2str(round(misorient.phi(i),1)) '\circ'],'Color','white','FontWeight','bold','FontSize',9)
end
ylabel('Angle [deg]')
xlabel('Time')
ylim([0 180])
xlim([dates(1) dates(end)])

fname=[sel_data.save_dir '/' sel_data.nc '.' sel_data.station '-Misorientation.png'];
print(fig,fname,'-dpng','-r600');
disp("Found misorientation")
close(fig);

end