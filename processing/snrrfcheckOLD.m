function [snrflag,rf,tt,trf,rft] = snrrfcheckOLD(phase,freqmax,thres)

dt = phase.dt;
offset2=round(length(phase.tr(1).data)/2)*dt;
slow = phase.pp./111.11;

for i = 1:3
phase.tr(i).data = detrend(phase.tr(i).data).*tukeywin(length(phase.tr(i).data),0.001)';
phase.tr(i).data = buttern_high(phase.tr(i).data,2,0.1,phase.dt);
end
[phase.tr(2).data,phase.tr(3).data] = rad_tra(phase.tr(2).data,phase.tr(3).data, phase.event.baz);
phase.comp = 'zrt';

nlayers = 2;
ldepth = [20 35 300];
velp = [5.8 6.5 8];
vels = [3.4 3.75 4.5];
ddepth = zeros(length(ldepth),1);
tau = ddepth;
ddepth(1) = ldepth(1);
sdelay = 0;
pdelay = 0;
for ii = 1:nlayers
    ddepth(ii+1) = ldepth(ii+1)-ldepth(ii);
    sdelay = sdelay+ddepth(ii)./vels(ii);
    pdelay = pdelay+ddepth(ii)./velp(ii);
    tau(ii) = sdelay-pdelay;
end
targkm = 0;
mlay = 1;
resz = targkm;
vdelay = 0;
for i = 1:length(ldepth)
    depth = ldepth(i);
    if targkm > depth 
        vdelay = tau(mlay);
        mlay = mlay+1;
        resz = targkm-depth;
    end
end
vdelay = vdelay+(resz./vels(mlay)-resz./velp(mlay));

tdur = 80;
dpts = round(tdur/dt);
dln2 = log(tdur/dt)/log(2);
dpad = round(2.^ceil(dln2));
tpad = dpad*dt;
tprior = 20;
m1 = round((offset2-tprior)/dt)-1;
m0 = m1 - dpts;
m2 = m1 + dpts;
npre = round(20/dt);
npost = round(60/dt);
ntot = npre+npost;
tt = dt*((1:ntot)-npre);

pnw = 2.5;
nwin = 3;
fmax = 6.0;
df = 1/tpad;
nf = round(fmax/df)+1;
% win_rfs = dpss(dpts,pnw,nwin);
load('windpss_3200_25_3.mat');

hr = zeros(nf,1);
ht = hr;
trt = phase.tr(strfind(phase.comp,'t')).data;
trr = phase.tr(strfind(phase.comp,'r')).data;
trz = phase.tr(strfind(phase.comp,'z')).data;
hr3 = [hr, hr, hr];
pfftr = hr3;
pfftt = hr3;
pfftz = hr3;
fftr = hr3;
fftt = hr3;
fftz = hr3;

tshft = 0;
if mlay>1
    for k = 1:mlay
        coss = sqrt(1-vels(k)*vels(k)*slow*slow);
        cosp = sqrt(1-velp(k)*velp(k)*slow*slow);
        tshft = tshft + ddepth(k)*(coss./vels(k)-cosp./velp(k));
    end
end
coss = sqrt(1-vels(mlay)*vels(mlay)*slow*slow);
cosp = sqrt(1-velp(mlay)*velp(mlay)*slow*slow);
tshft = tshft + resz.*(coss./vels(mlay)-cosp./velp(mlay));
tmshft = tshft/dt;
mshft = round(tmshft+0.5)-1;

% eigenspectra of the pre-event noise
for i = 1:length(phase.tr)
    % data = phase.tr(i).data-mean(phase.tr(i).data);
    data = detrend(phase.tr(i).data);
    if strcmp(phase.comp(i),'t')
        data = -data;
    end
    for kk = 1:nwin
        datapad = zeros(dpad,1);
        datapad(1:dpts) = data(m0:m1-1).*win_rfs(:,kk);
        fftdum = fft(datapad);
        if strcmp(phase.comp(i),'z')
            pfftz(1:nf,kk) = fftdum(1:nf);
        elseif strcmp(phase.comp(i),'r')
            pfftr(1:nf,kk) = fftdum(1:nf);
        elseif strcmp(phase.comp(i),'t')
            pfftt(1:nf,kk) = fftdum(1:nf);
        end
    end
end
% eigenspectra of the P-Coda
for i = 1:length(phase.tr)
    data = phase.tr(i).data-mean(phase.tr(i).data);
    if strcmp(phase.comp(i),'t')
        data = -data;
    end
    m11 = m1;
    m22 = m2;
    if ~strcmp(phase.comp(i),'z')
        m11 = m11+mshft;
        m22 = m22+mshft;
    end
    for kk = 1:nwin
        datapad = zeros(dpad,1);
        datapad(1:dpts) = data(m11:m22-1).*win_rfs(:,kk);
        fftdum = fft(datapad);
        if strcmp(phase.comp(i),'z')
            fftz(1:nf,kk) = fftdum(1:nf);
        elseif strcmp(phase.comp(i),'r')
            fftr(1:nf,kk) = fftdum(1:nf);
        elseif strcmp(phase.comp(i),'t')
            fftt(1:nf,kk) = fftdum(1:nf);
        end
    end
end

pzspec = zeros(nf,1);
zspec = zeros(nf,1);
tspec = zeros(nf,1);
rspec = zeros(nf,1);
zrcor = zeros(nf,1);
ztcor = zeros(nf,1);
pfft = pfftr;
for i = 1:nwin
    pzabs = real(conj(pfft(1:nf,i)).*pfft(1:nf,i));
    zabs = real(conj(fftz(1:nf,i)).*fftz(1:nf,i));
    tabs = real(conj(fftt(1:nf,i)).*fftt(1:nf,i));
    rabs = real(conj(fftr(1:nf,i)).*fftr(1:nf,i));
    zrdot = (conj(fftz(1:nf,i)).*fftr(1:nf,i));
    ztdot = (conj(fftz(1:nf,i)).*fftt(1:nf,i));
    zrcor = zrcor + zrdot;
    ztcor = ztcor + ztdot;
    zspec = zspec + zabs;
    rspec = rspec + rabs;
    tspec = tspec + tabs;
    pzspec = pzspec + pzabs;
end
denom = pzspec + zspec;
hr = zrcor./denom;
ht = ztcor./denom;
top = real(conj(zrcor).*zrcor);
denom = zspec.*rspec;
crfr = top./denom;
top = real(conj(ztcor).*ztcor);
denom = zspec.*tspec;
crft = top./denom;
top = 1-crfr;
denom = crfr*(nwin-1);
stdevr = sqrt(top./denom).*hr;
top = 1-crft;
denom = crft*(nwin-1);
stdevt = sqrt(top./denom).*ht;

% Finalize radial RF and find SNR
rff = hr;
fmax = freqmax;
nf2Hz = round(nf/3*fmax);
rramp = pi*((1:nf2Hz)-1)./nf2Hz;
ctap = zeros(nf,1);
ctaper = cos(rramp)+1;
ctap(1:nf2Hz) = ctaper;
rff = ctap.*rff;
rftemp = real(ifft(rff,dpad));
rf = [rftemp(end-npre+1:dpad)' rftemp(1:npost)']';

% Finalize transverse RF and find SNR
rff = ht;
rff = ctap.*rff;
rftemp = real(ifft(rff,dpad));
trf = [rftemp(end-npre+1:dpad)' rftemp(1:npost)']';

% rft = sqrt(rf.^2+trf.^2);
rft = rf;

ppeak = rft(tt>-0.5&tt<0.5);
pre = rft(tt>-10&tt<-1);
p = rft(tt>1&tt<10);
ppeak_max = mean(abs(ppeak));
p_max = max(p);
pre_max = max(abs(pre));
p_max2 = mean(abs(p));
pre_max2 = mean(abs(pre));
power_ratio = 10*log10((p_max/pre_max)^2);
power_ratio2 = 10*log10((p_max2/pre_max2)^2);
% disp(num2str(power_ratio))
% disp(num2str(power_ratio2))

snrflag = 0;
if (power_ratio>thres) && (power_ratio2>thres/2)
    snrflag = 1;
end

end