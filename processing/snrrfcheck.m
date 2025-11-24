function [snrflag,power_ratio,power_ratio2,rf,tt,trf,rft,test] = snrrfcheck(phase,freqmax,thres)

dt = phase.dt;
offset2 = seconds(phase.tt_abs-phase.time(1));

for i = 1:3
% phase.tr(i).data = detrend(phase.tr(i).data).*tukeywin(length(phase.tr(i).data),0.001)';
% phase.tr(i).data = buttern_high(phase.tr(i).data,2,0.1,phase.dt);
phase.tr(i).data = phase.tr(i).data';
end
% [phase.tr(2).data,phase.tr(3).data] = rad_tra(phase.tr(2).data,phase.tr(3).data, phase.event.baz);
[phase.tr(1).data,phase.tr(2).data,phase.tr(3).data]=rotVNE2VRT(phase.tr(1).data',phase.tr(2).data',phase.tr(3).data',phase.event.baz);
phase.comp = 'zrt';
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

fmax = 6.0;
df = 1/tpad;
nf = round(fmax/df)+1;
nwin = 3;
% win_rfs = dpss(dpts,pnw,nwin);
load('windpss_3200_25_3.mat');
win2 = tukeywin(3200,0.1);


hr = zeros(nf,1);
ht = hr;
hr3 = [hr, hr, hr];
pfftr = hr3;
pfftt = hr3;
pfftz = hr3;
fftr = hr3;
fftt = hr3;
fftz = hr3;

% eigenspectra of the pre-event noise
for i = 1:length(phase.tr)
    % data = phase.tr(i).data-mean(phase.tr(i).data);
    data = detrend(phase.tr(i).data);
    if strcmp(phase.comp(i),'t')
        data = -data;
    end
    for kk = 1:nwin
        datapad = zeros(dpad,1);
        datapad(1:dpts) = data(m0:m1-1).*win_rfs(:,kk);%.*win2;
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
    % data = phase.tr(i).data-mean(phase.tr(i).data);
    data = detrend(phase.tr(i).data);
    if strcmp(phase.comp(i),'t')
        data = -data;
    end
    m11 = m1;
    m22 = m2;
    for kk = 1:nwin
        datapad = zeros(dpad,1);
        datapad(1:dpts) = data(m11:m22-1).*win_rfs(:,kk);%.*win2;
        fftdum = fft(datapad);
        if strcmp(phase.comp(i),'z')
            fftz(1:nf,kk) = fftdum(1:nf);
            % if kk == 1
            % dataz = datapad(1:dpts);
            % end
        elseif strcmp(phase.comp(i),'r')
            fftr(1:nf,kk) = fftdum(1:nf);
            % if kk == 1
            % datar = datapad(1:dpts);
            % end
        elseif strcmp(phase.comp(i),'t')
            fftt(1:nf,kk) = fftdum(1:nf);
            % if kk == 1
            % datat = datapad(1:dpts);
            % end
        end
    end
end

pzspec = zeros(nf,1);
zspec = zeros(nf,1);
zrcor = zeros(nf,1);
ztcor = zeros(nf,1);
pfft = pfftr;
for i = 1:nwin
    pzabs = real(conj(pfft(1:nf,i)).*pfft(1:nf,i));
    zabs = real(conj(fftz(1:nf,i)).*fftz(1:nf,i));
    zrdot = (conj(fftz(1:nf,i)).*fftr(1:nf,i));
    ztdot = (conj(fftz(1:nf,i)).*fftt(1:nf,i));
    zrcor = zrcor + zrdot;
    ztcor = ztcor + ztdot;
    zspec = zspec + zabs;
    pzspec = pzspec + pzabs;
end
denom = pzspec + zspec;
hr = zrcor./denom;
ht = ztcor./denom;

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

ab = movmean(abs(rft),[round(0.75/dt) round(0.75/dt)]);
cd = movmean(abs(rft),[round((12-0.75)/dt) round(0.75/dt)]);
test = ab./cd;

% ppeak = rft(tt>-0.5&tt<0.5);
pre2 = rft(tt>-5&tt<-1);
p2 = rft(tt>1&tt<5);
% ppeak_max = mean(abs(ppeak));
% p_max = max(p);
% pre_max = max(abs(pre));
% p_max2 = mean(abs(p));
% pre_max2 = mean(abs(pre));
% power_ratio = 10*log10((p_max/pre_max)^2);
% power_ratio2 = 10*log10((p_max2/pre_max2)^2);
% ppeak = rft(tt>-1&tt<1);
pre = rft(tt>-10&tt<-0.5);
p = rft(tt>-0.1&tt<1);
% ppeak_max = mean(abs(ppeak));
p_max = max(p);
pre_max = max(abs(pre));
p_max2 = mean(abs(p2));
pre_max2 = mean(abs(pre2));
power_ratio = 10*log10((p_max/pre_max)^2);
% power_ratio2 = power_ratio;
power_ratio2 = 10*log10((p_max2/pre_max2)^2);
% disp(num2str(power_ratio))
% disp(num2str(power_ratio2))

snrflag = 0;
if ((power_ratio>thres) || (power_ratio2>(thres*2/3))) && test(abs(tt)==min(abs(tt))) > 2.1
    %((power_ratio>thres) || (power_ratio2>thres))
    snrflag = 1;
end

end