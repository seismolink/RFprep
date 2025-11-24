function [snrflag,power_ratio,power_ratio2,rf,tt,test] = snrsrfcheck(phase,freqmax,thres)

dt = phase.dt;
offset2 = seconds(phase.tt_abs-phase.time(1));

tdur = 160;
dpts = round(tdur/dt);
dln2 = log(tdur/dt)/log(2);
dpad = round(2.^ceil(dln2));
tpad = dpad*dt;
tprior = 20;
m1 = round((offset2-tprior)/dt)-1;
% m0 = m1 - dpts;
m2 = m1 + dpts;
npre = round(20/dt);
npost = round(140/dt);
ntot = npre+npost;
tt = - flip(dt*((1:ntot)-npre));

fmax = 6.0;
df = 1/tpad;
nf = round(fmax/df)+1;
nwin = 3;
% win_rfs = dpss(dpts,pnw,nwin);
load('windpss_3200_25_3.mat');
win2 = tukeywin(3200,0.1);


hz = zeros(nf,1);
hz3 = [hz, hz, hz];
pfftr = hz3;
pfftz = hz3;
fftr = hz3;
fftz = hz3;

% eigenspectra of the pre-event noise
for i = 1:length(phase.tr)
    % data = phase.tr(i).data-mean(phase.tr(i).data);
    data = detrend(phase.tr(i).data);
    for kk = 1:nwin
        datapad = zeros(dpad,1);
        datapad(1:dpts) = data(1:dpts).*win_rfs(:,kk).*win2;
        fftdum = fft(datapad);
        if strcmp(phase.comp(i),'z')
            pfftz(1:nf,kk) = fftdum(1:nf);
        elseif strcmp(phase.comp(i),'r')
            pfftr(1:nf,kk) = fftdum(1:nf);
        end
    end
end
% eigenspectra of the S-Coda
for i = 1:length(phase.tr)
    data = phase.tr(i).data-mean(phase.tr(i).data);
    m11 = m1;
    m22 = m2;
    for kk = 1:nwin
        datapad = zeros(dpad,1);
        datapad(1:dpts) = data(m11:m22-1).*win_rfs(:,kk).*win2;
        fftdum = fft(datapad);
        if strcmp(phase.comp(i),'z')
            fftz(1:nf,kk) = fftdum(1:nf);
            if kk == 1
            dataz = data(m11:m22-1).*win_rfs(:,kk).*win2;
            dataza = data;
            end
        elseif strcmp(phase.comp(i),'r')
            fftr(1:nf,kk) = fftdum(1:nf);
            if kk == 1
            datar = data(m11:m22-1).*win_rfs(:,kk).*win2;
            datara = data;
            end
        elseif strcmp(phase.comp(i),'t')
            if kk == 1
            datat = data(m11:m22-1).*win_rfs(:,kk).*win2;
            datata = data;
            end
        end
    end
end

prspec = zeros(nf,1);
% zspec = zeros(nf,1);
zrcor = zeros(nf,1);
% ztcor = zeros(nf,1);
pfft = pfftr;
rspec = zeros(nf,1);
for i = 1:nwin
    prabs = real(conj(pfft(1:nf,i)).*pfft(1:nf,i));
    % zabs = real(conj(fftz(1:nf,i)).*fftz(1:nf,i));
    rabs = real(conj(fftr(1:nf,i)).*fftr(1:nf,i));
    zrdot = (conj(fftr(1:nf,i)).*fftz(1:nf,i));
    % ztdot = (conj(fftz(1:nf,i)).*fftt(1:nf,i));
    zrcor = zrcor + zrdot;
    % ztcor = ztcor + ztdot;
    % zspec = zspec + zabs;
    prspec = prspec + prabs;
    rspec = rspec + rabs;
end
denom = prspec + rspec;
% denom = pzspec + zspec;
hz = zrcor./denom;
% ht = ztcor./denom;

% Finalize SRF and find SNR
rff = hz;
fmax = freqmax;
nf2Hz = round(nf/3*fmax);
rramp = pi*((1:nf2Hz)-1)./nf2Hz;
ctap = zeros(nf,1);
ctaper = cos(rramp)+1;
ctap(1:nf2Hz) = ctaper;
rff = ctap.*rff;
rftemp = real(ifft(rff,dpad));
rf = flip([rftemp(end-npre+1:dpad)' rftemp(1:npost)']');

% % Finalize transverse RF and find SNR
% rff = ht;
% rff = ctap.*rff;
% rftemp = real(ifft(rff,dpad));
% trf = [rftemp(end-npre+1:dpad)' rftemp(1:npost)']';
% 
% rft = sqrt(rf.^2+trf.^2);
% rft = rf;

rft = rf;

ab = movmean(abs(rft),[round(0.75/dt) round(0.75/dt)]);
cd = movmean(abs(rft),[round((12-0.75)/dt) round(0.75/dt)]);
test = ab./cd;

ppeak = rft(tt>-0.5&tt<0.5);
pre = rft(tt>-10&tt<-2);
p = rft(tt>1&tt<10);
ppeak_max = mean(abs(ppeak));
p_max = max(abs(p));
pre_max = max(abs(pre));
p_max2 = mean(abs(p));
pre_max2 = mean(abs(pre));
power_ratio = 10*log10((p_max/pre_max)^2);
power_ratio2 = 10*log10((p_max2/pre_max2)^2);
% disp(num2str(power_ratio))
% disp(num2str(power_ratio2))

snrflag = 0;
if ((power_ratio>thres) || (power_ratio2>thres)) %&& test(abs(tt)==min(abs(tt))) > 2.1
    %((power_ratio>thres) || (power_ratio2>thres))
    snrflag = 1;
end

end