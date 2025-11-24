function [hr,ht,stdevr,stdevt] = MTC2(phase,tshift,sel_data)

dt = phase.dt;
p_time=seconds(phase.tt_abs-phase.time(1));

dpts = sel_data.dpts;
dpad = sel_data.dpad;
tprior = sel_data.tprior;
m1 = round((p_time-tprior)/dt)-1;
m0 = m1 - dpts;
m2 = m1 + dpts;
pnw = 2.5;
nwin = 3;
nf = sel_data.nf;
% win_rfs = dpss(dpts,pnw,nwin);
load('windpss_3200_25_3.mat');
win2 = tukeywin(3200,0.1);

tmshft = tshift/dt;
mshft = round(tmshft+0.5);
comp = 'zrt';
if strcmp(phase.comp,'lqt')
    comp = 'lqt';
end

hr = zeros(nf,1);
ht = hr;
hr3 = [hr, hr, hr];
pfftr = hr3;
pfftt = hr3;
pfftz = hr3;
fftr = hr3;
fftt = hr3;
fftz = hr3;

% [phase.tr(2).data, phase.tr(3).data] =  rad_tra (phase.tr(2).data, phase.tr(3).data, 45);
if sel_data.srfflag
    if strcmp(phase.comp,'zrt')
        comp = 'rzt';
    elseif strcmp(phase.comp,'lqt')
        comp = 'qlt';
    end
end

% eigenspectra of the pre-event noise
for i = 1:length(phase.tr)
    % data = phase.tr(i).data-mean(phase.tr(i).data);
    data = detrend(phase.tr(i).data);
    if strcmp(phase.comp(i),comp(3))
        data = -data;
    end
    for kk = 1:1%nwin
        datapad = zeros(dpad,1);
        if sel_data.srfflag
            datapad(1:dpts) = data(1:dpts).*win_rfs(:,kk).*win2;
        else
            datapad(1:dpts) = data(m0:m1-1).*win_rfs(:,kk).*win2;
        end
        fftdum = fft(datapad);
        if strcmp(phase.comp(i),comp(1))
            pfftz(1:nf,kk) = fftdum(1:nf);
        elseif strcmp(phase.comp(i),comp(2))
            pfftr(1:nf,kk) = fftdum(1:nf);
        elseif strcmp(phase.comp(i),comp(3))
            pfftt(1:nf,kk) = fftdum(1:nf);
        end
    end
end
% eigenspectra of the P-Coda
for i = 1:length(phase.tr)
    data = phase.tr(i).data-mean(phase.tr(i).data);
    if strcmp(phase.comp(i),comp(3))
        data = -data;
    end
    m11 = m1;
    m22 = m2;
    if ~strcmp(phase.comp(i),comp(1))
        m11 = m11+mshft;
        m22 = m22+mshft;
    end
    for kk = 1:1%nwin
        datapad = zeros(dpad,1);
        datapad(1:dpts) = data(m11:m22-1).*win_rfs(:,kk).*win2;
        fftdum = fft(datapad);
        if strcmp(phase.comp(i),comp(1))
            fftz(1:nf,kk) = fftdum(1:nf);
        elseif strcmp(phase.comp(i),comp(2))
            fftr(1:nf,kk) = fftdum(1:nf);
        elseif strcmp(phase.comp(i),comp(3))
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
for i = 1:1%nwin
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
denom = crfr;%*(nwin-1);
stdevr = sqrt(top./denom).*abs(hr);
top = 1-crft;
denom = crft;%*(nwin-1);
stdevt = sqrt(top./denom).*abs(ht);

% [hr, ht] =  rad_tra (hr, ht, -45);
% [stdevr, stdevt] =  rad_tra (stdevr, stdevt, -45);
% ht = -ht;


end