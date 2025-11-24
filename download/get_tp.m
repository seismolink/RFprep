function [phases,tt,pp,nom]=get_tp(eq_dist,eq_depth,P,Pdiff,PKIKP,PKS,PKIKS,S,Sdiff,ScS,SKS,SKKS,SKIKS,PcS,pPcS,PS,PPS)

% table look up of travel times
% usage:
% eq_dist: distance earthquake - station
% eq_depth: event depth

% Copyright 2016 M.Reiss and G.Rümpker
% altered 2024 F.Link, J.Wolf and M.Reiss

% % load travel time files
% P = load('P.mat');
% Pdiff = load('Pdiff.mat');
% PKIKP = load('PKIKP.mat');
% PKS = load('PKS.mat');
% PKIKS = load('PKIKS.mat');
% S = load('S.mat');
% Sdiff = load('Sdiff.mat');
% ScS = load('ScS.mat');
% SKS = load('SKS.mat');
% SKKS = load('SKKS.mat');
% SKIKS = load('SKIKS.mat');
% PcS = load('PcS.mat');
% pPcS = load('pPcS.mat');
% PS = load('PS.mat');
% PPS = load('PPS.mat');

% table look up with interpolation
[t.P,p.P]=interp_tp(P,eq_dist,eq_depth);
[t.Pdiff,p.Pdiff]=interp_tp(Pdiff,eq_dist,eq_depth);
[t.PKIKP,p.PKIKP]=interp_tp(PKIKP,eq_dist,eq_depth);
[t.PKS,p.PKS]=interp_tp(PKS,eq_dist,eq_depth);
[t.PKIKS,p.PKIKS]=interp_tp(PKIKS,eq_dist,eq_depth);

[t.S,p.S]=interp_tp(S,eq_dist,eq_depth);
[t.Sdiff,p.Sdiff]=interp_tp(Sdiff,eq_dist,eq_depth);
[t.SKS,p.SKS]=interp_tp(SKS,eq_dist,eq_depth);
[t.SKKS,p.SKKS]=interp_tp(SKKS,eq_dist,eq_depth);
[t.SKIKS,p.SKIKS]=interp_tp(SKIKS,eq_dist,eq_depth);
[t.ScS,p.ScS]=interp_tp(ScS,eq_dist,eq_depth);

[t.PcS,p.PcS]=interp_tp(PcS,eq_dist,eq_depth);
[t.pPcS,p.pPcS]=interp_tp(pPcS,eq_dist,eq_depth);

[t.PS,p.PS]=interp_tp(PS,eq_dist,eq_depth);
[t.PPS,p.PPS]=interp_tp(PPS,eq_dist,eq_depth);

fn = fieldnames(t);

% remove non existent phases
for iF = 1:length(fn)
    af = fn{iF};
    if isnan(t.(af)) || isnan(p.(af))
        t = rmfield(t,af);
        p = rmfield(p,af);
    end
end


fn2 = fieldnames(t);
% rewrite fields
for iF2 = 1:length(fn2)
    af = fn2{iF2};
    phases(iF2).name = af;
    phases(iF2).tt = t.(af);
    phases(iF2).pp = p.(af);
    tt(iF2) = t.(af);
    pp(iF2) = p.(af);
    nom{iF2} = af;
end 

if exist('phases','var') == 0
    phases = 0;
end

end