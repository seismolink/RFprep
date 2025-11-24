function [P,Sv,Sh]=rotVRT2PSvSh(V,R,T,incident)
% [P,Sv,Sh]=rotVLT2PSvSh(V,L,T,i)
% Rotation from V-L-T to P-Sv-Sh
% The signals are in format SAC
%  incident: incidence angle (rad) of the arriving wave with respect to

M=[cosd(incident)   sind(incident)
    -sind(incident)  cosd(incident)
    ];

S=M*[V';R'];

P=S(1,:)';
Sv=S(2,:)';
Sh = T;
tr=P;tr=detrend(tr);tr=tr-mean(tr);P=tr;
tr=Sv;tr=detrend(tr);tr=tr-mean(tr);Sv=tr;
tr=Sh;tr=detrend(tr);tr=tr-mean(tr);Sh=tr;

end

