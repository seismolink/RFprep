function [V2,R,T]=rotVNE2VRT(V1,N,E,baz)

M=[1         0         0;
    0 -cosd(baz) -sind(baz);
    0  sind(baz) -cosd(baz)];
% M=[1         0         0;
%     0 cosd(baz) -sind(baz);
%     0  sind(baz) cosd(baz)];


S=M*[V1;N;E];

V2=S(1,:)';
R=S(2,:)';
T=S(3,:)';
tr=V2;tr=detrend(tr);tr=tr-mean(tr);V2=tr;
tr=R;tr=detrend(tr);tr=tr-mean(tr);R=tr;
tr=T;tr=detrend(tr);tr=tr-mean(tr);T=tr;

end