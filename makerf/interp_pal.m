function Cout=interp_pal(Cin,n)
% Cout=interp_pal(Cin,n)
% interpolation d'une palette de couleur Cin
% a un facteur n

xCin=[0:size(Cin,1)-1]*n;
xCout=[0:max(xCin)];

Cout1=interp1(xCin,Cin(:,1),xCout);
Cout2=interp1(xCin,Cin(:,2),xCout);
Cout3=interp1(xCin,Cin(:,3),xCout);
Cout=[Cout1;Cout2;Cout3];
Cout=Cout';
