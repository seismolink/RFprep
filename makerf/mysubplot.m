function [object] = mysubplot(H,posx,posy,width,hight)
pos = get(H,'Position');
un = get(H,'units');
if strcmp(un,'pixels')
    posxp = pos(3)*posx+pos(1);
    posyp = pos(4)*posy+pos(2);
    wp = pos(3)*width;
    hp = pos(4)*hight;
    dwp = wp/10;
    dhp = hp/10;
    dd = min(dwp,dhp);
    dw = dd/pos(3);
    dh = dd/pos(4);
else
    dw = width/10;
    dh = hight/10;
end
set(0, 'CurrentFigure', H)
a=subplot('Position',[posx+dw posy+dh width-3.5*dw hight-3.5*dh]);
% set(a,'YTick',[],'XTick',[]);
set(a,'XGrid','on','YGrid','on');
% text(0.5,0.5,'this is no. 1')
hold on
b=subplot('Position',[posx+dw posy+hight-2.5*dh width-3.5*dw 1.5*dh]);
set(b,'XTickLabel',[],'YTickLabel',[]);
set(b,'XGrid','on');%'YGrid','on');
% text(0.5,0.5,'this is no. 2')
hold on;
c=subplot('Position',[posx+width-2.5*dw posy+dh 1.5*dw hight-3.5*dh]);
set(c,'YTickLabel',[],'XTickLabel',[]);
set(c,'YGrid','on');%,'YGrid','on');
% text(0.5,0.5,'this is no. 3')

object.a = a;
object.b = b;
object.c = c;