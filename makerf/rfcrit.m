function [rfflag,r1,r2,r3,r4,r5] = rfcrit(temp,temp2)
% temp = radial
% temp2 = transverse
rlim = 1.5;
rlim2 = 1;
r1 = max(abs(temp(40:200)))./max(abs(temp(end-44:end-5)));
r2 = max(abs(temp(1:100)))./max(abs(temp2(1:100)));
r3 = max(abs(temp2(40:200)))./max(abs(temp2(end-44:end-5)));
% r1 = sum(abs(temp(1:20)))/sum(abs(temp(end-29:end-10)));
r2 = sum(abs(temp(1:100)))/sum(abs(temp(end-109:end-10)));
r3 = sum(abs(temp2(1:100)))/sum(abs(temp2(end-99:end)));
r4 = max(abs(temp2(1:300)))/max(abs(temp2(end-29:end-1)));
r5 = sum(abs(temp(1:100)))/sum(abs(temp2(1:100)));
% if r4 < rlim || r2 < rlim || r1 < rlim || r3 < rlim || r5 < rlim
if r1 < rlim || r2 < rlim2 || r3 < rlim2
    rfflag = 0;
else
    rfflag = 1;
end
end