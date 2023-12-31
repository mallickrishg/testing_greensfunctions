function [s12,s13] = calcstress_antiplanefault(G,x2,x3,y2,y3,Wf,dip)

s12= G*( ...
    -(x3-y3)./((x2-y2).^2+(x3-y3).^2) ...
    +(x3-y3-Wf*sind(dip))./((x2-y2-Wf*cosd(dip)).^2+(x3-y3-Wf*sind(dip)).^2) ...
    )/2/pi;

s13= G*( ...
    (x2-y2)./((x2-y2).^2+(x3-y3).^2) ...
    -(x2-y2-Wf*cosd(dip))./((x2-y2-Wf*cosd(dip)).^2+(x3-y3-Wf*sind(dip)).^2) ...
    )/2/pi;

end