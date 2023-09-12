function u1 = calcdisp_antiplanefault(x2,x3,y2,y3,W,dip)
% calculate displacements for a dipping anti-plane finite fault 
% x2,x3 - observation points
% y2,y3 - source location
% W,dip - width and dip of source patch
% Rishav Mallick, 2023, Caltech Seismolab

u1 = (atan2((x3-y3),(x2-y2)) - atan2((x3-y3-W*sind(dip)),(x2-y2-W*cosd(dip))))/2/pi;
index = abs(u1) >= 0.5;
u1(index) = u1(index) - sign(u1(index)).*1;

end