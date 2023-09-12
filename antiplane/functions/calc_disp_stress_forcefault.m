function [u,s12,s13] = calc_disp_stress_forcefault(X,Y,a,dip)
% the elementary solutions are for a line centered at (0,0) 
% extending from [-a <= y <= a]

beta = dip-90;

% first rotate from [X,Y] to [x,y] by dip angle
R = [cosd(beta),-sind(beta);...
     sind(beta), cosd(beta)];
rot_coords = [X,Y]*R;
x = rot_coords(:,1);
y = rot_coords(:,2);

u = (-4*a + 2.*x.*(atan2(2*a.*x,x.^2+y.^2-a^2)) +...
     (a - y).*log(x.^2 + (a - y).^2) + ...
     (a + y).*log(x.^2 + (a + y).^2))./(4*pi);

s12_o = (atan2(2*a.*x,x.^2+y.^2-a^2))./(2*pi);

s13_o = (log(x.^2 + (a + y).^2) - log(x.^2 + (a - y).^2))./(4*pi);

stress = [s12_o,s13_o]*R';
s12 = stress(:,1);
s13 = stress(:,2);

end