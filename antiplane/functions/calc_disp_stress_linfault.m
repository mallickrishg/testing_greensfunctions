function [u,s12,s13] = calc_disp_stress_linfault(x,y,w,a0,a1)

u = [a0/2/pi.*atan2(2*w*x,x.^2 + y.^2 - w.^2),... 
    a1.*(y./2/pi.*atan2(2*w.*x,x.^2 + y.^2 - w.^2) + x./4/pi.*log((x.^2 + (y-w).^2)./(x.^2 + (y+w).^2))) ];

end