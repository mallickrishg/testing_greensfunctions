function [u,s12,s13] = calc_disp_stress_linfault(x,y,w,a0,a1)

% TODO - need to add rotation for dip

u = [a0/2/pi.*atan2(2*w*x,x.^2 + y.^2 - w.^2),... 
    a1.*(y./2/pi.*atan2(2*w.*x,x.^2 + y.^2 - w.^2) + x./4/pi.*log((x.^2 + (y-w).^2)./(x.^2 + (y+w).^2))) ];

s12 = [-a0*w/pi.*(x.^2+w.^2-y.^2)./(w^4 + 2*w^2.*(x.^2-y.^2)+(x.^2+y.^2).^2),...
    a1.*(y./pi.*(x.^2+y.^2-w.^2)./(w^4 + 2*w^2.*(x.^2-y.^2)+(x.^2+y.^2).^2)+1/4/pi.*log((x.^2 + (y-w).^2)./(x.^2 + (y+w).^2)))];

s13 = [-2*a0.*w.*x.*y./pi./(w^4 + 2*w^2.*(x.^2-y.^2)+(x.^2+y.^2).^2),...
    a1.*(-w.*x./pi.*(w^2+x.^2+y.^2)./(w^4 + 2*w^2.*(x.^2-y.^2)+(x.^2+y.^2).^2) +1/2/pi.*atan2(2*w.*x,x.^2 + y.^2 - w.^2))];

end