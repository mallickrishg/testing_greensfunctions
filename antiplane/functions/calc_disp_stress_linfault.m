function [u,s12,s13] = calc_disp_stress_linfault(xo,yo,xloc,yloc,w,alpha)
% compute displacement and stress in response to a spatially variable 
% fault source with slip describes as s(ζ) = α_0 + α_1*ζ
% as [N x 2] matrices
% INPUTS:
% xo,yo - observation points: each variable as a (N x 1) vector
% xloc,yloc - location of fault center
% w - fault half-width
% alpha - source terms, where alpha(1) = alpha_0, alpha(2) = alpha_1
% OUTPUTS:
% u1 - displacement kernel as [u(α_0) , u(α_1)] matrix
% s12,s13 - stress kernels as [N x 2] matrices each
% 
% Rishav Mallick, 2023, Caltech Seismolab

% TODO - need to add rotation for dip
x = xo-xloc;
y = yo-yloc;
a0 = alpha(1);
a1 = alpha(2);

u = [a0/2/pi.*atan2(2*w*x,x.^2 + y.^2 - w.^2),... 
    a1.*(y./2/pi.*atan2(2*w.*x,x.^2 + y.^2 - w.^2) + x./4/pi.*log((x.^2 + (y-w).^2)./(x.^2 + (y+w).^2))) ];

s12 = [-a0*w/pi.*(x.^2+w.^2-y.^2)./(w^4 + 2*w^2.*(x.^2-y.^2)+(x.^2+y.^2).^2),...
    a1.*(y./pi.*(x.^2+y.^2-w.^2)./(w^4 + 2*w^2.*(x.^2-y.^2)+(x.^2+y.^2).^2)+1/4/pi.*log((x.^2 + (y-w).^2)./(x.^2 + (y+w).^2)))];

s13 = [-2*a0.*w.*x.*y./pi./(w^4 + 2*w^2.*(x.^2-y.^2)+(x.^2+y.^2).^2),...
    a1.*(-w.*x./pi.*(w^2+x.^2+y.^2)./(w^4 + 2*w^2.*(x.^2-y.^2)+(x.^2+y.^2).^2) +1/2/pi.*atan2(2*w.*x,x.^2 + y.^2 - w.^2))];

end