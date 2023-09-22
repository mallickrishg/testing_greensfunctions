function [u,s12,s13] = calc_disp_stress_quadfault(xo,yo,xloc,yloc,w,alpha,dip)
% compute displacement and stress in response to a spatially variable 
% fault source with slip describes as s(ζ) = α_0 + α_1*ζ + α_2*ζ^2
% as [N x 3] matrices
% INPUTS:
% xo,yo - observation points: each variable as a (N x 1) vector
% xloc,yloc - location of fault center
% w - fault half-width
% alpha - source terms, where alpha = [α_0,α_1,α_2]
% dip - fault dip angle (in degrees)
% OUTPUTS:
% u1 - displacement kernel as [u(α_0), u(α_1), u(α_2)] matrix
% s12,s13 - stress kernels as [N x 3] matrices each
% 
% Rishav Mallick, 2023, Caltech Seismolab

x = xo-xloc;
y = yo-yloc;
a0 = alpha(1);
a1 = alpha(2);
a2 = alpha(3);

beta = 90-dip;

% first rotate from [X,Y] to [x,y] by dip angle
R = [cosd(beta),-sind(beta);...
     sind(beta), cosd(beta)];
rot_coords = [x,y]*R;
x = rot_coords(:,1);
y = rot_coords(:,2);

% displacement
u_0 = 1/2/pi.*atan2(2*w*x,x.^2 + y.^2 - w.^2);
u_1 = (y./(2*pi).*atan2(2*w.*x,x.^2 + y.^2 - w.^2) +...
    x./4/pi.*log((x.^2 + (y-w).^2)./(x.^2 + (y+w).^2)));
u_2 = 1/2/pi.*(2*w.*x + ...
    (y.^2-x.^2).*atan2(2*w.*x,x.^2 + y.^2 - w.^2) +... 
    x.*y.*log((x.^2 + (y-w).^2)./(x.^2 + (y+w).^2)));

u = [a0.*u_0, a1.*u_1, a2.*u_2];

% stress components
sx1 = -w/pi.*(x.^2+w.^2-y.^2)./(w^4 + 2*w^2.*(x.^2-y.^2)+(x.^2+y.^2).^2);
sx2 = (w.*y./pi.*(x.^2+y.^2-w.^2)./(w^4 + 2*w^2.*(x.^2-y.^2)+(x.^2+y.^2).^2)+...
    1/4/pi.*log((x.^2 + (y-w).^2)./(x.^2 + (y+w).^2)));
sx3 = 1/2/pi.*(w.*(4 + w.*(y-w)./(x.^2 + (w-y).^2) - w.*(y+w)./(x.^2 + (w+y).^2))...
    -2*x.*atan2(2*w.*x,x.^2 + y.^2 - w.^2) +... 
    y.*log((x.^2 + (y-w).^2)./(x.^2 + (y+w).^2)));

s12 = [a0.*sx1, a1.*sx2, a2.*sx3];

sy1 = -2*w.*x.*y./pi./(w^4 + 2*w^2.*(x.^2-y.^2)+(x.^2+y.^2).^2);
sy2 = (-w.*x./pi.*(w^2+x.^2+y.^2)./(w^4 + 2*w^2.*(x.^2-y.^2)+(x.^2+y.^2).^2) +...
    1/2/pi.*atan2(2*w.*x,x.^2 + y.^2 - w.^2));
sy3 = 1/(2*pi).*(2.*y.*atan2(2*w.*x,x.^2 + y.^2 - w.^2) +...
    -4*(w^3).*x.*y./(w^4 + 2*w^2.*(x.^2-y.^2)+(x.^2+y.^2).^2)+...
    x.*log((x.^2 + (y-w).^2)./(x.^2 + (y+w).^2)));

s13 = [a0.*sy1, a1.*sy2, a2.*sy3];

% rotate stress components back to original coordinates
for i = 1:3
    stress = [s12(:,i),s13(:,i)]*R';
    s12(:,i) = stress(:,1);
    s13(:,i) = stress(:,2);
end

end

















