function [Disp,Stress] = LinForceKernelFS(xo,yo,xf,yf,w,mu)
% compute displacement and stress kernels for a linearly varying force on a
% horizontal source element (-w <= x <= w, y = 0)
% INPUTS
% x,y - observation locations provided as individual vectors [Nobs x 1]
% xf,yf - source element center location (scalars)
% w - source element half-length
% fx,fy - force components (scalars)
% nu,mu - Elastic parameters
% OUTPUTS
% Disp - 4-d displacement kernels [Nobs x (ux or uy) x (fx or fy) x 2 basis functions]
% Stress - 4-d stress_kernels     [Nobs x (sxx,sxy,syy) x (fx or fy) x 2 basis functions]
% 
% AUTHORS
% Rishav Mallick, JPL, 2024

Nobs = length(xo);

x = xo - xf;
y = yo - yf;

u_1 = @(w) 0.5.*((1/8).*pi.^(-1).*w.^(-1).*((-4).*w.*(2.*w+x)+4.*(w+x).*y.*atan2(w-x,y)+...
    4.*(w+x).*y.*atan2(w+x,y)+((w+(-1).*x) ...
  .*(3.*w+x)+y.^2).*log((w+(-1).*x).^2+y.^2)+(w+x+(-1).*y).*(w+x+y) ...
  .*log((w+x).^2+y.^2)));

u_2 = @(w) 0.5.*((1/8).*pi.^(-1).*w.^(-1).*(4.*w.*((-2).*w+x)+4.*(w+(-1).*x).*y.* ...
  atan2(w-x,y)+4.*(w+(-1).*x).*y.*atan2(w+x,y) ...
  +(w+(-1).*x+(-1).*y).*(w+(-1).*x+y).*log((w+(-1).*x).^2+y.^2)+(( ...
  3.*w+(-1).*x).*(w+x)+y.^2).*log((w+x).^2+y.^2)));

ex_1 = @(w) (1/8).*pi.^(-1).*w.^(-1).*((-4).*w+2.*y.*(atan2(w-x,y)+atan2(w+x,y))+...
    (-1).*(w+x).*log((w+(-1).*x).^2+y.^2)+(w+x).*log((w+x).^2+y.^2));

ex_2 = @(w) (-1/8).*pi.^(-1).*w.^(-1).*((-4).*w+2.*y.*(atan2(w-x,y)+atan2(w+x,y))+(w+(-1).*x).*(log((w+(-1).*x).^2+y.^2)+ ...
  (-1).*log((w+x).^2+y.^2)));

ey_1 = @(w) (1/8).*pi.^(-1).*w.^(-1).*(2.*(w+x).*(atan2(w-x,y)+ ...
  atan2(w+x,y))+y.*(log((w+(-1).*x).^2+y.^2)+(-1).*log((w+x) ...
  .^2+y.^2)));

ey_2 = @(w) (1/8).*pi.^(-1).*w.^(-1).*(2.*(w+(-1).*x).*(atan2(w-x,y)+atan2(w+x,y))+...
    y.*((-1).*log((w+(-1).*x).^2+y.^2)+log((w+x).^2+y.^2)));

% Store displacement kernels (2-d matrix)
% Disp_kernels - [Nobs x 2 basis functions]
Disp = zeros(Nobs,2);
Disp(:,:) = [u_1(w),u_2(w)];


% Store stress kernels (3-d matrix)
% Stress_kernels - [Nobs x (sx or sy) x 2 basis functions]
Stress = zeros(Nobs,2,2);

Stress(:,1,:) = mu.*[ex_1(w),ex_2(w)];
Stress(:,2,:) = mu.*[ey_1(w),ey_2(w)];


end