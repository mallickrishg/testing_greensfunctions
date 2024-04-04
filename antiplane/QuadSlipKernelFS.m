function [Disp,Stress] = QuadSlipKernelFS(xo,yo,xf,yf,w,mu)
% compute displacement and stress kernels for a quadratically varying slip 
% on a horizontal source element (-w <= x <= w, y = 0)
% INPUTS
% x,y - observation locations provided as individual vectors [Nobs x 1]
% xf,yf - source element center location (scalars)
% w - source element half-length
% mu - Shear modulus
% OUTPUTS
% Disp - displacement kernels [Nobs x 3 basis functions]
% Stress - 3-d stress_kernels     [Nobs x (sx or sy) x 3 basis functions]
% 
% AUTHORS
% Rishav Mallick, JPL, 2024

x = xo - xf;
y = yo - yf;
Nobs = length(x(:,1));

u1 = (3/16).*w.^(-2).*pi.^(-1).*(6.*w.*y+((-2).*w.*x+3.*(x+(-1).*y).*( ...
  x+y)).*atan((w+(-1).*x).*y.^(-1))+((-2).*w.*x+3.*(x+(-1).*y).*(x+ ...
  y)).*atan((w+x).*y.^(-1))+(w+(-3).*x).*y.*((-1).*log((w+(-1).*x) ...
  .^2+y.^2)+log((w+x).^2+y.^2)));

u2 = (1/8).*w.^(-2).*pi.^(-1).*((-18).*w.*y+(4.*w.^2+9.*y.^2).*atan((w+ ...
  (-1).*x).*y.^(-1))+9.*x.^2.*atan(((-1).*w+x).*y.^(-1))+(4.*w.^2+( ...
  -9).*x.^2+9.*y.^2).*atan((w+x).*y.^(-1))+9.*x.*y.*((-1).*log((w+( ...
  -1).*x).^2+y.^2)+log((w+x).^2+y.^2)));

u3 = (3/16).*w.^(-2).*pi.^(-1).*(6.*w.*y+(2.*w.*x+3.*(x+(-1).*y).*(x+y) ...
  ).*atan((w+(-1).*x).*y.^(-1))+(2.*w.*x+3.*(x+(-1).*y).*(x+y)).* ...
  atan((w+x).*y.^(-1))+(w+3.*x).*y.*(log((w+(-1).*x).^2+y.^2)+(-1).* ...
  log((w+x).^2+y.^2)));

% Store displacement kernels (2-d matrix)
% Disp_kernels - [Nobs x 3 basis functions]
Disp = [u1,u2,u3];

ex1 = (3/32).*w.^(-2).*pi.^(-1).*((-2).*(w+(-3).*x).*(atan((w+(-1).*x).* ...
  y.^(-1))+atan((w+x).*y.^(-1)))+y.*(w.^2.*((-1).*((w+(-1).*x).^2+ ...
  y.^2).^(-1)+5.*((w+x).^2+y.^2).^(-1))+3.*log((w+(-1).*x).^2+y.^2)+ ...
  (-3).*log((w+x).^2+y.^2)));

ex2 = (1/16).*w.^(-2).*pi.^(-1).*((-18).*x.*(atan((w+(-1).*x).*y.^(-1))+ ...
  atan((w+x).*y.^(-1)))+y.*(20.*w.^3.*x.*(w.^4+2.*w.^2.*((-1).*x.^2+ ...
  y.^2)+(x.^2+y.^2).^2).^(-1)+(-9).*log((w+(-1).*x).^2+y.^2)+9.*log( ...
  (w+x).^2+y.^2)));

ex3 = (3/32).*w.^(-2).*pi.^(-1).*(2.*(w+3.*x).*(atan((w+(-1).*x).*y.^( ...
  -1))+atan((w+x).*y.^(-1)))+y.*(w.^2.*((-5).*((w+(-1).*x).^2+y.^2) ...
  .^(-1)+((w+x).^2+y.^2).^(-1))+3.*log((w+(-1).*x).^2+y.^2)+(-3).* ...
  log((w+x).^2+y.^2)));

sx = mu.*[ex1,ex2,ex3];

ey1 = (3/32).*w.^(-2).*pi.^(-1).*(w.*(12+w.*((-1).*w+x).*((w+(-1).*x) ...
  .^2+y.^2).^(-1)+(-5).*w.*(w+x).*((w+x).^2+y.^2).^(-1))+(-6).*y.*( ...
  atan((w+(-1).*x).*y.^(-1))+atan((w+x).*y.^(-1)))+(-1).*(w+(-3).*x) ...
  .*(log((w+(-1).*x).^2+y.^2)+(-1).*log((w+x).^2+y.^2)));

ey2 = (-1/16).*w.^(-2).*pi.^(-1).*(w.*(36+5.*w.*((-1).*w+x).*((w+(-1).* ...
  x).^2+y.^2).^(-1)+(-5).*w.*(w+x).*((w+x).^2+y.^2).^(-1))+(-18).* ...
  y.*(atan((w+(-1).*x).*y.^(-1))+atan((w+x).*y.^(-1)))+9.*x.*log((w+ ...
  (-1).*x).^2+y.^2)+(-9).*x.*log((w+x).^2+y.^2));

ey3 = (3/32).*w.^(-2).*pi.^(-1).*(w.*(12+5.*w.*((-1).*w+x).*((w+(-1).*x) ...
  .^2+y.^2).^(-1)+(-1).*w.*(w+x).*((w+x).^2+y.^2).^(-1))+(-6).*y.*( ...
  atan((w+(-1).*x).*y.^(-1))+atan((w+x).*y.^(-1)))+(w+3.*x).*(log(( ...
  w+(-1).*x).^2+y.^2)+(-1).*log((w+x).^2+y.^2)));

sy = mu.*[ey1,ey2,ey3];

% Store stress kernels (3-d matrix)
% Stress_kernels - [Nobs x (sx or sy) x 3 basis functions]
Stress = zeros(Nobs,2,3);
Stress(:,1,:) = sx;
Stress(:,2,:) = sy;

end