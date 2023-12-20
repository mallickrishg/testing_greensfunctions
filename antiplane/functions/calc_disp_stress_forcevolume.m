function [u1,s12,s13] = calc_disp_stress_forcevolume(X,Y,Lx,Ly)
% compute displacements and gradients from a rectangular source computed
% using numerical integration of point source greens functions
% INPUTS
% X,Y - observation points relative to the center of the source (X-x0),(Y-y0)
% Lx,Ly - dimensions of the source
% OUTPUTS
% u1 - displacements
% s12,s13 - x,y displacement gradients
% AUTHOR:
% Rishav Mallick, JPL, 2023

u1 = zeros(size(X));
s12 = zeros(size(X));
s13 = zeros(size(X));

parfor k = 1:numel(x_mat)

    fun_u = @(x0,y0) pointforce_disp(X(k),Y(k),x0,y0);
    fun_s12 = @(x0,y0) pointforce_s12(X(k),Y(k),x0,y0);
    fun_s13 = @(x0,y0) pointforce_s13(X(k),Y(k),x0,y0);
    
    % rectangle source
    u1(k) = integral2(fun_u,-Lx/2,Lx/2,-Ly/2,Ly/2);
    s12(k) = integral2(fun_s12,-Lx/2,Lx/2,-Ly/2,Ly/2);
    s13(k) = integral2(fun_s13,-Lx/2,Lx/2,-Ly/2,Ly/2);
  
end

end

function u = pointforce_disp(x,y,xs,ys)
r = sqrt((x-xs).^2 + (y-ys).^2);
u = 1/2/pi.*log(r);
end

function s12 = pointforce_s12(x,y,xs,ys)
r = sqrt((x-xs).^2 + (y-ys).^2);
s12 = 1/2/pi.*(x-xs)./(r.^2);
end

function s13 = pointforce_s13(x,y,xs,ys)
r = sqrt((x-xs).^2 + (y-ys).^2);
s13 = 1/2/pi.*(y-ys)./(r.^2);
end