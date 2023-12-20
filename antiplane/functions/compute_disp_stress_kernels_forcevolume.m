function [Disp,Dgradients] = compute_disp_stress_kernels_forcevolume(rcv,obs)
% use 2-d source solutions to construct disp and disp-gradient kernels
% INPUTS
% rcv - data structure containing source geometry
% obs - Nobs x 2 matrix containing (x,z) locations of observation points
% OUTPUTS
% Disp - Nobs x rcv.N displacement kernels at obs
% Dgradients - Nobs x rcv.N x 2 kernels for du/dx & du/dz
% AUTHOR:
% Rishav Mallick, JPL, 2023

x = obs(:,1);
y = obs(:,2);
Nobs = length(x);

Disp = zeros(Nobs,rcv.N);
Dgradients = zeros(Nobs,rcv.N,2);

% construct kernels
for i = 1:rcv.N
    [u1,s12,s13] = calc_disp_stress_forcevolume(x-rcv.xc(i,1),y-rcv.xc(i,2),rcv.Lx(i),rcv.Ly(i));
        
    Disp(:,i) = u1;
    Dgradients(:,i,1) = s12;
    Dgradients(:,i,2) = s13;
end

end