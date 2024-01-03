function [Disp,Strain] = compute_disp_stress_kernels_force(rcv,obs)
% compute displacement & displacement gradient kernels for a discretized
% line mesh using uniform force greens function solutions
% AUTHOR
% Rishav Mallick, JPL, 2023

x = obs(:,1);
y = obs(:,2);

Disp = zeros(length(x),rcv.N);
Strain = zeros(length(x),rcv.N,2);

dr = -1e-12;

for i = 1:rcv.N
    a = rcv.W(i)/2;
    [u1,~,~] = calc_disp_strain_forcefault(x-rcv.x2c(i),y-rcv.x3c(i),a,rcv.dip(i));
    
    % shift coincident strain kernel by epsilon
    r = sqrt((x-rcv.x2c(i)).^2 + (y-rcv.x3c(i)).^2);
    index = r == 0;
    xmod = x-rcv.x2c(i);
    ymod = y-rcv.x3c(i);
    xmod(index) = xmod(index) + dr*rcv.nvec(i,1);
    ymod(index) = ymod(index) + dr*rcv.nvec(i,2);
    [~,s12,s13] = calc_disp_strain_forcefault(xmod,ymod,a,rcv.dip(i));

    Disp(:,i) = u1;
    Strain(:,i,1) = s12;
    Strain(:,i,2) = s13;
end

end