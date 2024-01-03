function [Disp,Dgradient] = compute_disp_stress_kernels_fault(rcv,obs)
% compute displacement & displacement gradient kernels for a discretized
% fault mesh using uniform slip greens function solutions
% AUTHOR
% Rishav Mallick, JPL, 2023

x = obs(:,1);
y = obs(:,2);

Disp = zeros(length(x),rcv.N);
Dgradient = zeros(length(x),rcv.N,2);

dr = -1e-12;

for i = 1:rcv.N
    [s12,s13] = calcstress_antiplanefault(1,x,y,rcv.x2(i),rcv.x3(i),rcv.W(i),rcv.dip(i));
    
    % shift coincident displacement kernel by epsilon
    r = sqrt((x-rcv.x2c(i)).^2 + (y-rcv.x3c(i)).^2);
    index = r == 0;
    xmod = x;
    ymod = y;
    xmod(index) = xmod(index) + dr*rcv.nvec(i,1);
    ymod(index) = ymod(index) + dr*rcv.nvec(i,2);
    u1 = calcdisp_antiplanefault(xmod,ymod,rcv.x2(i),rcv.x3(i),rcv.W(i),rcv.dip(i));

    Disp(:,i) = u1;
    Dgradient(:,i,1) = s12;
    Dgradient(:,i,2) = s13;
end

end