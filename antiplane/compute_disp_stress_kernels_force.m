function [Disp,Stress] = compute_disp_stress_kernels_force(G,rcv,obs)

x = obs(:,1);
y = obs(:,2);

Disp = zeros(length(x),rcv.N);
Stress = zeros(length(x),rcv.N,2);

dr = -1e-12;

for i = 1:rcv.N
    a = rcv.W(i)/2;
    [u1,~,~] = calc_disp_stress_forcefault(x-rcv.xc(i,1),y-rcv.xc(i,2),a,rcv.dip(i));
    
    % shift coincident displacement kernel by epsilon
    r = sqrt((x-rcv.xc(i,1)).^2 + (y-rcv.xc(i,2)).^2);
    index = r == 0;
    xmod = x-rcv.xc(i,1);
    ymod = y-rcv.xc(i,2);
    xmod(index) = xmod(index) + dr*rcv.nv(i,1);
    ymod(index) = ymod(index) + dr*rcv.nv(i,2);
    [~,s12,s13] = calc_disp_stress_forcefault(xmod,ymod,a,rcv.dip(i));

    Disp(:,i) = u1./G;
    Stress(:,i,1) = s12;
    Stress(:,i,2) = s13;
end

end