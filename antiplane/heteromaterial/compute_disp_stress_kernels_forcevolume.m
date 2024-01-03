function [Disp,Strain] = compute_disp_stress_kernels_forcevolume(rcv,obs)
% compute displacement & displacement gradient kernels for a discretized
% volume mesh using uniform force greens function solutions
% AUTHOR
% Rishav Mallick, JPL, 2024

x = obs(:,1);
y = obs(:,2);

Disp = zeros(length(x),rcv.N);
Strain = zeros(length(x),rcv.N,2);

for i = 1:rcv.N

    x_mat = x-rcv.xc(i,1);
    y_mat = y-rcv.xc(i,2);    
    
    Rx = rcv.Lx(i);
    Rz = rcv.Lz(i);

    parfor k = 1:numel(x_mat)
        fun_u = @(x0,y0) pointforce_disp(x_mat(k),y_mat(k),x0,y0);
        fun_s12 = @(x0,y0) pointforce_s12(x_mat(k),y_mat(k),x0,y0);
        fun_s13 = @(x0,y0) pointforce_s13(x_mat(k),y_mat(k),x0,y0);

        % rectangle source
        u1(k) = integral2(fun_u,-Rx/2,Rx/2,-Rz/2,Rz/2);
        s12(k) = integral2(fun_s12,-Rx/2,Rx/2,-Rz/2,Rz/2);
        s13(k) = integral2(fun_s13,-Rx/2,Rx/2,-Rz/2,Rz/2);

    end

    Disp(:,i) = u1;
    Strain(:,i,1) = s12;
    Strain(:,i,2) = s13;
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