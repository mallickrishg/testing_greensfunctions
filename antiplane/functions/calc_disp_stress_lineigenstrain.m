function [u1,s12,s13] = calc_disp_stress_lineigenstrain(x_mat,y_mat,xloc,yloc,Rx,Ry,alpha)
% compute displacement and stress in response to a spatially variable eigen
% strain source (alpha_0 + alpha_1*x).

% source strength terms (alpha_0 + alpha_1*x)
alpha_0 = alpha(1);
alpha_1 = alpha(2);

%% numerical integration (with matlab integral)
Lu1 = zeros(size(x_mat));
Ls12 = zeros(size(x_mat));
Ls13 = zeros(size(x_mat));

tic    
% volume source
parfor k=1:numel(x_mat)
    fun_u = @(x0,y0) pointforce_disp(x_mat(k)-xloc,y_mat(k)-yloc,x0,y0);
    fun_s12 = @(x0,y0) pointforce_s12(x_mat(k)-xloc,y_mat(k)-yloc,x0,y0);
    fun_s13 = @(x0,y0) pointforce_s13(x_mat(k)-xloc,y_mat(k)-yloc,x0,y0);

    % rectangle source
    Lu1(k) = integral2(fun_u,-Rx,Rx,-Ry,Ry);
    Ls12(k) = integral2(fun_s12,-Rx,Rx,-Ry,Ry);
    Ls13(k) = integral2(fun_s13,-Rx,Rx,-Ry,Ry);
end  
toc

% vertical line source at x = ±Rx
[Klu1,Kls12,Kls13] = calc_disp_stress_forcefault(x_mat(:)-(xloc-Rx),y_mat(:)-yloc,Ry,90);
[Kru1,Krs12,Krs13] = calc_disp_stress_forcefault(x_mat(:)-(xloc+Rx),y_mat(:)-yloc,Ry,90);

%% compute resulting displacement and stresses
u1 = Lu1(:).*alpha_1 - Kru1.*(alpha_0 + alpha_1*Rx) + Klu1.*(alpha_0-alpha_1*Rx);
s12 = Ls12(:).*alpha_1 - Krs12.*(alpha_0 + alpha_1*Rx) + Kls12.*(alpha_0-alpha_1*Rx);
s13 = Ls13(:).*alpha_1 - Krs13.*(alpha_0 + alpha_1*Rx) + Kls13.*(alpha_0-alpha_1*Rx);

end

%% define function
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