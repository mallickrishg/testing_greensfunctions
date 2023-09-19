% script to calculate the displacement and stresses from 
% an eigen strain source whose strength is spatially varying as
% edot(x) = α + βx

clear
addpath functions/
addpath ~/Dropbox/scripts/utils/

% Set model parameters 
% Elasticity parameters
mu = 1;

% source strength terms (alpha_0 + alpha_1*x)
alpha_0 = 1;
alpha_1 = 0;

% dimensions of source
% specify location and dimensions of cources
xloc = 0;
yloc = 0;
Rx = 1;
Ry = 1.5;

%% discretize evaluation points
nx = 100;
ny = nx;

x_vec = linspace(-4, 4, nx);
y_vec = linspace(-4, 4, ny);
[x_mat, y_mat] = meshgrid(x_vec, y_vec);

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
    Ls12(k) = mu.*integral2(fun_s12,-Rx,Rx,-Ry,Ry);
    Ls13(k) = mu.*integral2(fun_s13,-Rx,Rx,-Ry,Ry);
end  
toc

% line source at x = ±Rx
[Klu1,Kls12,Kls13] = calc_disp_stress_forcefault(x_mat(:)-(xloc-Rx),y_mat(:)-yloc,Ry,90);
[Kru1,Krs12,Krs13] = calc_disp_stress_forcefault(x_mat(:)-(xloc+Rx),y_mat(:)-yloc,Ry,90);

%% compute resulting displacement and stresses
u1 = Lu1(:).*alpha_1 - Kru1.*(alpha_0 + alpha_1*Rx) + Klu1.*(alpha_0-alpha_1*Rx);
s12 = Ls12(:).*alpha_1 - Krs12.*(alpha_0 + alpha_1*Rx) + Kls12.*(alpha_0-alpha_1*Rx);
s13 = Ls13(:).*alpha_1 - Krs13.*(alpha_0 + alpha_1*Rx) + Kls13.*(alpha_0-alpha_1*Rx);

figure(1),clf
pcolor(x_vec,y_vec,reshape(u1,ny,nx)), shading interp, hold on
contour(x_vec,y_vec,reshape(u1,ny,nx),linspace(-1,1,21).*max(abs(u1(:))),'k-')
axis tight equal
colorbar
clim([-1 1].*max(abs(u1(:))))
colormap bluewhitered(20)

figure(2),clf
subplot(2,1,1)
pcolor(x_vec,y_vec,reshape(s12,ny,nx)), shading interp, hold on
contour(x_vec,y_vec,reshape(s12,ny,nx),linspace(-1,1,21).*max(abs(s12(:))),'k-')
axis tight equal
clim([-1 1].*max(abs(s12(:))))
colorbar

subplot(2,1,2)
pcolor(x_vec,y_vec,reshape(s13,ny,nx)), shading interp, hold on
contour(x_vec,y_vec,reshape(s13,ny,nx),linspace(-1,1,21).*max(abs(s13(:))),'k-')
axis tight equal
clim([-1 1].*max(abs(s13(:))))
colorbar

colormap bluewhitered(20)

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