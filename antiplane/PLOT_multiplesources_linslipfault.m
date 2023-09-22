
clear
addpath functions/

% Set model parameters 
% Elasticity parameters
mu = 1;

% source strength terms (alpha_0 + alpha_1*x)
Nsources = 3;
alpha_0 = 1;
alpha_1 = 1;

% dimensions of source
% specify location and dimensions of cources
xloc = ones(Nsources,1).*0;
yloc = [2,0,-2];
w = ones(Nsources,1).*1;

% discretize evaluation points
nx = 500;
ny = nx;

x_vec = linspace(-4, 4, nx);
y_vec = linspace(-4, 4, ny);
[x_mat, y_mat] = meshgrid(x_vec, y_vec);

%% compute displacement and stress
Ku = zeros(nx*ny,2,Nsources);
Ks12 = zeros(nx*ny,2,Nsources);
Ks13 = zeros(nx*ny,2,Nsources);

for i = 1:Nsources
    [u1,s12,s13] = calc_disp_stress_linfault(x_mat(:),y_mat(:),...
        xloc(i),yloc(i),w(i),[alpha_0,alpha_1]);
    Ku(:,:,i) = u1;
    Ks12(:,:,i) = s12;
    Ks13(:,:,i) = s13;
end

%% plot displacement and stresses

alpha0_vec = [1,1.75,0.75];
alpha1_vec = [-1,0.25,0.75];
sources = [alpha0_vec;alpha1_vec];

% use tensor products to contract 3-d matrices
u1 = tensorprod(Ku,sources,[2 3],[1 2]);
s12 = tensorprod(Ks12,sources,[2 3],[1 2]);
s13 = tensorprod(Ks13,sources,[2 3],[1 2]);

figure(1),clf
pcolor(x_vec,y_vec,reshape(u1,ny,nx)), shading interp, hold on
contour(x_vec,y_vec,reshape(u1,ny,nx),linspace(-1,1,11).*max(abs(u1(:))),'k-')
axis tight equal
cb=colorbar; cb.Location='northoutside';
clim([-1 1].*max(abs(u1(:))))
colormap bluewhitered(1000)

figure(2),clf
subplot(1,2,1)
pcolor(x_vec,y_vec,reshape(s12,ny,nx)), shading interp, hold on
contour(x_vec,y_vec,reshape(s12,ny,nx),linspace(-1,1,21).*max(abs(s12(:))),'k-')
axis tight equal
% clim([-1 1].*max(abs(s12(:))))
clim([-1 1].*1)
colorbar

subplot(1,2,2)
pcolor(x_vec,y_vec,reshape(s13,ny,nx)), shading interp, hold on
contour(x_vec,y_vec,reshape(s13,ny,nx),linspace(-1,1,21).*max(abs(s13(:))),'k-')
axis tight equal
% clim([-1 1].*max(abs(s13(:))))
clim([-1 1].*1)
colorbar

colormap bluewhitered(1000)

figure(3),clf
subplot(2,1,1)
toplot = reshape(u1,ny,nx);
plot(-toplot(:,nx/2)+toplot(:,nx/2+1),y_vec,'-','LineWidth',2)
axis tight, grid on, box on
subplot(2,1,2)
toplot = reshape(s12,ny,nx);
plot(toplot(:,nx/2),y_vec,'-','LineWidth',2), hold on
toplot = reshape(s13,ny,nx);
plot(toplot(:,nx/2),y_vec,'-','LineWidth',2)
axis tight
grid on, box on

