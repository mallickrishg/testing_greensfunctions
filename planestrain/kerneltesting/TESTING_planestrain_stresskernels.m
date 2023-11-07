% script to test how stress singularities manifest in the deviatoric stress
% kernels and the scale at which they are most important 
% for eq-cycle calculations
% 
% Authors:
% Rishav Mallick, JPL/Caltech

clear

addpath ~/Documents/GitHub/utils/

% construct mesh
x2extent = 10e3;
x3extent = 20e3;

Nx2 = 14;
Nx3 = x3extent*(Nx2)/2/x2extent;

x3shift = 100e3;

shz = create_shzmesh(x2extent,Nx2,x3extent,Nx3,x3shift);

%% load kernels for displacement and stress
tic
evl = compute_shzkernels_planestrain(shz);
toc
%% construct deviatoric stress kernels and remove positive eigen values

L2222 = (evl.LL2222 - evl.LL3322 - evl.LL2233 + evl.LL3333)./2;
L2322 = (evl.LL2322 - evl.LL2333)./2;
L2223 = (evl.LL2223 - evl.LL3323);
L2323 = evl.LL2323;

%  combine all individual kernels
bigL = [L2222 L2322;L2223 L2323];

% Lorig = [evl.LL2222, 0.*evl.LL2322, 0.*evl.LL3322;...
%          0.*evl.LL2223, evl.LL2323, 0.*evl.LL3323;...
%          0.*evl.LL2233, 0.*evl.LL2333, evl.LL3333];
[Evec,eig_ps] = eig(bigL);

% remove eigen values that cause instabilities 
eig_ps_c = diag(eig_ps);
eig_ps_c(real(eig_ps_c)>0) = 0;
bigLmod = Evec*diag(eig_ps_c)/Evec;

% store modified kernels
L2222m = bigLmod(1:end/2,1:end/2);
L2322m = bigLmod(1:end/2,end/2+1:end);
L2223m = bigLmod(end/2+1:end,1:end/2);
L2323m = bigLmod(end/2+1:end,end/2+1:end);

% plot eigen values
figure(1),clf
subplot(211)
stem(sort(diag(eig_ps)),'k.'), hold on
stem(sort(eig_ps_c),'r.')
axis tight, grid on

subplot(212)
plot(real(eig_ps),imag(eig_ps),'ko'), hold on
plot(real(eig_ps_c),imag(eig_ps_c),'r.')
axis tight, grid on
xlim([-1 1].*max(abs(get(gca,'XLim'))))
set(findobj(gcf,'type','axes'),'Fontsize',15,'Linewidth',1)
%% prescribe spatial distribution of strain
x20 = 0;
x30 = -(x3shift + x3extent/2);
r = 3e3;
dummy = zeros(shz.N,1);

shzindex = sqrt((shz.x3-x30).^2 + (shz.x2-x20).^2)<=r;
dummy(shzindex) = 0.2/shz.W(1);

% plot stress components for a given spatial strain distribution 
% (original kernels)
figure(100),clf
set(gcf,'Name','Original Kernels')
subplot(2,2,1)
toplot = L2222*dummy;
plotshz(shz,toplot,1), shading flat
axis tight equal
cb = colorbar;
clim([-1 1].*2)
colormap turbo(20)
title('\sigma_{22}^{d} from \epsilon_{22}^d')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,2)
toplot = L2322*dummy;
plotshz(shz,toplot,1), shading flat
axis tight equal
cb = colorbar;
clim([-1 1].*2)
title('\sigma_{22}^d from \epsilon_{23}')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,3)
toplot = L2223*dummy;
plotshz(shz,toplot,1), shading flat
axis tight equal
cb = colorbar;
clim([-1 1].*2)
title('\sigma_{23} from \epsilon_{22}^d')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,4)
toplot = L2323*dummy;
plotshz(shz,toplot,1), shading flat
axis tight equal
cb = colorbar;
clim([-1 1].*2)
title('\sigma_{23} from \epsilon_{23}')
set(gca,'FontSize',12,'Linewidth',2)

% plot stress components for a given spatial strain distribution 
% (modified kernels)
figure(101),clf
set(gcf,'Name','Modified Kernels')
subplot(2,2,1)
toplot = L2222m*dummy;
plotshz(shz,toplot,1), shading flat
axis tight equal
cb = colorbar;
clim([-1 1].*2)
colormap turbo(20)
title('\sigma_{22}^{d} from \epsilon_{22}^d')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,2)
toplot = L2322m*dummy;
plotshz(shz,toplot,1), shading flat
axis tight equal
cb = colorbar;
clim([-1 1].*2)
title('\sigma_{22}^d from \epsilon_{23}')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,3)
toplot = L2223m*dummy;
plotshz(shz,toplot,1), shading flat
axis tight equal
cb = colorbar;
clim([-1 1].*2)
title('\sigma_{23} from \epsilon_{22}^d')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,4)
toplot = L2323m*dummy;
plotshz(shz,toplot,1), shading flat
axis tight equal
cb = colorbar;
clim([-1 1].*2)
title('\sigma_{23} from \epsilon_{23}')
set(gca,'FontSize',12,'Linewidth',2)

% plot differences in stress components
% original - modified
figure(102),clf
set(gcf,'Name','Original - Modified Kernels')
subplot(2,2,1)
toplot = (L2222-L2222m)*dummy;
plotshz(shz,toplot,1), shading flat
axis tight equal
cb = colorbar;
clim([-1 1].*2)
colormap bluewhitered(1000)
title('\sigma_{22}^{d} from \epsilon_{22}^d')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,2)
toplot = (L2322-L2322m)*dummy;
plotshz(shz,toplot,1), shading flat
axis tight equal
cb = colorbar;
clim([-1 1].*2)
title('\sigma_{22}^d from \epsilon_{23}')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,3)
toplot = (L2223-L2223m)*dummy;
plotshz(shz,toplot,1), shading flat
axis tight equal
cb = colorbar;
clim([-1 1].*2)
title('\sigma_{23} from \epsilon_{22}^d')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,4)
toplot = (L2323-L2323m)*dummy;
plotshz(shz,toplot,1), shading flat
axis tight equal
cb = colorbar;
clim([-1 1].*.2)
title('\sigma_{23} from \epsilon_{23}')
set(gca,'FontSize',12,'Linewidth',2)
set(findobj(gcf,'type','axes'),'CLim',[-1 1]*0.2)

%% show how the kernels were altered
figure(105),clf
subplot(1,3,1)
imagesc(bigL)
colorbar
clim([-1 1].*1000)
axis tight equal

subplot(1,3,2)
imagesc(bigLmod)
colorbar
clim([-1 1].*1000)
axis tight equal

subplot(1,3,3)
imagesc(bigL-bigLmod)
colorbar
clim([-1 1].*100)
axis tight equal
colormap bluewhitered(40)
