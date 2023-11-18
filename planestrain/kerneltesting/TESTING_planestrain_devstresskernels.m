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

Nx2 = 20;
Nx3 = x3extent*(Nx2)/2/x2extent;

x3shift = 100e3;

shz = create_shzmesh(x2extent,Nx2,x3extent,Nx3,x3shift);

nobs = 20;
xobs = linspace(-x2extent,x2extent,nobs);
zobs = linspace(-(x3shift+x3extent),-x3shift,nobs);
[xobs,zobs] = meshgrid(xobs,zobs);
obs = [xobs(:), zobs(:)];

%% load kernels for displacement and stress
tic
evl_orig = compute_shzstresskernels_planestrain(shz,1);
% stresskernel = evl_orig;
dispkernel = compute_shzdispkernels_planestrain(shz,obs);
toc
%% construct deviatoric stress kernels and remove positive eigen values

% L2222 = (stresskernel.LL2222 - stresskernel.LL3322 - stresskernel.LL2233 + stresskernel.LL3333)./2;
% L2322 = (stresskernel.LL2322 - stresskernel.LL2333)./2;
% L2223 = (stresskernel.LL2223 - stresskernel.LL3323);
% L2323 = stresskernel.LL2323;

% impose viscosity shape function
x20 = 0;
x30 = -(x3shift + x3extent/2);
% r0 = 5e3;
% r = sqrt((shz.x2-x20).^2 + (shz.x3-x30).^2 );
% eta_f = 1 - 0.*exp(-(r.^2)./(2*r0^2));
eta_matrix = 1;%repmat(eta_f,1,shz.N);

% store original deviatoric stress kernels
L2222o = (evl_orig.LL2222 - evl_orig.LL3322 - evl_orig.LL2233 + evl_orig.LL3333)./(2.*eta_matrix);
L2322o = (evl_orig.LL2322 - evl_orig.LL2333)./(2.*eta_matrix);
L2223o = (evl_orig.LL2223 - evl_orig.LL3323)./eta_matrix;
L2323o = evl_orig.LL2323./eta_matrix;

%  combine all individual deviatorickernels
% bigL = [L2222./eta_matrix L2322./eta_matrix;...
        % L2223./eta_matrix L2323./eta_matrix];
bigL = [L2222o L2322o;...
        L2223o L2323o];

[Evec,eig_ps] = eig(bigL);

% remove eigen values that cause instabilities 
eig_ps_c = diag(eig_ps);
eig_ps_c(real(eig_ps_c)>0) = 0;
bigLmod = real(Evec*diag(eig_ps_c)/Evec);

% store modified kernels
L2222m = bigLmod(1:end/2,1:end/2);
L2322m = bigLmod(1:end/2,end/2+1:end);
L2223m = bigLmod(end/2+1:end,1:end/2);
L2323m = bigLmod(end/2+1:end,end/2+1:end);

% plot eigen values
figure(1),clf
subplot(211)
stem(sort(real(diag(eig_ps))),'k.'), hold on
stem(sort(real(eig_ps_c)),'r.')
axis tight, grid on

subplot(212)
plot(real(eig_ps),imag(eig_ps),'ko'), hold on
plot(real(eig_ps_c),imag(eig_ps_c),'r.')
axis tight, grid on
xlim([-1 1].*max(abs(get(gca,'XLim'))))
set(findobj(gcf,'type','axes'),'Fontsize',15,'Linewidth',1)

%% %%%%%%%%%%% prescribe spatial distribution of strain %%%%%%%%%%%%%
% x20 = 0;
% x30 = -(x3shift + x3extent/2);
r0 = 4.1e3;
strainmax = 1e-3;%1000/(r0^2);
r = sqrt((shz.x3-x30).^2 + (shz.x2-x20).^2);
rmax = 1.3*r0;
% rmax = r0 + 3e3;

strain_source = zeros(shz.N,1);
shzindex = r < r0;
strain_source(shzindex) = strainmax;
shzindex = r >= r0 & r <= rmax;
strain_source(shzindex) = strainmax.*abs((r(shzindex)-rmax)./(rmax-r0));

% dummy = 3e-4.*(exp(-(r.^2)./2/r0^2)).^2;
% circle source
theta = linspace(0,360,100);
xcircle = r0.*cosd(theta)./1e3;
ycircle = (r0.*sind(theta)+x30)./1e3;

figure(2),clf
subplot(2,1,1)
plotshz(shz,strain_source./strainmax,1), axis tight equal
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
colorbar
colormap(parula(100))
subplot(2,1,2)
plot(shz.x2(:)./1e3,strain_source./strainmax,'o','LineWidth',2), hold on
plot([r0 rmax]./1e3,[1,0],'k-','LineWidth',2)
plot(-[r0 rmax]./1e3,[1,0],'k-','LineWidth',2)
plot(r0./1e3*[1 -1; 1 -1],([1;1]*get(gca,'YLim'))','k--')
axis tight

% plot stress components for a given spatial strain distribution 
% (original kernels)
figure(100),clf
set(gcf,'Name','Original Kernels')
subplot(2,2,1)
toplot = L2222o*strain_source;
maxplotval = max(abs(toplot));
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
% quiver(obs(:,1)./1e3,obs(:,2)./1e3,(dispkernel.Gx22-dispkernel.Gx33)*strain_source,(dispkernel.Gz22-dispkernel.Gz33)*strain_source,'r')
axis tight equal
cb = colorbar;
clim([-1 1].*maxplotval)
colormap turbo(20)
title('\sigma_{22}^{d} from \epsilon_{22}^d')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,2)
toplot = L2322o*strain_source;
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
% quiver(obs(:,1)./1e3,obs(:,2)./1e3,dispkernel.Gx23*strain_source,dispkernel.Gz23*strain_source,'r')
axis tight equal
cb = colorbar;
clim([-1 1].*maxplotval)
title('\sigma_{22}^d from \epsilon_{23}')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,3)
toplot = L2223o*strain_source;
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
cb = colorbar;
clim([-1 1].*maxplotval)
title('\sigma_{23} from \epsilon_{22}^d')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,4)
toplot = L2323o*strain_source;
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
cb = colorbar;
clim([-1 1].*maxplotval)
title('\sigma_{23} from \epsilon_{23}')
set(gca,'FontSize',12,'Linewidth',2)

% plot stress components for a given spatial strain distribution 
% (modified kernels)
figure(101),clf
set(gcf,'Name','Modified Kernels')
subplot(2,2,1)
toplot = L2222m*strain_source;
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
cb = colorbar;
clim([-1 1].*maxplotval)
colormap turbo(20)
title('\sigma_{22}^{d} from \epsilon_{22}^d')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,2)
toplot = L2322m*strain_source;
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
cb = colorbar;
clim([-1 1].*maxplotval)
title('\sigma_{22}^d from \epsilon_{23}')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,3)
toplot = L2223m*strain_source;
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
cb = colorbar;
clim([-1 1].*maxplotval)
title('\sigma_{23} from \epsilon_{22}^d')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,4)
toplot = L2323m*strain_source;
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
cb = colorbar;
clim([-1 1].*maxplotval)
title('\sigma_{23} from \epsilon_{23}')
set(gca,'FontSize',12,'Linewidth',2)

% plot differences in stress components
% original - modified
figure(102),clf
set(gcf,'Name','Original - Modified Kernels')
subplot(2,2,1)
toplot = (L2222o-L2222m)*strain_source;
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
cb = colorbar;
clim([-1 1])
colormap bluewhitered(1000)
title('\sigma_{22}^{d} from \epsilon_{22}^d')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,2)
toplot = (L2322o-L2322m)*strain_source;
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
cb = colorbar;
clim([-1 1])
title('\sigma_{22}^d from \epsilon_{23}')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,3)
toplot = (L2223o-L2223m)*strain_source;
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
cb = colorbar;
clim([-1 1])
title('\sigma_{23} from \epsilon_{22}^d')
set(gca,'FontSize',12,'Linewidth',2)

subplot(2,2,4)
toplot = (L2323o-L2323m)*strain_source;
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
cb = colorbar;
clim([-1 1])
title('\sigma_{23} from \epsilon_{23}')
set(gca,'FontSize',12,'Linewidth',2)
set(findobj(gcf,'type','axes'),'CLim',[-1 1]*maxplotval/10)

% plot differences in stress components as percentage
% original - modified
figure(103),clf
set(gcf,'Name','Percentage change in kernels')
subplot(2,2,1)
toplot = 100.*abs(((L2222o-L2222m)*strain_source)./(L2222o*strain_source));
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
cb = colorbar;cb.Label.String = '% change';
% clim([-1 1])
colormap(flipud(hot(100)))
title('\sigma_{22}^{d} from \epsilon_{22}^d')

subplot(2,2,2)
toplot = 100.*abs(((L2322o-L2322m)*strain_source)./(L2322o*strain_source));
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
cb = colorbar;cb.Label.String = '% change';
title('\sigma_{22}^d from \epsilon_{23}')

subplot(2,2,3)
toplot = 100.*abs(((L2223o-L2223m)*strain_source)./(L2223o*strain_source));
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
cb = colorbar;cb.Label.String = '% change';
title('\sigma_{23} from \epsilon_{22}^d')

subplot(2,2,4)
toplot = 100.*abs(((L2323o-L2323m)*strain_source)./(L2323o*strain_source));
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
cb = colorbar;cb.Label.String = '% change';
title('\sigma_{23} from \epsilon_{23}')
set(findobj(gcf,'type','axes'),'CLim',10.^[0 2],'ColorScale','log','FontSize',20,'Linewidth',2,'TickDir','both')

% plot σ22 component as a function of distance from center of source
figure(104),clf
set(gcf,'Color','w')
subplot(2,1,1)
toplot = (L2222o*strain_source);
plot(r./1e3,toplot,'ko','Linewidth',1,'MarkerFaceColor','r','MarkerSize',10), hold on
toplot = (L2222m*strain_source);
plot(r./1e3,toplot,'bo','Linewidth',1.5)
axis tight
xlabel('r (km)'), ylabel('\Delta\tau (MPa)')
ylim([-1 1]*maxplotval)
xlim([0 x2extent*sqrt(2)/1e3])
grid on
set(gca,'FontSize',25,'Linewidth',1.5)
yyaxis right
plot(r./1e3,strain_source./strainmax,'o','LineWidth',2), hold on
plot([r0 rmax]./1e3,[1,0],'-','LineWidth',2)
ylabel('normalized \epsilon_{source}')
legend('Original','regularized','Strain','box','off')

subplot(2,1,2)
toplot = 100.*abs(((L2222o-L2222m)*strain_source)./(abs(L2222o*strain_source)+1));
semilogy(r./1e3,toplot,'ko','Linewidth',1,'MarkerFaceColor','r','MarkerSize',10)
axis tight, grid on
xlabel('r (km)'), ylabel('% error in \sigma_{22}')
ylim(10.^[-2 2])
xlim([0 x2extent*sqrt(2)/1e3])
set(gca,'FontSize',25,'Linewidth',1.5)
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
