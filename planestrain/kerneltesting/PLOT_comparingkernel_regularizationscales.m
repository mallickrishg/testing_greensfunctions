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

Nx2 = 60;
Nx3 = x3extent*(Nx2)/2/x2extent;

x3shift = 100e3;

shz = create_shzmesh(x2extent,Nx2,x3extent,Nx3,x3shift);

%% load kernels for displacement and stress
tic
evl_orig = compute_shzstresskernels_planestrain(shz,1);
toc
%% construct deviatoric stress kernels and remove positive eigen values

eta_matrix = 1;

% store original deviatoric stress kernels
L2222o = (evl_orig.LL2222 - evl_orig.LL3322 - evl_orig.LL2233 + evl_orig.LL3333)./(2.*eta_matrix);
L2322o = (evl_orig.LL2322 - evl_orig.LL2333)./(2.*eta_matrix);
L2223o = (evl_orig.LL2223 - evl_orig.LL3323)./eta_matrix;
L2323o = evl_orig.LL2323./eta_matrix;

%  combine all individual deviatorickernels
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
x20 = 0;
x30 = -(x3shift + x3extent/2);
r0 = 4e3;
strainmax = 1e-3;
r = sqrt((shz.x3-x30).^2 + (shz.x2-x20).^2);
rmax = 1*r0;

% Circular source (with linear taper)
strain_source = zeros(shz.N,1);
shzindex = r < r0;
strain_source(shzindex) = strainmax;
shzindex = r >= r0 & r <= rmax;
strain_source(shzindex) = strainmax.*abs((r(shzindex)-rmax)./(rmax-r0));

% Gaussian source
% strain_source = strainmax.*exp(-r.^2/(2*rmax^2));

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
set(gcf,'Name','Percentage change in kernels','Position',[0 0 3 2]*400)
subplot(2,2,1)
toplot = 100.*abs(((L2222o-L2222m)*strain_source)./(L2222o*strain_source));
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
% cb = colorbar;cb.Label.String = '% change';
xlabel('x (km)'), ylabel('Depth (km)')
colormap(flipud(hot(1000)))
title('\sigma_{xx}^{d} from \epsilon_{xx}^d')

subplot(2,2,2)
toplot = 100.*abs(((L2322o-L2322m)*strain_source)./(L2322o*strain_source));
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
% cb = colorbar;cb.Label.String = '% change';
xlabel('x (km)'), ylabel('Depth (km)')
title('\sigma_{xx}^d from \epsilon_{xz}')

subplot(2,2,3)
toplot = 100.*abs(((L2223o-L2223m)*strain_source)./(L2223o*strain_source));
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
% cb = colorbar;cb.Label.String = '% change';
xlabel('x (km)'), ylabel('Depth (km)')
title('\sigma_{xz} from \epsilon_{xx}^d')

subplot(2,2,4)
toplot = 100.*abs(((L2323o-L2323m)*strain_source)./(L2323o*strain_source));
plotshz(shz,toplot,1), shading flat
hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal
% cb = colorbar;cb.Label.String = '% change';
title('\sigma_{xz} from \epsilon_{xz}')
xlabel('x (km)'), ylabel('Depth (km)')
set(findobj(gcf,'type','axes'),'CLim',10.^[0 2],'ColorScale','log','FontSize',20,'Linewidth',2,'TickDir','both','box','on')
print('kerneldifference_image','-djpeg','-r200')

figure(201),clf
cb = colorbar;cb.Label.String = '% change';cb.LineWidth=1.5;cb.TickDirection='both';
colormap(flipud(hot(1000)))
set(gca,'FontSize',20,'CLim',10.^[0 2],'ColorScale','log')
% print('kerneldifference_image_colorbar','-djpeg','-r200')

%% %%%%%%%% plot σ22 component as a function of distance from center of source %%%%%%%%%%
multiplier_vec = [1,1.25,1.5,1.75,2.0];
cspec = cool(length(multiplier_vec));

figure(104),clf
set(gcf,'Color','w','Position',[0 0 3 2]*400)
for i = 1:length(multiplier_vec)
    rmax = multiplier_vec(i)*r0;

    % Circular source (with linear taper)
    strain_source = zeros(shz.N,1);
    shzindex = r < r0;
    strain_source(shzindex) = strainmax;
    shzindex = r >= r0 & r <= rmax;
    strain_source(shzindex) = strainmax.*abs((r(shzindex)-rmax)./(rmax-r0));

    index = abs(shz.x3 + 110e3) == min(abs(shz.x3 + 110e3));

    figure(104)
    subplot(2,2,1)
    toplot = (L2222o*strain_source);
    plot(r(index)./1e3,toplot(index),'-','Linewidth',2,'Color',cspec(i,:)), hold on
    toplot = (L2222m*strain_source);
    plot(r(index)./1e3,toplot(index),'ko','MarkerFaceColor',cspec(i,:))
    axis tight
    xlabel('r (km)'), ylabel('\Delta\sigma (MPa)')
    ylim([-22 15])
    xlim([0 x2extent/1e3])
    % grid on
    legend('original','regularized','box','off')
    set(gca,'FontSize',20,'Linewidth',1.5)
    
    subplot(2,2,3)
    plot(r(index)./1e3,strain_source(index)./strainmax,'o-','LineWidth',2,'Color',cspec(i,:)), hold on
    % plot(r./1e3,strain_source./strainmax,'.','LineWidth',2)
    ylim([0 1])
    xlim([0 x2extent/1e3])
    ylabel('normalized \epsilon_{source}'), xlabel('r (km)')    
    % title('horizontal transect')
    set(gca,'FontSize',20,'Linewidth',1.5)

    % subplot(2,2,3)
    % toplot = 100.*abs(((L2222o-L2222m)*strain_source)./(abs(L2222o*strain_source)+1));
    % semilogy(r(index)./1e3,toplot(index),'ko-','Linewidth',1,'MarkerFaceColor','r','MarkerSize',10), hold on
    % axis tight, grid on
    % xlabel('r (km)'), ylabel('% error')
    % ylim(10.^[-2 2])
    % xlim([0 x2extent/1e3])
    % set(gca,'FontSize',20,'Linewidth',1.5)


    % diagonal transect
    % index = abs(shz.x2 - 0e3) == min(abs(shz.x2 - 0e3));
    index_diagonal = abs(1*shz.x2 + 1*shz.x3 + 110e3)/sqrt(1^2 + 1^2) <= 0.2e3;

    subplot(2,2,2)
    toplot = (L2222o*strain_source);
    plot(r(index_diagonal)./1e3,toplot(index_diagonal),'-','Linewidth',2,'Color',cspec(i,:)), hold on
    toplot = (L2222m*strain_source);
    plot(r(index_diagonal)./1e3,toplot(index_diagonal),'ko','MarkerFaceColor',cspec(i,:))
    axis tight
    xlabel('r (km)'), ylabel('\Delta\sigma (MPa)')
    ylim([-22 10])
    xlim([0 x2extent/1e3])
    grid off
    set(gca,'FontSize',20,'Linewidth',1.5)
    
    subplot(2,2,4)
    plot(r(index_diagonal)./1e3,strain_source(index_diagonal)./strainmax,'-','LineWidth',2,'Color',cspec(i,:)), hold on
    % plot(r./1e3,strain_source./strainmax,'.','LineWidth',2)
    ylabel('normalized \epsilon_{source}'), xlabel('r (km)')
    ylim([0 1])
    xlim([0 x2extent/1e3])
    set(gca,'FontSize',20,'Linewidth',1.5)
    % title('diagonal transect')

    % subplot(2,2,4)
    % toplot = 100.*abs(((L2222o-L2222m)*strain_source)./(abs(L2222o*strain_source)+1));
    % semilogy(r(index_diagonal)./1e3,toplot(index_diagonal),'ko-','Linewidth',1,'MarkerFaceColor','r','MarkerSize',10)
    % axis tight, grid on
    % xlabel('r (km)'), ylabel('% error')
    % ylim(10.^[-2 2])
    % xlim([0 x2extent/1e3])
    % set(gca,'FontSize',20,'Linewidth',1.5)
end
print('kerneldifference_transects','-djpeg','-r200')

% figure(200),clf
% plotshz(shz,(index|index_diagonal).*1,1), shading flat
% hold on
% plot(xcircle,ycircle,'k-','LineWidth',2)
% axis tight equal
% cb = colorbar;
% set(gca,'FontSize',12,'Linewidth',2)