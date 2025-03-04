% script to convert Mathematica expressions to MATLAB
% for Kelvin problem with linearly varying body forces along a line segment
% 
% AUTHORS:
% Rishav Mallick, JPL, 2024

clear

mu = 1;
nu = 0.25;

Nobs = 1000;
xo = linspace(-2,2,Nobs)';
yo = zeros(Nobs,1);

% Displacement kernels copied over from Mathematica for a horizontal line
% element present at y = 0, -w <= x <= w
% 
% ux,uy kernels have two linear basis functions
% f1: goes from (-w,0) to (w,1)
% f2: goes from (-w,1) to (w,0)

% ux kernels (is non-zero at yo = 0 only for fx)
ux_1 = @(fx,fy,w) ((1/8).*fx.*w.^(-1).*(w+(-1).*xo).*(3.*w+xo).*mu.^(-1).*nu.*(pi+(-1) ...
  .*pi.*nu).^(-1).*log((w+(-1).*xo).^2+yo.^2)+(1/32).*fx.*pi.^(-1).* ...
  w.^(-1).*mu.^(-1).*((-1)+nu).^(-1).*(4.*w.*(8.*w.*((-1)+nu)+xo.*((-3) ...
  +4.*nu))+(-16).*(w+xo).*yo.*((-1)+nu).*atan2((w-xo),yo)+ ...
  (-16).*(w+xo).*yo.*((-1)+nu).*atan2((w+xo),yo)+(3.*(w+(-1).* ...
  xo).*(3.*w+xo)+yo.^2.*(5+(-4).*nu)).*log((w+(-1).*xo).^2+yo.^2)+( ...
  3.*(w+xo).^2+(-5).*yo.^2+(-4).*(w+xo+(-1).*yo).*(w+xo+yo).*nu).* ...
  log((w+xo).^2+yo.^2)) + ...
  (1/16).*fy.*pi.^(-1).*w.^(-1).*yo.*mu.^(-1).*((-1)+nu).^(-1).*(4.*w+ ...
  (-2).*yo.*(atan2((w-xo),yo)+atan2((w+xo),yo))+( ...
  w+xo).*(log((w+(-1).*xo).^2+yo.^2)+(-1).*log((w+xo).^2+yo.^2))))./2;

ux_2 = @(fx,fy,w) (1/64).*fx.*pi.^(-1).*w.^(-1).*(w+(-1).*xo).^2.*mu.^(-1).*(3+(-4).* ...
  nu).*((-1)+nu).^(-1).*log((w+(-1).*xo).^2+yo.^2)+(1/32).*pi.^(-1).* ...
  w.^(-1).*(w+(-1).*xo).*mu.^(-1).*((-1)+nu).^(-1).*((-8).*fx.*yo.*(( ...
  -1)+nu).*(atan2((w-xo),yo)+atan2((w+xo),yo))+fy.* ...
  yo.*log((w+(-1).*xo).^2+yo.^2))+(1/64).*pi.^(-1).*w.^(-1).*mu.^(-1) ...
  .*((-1)+nu).^(-1).*(4.*w.*((-2).*fy.*yo+fx.*xo.*(3+(-4).*nu)+8.*fx.* ...
  w.*((-1)+nu))+yo.^2.*(4.*fy.*(atan2((w-xo),yo)+atan2((w+xo),yo))+...
  fx.*((-5)+4.*nu).*log((w+(-1).*xo).^2+yo.^2))+(2.* ...
  fy.*((-1).*w+xo).*yo+fx.*(3.*(3.*w+(-1).*xo).*(w+xo)+5.*yo.^2+(-4) ...
  .*((3.*w+(-1).*xo).*(w+xo)+yo.^2).*nu)).*log((w+xo).^2+yo.^2));

% uy kernels (is non-zero at yo = 0 only for fy)
uy_1 = @(fx,fy,w) (1/32).*fx.*pi.^(-1).*w.^(-1).*yo.*mu.^(-1).*((-1)+nu).^(-1).*(4.*w+ ...
  (-2).*yo.*(atan2(w-xo,yo)+atan2(w+xo,yo))+( ...
  w+xo).*(log((w+(-1).*xo).^2+yo.^2)+(-1).*log((w+xo).^2+yo.^2)))+( ...
  1/64).*fy.*pi.^(-1).*w.^(-1).*mu.^(-1).*((-1)+nu).^(-1).*(4.*w.*(2.* ...
  w+xo).*((-3)+4.*nu)+8.*(w+xo).*yo.*((-1)+2.*nu).*atan2((-w+xo),yo)+...
  (-8).*(w+xo).*yo.*((-1)+2.*nu).*atan2(w+xo,yo)+( ...
  w.^2.*(9+(-12).*nu)+yo.^2.*(1+(-4).*nu)+2.*w.*xo.*((-3)+4.*nu)+ ...
  xo.^2.*((-3)+4.*nu)).*log((w+(-1).*xo).^2+yo.^2)+(-1).*((-3).*(w+ ...
  xo).^2+yo.^2+4.*(w+xo+(-1).*yo).*(w+xo+yo).*nu).*log((w+xo).^2+ ...
  yo.^2));

uy_2 = @(fx,fy,w) fx.*((-1/32).*pi.^(-1).*w.^(-1).*yo.*mu.^(-1).*((-1)+nu).^(-1).*(4.* ...
  w+(-2).*yo.*(atan2(w-xo,yo)+atan2(w+xo,yo))+ ...
  ((-1).*w+xo).*log((w+(-1).*xo).^2+yo.^2))+(w+(-1).*xo).*yo.*(32.* ...
  pi.*w.*mu+(-32).*pi.*w.*mu.*nu).^(-1).*log(w.^2+2.*w.*xo+xo.^2+yo.^2) ...
  )+fy.*((-1/8).*pi.^(-1).*w.^(-1).*(w+(-1).*xo).*yo.*mu.^(-1).*((-1) ...
  +nu).^(-1).*((-1)+2.*nu).*(atan2(w-xo,yo)+atan2((w+xo),yo))+...
  (1/64).*pi.^(-1).*w.^(-1).*(w+(-1).*xo).^2.*mu.^(-1).* ...
  (3+(-4).*nu).*((-1)+nu).^(-1).*log((w+(-1).*xo).^2+yo.^2)+(1/64).* ...
  pi.^(-1).*w.^(-1).*mu.^(-1).*((-1)+nu).^(-1).*(4.*w.*(2.*w+(-1).*xo) ...
  .*((-3)+4.*nu)+yo.^2.*((-1)+4.*nu).*log((w+(-1).*xo).^2+yo.^2)+(3.*( ...
  3.*w+(-1).*xo).*(w+xo)+yo.^2+(-4).*((3.*w+(-1).*xo).*(w+xo)+yo.^2) ...
  .*nu).*log((w+xo).^2+yo.^2)));

% sxy kernels
sxy_1 = @(fx,fy,w) (-1/8).*fx.*pi.^(-1).*w.^(-1).*((-1)+nu).^(-1).*(2.*w.*(w+(-1).*xo) ...
  .*yo.*((w+(-1).*xo).^2+yo.^2).^(-1)+2.*(w+xo).*((-1)+nu).*atan2((w-xo),yo)+...
  2.*(w+xo).*((-1)+nu).*atan2((w+xo),yo)+( ...
  1/2).*yo.*((-3)+2.*nu).*(log((w+(-1).*xo).^2+yo.^2)+(-1).*log((w+ ...
  xo).^2+yo.^2)))+(1/8).*fy.*pi.^(-1).*w.^(-1).*((-1)+nu).^(-1).*(( ...
  -2).*w.*(w+(-1).*xo).^2.*((w+(-1).*xo).^2+yo.^2).^(-1)+4.*w.*nu+2.* ...
  yo.*nu.*atan2((-w+xo),yo)+(-2).*yo.*nu.*atan2((w+xo),yo)+...
  (1/2).*(w+xo).*((-1)+2.*nu).*(log((w+(-1).*xo).^2+yo.^2)+ ...
  (-1).*log((w+xo).^2+yo.^2)));

sxy_2 = @(fx,fy,w) (-1/16).*fx.*pi.^(-1).*w.^(-1).*((-1)+nu).^(-1).*(4.*w.*(w+xo).* ...
  yo.*((w+xo).^2+yo.^2).^(-1)+4.*(w+(-1).*xo).*((-1)+nu).*atan2((w-xo),yo)+4.*(w+(-1).*xo).*((-1)+nu).*atan2(w+xo,yo)+...
  (-1).*yo.*((-3)+2.*nu).*(log((w+(-1).*xo).^2+yo.^2)+(-1).*log( ...
  (w+xo).^2+yo.^2)))+(-1/16).*fy.*pi.^(-1).*w.^(-1).*((-1)+nu).^(-1) ...
  .*((-4).*w.*(w+xo).^2.*((w+xo).^2+yo.^2).^(-1)+8.*w.*nu+(-4).*yo.* ...
  nu.*atan2(w-xo,yo)+(-4).*yo.*nu.*atan2((w+xo),yo)+...
  ((-1).*w+xo).*((-1)+2.*nu).*(log((w+(-1).*xo).^2+yo.^2)+(-1).* ...
  log((w+xo).^2+yo.^2)));

% sxx kernels
sxx_1 = @(fx,fy,w) (1/8).*fx.*pi.^(-1).*w.^(-1).*((-1)+nu).^(-1).*(((w+(-1).*xo).^2+ ...
  yo.^2).^(-1).*((-6).*w.*(w+(-1).*xo).^2+(-8).*w.*yo.^2)+4.*w.*nu+( ...
  -2).*yo.*((-2)+nu).*atan2((w-xo),yo)+(-2).*yo.*((-2)+nu) ...
  .*atan2((w+xo),yo)+(1/2).*(w+xo).*((-3)+2.*nu).*(log((w+(-1) ...
  .*xo).^2+yo.^2)+(-1).*log((w+xo).^2+yo.^2)))+(1/8).*fy.*pi.^(-1).* ...
  w.^(-1).*((-1)+nu).^(-1).*(2.*w.*((-1).*w+xo).*yo.*((w+(-1).*xo) ...
  .^2+yo.^2).^(-1)+2.*(w+xo).*nu.*atan2((w-xo),yo)+2.*(w+ ...
  xo).*nu.*atan2((w+xo),yo)+(1/2).*yo.*(1+2.*nu).*(log((w+(-1).* ...
  xo).^2+yo.^2)+(-1).*log((w+xo).^2+yo.^2)));

sxx_2 = @(fx,fy,w) (1/16).*fx.*pi.^(-1).*w.^(-1).*((-1)+nu).^(-1).*(4.*w.*((w+xo).^2+ ...
  yo.^2).^(-1).*(3.*(w+xo).^2+4.*yo.^2+(-2).*((w+xo).^2+yo.^2).*nu)+ ...
  4.*yo.*((-2)+nu).*atan2(w-xo,yo)+4.*yo.*((-2)+nu).* ...
  atan2(w+xo,yo)+(w+(-1).*xo).*((-3)+2.*nu).*(log((w+(-1).* ...
  xo).^2+yo.^2)+(-1).*log((w+xo).^2+yo.^2)))+(1/16).*fy.*pi.^(-1).* ...
  w.^(-1).*((-1)+nu).^(-1).*((-4).*w.*(w+xo).*yo.*((w+xo).^2+yo.^2) ...
  .^(-1)+4.*(w+(-1).*xo).*nu.*atan2(w-xo,yo)+4.*(w+(-1) ...
  .*xo).*nu.*atan2(w+xo,yo)+(-1).*yo.*(1+2.*nu).*(log((w+(-1) ...
  .*xo).^2+yo.^2)+(-1).*log((w+xo).^2+yo.^2)));

% syy kernels
syy_1 = @(fx,fy,w) (-1/16).*fx.*pi.^(-1).*w.^(-1).*((-1)+nu).^(-1).*(((w+(-1).*xo).^2+ ...
  yo.^2).^(-1).*((-4).*w.*(w+(-1).*xo).^2+(-8).*w.*yo.^2)+8.*w.*nu+( ...
  -4).*yo.*((-1)+nu).*atan2((w-xo),yo)+(-4).*yo.*((-1)+nu) ...
  .*atan2((w+xo),yo)+(w+xo).*((-1)+2.*nu).*(log((w+(-1).*xo) ...
  .^2+yo.^2)+(-1).*log((w+xo).^2+yo.^2)))+(-1/8).*fy.*pi.^(-1).*w.^( ...
  -1).*((-1)+nu).^(-1).*(2.*w.*((-1).*w+xo).*yo.*((w+(-1).*xo).^2+ ...
  yo.^2).^(-1)+2.*(w+xo).*((-1)+nu).*atan2((w-xo),yo)+2.* ...
  (w+xo).*((-1)+nu).*atan2((w+xo),yo)+(1/2).*yo.*((-1)+2.*nu).*( ...
  log((w+(-1).*xo).^2+yo.^2)+(-1).*log((w+xo).^2+yo.^2)));

syy_2 = @(fx,fy,w) (-1/8).*fx.*pi.^(-1).*w.^(-1).*((-1)+nu).^(-1).*(2.*w.*((w+xo).^2+ ...
  yo.^2).^(-1).*((w+xo).^2+2.*yo.^2+(-2).*((w+xo).^2+yo.^2).*nu)+2.* ...
  yo.*((-1)+nu).*atan2(w-xo,yo)+2.*yo.*((-1)+nu).*atan2((w+xo),yo)+(1/2).*(w+(-1).*xo).*((-1)+2.*nu).*(log((w+(-1).* ...
  xo).^2+yo.^2)+(-1).*log((w+xo).^2+yo.^2)))+(-1/8).*fy.*pi.^(-1).* ...
  w.^(-1).*((-1)+nu).^(-1).*((-2).*w.*(w+xo).*yo.*((w+xo).^2+yo.^2) ...
  .^(-1)+2.*(w+(-1).*xo).*((-1)+nu).*atan2(w-xo,yo)+2.* ...
  (w+(-1).*xo).*((-1)+nu).*atan2(w+xo,yo)+(-1/2).*yo.*((-1)+ ...
  2.*nu).*(log((w+(-1).*xo).^2+yo.^2)+(-1).*log((w+xo).^2+yo.^2)));

% set element length
w = 1;

% Displacements
figure(1),clf
subplot(2,2,1)
plot(xo,ux_1(1,0,w),'Linewidth',2), hold on
plot(xo,ux_2(1,0,w),'Linewidth',2)
xlabel('x')
title('u_x (f_x)','FontWeight','normal')

subplot(2,2,3)
plot(xo,ux_1(0,1,w),'Linewidth',2), hold on
plot(xo,ux_2(0,1,w),'Linewidth',2)
xlabel('x')
title('u_x (f_y)','FontWeight','normal')

subplot(2,2,2)
plot(xo,uy_1(1,0,w),'Linewidth',2), hold on
plot(xo,uy_2(1,0,w),'Linewidth',2)
xlabel('x')
title('u_y (f_x)','FontWeight','normal')

subplot(2,2,4)
plot(xo,uy_1(0,1,w),'Linewidth',2), hold on
plot(xo,uy_2(0,1,w),'Linewidth',2)
xlabel('x')
title('u_y (f_y)','FontWeight','normal')
set(findobj(gcf,'type','axes'),'FontSize',20,'LineWidth',1,'TickDir','both','YLim',[-1 1]*0.2);
for i = 1:4
    subplot(2,2,i)
    plot([1 1]*w,get(gca,'YLim'),'k--')
    plot([1 1]*-w,get(gca,'YLim'),'k--')
end

% Stresses
% sxy
figure(11),clf
subplot(2,3,1)
plot(xo,sxy_1(1,0,w),'Linewidth',2), hold on
plot(xo,sxy_2(1,0,w),'Linewidth',2)
plot(xo,sxy_1(1,0,w)+sxy_2(1,0,w),'k--','Linewidth',2)
xlabel('x')
title('\sigma_{xy} (f_x)','FontWeight','normal')

subplot(2,3,4)
plot(xo,sxy_1(0,1,w),'Linewidth',2), hold on
plot(xo,sxy_2(0,1,w),'Linewidth',2)
plot(xo,sxy_1(0,1,w)+sxy_2(0,1,w),'k--','Linewidth',2)
xlabel('x')
title('\sigma_{xy} (f_y)','FontWeight','normal')

% sxx
subplot(2,3,2)
plot(xo,sxx_1(1,0,w),'Linewidth',2), hold on
plot(xo,sxx_2(1,0,w),'Linewidth',2)
plot(xo,sxx_1(1,0,w)+sxx_2(1,0,w),'k--','Linewidth',2)
xlabel('x')
title('\sigma_{xx} (f_x)','FontWeight','normal')

subplot(2,3,5)
plot(xo,sxx_1(0,1,w),'Linewidth',2), hold on
plot(xo,sxx_2(0,1,w),'Linewidth',2)
plot(xo,sxx_1(0,1,w)+sxx_2(0,1,w),'k--','Linewidth',2)
xlabel('x')
title('\sigma_{xx} (f_y)','FontWeight','normal')

% syy
subplot(2,3,3)
plot(xo,syy_1(1,0,w),'Linewidth',2), hold on
plot(xo,syy_2(1,0,w),'Linewidth',2)
plot(xo,syy_1(1,0,w)+syy_2(1,0,w),'k--','Linewidth',2)
xlabel('x')
title('\sigma_{yy} (f_x)','FontWeight','normal')

subplot(2,3,6)
plot(xo,syy_1(0,1,w),'Linewidth',2), hold on
plot(xo,syy_2(0,1,w),'Linewidth',2)
plot(xo,syy_1(0,1,w)+syy_2(0,1,w),'k--','Linewidth',2)
xlabel('x')
title('\sigma_{yy} (f_y)','FontWeight','normal')

set(findobj(gcf,'type','axes'),'FontSize',20,'LineWidth',1,'TickDir','both','YLim',[-1 1]*0.5);
for i = 1:6
    subplot(2,3,i)
    plot([1 1]*w,get(gca,'YLim'),'k--')
    plot([1 1]*-w,get(gca,'YLim'),'k--')
end


%% plot u,σ in the bulk
Nx = 100;
Ny = 100;
xg = linspace(-2,2,Nx);
yg = linspace(-2,2,Ny);
[xo,yo] = meshgrid(xg,yg);

[Disp,Stress,Strain] = LinForceKernelFS(xo(:),yo(:),0,0,w,nu,mu);

% Disp_kernels - [Nobs x (ux or uy) x (fx or fy) x 2 basis functions]
figure(100),clf
for i = 1:2 % ux or uy
    for j = 1:2 % fx or fy
        plotval = (i-1)*2 + j;
        subplot(2,2,plotval)
        for k = 1 % basis functions
            toplot = Disp(:,i,j,k);
            pcolor(xg,yg,reshape(toplot,Ny,Nx)), shading interp
            colorbar
        end
        if j == 1
            title('f_x')
        else
            title('f_y')
        end
        axis tight equal
    end
end

% Stress or strain kernels - [Nobs x (sxx,sxy,syy) x (fx or fy) x 2 basis functions]
figure(101),clf
for i = 1:3 % sxx, sxy or syy
    for j = 1:2 % fx or fy
        plotval = (i-1)*2 + j;
        subplot(3,2,plotval)
        for k = 1 % basis functions
            toplot = Stress(:,i,j,k);          
            pcolor(xg,yg,reshape(toplot,Ny,Nx)), shading interp
            colorbar, clim([-1 1]*0.5)
            colormap(turbo(20))
        end
        if j == 1
            title('f_x')
        else
            title('f_y')
        end
        axis tight equal
        set(gca,'FontSize',15,'LineWidth',1.5,'TickDir','out')
    end
end

%% test strain using FD displacement gradient calculations
[duxdx,duxdy] = gradient(reshape(Disp(:,1,2,1),Ny,Nx),xg,yg);
[duydx,duydy] = gradient(reshape(Disp(:,2,2,1),Ny,Nx),xg,yg);
toplot = (duxdy + duydx)./2;% exy
% toplot = duxdx;% exx
% toplot = duydy;% eyy

figure(200),clf
subplot(2,1,1)
pcolor(xg,yg,toplot), shading interp
colorbar, 
clim([-1 1]*0.25)
colormap(turbo(20))
axis tight equal

subplot(2,1,2)
toplot = reshape(Strain(:,2,2,1),Ny,Nx);
pcolor(xg,yg,toplot), shading interp
colorbar, 
clim([-1 1]*0.25)
colormap(turbo(20))
axis tight equal



