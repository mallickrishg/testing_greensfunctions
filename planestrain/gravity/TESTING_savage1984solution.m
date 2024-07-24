
clear

nu = 0.25;

tectx = 2;

rhog=1*9.81;

a = 2;
b = -1;
n_u = 50;
n_v = 50;

u = linspace(0,4,n_u);
v = linspace(0,-4,n_v);
[u,v] = meshgrid(u,v);

[Sxx,Syy,Sxy,X,Y] =...
    Savage1984_GravityValleyStress_CoordsIn(tectx,rhog,nu,a,b,u,v);

figure(1),clf
subplot(1,3,1)
pcolor(X,Y,Sxx), shading interp
colorbar
axis tight equal
clim([-1 1].*rhog*max(abs(v(:)))*nu/(1-nu))
subplot(1,3,2)
pcolor(X,Y,Syy), shading interp
colorbar
axis tight equal
clim([-1 1].*rhog*max(abs(v(:))))
subplot(1,3,3)
pcolor(X,Y,Sxy), shading interp
colorbar
axis tight equal
clim([-1 1]*5)
colormap(bluewhitered(20))

set(findall(gcf, 'Type', 'axes'), 'Fontsize', 20,'Linewidth',2);