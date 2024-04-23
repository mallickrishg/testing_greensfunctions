
clear

% construct source
source = create_horizontalinterfaces(-1,1,2,0);

nx = 50;
nz = 50;
ox = linspace(-2,2,nx);
oz = linspace(-2,2,nz);
[xg,zg] = meshgrid(ox,oz);
obs = [xg(:),zg(:)];

[Disp,Dgradient] = compute_disp_stress_kernels_fault(source,[obs(:,1),obs(:,2)]);
slip = -1;

u0 = Disp*slip;
ux0 = Dgradient(:,:,1)*slip;
uz0 = Dgradient(:,:,2)*slip;

figure(1),clf
pcolor(ox,oz,reshape(u0,nz,nx)), shading interp, hold on
contour(ox,oz,reshape(u0,nz,nx),[-1:0.2:1]*0.5,'k-')
axis tight equal
clim([-1 1]*0.5)
cb=colorbar;cb.Label.String = 'u';
ylabel('z'),xlabel('x')
colormap(ttscm('turku',10))
set(gca,'FontSize',15,'Linewidth',1,'TickDir','both')

figure(2),clf
subplot(2,1,1)
pcolor(ox,oz,reshape(ux0,nz,nx)), shading interp, hold on
contour(ox,oz,reshape(ux0,nz,nx),[-1:0.1:1],'k-')
axis tight equal
clim([-1 1])
cb=colorbar;cb.Label.String = 'u_{,x}';
ylabel('z'),xlabel('x')
set(gca,'FontSize',15,'Linewidth',1,'TickDir','both')

subplot(2,1,2)
pcolor(ox,oz,reshape(uz0,nz,nx)), shading interp, hold on
contour(ox,oz,reshape(uz0,nz,nx),[-1:0.1:1],'k-')
axis tight equal
clim([-1 1])
colormap(ttscm('vik',20))
cb=colorbar;cb.Label.String = 'u_{,z}';
ylabel('z'),xlabel('x')
set(gca,'FontSize',15,'Linewidth',1,'TickDir','both')