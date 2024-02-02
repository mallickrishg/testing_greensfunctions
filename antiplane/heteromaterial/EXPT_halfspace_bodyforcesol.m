% script to solve quasi-static half-space elasticity problem
% with equivalent body force solution for horizontal interfaces
% 
% AUTHOR:
% Rishav Mallick, JPL, 2024

clear

% topo mesh
Ltopo = 50;
Ntopo = Ltopo*25;

% construct source
rcv = create_verticalinterfaces(-5,-0,201,0);
% specify slip taper length-scale
Ltaper = 1;

% construct free surface
topo = create_horizontalinterfaces(-Ltopo/2,Ltopo/2,Ntopo+1,0);

%% compute tapered source term at the interfaces
source = ones(rcv.N,1);

source(rcv.x3c > -0) = 0;
source(rcv.x3c < -Ltaper) = 1;
index = (rcv.x3c<=-0) & (rcv.x3c>=-Ltaper);
source(index) = linspace(1,0,length(find(index)));

deltashift = -1e-9;

[Disp,Dgradient] = compute_disp_stress_kernels_fault(rcv,[topo.x2c + deltashift.*topo.nvec(:,1),...
                                                          topo.x3c + deltashift.*topo.nvec(:,2)]);

u0 = Disp*source;
ux0 = Dgradient(:,:,1)*source;
uz0 = Dgradient(:,:,2)*source;

[Ku,Dgradient] = compute_disp_stress_kernels_force(topo,...
                                                         [topo.x2c + deltashift.*topo.nvec(:,1),...
                                                          topo.x3c + deltashift.*topo.nvec(:,2)]);
Kux = Dgradient(:,:,1);
Kuz = Dgradient(:,:,2);

% solve matrix equations
bemsol = -Kuz\uz0;

%% surface and bulk displacements

Nobs = 100;
obs = [linspace(-2,2,Nobs)',zeros(Nobs,1)];

% source contribution
[Disp,~] = compute_disp_stress_kernels_fault(rcv,obs);
uplot_0 = Disp*source;

% compute solution from body force equivalent
[Ku,~] = compute_disp_stress_kernels_force(topo,obs);
uplot_bem = Ku*bemsol;

figure(2),clf
plot(obs(:,1),uplot_0,'-','LineWidth',4), hold on
% plot(obs(:,1),uplot_bem,'-','LineWidth',3)
% plot(obs(:,1),uplot_0 + uplot_bem,'k-','LineWidth',3)
axis tight
% legend('full space','topography','BEM','location','best')
xlabel('x'), ylabel('displacement')
ylim([-1 1]*0.5)
set(gca,'FontSize',15,'Linewidth',1)
% print('disp_surf_source','-djpeg','-r200')

figure(3),clf
plot(obs(:,1),uplot_0 + uplot_bem,'k-','LineWidth',3)
xlabel('x'), ylabel('displacement')
set(gca,'FontSize',15,'Linewidth',1)
ylim([-1 1]*0.5)
% print('disp_surf_bem','-djpeg','-r200')

%% plot displacements in the bulk
nx = 50;
nz = 100;
ox = linspace(-2,2,nx);
oz = linspace(-2,0,nz);
[xg,zg] = meshgrid(ox,oz);
obs = [xg(:),zg(:)];

% source contribution
[Disp,Dgradient] = compute_disp_stress_kernels_fault(rcv,obs);
uplot_0 = Disp*source;
e12_0 = Dgradient(:,:,1)*source;
e13_0 = Dgradient(:,:,2)*source;

% compute solution from body force equivalent
[Ku,Dgradient] = compute_disp_stress_kernels_force(topo,obs);
uplot_bem = Ku*bemsol;
e12_bem = Dgradient(:,:,1)*bemsol;
e13_bem = Dgradient(:,:,2)*bemsol;

figure(4),clf
pcolor(ox,oz,reshape(uplot_0,nz,nx)), shading interp, hold on
contour(ox,oz,reshape(uplot_0,nz,nx),[-1:0.1:1]*0.5,'k-')
axis tight equal
clim([-1 1]*0.5)
cb=colorbar;cb.Label.String = 'u (source)';
ylabel('z'),xlabel('x')
set(gca,'FontSize',15,'Linewidth',1,'TickDir','both')
colormap(ttscm('vik',1000))
% print('disp_bulk_source','-djpeg','-r200')

figure(5),clf
pcolor(ox,oz,reshape(uplot_0+uplot_bem,nz,nx)), shading interp, hold on
contour(ox,oz,reshape(uplot_0+uplot_bem,nz,nx),[-1:0.1:1]*0.5,'k-')
axis tight equal
clim([-1 1]*0.5)
cb=colorbar;cb.Label.String = 'u (BEM)';
ylabel('z'),xlabel('x')
set(gca,'FontSize',15,'Linewidth',1,'TickDir','both')
colormap(ttscm('vik',1000))
% print('disp_bulk_bem','-djpeg','-r200')

figure(10),clf
subplot(2,1,1)
pcolor(ox,oz,reshape(e12_0+e12_bem,nz,nx)), shading interp, hold on
contour(ox,oz,reshape(e12_0+e12_bem,nz,nx),[-1:0.1:1]*0.5,'k-')
axis tight equal
clim([-1 1]*0.5)
cb=colorbar;cb.Label.String = '\epsilon_x (BEM)';
ylabel('z'),xlabel('x')
set(gca,'FontSize',15,'Linewidth',1,'TickDir','both')

subplot(2,1,2)
pcolor(ox,oz,reshape(e13_0+e13_bem,nz,nx)), shading interp, hold on
contour(ox,oz,reshape(e13_0+e13_bem,nz,nx),[-1:0.1:1]*0.5,'k-')
axis tight equal
clim([-1 1]*0.5)
cb=colorbar;cb.Label.String = '\epsilon_z (BEM)';
ylabel('z'),xlabel('x')
set(gca,'FontSize',15,'Linewidth',1,'TickDir','both')
colormap(ttscm('vik',20))

figure(11),clf
toplot = sqrt((e12_0+e12_bem).^2 + (e13_0+e13_bem).^2);
pcolor(ox,oz,reshape(toplot,nz,nx)), shading interp, hold on
contour(ox,oz,reshape(toplot,nz,nx),[0:0.1:1],'w-')
plot(ox,-Ltaper.*ones(nx,1),'r--')
axis tight equal
clim([0 1])
cb=colorbar;cb.Label.String = '|\epsilon| (BEM)';
ylabel('z'),xlabel('x')
set(gca,'FontSize',15,'Linewidth',1,'TickDir','both')
colormap(magmacolor(1000))

figure(12),clf
plot(reshape(toplot,nz,nx),oz,'.-'), hold on
plot(abs(oz),oz,'k-','LineWidth',1)
plot(abs(oz),-Ltaper.*ones(nz,1),'r-')
xlabel('strain'), ylabel('z')
set(gca,'FontSize',15,'Linewidth',1)
