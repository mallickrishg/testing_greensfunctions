% script to solve quasi-static half-space elasticity problem
% with equivalent body force solution for horizontal interfaces
% 
% AUTHOR:
% Rishav Mallick, JPL, 2024

clear

% topo mesh
Ltopo = 50;
Ntopo = Ltopo*100;

% construct source
rcv = create_verticalinterfaces(-1.5,-0.5,2,0);

% construct free surface
topo = create_horizontalinterfaces(-Ltopo/2,Ltopo/2,Ntopo+1,0);

%% compute source term at the interfaces
source = ones(rcv.N,1);
deltashift = -1e-5;

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
plot(obs(:,1),uplot_bem,'-','LineWidth',3)
plot(obs(:,1),uplot_0 + uplot_bem,'k-','LineWidth',3)
axis tight, grid on
legend('full space','topography','BEM','location','best')
xlabel('x'), ylabel('displacement')
% ylim([-1 1]*0.45)
set(gca,'FontSize',15,'Linewidth',1)
print('disp_surf_source','-djpeg','-r200')

figure(3),clf
plot(obs(:,1),uplot_0 + uplot_bem,'k-','LineWidth',3)
xlabel('x'), ylabel('displacement')
set(gca,'FontSize',15,'Linewidth',1)
print('disp_surf_bem','-djpeg','-r200')

%% plot displacements in the bulk
nx = 50;
nz = 51;
ox = linspace(-2,2,nx);
oz = linspace(-2,0,nz);
[xg,zg] = meshgrid(ox,oz);
obs = [xg(:),zg(:)];

% source contribution
[Disp,~] = compute_disp_stress_kernels_fault(rcv,obs);
uplot_0 = Disp*source;

% compute solution from body force equivalent
[Ku,~] = compute_disp_stress_kernels_force(topo,obs);
uplot_bem = Ku*bemsol;

figure(4),clf
pcolor(ox,oz,reshape(uplot_0,nz,nx)), shading interp, hold on
contour(ox,oz,reshape(uplot_0,nz,nx),[-1:0.1:1]*0.5,'k-')
axis tight equal
clim([-1 1]*0.5)
cb=colorbar;cb.Label.String = 'u (source)';
ylabel('z'),xlabel('x')
set(gca,'FontSize',20,'Linewidth',1,'TickDir','both')
colormap(ttscm('vik',1000))
print('disp_bulk_source','-djpeg','-r200')

figure(5),clf
pcolor(ox,oz,reshape(uplot_0+uplot_bem,nz,nx)), shading interp, hold on
contour(ox,oz,reshape(uplot_0+uplot_bem,nz,nx),[-1:0.1:1]*0.5,'k-')
axis tight equal
clim([-1 1]*0.5)
cb=colorbar;cb.Label.String = 'u (BEM)';
ylabel('z'),xlabel('x')
set(gca,'FontSize',20,'Linewidth',1,'TickDir','both')
colormap(ttscm('vik',1000))
print('disp_bulk_bem','-djpeg','-r200')
