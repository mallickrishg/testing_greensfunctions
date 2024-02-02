% script to solve quasi-static elasticity problem
% with equivalent body force solution for a fault in a plate
% 
% AUTHOR:
% Rishav Mallick, JPL, 2024

clear

% topo mesh
Ltopo = 40;
Ntopo = Ltopo*37;

% construct source
rcv = create_verticalinterfaces(-2,-0.5,2,0);

% construct free surface
topo = create_horizontalinterfaces(-Ltopo/2,Ltopo/2,Ntopo+1,[0,-2]);

%% compute source term at the interfaces
source = ones(rcv.N,1);
deltashift = 1e-5.*[-ones(topo.N/2,1);ones(topo.N/2,1)];

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
subplot(2,1,1)
plot(obs(:,1),uplot_0,'-','LineWidth',4), hold on
axis tight
xlabel('x'), ylabel('displacement')
ylim([-1 1]*0.5)
set(gca,'FontSize',15,'Linewidth',1)

subplot(2,1,2)
plot(obs(:,1),uplot_0 + uplot_bem,'k-','LineWidth',3)
xlabel('x'), ylabel('displacement')
set(gca,'FontSize',15,'Linewidth',1)
ylim([-1 1]*0.5)

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
subplot(2,1,1)
pcolor(ox,oz,reshape(uplot_0,nz,nx)), shading interp, hold on
contour(ox,oz,reshape(uplot_0,nz,nx),[-1:0.1:1]*0.5,'k-')
axis tight equal
clim([-1 1]*0.5)
cb=colorbar;cb.Label.String = 'u (source)';
ylabel('z'),xlabel('x')
set(gca,'FontSize',15,'Linewidth',1,'TickDir','both')
colormap(ttscm('vik',1000))

subplot(2,1,2)
pcolor(ox,oz,reshape(uplot_0+uplot_bem,nz,nx)), shading interp, hold on
contour(ox,oz,reshape(uplot_0+uplot_bem,nz,nx),[-1:0.1:1]*0.5,'k-')
axis tight equal
clim([-1 1]*0.5)
cb=colorbar;cb.Label.String = 'u (BEM)';
ylabel('z'),xlabel('x')
set(gca,'FontSize',15,'Linewidth',1,'TickDir','both')
colormap(ttscm('vik',1000))
