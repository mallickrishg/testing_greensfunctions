% script to solve heterogeneous quasi-static half-space elasticity problem
% with equivalent body force solution for distributed body forces
% 
% AUTHOR:
% Rishav Mallick, JPL, 2024

clear

% bimaterial mesh parameters
Lmesh = 10;
Nmesh = Lmesh*1;

% construct source
rcv = create_verticalinterfaces(-1,1,2,0);

bimat = construct_box(Lmesh, Lmesh, Nmesh, Nmesh);

% shear moduli structure
alpha_x = 0.*ones(bimat.N,1);
alpha_z = ones(bimat.N,1)*0;
alpha_z(abs(bimat.xc(:,2)) <= 1) = 1;

%% compute source term at the interfaces
source = ones(rcv.N,1);

[Disp,Dgradient] = compute_disp_stress_kernels_fault(rcv,bimat.xc);

u0 = Disp*source;
ux0 = Dgradient(:,:,1)*source;
uz0 = Dgradient(:,:,2)*source;

%% solve same problem using body force kernels

[Ku,Dgradient] = compute_disp_stress_kernels_forcevolume(bimat,bimat.xc);
Kux = Dgradient(:,:,1);
Kuz = Dgradient(:,:,2);

%% Body force equation is:
% f_b = α_x(du0/dx + Kx.f_b) + α_z(du0/dz + Kz.f_b)

bigmat = (eye(bimat.N) - repmat(alpha_x,1,bimat.N).*Kux - repmat(alpha_z,1,bimat.N).*Kuz);
rhsvec = alpha_x.*ux0 + alpha_z.*uz0;

bemsol = bigmat\rhsvec;

%% plot body force
figure(1),clf
scatter(bimat.xc(:,1),bimat.xc(:,2),150,bemsol,'filled','MarkerEdgeColor','k')
axis tight equal, grid on, box on
colormap bluewhitered
colorbar

%% compute displacements in the medium
xo = linspace(-5,5,40);
zo = linspace(-5,5,40);
[xg,zg] = meshgrid(xo,zo);
obs = [xg(:),zg(:)];
% source contribution
[Disp,Dgradient0] = compute_disp_stress_kernels_fault(rcv,obs);
uplot_0 = Disp*source;

% compute solution from body force equivalent
[Ku,Dgradient] = compute_disp_stress_kernels_forcevolume(bimat,obs);

%% plot BEM solution
uplot_bem = uplot_0 + Ku*bemsol;
e12_bem = Dgradient0(:,:,1)*source + Dgradient(:,:,1)*bemsol;
e13_bem = Dgradient0(:,:,2)*source + Dgradient(:,:,2)*bemsol;


figure(2),clf
pcolor(xo,zo,reshape(uplot_bem,length(zo),length(xo))), shading interp, hold on
contour(xo,zo,reshape(uplot_bem,length(zo),length(xo)),[-1:0.1:1]*0.5,'k')
contour(xo,zo,reshape(uplot_0,length(zo),length(xo)),[-1:0.1:1]*0.5,'k--')
axis tight equal, grid on, box on
colormap bluewhitered(20)
colorbar
clim([-1 1]*0.5)
set(gca,'Fontsize',20,'LineWidth',2,'TickDir','both')












