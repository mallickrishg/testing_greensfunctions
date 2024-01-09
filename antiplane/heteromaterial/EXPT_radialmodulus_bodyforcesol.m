% script to solve half-space quasi-static elasticity problem
% with equivalent body force solution for
% depth-dependent shear modulus structure
% 
% AUTHOR:
% Rishav Mallick, 2024, JPL

clear

% bimaterial mesh parameters
Lmesh = 20;
Nmesh = 200;
bimaterial_z = linspace(-10,-0.1,10);%[-3,-2,-1];

% shear moduli structure (must have dimensions = len(bimaterial) + 1)
mu_structure = linspace(1,1,length(bimaterial_z)+1);

% construct source
rcv = create_verticalinterfaces(-1e5,-0.5,2,0);

% construct bimaterial interface
bimat = create_horizontalinterfaces(-Lmesh/2,Lmesh/2,Nmesh+1,bimaterial_z);

% construct free surface
topo = create_horizontalinterfaces(-Lmesh/2,Lmesh/2,Nmesh+1,0);

% join bimat+topo
indextopo = [false(bimat.N,1);true(topo.N,1)];
combostructure = append_structures(bimat,topo);
%% compute source term at the interface
source = ones(rcv.N,1);
deltashift = -1e-5;

[Disp,Dgradient] = compute_disp_stress_kernels_fault(rcv,[combostructure.x2c + deltashift.*indextopo.*combostructure.nvec(:,1),...
                                                          combostructure.x3c + deltashift.*indextopo.*combostructure.nvec(:,2)]);

u0 = Disp*source;
ux0 = Dgradient(:,:,1)*source;
uz0 = Dgradient(:,:,2)*source;

%% solve problem using body force kernels
[Ku,Dgradient] = compute_disp_stress_kernels_force(combostructure,...
                                                         [combostructure.x2c + deltashift.*indextopo.*combostructure.nvec(:,1),...
                                                          combostructure.x3c + deltashift.*indextopo.*combostructure.nvec(:,2)]);

% Kux = Dgradient(:,:,1);
Kuz = Dgradient(:,:,2);

% the body force equivalent due to a change in shear modulus is only at the
% interface and is a term fb = (dlogμ/dz)(du/dz)
% du/dz = du0/dz + Kuz.fb
% we can solve for the unknown fb as 
% fb = α(du0/dz) + α(Kuz.fb)
% simplifying: 
% fb = [I - αKuz]\(α)(du0/dz)
% we also need to solve for fictitious forces on the topo interface

alpha = ones(bimat.N,1);
for i = 1:length(bimaterial_z) 
    index = (bimat.x3c == bimaterial_z(i));
    alpha(index) = (mu_structure(i)-mu_structure(i+1))/(mu_structure(i+1));
end

% alpha = ones(bimat.N,1);
% index = bimat.x3c>= 0;
% alpha(index) = (1-mu_ratio)/(mu_ratio);
% alpha(~index) = -(1-mu_ratio);

bigmat = [(eye(bimat.N) - repmat(alpha,1,bimat.N).*Kuz(~indextopo,~indextopo)), -repmat(alpha,1,topo.N).*Kuz(~indextopo,indextopo);...
           Kuz(indextopo,~indextopo)                                          , Kuz(indextopo,indextopo)];
rhsvec = [(alpha.*uz0(~indextopo));-uz0(indextopo)];

% solve matrix equations
bemsol = bigmat\rhsvec;
fb = bemsol(~indextopo);
ft = bemsol(indextopo);

% fb = (eye(bimat.N) - repmat(alpha,1,bimat.N).*Kuz)\(alpha.*uz0);

figure(1),clf
% stem(bimat.x2c,fb,'r-','Linewidth',1,'MarkerFaceColor','r'), hold on
for i = 1:bimat.N
    toplotx = bimat.x2(i)+bimat.W(i)*cosd(bimat.dip(i));
    plot([bimat.x2(i),bimat.x2(i),toplotx,toplotx],fb(i).*[0,1,1,0],'k-','Linewidth',1), hold on
end
for i = 1:topo.N
    toplotx = topo.x2(i)+topo.W(i)*cosd(topo.dip(i));
    plot([topo.x2(i),topo.x2(i),toplotx,toplotx],ft(i).*[0,1,1,0],'r-','Linewidth',.1)
end
axis tight, grid on
xlim([-1 1]*5)
xlabel('x'), ylabel('body force strength')
set(gca,'FontSize',15,'Linewidth',1)

%% compute u and du/dx,du/dz in the medium
Nobs = 1000;
% obs = [zeros(Nobs,1),linspace(-3,3,Nobs)'];
obs = [linspace(-2,2,Nobs)',zeros(Nobs,1)-0.05];

% provide shear modulus structure
% muplot = ones(Nobs,1).*mu_ratio;
% muplot(obs(:,2)>bimaterial_z(2) | obs(:,2)<bimaterial_z(1)) = 1;

% source contribution
[Disp,Dgradient0] = compute_disp_stress_kernels_fault(rcv,obs);
uplot_0 = Disp*source;

% compute solution from body force equivalent
[Ku,Dgradient] = compute_disp_stress_kernels_force(combostructure,obs);
uplotf = uplot_0 + Ku*bemsol;
e12f = Dgradient0(:,:,1)*source + Dgradient(:,:,1)*bemsol;
e13f = Dgradient0(:,:,2)*source + Dgradient(:,:,2)*bemsol;
% s13f = muplot.*(e13f);

figure(2),clf
plot(obs(:,1),e12f,'-','LineWidth',2), hold on
plot(obs(:,1),e13f,'r-','LineWidth',2)
grid on, axis tight
ylabel('Strain'),xlabel('x')
ylim([-0.2 1.5])
set(gca,'FontSize',15,'Linewidth',1)

figure(3),clf
subplot(2,1,1)
plot(obs(:,1),uplot_0,'Linewidth',2), hold on
plot(obs(:,1),Ku*bemsol,'r-','LineWidth',2)
axis tight, grid on
ylabel('u (individual contribution)')
ylim([-1 1]*0.5)
set(gca,'FontSize',15,'Linewidth',1)
subplot(2,1,2)
plot(obs(:,1),uplotf,'k-','LineWidth',3)
axis tight, grid on
ylim([-1 1]*0.5)
ylabel('u_{total}'),xlabel('x')
set(gca,'FontSize',15,'Linewidth',1)

figure(10),clf
plot(obs(:,1),uplotf,'k-','LineWidth',3), hold on
axis tight, grid on
ylim([-1 1]*0.5)
ylabel('u_{total}'),xlabel('x')
set(gca,'FontSize',15,'Linewidth',1)

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
[Ku,~] = compute_disp_stress_kernels_force(combostructure,obs);
uplotf = uplot_0 + Ku*bemsol;

figure(4),clf
subplot(2,1,1)
pcolor(ox,oz,reshape(uplot_0,nz,nx)), shading interp, hold on
contour(ox,oz,reshape(uplot_0,nz,nx),[-1:0.1:1]*0.5,'k-')
axis tight equal
clim([-1 1]*0.5)
cb=colorbar;cb.Label.String = 'u (source)';
ylabel('z'),xlabel('x')
set(gca,'FontSize',15,'Linewidth',1,'TickDir','both')

subplot(2,1,2)
pcolor(ox,oz,reshape(uplotf,nz,nx)), shading interp, hold on
contour(ox,oz,reshape(uplotf,nz,nx),[-1:0.1:1]*0.5,'k-')
axis tight equal
clim([-1 1]*0.5)
cb=colorbar;cb.Label.String = 'u (BEM)';
colormap bluewhitered(20)
ylabel('z'),xlabel('x')
set(gca,'FontSize',15,'Linewidth',1,'TickDir','both')
