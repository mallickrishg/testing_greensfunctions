% script to solve heterogeneous material quasi-static half-space elasticity problem
% with equivalent body force solution for vertical horizontal interfaces
% 
% AUTHOR:
% Rishav Mallick, JPL, 2024

clear

% bimaterial mesh parameters
Lmesh = 20;
Nmesh = Lmesh*4;
bimaterial_x = 0.5*[-1,-0.5,0.5,1];
% topo mesh
Ltopo = 20;
Ntopo = Ltopo*64;

% shear moduli structure
mu_structure = [1,1/2,1/4,1/2,1];

% construct source
rcv = create_verticalinterfaces(-20,-0.5,2,0);
% construct bimaterial interface
bimat = create_verticalinterfaces(-Lmesh,0,Nmesh+1,bimaterial_x);
% construct free surface
topo = create_horizontalinterfaces(-Ltopo/2,Ltopo/2,Ntopo+1,0);

% join bimat+topo
indextopo = [false(bimat.N,1);true(topo.N,1)];
combostructure = append_structures(bimat,topo);

%% compute source term at the interfaces
source = ones(rcv.N,1);
deltashift = -1e-5;

[Disp,Dgradient] = compute_disp_stress_kernels_fault(rcv,[combostructure.x2c + deltashift.*indextopo.*combostructure.nvec(:,1),...
                                                          combostructure.x3c + deltashift.*indextopo.*combostructure.nvec(:,2)]);

u0 = Disp*source;
ux0 = Dgradient(:,:,1)*source;
uz0 = Dgradient(:,:,2)*source;

%% solve same problem using body force kernels
% deltashift = -1e-5;
[Ku,Dgradient] = compute_disp_stress_kernels_force(combostructure,...
                                                         [combostructure.x2c + deltashift.*indextopo.*combostructure.nvec(:,1),...
                                                          combostructure.x3c + deltashift.*indextopo.*combostructure.nvec(:,2)]);
Kux = Dgradient(:,:,1);
Kuz = Dgradient(:,:,2);

% the body force equivalent due to a change in shear modulus is only at the
% interface and is a term fb = (dlogμ/dx)(du/dx)
% du/dx = du0/dx + Kx_bb.fb + Kx_bt.ft
% we can solve for the unknown fb as 
% fb = αx(du0/dx) + (αx)(Kx_bb.fb + Kx_bt.ft)
% we also need to solve for fictitious forces on the topo interface

alpha = ones(bimat.N,1);
for i = 1:length(bimaterial_x) 
    index = (bimat.x2c == bimaterial_x(i));
    alpha(index) = (mu_structure(i)-mu_structure(i+1))/(mu_structure(i+1));
end
% alpha(~index) = -(1-mu_ratio);

bigmat = [(eye(bimat.N) - repmat(alpha,1,bimat.N).*Kux(~indextopo,~indextopo)), -repmat(alpha,1,topo.N).*Kux(~indextopo,indextopo);...
           Kuz(indextopo,~indextopo)                                          , Kuz(indextopo,indextopo)];
rhsvec = [(alpha.*ux0(~indextopo));-uz0(indextopo)];
% solve matrix equations
bemsol = bigmat\rhsvec;
fb = bemsol(~indextopo);
ft = bemsol(indextopo);

%
figure(1),clf
subplot(2,1,1)
stem(bimat.x3c,fb,'-','Linewidth',1,'MarkerFaceColor','k')
axis tight, grid on
% xlim([-3 0])
xlabel('z'), ylabel('body force strength')
set(gca,'FontSize',15,'Linewidth',1)

subplot(2,1,2)
stem(topo.x2c,ft,'r-','Linewidth',1,'MarkerFaceColor','r')
axis tight, grid on
xlabel('x'), ylabel('topo force strength')
set(gca,'FontSize',15,'Linewidth',1)
xlim([-1 1]*2)
%% compute displacement at the free surface

Nobs = 1000;
obs = [linspace(-2,2,Nobs)',zeros(Nobs,1)-1e-1];
% provide shear modulus structure
muplot = ones(Nobs,1).*mu_structure;
muplot(obs(:,1)>bimaterial_x(2) | obs(:,1)<bimaterial_x(1)) = 1;

% source contribution
[Disp,Dgradient0] = compute_disp_stress_kernels_fault(rcv,obs);
uplot_0 = Disp*source;

% compute solution from body force equivalent
[Ku,Dgradient] = compute_disp_stress_kernels_force(combostructure,obs);
uplot_bem = uplot_0 + Ku*bemsol;
e12_bem = Dgradient0(:,:,1)*source + Dgradient(:,:,1)*bemsol;
e13_bem = Dgradient0(:,:,2)*source + Dgradient(:,:,2)*bemsol;

figure(2),clf
subplot(2,1,1)
plot(obs(:,1),uplot_0,'LineWidth',4), hold on
plot(obs(:,1),Ku(:,indextopo)*ft + uplot_0,'LineWidth',3)
plot(obs(:,1),Ku(:,~indextopo)*fb,'LineWidth',2)
plot(obs(:,1),uplot_bem,'k-','LineWidth',3)
axis tight, grid on
legend('full space','half-space','interface','hetero BEM')
xlabel('x'), ylabel('displacement')
% ylim([-1 1]*0.45)
set(gca,'FontSize',15,'Linewidth',1)

subplot(2,1,2)
plot(obs(:,1),uplot_bem,'k-','LineWidth',3),hold on
plot(obs(:,1),Ku(:,indextopo)*ft + uplot_0,'LineWidth',3)
axis tight, grid on
xlabel('x'), ylabel('displacement')
ylim([-1 1]*0.45)
set(gca,'FontSize',15,'Linewidth',1)

figure(3),clf
subplot(2,1,1)
plot(obs(:,1),uplot_bem,'k-','LineWidth',3)
axis tight, grid on
xlabel('x'), ylabel('displacement')
ylim([-1 1]*0.45)
set(gca,'FontSize',15,'Linewidth',1)

subplot(2,1,2)
plot(obs(:,1),e12_bem,'-','LineWidth',2),hold on
plot(obs(:,1),e13_bem,'-','LineWidth',2)
axis tight, grid on
xlabel('x'), ylabel('strain')
set(gca,'FontSize',15,'Linewidth',1)


