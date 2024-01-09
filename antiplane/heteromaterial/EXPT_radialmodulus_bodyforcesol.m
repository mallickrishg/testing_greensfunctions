% script to solve bimaterial quasi-static elasticity problem
% with equivalent body force solution for multiple horizontal interfaces

clear

% bimaterial mesh parameters
Lmesh = 20;
Nmesh = 320;
bimaterial_z = 0.25*[-1 1];
% ratio of shear moduli
mu_ratio = 3;

% construct source
rcv = create_horizontalinterfaces(-0.5,0.5,2,0);
% construct bimaterial interface
bimat = create_horizontalinterfaces(-Lmesh/2,Lmesh/2,Nmesh+1,bimaterial_z);

%% compute source term at the interface
source = ones(rcv.N,1);

[Disp,Dgradient] = compute_disp_stress_kernels_fault(rcv,[bimat.x2c,bimat.x3c]);

% we only need u and du/dz for continuity across the interface
u0 = Disp*source;
% ux0 = Dgradient(:,:,1)*source;
uz0 = Dgradient(:,:,2)*source;

%% solve same problem using body force kernels
[Ku,Dgradient] = compute_disp_stress_kernels_force(bimat,[bimat.x2c,bimat.x3c]);
Kuz = Dgradient(:,:,2);

% the body force equivalent due to a change in shear modulus is only at the
% interface and is a term fb = (dlogμ/dz)(du/dz)
% du/dz = du0/dz + Kuz.fb
% we can solve for the unknown fb as 
% fb = (dlogμ/dz)(du0/dz) + (dlogμ/dz)(Kuz.fb)
% simplifying: 
% fb = [I - (dlogμ/dz)Kuz]\(dlogμ/dz)(du0/dz)

alpha = ones(bimat.N,1);
index = bimat.x3c>= 0;
alpha(index) = (1-mu_ratio)/(mu_ratio);
alpha(~index) = -(1-mu_ratio);

fb = (eye(bimat.N) - repmat(alpha,1,bimat.N).*Kuz)\(alpha.*uz0);

figure(1),clf
stem(bimat.x2c,fb,'r-','Linewidth',1,'MarkerFaceColor','r'), hold on
for i = 1:bimat.N
    toplotx = bimat.x2(i)+bimat.W(i)*cosd(bimat.dip(i));
    plot([bimat.x2(i),bimat.x2(i),toplotx,toplotx],fb(i).*[0,1,1,0],'k-','Linewidth',1)
end
axis tight, grid on
xlim([-1 1]*5)
xlabel('x'), ylabel('body force strength')
set(gca,'FontSize',15,'Linewidth',1)

%% compute u and du/dx,du/dz in the medium
Nobs = 1000;
obs = [zeros(Nobs,1),linspace(-3,3,Nobs)'];
% provide shear modulus structure
muplot = ones(Nobs,1).*mu_ratio;
muplot(obs(:,2)>bimaterial_z(2) | obs(:,2)<bimaterial_z(1)) = 1;

% source contribution
[Disp,Dgradient0] = compute_disp_stress_kernels_fault(rcv,obs);
uplot_0 = Disp*source;

% compute solution from body force equivalent
[Ku,Dgradient] = compute_disp_stress_kernels_force(bimat,obs);
uplotf = uplot_0 + Ku*fb;
e13f = Dgradient0(:,:,2)*source + Dgradient(:,:,2)*fb;
s13f = muplot.*(e13f);

figure(2),clf
plot(obs(:,2),s13f,'r-','LineWidth',2)
grid on, axis tight
ylabel('\sigma_{13}'),xlabel('z')
set(gca,'FontSize',15,'Linewidth',1)

figure(3),clf
subplot(2,1,1)
plot(obs(:,2),uplot_0,'Linewidth',2), hold on
plot(obs(:,2),Ku*fb,'r-','LineWidth',2)
axis tight, grid on
ylabel('u (individual contribution)')
% ylim([-1 1]*0.5)
set(gca,'FontSize',15,'Linewidth',1)
subplot(2,1,2)
plot(obs(:,2),uplotf,'k-','LineWidth',2)
axis tight, grid on
ylabel('u_{total}'),xlabel('z')
set(gca,'FontSize',15,'Linewidth',1)
