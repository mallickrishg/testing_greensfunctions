% script to solve bimaterial quasi-static elasticity problem and compare
% with equivalent body force solution

clear

% bimaterial mesh parameters
Lmesh = 80;
Nmesh = 320*4;
bimaterial_z = 0.;
% ratio of shear moduli
mu_ratio = 3;

% construct source
rcv = create_horizontalinterfaces(-0.5,0.5,2,0);
% construct bimaterial interface
bimat = create_horizontalinterfaces(-Lmesh/2,Lmesh/2,Nmesh+1,bimaterial_z);

% index for overlap between rcv & bimat
% overlapindex = bimat.x2 >= rcv.x2 & bimat.x2+bimat.W <=rcv.x2+rcv.W & bimat.x3c == rcv.x3c;

%% compute source term at the interface
source = ones(rcv.N,1);

[Disp,Dgradient] = compute_disp_stress_kernels_fault(rcv,[bimat.x2c,bimat.x3c]);

% we only need u and du/dz for continuity across the interface
u0 = Disp*source;
% ux0 = Dgradient(:,:,1)*source;
uz0 = Dgradient(:,:,2)*source;

%% compute self interaction kernel for interface and solve matrix equations
[Ku,Dgradient] = compute_disp_stress_kernels_fault(bimat,[bimat.x2c,bimat.x3c]);
Kuz = Dgradient(:,:,2);

% enforce continuity of 'u' and 'du/dz' as matrix equations
% u0_l + Ku_l*s_l = u0_r + Ku_r*s_r
% u0 - Ku*s_l = u0 + Ku_r*s_r (-ve sign for Ku_l = -Ku because of opposite direction evaluation)
% this condition typically implies s_l = -s_r
% 
% G_l*(uz0_l + Kuz_l*s_l) = G_r*(uz0_r + Kuz_r*s_r)
% uz0 + Kuz*s_l = G_r/G_l*uz0 + G_r/G_l*Kuz*s_r
% 
% rewrite as:
% [Ku, Ku].[s_l;s_r] = 0 
% [Kuz, -G_r/G_l*Kuz].[s_l;s_r] = (G_r/G_l - 1)*uz0

bigK = [Ku, Ku;...
        Kuz, -mu_ratio*Kuz];

rhs = [zeros(bimat.N,1);(mu_ratio-1).*uz0];

% solve for fictitious slips
bemsol = bigK\rhs;
bemsol_l = bemsol(1:end/2);
bemsol_r = bemsol(1+end/2:end);

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

alpha = (1-mu_ratio)/(mu_ratio)*ones(bimat.N,1);
fb = (eye(bimat.N) - repmat(alpha,1,bimat.N).*Kuz)\(alpha.*uz0);

figure(1),clf
stem(bimat.x2c,bemsol_l,'b-','Linewidth',2), hold on
stem(bimat.x2c,bemsol_r,'-','Linewidth',1)
stem(bimat.x2c,fb,'k-','Linewidth',1)
axis tight, grid on
xlim([-1 1]*5)

%% compute u and du/dx,du/dz in the medium
Nobs = 1000;
obs = [zeros(Nobs,1)+min(abs(bimat.x2c)),linspace(-3,3,Nobs)'];

% source contribution
[Disp,Dgradient0] = compute_disp_stress_kernels_fault(rcv,obs);
uplot_0 = Disp*source;
s13plot = zeros(Nobs,1);

% interface contributions
[Disp,Dgradient] = compute_disp_stress_kernels_fault(bimat,obs);

index = bimaterial_z <= obs(:,2);
uplot = uplot_0.*0;
uplot(index) = Disp(index,:)*bemsol_r;
s13plot(index) = mu_ratio.*(Dgradient0(index,:,2)*source + Dgradient(index,:,2)*bemsol_r);

index = bimaterial_z > obs(:,2);
uplot(index) = Disp(index,:)*bemsol_l;
s13plot(index) = (Dgradient0(index,:,2)*source + Dgradient(index,:,2)*bemsol_l);

% total displacement (source + fictitious slip)
toplot = (uplot + uplot_0);

figure(2),clf
subplot(2,1,1)
plot(obs(:,2),uplot_0,'Linewidth',2), hold on
plot(obs(:,2),uplot,'Linewidth',2)
grid on, axis tight
legend('source','interface')
xlabel('z'),ylabel('u')
set(gca,'FontSize',15,'Linewidth',1)
subplot(2,1,2)
plot(obs(:,2),uplot_0,'Linewidth',1), hold on
plot(obs(:,2),toplot,'r-','LineWidth',2)
axis tight, grid on
legend('source','total')
xlabel('z'),ylabel('u')
set(gca,'FontSize',15,'Linewidth',1)

% compute solution from body force equivalent
[Ku,Dgradient] = compute_disp_stress_kernels_force(bimat,obs);
uplotf = uplot_0 + Ku*fb;
s13f = zeros(Nobs,1);
e13f = Dgradient0(:,:,2)*source + Dgradient(:,:,2)*fb;
index = obs(:,2) >= bimaterial_z;
s13f(index) = mu_ratio.*(e13f(index));
s13f(~index) = e13f(~index);

figure(3),clf
plot(obs(:,2),s13plot,'b-','Linewidth',2), hold on
plot(obs(:,2),s13f,'r-','LineWidth',2)
grid on, axis tight
legend('slip','force')
ylabel('\sigma_{13}'),xlabel('z')
set(gca,'FontSize',15,'Linewidth',1)
% ylim([0 1])


figure(4),clf
subplot(2,1,1)
plot(obs(:,2),uplot,'r-','Linewidth',2), hold on
plot(obs(:,2),Ku*fb,'k-','LineWidth',2)
axis tight, grid on
legend('slip','force')
ylabel('u (interface contribution)')
ylim([0 0.5])
set(gca,'FontSize',15,'Linewidth',1)
subplot(2,1,2)
plot(obs(:,2),toplot,'r-','LineWidth',2), hold on
plot(obs(:,2),uplotf,'k-','LineWidth',2)
axis tight, grid on
ylabel('u_{total}'),xlabel('z')
set(gca,'FontSize',15,'Linewidth',1)

return

%% plot 2d bulk deformation
nx = 70;
nz = 70;
xplot = linspace(-2,2,nx);
zplot = linspace(-2,2,nz);
[X,Z] = meshgrid(xplot,zplot);
obs = [X(:),Z(:)];

% source contribution
[Disp,Dgradient0] = compute_disp_stress_kernels_fault(rcv,obs);
uplot_0 = Disp*source;
s13plot = zeros(nx*nz,1);

% interface contributions
[Disp,Dgradient] = compute_disp_stress_kernels_fault(bimat,obs);

index = bimaterial_z <= obs(:,2);
uplot = uplot_0.*0;
uplot(index) = Disp(index,:)*bemsol_r;
s13plot(index) = mu_ratio.*(Dgradient0(index,:,2)*source + Dgradient(index,:,2)*bemsol_r);

index = bimaterial_z > obs(:,2);
uplot(index) = Disp(index,:)*bemsol_l;
s13plot(index) = (Dgradient0(index,:,2)*source + Dgradient(index,:,2)*bemsol_l);

% plot displacement field
% total displacement (source + fictitious slip)
toplot = reshape(uplot + uplot_0,nz,nx);

figure(4),clf
pcolor(xplot,zplot,toplot), shading interp
hold on
contour(xplot,zplot,toplot,20,'k-')
axis tight equal
colorbar
clim([-1 1]*1)
colormap(bluewhitered(1000))

% plot stress field field
figure(5),clf
toplot = reshape(s13plot,nz,nx);
pcolor(xplot,zplot,toplot), shading interp
hold on
contour(xplot,zplot,toplot,20,'k-')
axis tight equal
colorbar
clim([-1 1]*1)
colormap(turbo(100))



