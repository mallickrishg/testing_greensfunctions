clear

addpath functions/

% Elastic parameters
G = 1;
% construct mesh for triangle source
A =[0,-0.5];
B = [-1,-2];
C = [2,-2.5];

% dimensions of individual elements
dx = 0.5;
Ntopo = 200;
rcv = construct_triangle(A,B,C,dx,Ntopo);


%% use cuboid kernels to calculate displacement and stress

obs = rcv.xc - rcv.nv.*1e-12;% -ve normal vector means outside the domain

source_strain = 0.5;

% e12 component
[s12,s13]=computeStressAntiplaneTriangleShearZone(obs(:,1),-obs(:,2),A,B,C,1,0,G);
Stress_12 = [s12,s13];

% % e13 component
% [s12,s13]=computeStressAntiplaneTriangleShearZone(obs(:,1),-obs(:,2),A,B,C,0,1,G);
% Stress_13 = [s12,s13];

s12_shz = Stress_12(:,1).*source_strain;
s13_shz = Stress_12(:,2).*source_strain;

% calculate displacements
u_12=computeDisplacementAntiplaneTriangleShearZone(obs(:,1),-obs(:,2),A,B,C,1,0);
u_13=computeDisplacementAntiplaneTriangleShearZone(obs(:,1),-obs(:,2),A,B,C,0,1);
Disp_shz = [u_12,u_13];

u1_shz = Disp_shz(:,1).*source_strain;

%% use BEM to solve the same problem
% observation points
% obs = rcv.xc - rcv.nv.*1e-6;

% compute displacement and stress kernels (for a slip source)
% Disp - [Nobs x Nsources]
% Stress - [Nobs x Nsources x 2], for σ12, σ13
% [Disp_f,Stress_f] = compute_disp_stress_kernels_fault(G,rcv,obs);

% compute displacement and stress kernels (for a traction source)
[Disp_f,Stress_f] = compute_disp_stress_kernels_force(G,rcv,obs);

Tau_f = Stress_f(:,:,1).*rcv.nv(:,1) + Stress_f(:,:,2).*rcv.nv(:,2);

% solve BVP using BEM (boundary element interaction kernel)
K = Tau_f;

% specify right-hand side of BVP
BC = s12_shz.*rcv.nv(:,1) + s13_shz.*-rcv.nv(:,2);
BC(end-Ntopo+1:end) = 0;

% solve linear BEM problem
source_bem = K\BC;

% plot BEM slip distribution
figure(5),clf
subplot(211)
plot3(rcv.xc(1:rcv.N-Ntopo,1),rcv.xc(1:rcv.N-Ntopo,2),source_bem(1:rcv.N-Ntopo),'ko-','LineWidth',2)
box on, grid on
axis tight 
ylabel('y'), xlabel('x'), zlabel('slip')

subplot(212)
plot(source_bem(1:rcv.N-Ntopo),'.-')
ylabel('source strength'), xlabel('node')

%% plot results in the bulk
% calculate displacements and stresses in the medium
nx = 100;
ny = nx/2;

x = linspace(-4,4,nx);
y = linspace(-3,0,ny);
[X,Y] = meshgrid(x,y);
obs_plot = [X(:), Y(:)];

% compute u,σ using BEM
[Disp_f,Stress_f] = compute_disp_stress_kernels_force(G,rcv,obs_plot);
u1_bem_bulk = Disp_f*source_bem;
s12_bem_bulk = Stress_f(:,:,1)*source_bem;
s13_bem_bulk = Stress_f(:,:,2)*source_bem;

% % make s12 discontinuous
% in = inpolygon(obs_plot(:,1),obs_plot(:,2),rcv.xc(:,1),rcv.xc(:,2));
% s12_bem_bulk(in) = s12_bem_bulk(in) - 2*G*source_strain;

% compute u,σ using Barbot 2018 solutions
[s12,s13]=computeStressAntiplaneTriangleShearZone(obs_plot(:,1),-obs_plot(:,2),A,B,C,1,0,G);
Stress_12 = [s12,s13];
s12_shz_bulk = Stress_12(:,1).*source_strain;
s13_shz_bulk = Stress_12(:,2).*source_strain;

% calculate displacements
u_12=computeDisplacementAntiplaneTriangleShearZone(obs_plot(:,1),-obs_plot(:,2),A,B,C,1,0);
u_13=computeDisplacementAntiplaneTriangleShearZone(obs_plot(:,1),-obs_plot(:,2),A,B,C,0,1);
Disp_shz = [u_12,u_13];
u1_shz_bulk = Disp_shz(:,1).*source_strain;

figure(100),clf
subplot(3,2,1)
pcolor(x,y,reshape(u1_bem_bulk,ny,nx)), shading interp
hold on
plot([rcv.x1,rcv.x2]',[rcv.y1,rcv.y2]','k.-')
axis tight equal
colorbar
clim([-1 1])
xlabel('x'), ylabel('y')
title('u_1 (BEM)')

subplot(3,2,3)
pcolor(x,y,reshape(s12_bem_bulk,ny,nx)), shading interp
hold on
plot([rcv.x1,rcv.x2]',[rcv.y1,rcv.y2]','k.-')
axis tight equal
clim([-1 1])
colorbar
xlabel('x'), ylabel('y')
title('\sigma_{12} (BEM)')

subplot(3,2,5)
pcolor(x,y,-reshape(s13_bem_bulk,ny,nx)), shading interp
hold on
plot([rcv.x1,rcv.x2]',[rcv.y1,rcv.y2]','k.-')
axis tight equal
clim([-1 1])
colorbar
xlabel('x'), ylabel('y')
title('\sigma_{13} (BEM)')

subplot(3,2,2)
pcolor(x,y,reshape(u1_shz_bulk,ny,nx)), shading interp
hold on
plot([rcv.x1,rcv.x2]',[rcv.y1,rcv.y2]','k.-')
axis tight equal
colorbar
clim([-1 1])
xlabel('x'), ylabel('y')
title('u_1 (L&B16)')

subplot(3,2,4)
pcolor(x,y,reshape(s12_shz_bulk,ny,nx)), shading interp
hold on
plot([rcv.x1,rcv.x2]',[rcv.y1,rcv.y2]','k.-')
axis tight equal
clim([-1 1])
colorbar
xlabel('x'), ylabel('y')
title('\sigma_{12} (L&B16)')

subplot(3,2,6)
pcolor(x,y,reshape(s13_shz_bulk,ny,nx)), shading interp
hold on
plot([rcv.x1,rcv.x2]',[rcv.y1,rcv.y2]','k.-')
axis tight equal
clim([-1 1])
colorbar
xlabel('x'), ylabel('y')
title('\sigma_{13} (L&B16)')

colormap bluewhitered(20)
set(findobj(gcf,'type','axes'),'FontSize',15,'LineWidth', 1,'XLim',[-4,4]);




