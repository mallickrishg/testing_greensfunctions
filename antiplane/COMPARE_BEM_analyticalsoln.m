clear
addpath functions/

% Elastic parameters
G = 1;
% construct mesh for cuboid boundary
% bottom left corner
x0 = -0.5;
y0 = -1;

% dimensions of box
L_x = 2*abs(x0);
L_y = 2*abs(y0);
dx = 0.05;
rcv = construct_box(x0,y0,L_x,L_y,dx);
left = rcv.xc(:,1) == x0;
right = rcv.xc(:,1) == x0 + L_x;
bot = rcv.xc(:,2) == y0;
top = rcv.xc(:,2) == y0 + L_y;

%% use cuboid kernels to calculate displacement and stress
% calculate stresses from a shear zone
% obs = rcv.xc;
obs = rcv.xc - rcv.nv.*1e-9;% -ve normal vector means outside the domain

source_strain = 0.7;

[Stress_12,~] = calc_stressgreensfunctions_antiplaneshz(G,...
                        obs(:,1) - (x0 + L_x/2),obs(:,2),...
                        y0,L_x,L_y);
s12_shz = Stress_12(:,1).*source_strain;
s13_shz = Stress_12(:,2).*source_strain;

% calculate displacements
[Disp_shz] = calc_dispgreensfunctions_antiplaneshz(obs(:,1) - (x0 + L_x/2),obs(:,2),...
                                               y0,L_x,L_y);
u1_shz = Disp_shz(:,1).*source_strain;

figure(11),clf
scatter(obs(:,1),obs(:,2),100,u1_shz,'o','filled','MarkerEdgeColor','k')
axis equal
box on
colorbar
title('u_1 (L&B16)')
clim([-1 1]*max(abs(u1_shz(:))))
colormap bluewhitered(20)
set(gca,'Fontsize',20,'Linewidth',1)

figure(12),clf
subplot(2,1,1)
scatter(obs(:,1),obs(:,2),100,s12_shz,'o','filled','MarkerEdgeColor','k')
axis equal
box on
colorbar
clim([-1 1]*max(abs(s12_shz(:))))
title('\sigma_{12} (L&B16)')

subplot(2,1,2)
scatter(obs(:,1),obs(:,2),100,s13_shz,'o','filled','MarkerEdgeColor','k')
axis equal
box on
colorbar
clim([-1 1]*max(abs(s13_shz(:))))
title('\sigma_{13} (L&B16)')

colormap bluewhitered(20)
set(findobj(gcf,'type','axes'),'FontSize',20,'LineWidth', 1);

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

% solve BVP using BEM
K = zeros(rcv.N,rcv.N);

% assemble kernel
K(left,:) = Tau_f(left,:);
K(right,:) = Tau_f(right,:);
K(bot,:) = Tau_f(bot,:);
K(top,:) = Tau_f(top,:);

% specify right-hand side of BVP
BC = zeros(rcv.N,1);
BC_tau = s12_shz.*rcv.nv(:,1) + s13_shz.*rcv.nv(:,2);
BC_disp = u1_shz;
BC(left|right) = BC_tau(left|right);
BC(top|bot) = BC_tau(top|bot);

% solve linear BEM problem
slip_bem = pinv(K)*BC;

% plot slip distribution
figure(5),clf
subplot(211)
plot3(rcv.xc(:,1),rcv.xc(:,2),slip_bem,'k.-','LineWidth',2)
box on, grid on
axis tight 
ylabel('y'), xlabel('x'), zlabel('slip')

subplot(212)
plot(slip_bem,'.-')
ylabel('slip'), xlabel('node')

% plot displacements and stresses on the boundary
figure(1),clf
scatter(obs(:,1),obs(:,2),100,Disp_f*slip_bem,'o','filled','MarkerEdgeColor','k')
axis equal
box on
colorbar
title('u_1 (BEM)')
colormap bluewhitered(20)
set(gca,'Fontsize',20,'Linewidth',1)

figure(2),clf
subplot(2,1,1)
scatter(obs(:,1),obs(:,2),100,Stress_f(:,:,1)*slip_bem,'o','filled','MarkerEdgeColor','k')
axis equal
box on
colorbar
title('\sigma_{12}')
clim([-1 1])

subplot(2,1,2)
scatter(obs(:,1),obs(:,2),100,Stress_f(:,:,2)*slip_bem,'o','filled','MarkerEdgeColor','k')
axis equal
box on
colorbar
clim([-1 1])
title('\sigma_{13}')
colormap bluewhitered(20)
set(findobj(gcf,'type','axes'),'FontSize',20,'LineWidth', 1);

%% compare displacements & stresses at the boundary
figure(20),clf
subplot(311)
plot(u1_shz,'o-'), hold on
plot(Disp_f*slip_bem,'r.-')
legend('L&B16','bem')
ylabel('u_1'), xlabel('node')

subplot(312)
plot(s12_shz,'o-'), hold on
plot(Stress_f(:,:,1)*slip_bem,'r.-')
ylabel('\sigma_{12}'), xlabel('node')

subplot(313)
plot(s13_shz,'o-'), hold on
plot(Stress_f(:,:,2)*slip_bem,'r.-')
ylabel('\sigma_{13}'), xlabel('node')

%% plot results in the bulk
% calculate displacements and stresses in the medium
nx = 100;
ny = nx/2;

x = linspace(-6,6,nx);
y = linspace(-3,3,ny);
[X,Y] = meshgrid(x,y);
obs_plot = [X(:), Y(:)];

% compute u,σ using BEM
[Disp_f,Stress_f] = compute_disp_stress_kernels_force(G,rcv,obs_plot);
u1_bem_bulk = Disp_f*slip_bem;
s12_bem_bulk = Stress_f(:,:,1)*slip_bem;
s13_bem_bulk = Stress_f(:,:,2)*slip_bem;

% make s12 discontinuous
in = inpolygon(obs_plot(:,1),obs_plot(:,2),rcv.xc(:,1),rcv.xc(:,2));
s12_bem_bulk(in) = s12_bem_bulk(in) - 2*G*source_strain;

% compute u,σ using L&B 2016 solutions
[Stress_12,~] = calc_stressgreensfunctions_antiplaneshz(G,...
                        obs_plot(:,1) - (x0 + L_x/2),obs_plot(:,2),...
                        y0,L_x,L_y);
s12_shz_bulk = Stress_12(:,1).*source_strain;
s13_shz_bulk = Stress_12(:,2).*source_strain;
% calculate displacements
[Disp_f] = calc_dispgreensfunctions_antiplaneshz(obs_plot(:,1) - (x0 + L_x/2),obs_plot(:,2),...
                                               y0,L_x,L_y);
u1_shz_bulk = Disp_f(:,1).*source_strain;

figure(100),clf
subplot(3,2,1)
pcolor(x,y,reshape(u1_bem_bulk,ny,nx)), shading interp
hold on
plot([rcv.x1,rcv.x2],[rcv.y1,rcv.y2],'k.-')
axis tight equal
colorbar
clim([-1 1])
xlabel('x'), ylabel('y')
title('u_1 (BEM)')

subplot(3,2,3)
pcolor(x,y,reshape(s12_bem_bulk,ny,nx)), shading interp
hold on
plot([rcv.x1,rcv.x2],[rcv.y1,rcv.y2],'k.-')
axis tight equal
clim([-1 1])
colorbar
xlabel('x'), ylabel('y')
title('\sigma_{12} (BEM)')

subplot(3,2,5)
pcolor(x,y,reshape(s13_bem_bulk,ny,nx)), shading interp
hold on
plot([rcv.x1,rcv.x2],[rcv.y1,rcv.y2],'k.-')
axis tight equal
clim([-1 1])
colorbar
xlabel('x'), ylabel('y')
title('\sigma_{13} (BEM)')

subplot(3,2,2)
pcolor(x,y,reshape(u1_shz_bulk,ny,nx)), shading interp
hold on
plot([rcv.x1,rcv.x2],[rcv.y1,rcv.y2],'k.-')
axis tight equal
colorbar
clim([-1 1])
xlabel('x'), ylabel('y')
title('u_1 (L&B16)')

subplot(3,2,4)
pcolor(x,y,reshape(s12_shz_bulk,ny,nx)), shading interp
hold on
plot([rcv.x1,rcv.x2],[rcv.y1,rcv.y2],'k.-')
axis tight equal
clim([-1 1])
colorbar
xlabel('x'), ylabel('y')
title('\sigma_{12} (L&B16)')

subplot(3,2,6)
pcolor(x,y,reshape(s13_shz_bulk,ny,nx)), shading interp
hold on
plot([rcv.x1,rcv.x2],[rcv.y1,rcv.y2],'k.-')
axis tight equal
clim([-1 1])
colorbar
xlabel('x'), ylabel('y')
title('\sigma_{13} (L&B16)')

colormap bluewhitered(20)
set(findobj(gcf,'type','axes'),'FontSize',15,'LineWidth', 1);





