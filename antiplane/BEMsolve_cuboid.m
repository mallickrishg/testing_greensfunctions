clear
addpath functions/

% Elastic parameters
G = 1;
% construct mesh for cuboid boundary
% bottom left corner
x0 = -1;
y0 = -1.5;

% dimensions of box
L_x = 2;
L_y = 3;
dx = 0.1;
rcv = construct_box(x0,y0,L_x,L_y,dx);
left = rcv.xc(:,1) == x0;
right = rcv.xc(:,1) == x0 + L_x;
bot = rcv.xc(:,2) == y0;
top = rcv.xc(:,2) == y0 + L_y;

% observation points
obs = rcv.xc + rcv.nv*1e-12;

% compute displacement and stress kernels (for a slip source)
% Disp - [Nobs x Nsources]
% Stress - [Nobs x Nsources x 2], for σ12, σ13
[Disp_f,Stress_f] = compute_disp_stress_kernels_fault(G,rcv,obs);

Tau_f = Stress_f(:,:,1).*rcv.nv(:,1) + Stress_f(:,:,2).*rcv.nv(:,2);

% compute displacement and stress kernels (for a traction source)
[Disp_t,Stress_t] = compute_disp_stress_kernels_force(G,rcv,obs);

Tau_t = Stress_t(:,:,1).*rcv.nv(:,1) + Stress_t(:,:,2).*rcv.nv(:,2);
% Tau_t = Stress_t(:,:,1);

% solve BVP using BEM
K = zeros(rcv.N,rcv.N);

% assemble kernel
K(left,:) = Disp_t(left,:);
K(right,:) = Disp_t(right,:);
K(bot,:) = Tau_t(bot,:);
K(top,:) = Tau_t(top,:);

% specify right-hand side of BVP
BC = zeros(rcv.N,1);
% BC(left) = -(linspace(-1,0,length(find(left))).^2);
% BC(right) = linspace(0,1,length(find(right))).^2;
% linspace(-1,-1,length(find(left)))
BC(left) = -.5;
BC(right) = -.5;
BC(top) = 0;
BC(bot) = 0;

% solve linear problem
slip_bem = K\BC;
% slip_bem = ones(rcv.N,1);
% slip_bem(left|right|bot) = 0;

figure(100),clf
plot3(rcv.xc(:,1),rcv.xc(:,2),slip_bem,'k.-','LineWidth',2)
box on, grid on
axis tight equal

figure(1),clf
scatter(obs(:,1),obs(:,2),100,Disp_t*slip_bem,'o','filled','MarkerEdgeColor','k')
axis equal
box on
colorbar
title('u_1')
colormap bluewhitered

figure(2),clf
subplot(2,1,1)
scatter(obs(:,1),obs(:,2),100,Stress_t(:,:,1)*slip_bem,'o','filled','MarkerEdgeColor','k')
axis equal
box on
colorbar
title('\sigma_{12}')
clim([-1 1])

subplot(2,1,2)
scatter(obs(:,1),obs(:,2),100,Stress_t(:,:,2)*slip_bem,'o','filled','MarkerEdgeColor','k')
axis equal
box on
colorbar
clim([-1 1])
title('\sigma_{13}')
colormap bluewhitered

% calculate displacements and stresses in the medium
nx = 100;
ny = nx/2;

x = linspace(-6,6,nx);
y = linspace(-3,3,ny);
[X,Y] = meshgrid(x,y);

obs_plot = [X(:), Y(:)];
[Disp_f,Stress_f] = compute_disp_stress_kernels_fault(G,rcv,obs_plot);
[Disp_t,Stress_t] = compute_disp_stress_kernels_force(G,rcv,obs_plot);
u1_bem = Disp_t*slip_bem;
s12_bem = Stress_t(:,:,1)*slip_bem;
s13_bem = Stress_t(:,:,2)*slip_bem;

figure(3),clf
subplot(3,1,1)
pcolor(x,y,reshape(u1_bem,ny,nx)), shading interp
hold on
contour(x,y,reshape(u1_bem,ny,nx),[-1:0.2:1],'k-','LineWidth',1)
plot([rcv.x1,rcv.x2],[rcv.y1,rcv.y2],'k.-')
axis tight equal
colorbar
clim([-1 1])
xlabel('x'), ylabel('y')
title('u_z (BEM)')

subplot(3,1,2)
pcolor(x,y,reshape(s12_bem,ny,nx)), shading interp
hold on
plot([rcv.x1,rcv.x2],[rcv.y1,rcv.y2],'k.-')
axis tight equal
clim([-1 1])
colorbar
xlabel('x'), ylabel('y')
title('\sigma_{xz}')

subplot(3,1,3)
pcolor(x,y,reshape(s13_bem,ny,nx)), shading interp
hold on
plot([rcv.x1,rcv.x2],[rcv.y1,rcv.y2],'k.-')
axis tight equal
clim([-1 1])
colorbar
xlabel('x'), ylabel('y')
title('\sigma_{yz}')

colormap bluewhitered(20)
set(findobj(gcf,'type','axes'),'FontSize',15,'LineWidth', 1);

% %% use cuboid kernels to calculate displacement and stress
% % calculate stresses from a shear zone
% source_strain = (L_x+L_y)/L_x/L_y;
% 
% [Stress_12,Stress_13] = calc_stressgreensfunctions_antiplaneshz(G,...
%                         obs_plot(:,1) - (x0 + L_x/2),obs_plot(:,2),...
%                         y0,L_x,L_y);
% s12_shz = Stress_12(:,1).*source_strain;
% s13_shz = Stress_12(:,2).*source_strain;
% 
% % calculate displacements
% [Disp_f] = calc_dispgreensfunctions_antiplaneshz(obs_plot(:,1) - (x0 + L_x/2),obs_plot(:,2),...
%                                                y0,L_x,L_y);
% u1_shz = Disp_f(:,1).*source_strain;
% 
% figure(11),clf
% % scatter(obs_plot(:,1),obs_plot(:,2),500,u1_shz,'o','filled','MarkerEdgeColor','k')
% % axis tight equal
% % box on
% % colorbar
% % colormap bluewhitered
% 
% subplot(3,1,1)
% pcolor(x,y,reshape(u1_shz,ny,nx)), shading interp
% hold on
% plot([rcv.x1,rcv.x2],[rcv.y1,rcv.y2],'k.-')
% axis tight equal
% colorbar
% clim([-1 1])
% xlabel('x'), ylabel('y')
% title('u_z (Lambert & Barbot, 2016)')
% 
% subplot(3,1,2)
% pcolor(x,y,reshape(s12_shz,ny,nx)), shading interp
% hold on
% plot([rcv.x1,rcv.x2],[rcv.y1,rcv.y2],'k.-')
% axis tight equal
% clim([-1 1])
% colorbar
% xlabel('x'), ylabel('y')
% title('\sigma_{xz}')
% 
% subplot(3,1,3)
% pcolor(x,y,reshape(s13_shz,ny,nx)), shading interp
% hold on
% plot([rcv.x1,rcv.x2],[rcv.y1,rcv.y2],'k.-')
% axis tight equal
% clim([-1 1])
% colorbar
% xlabel('x'), ylabel('y')
% title('\sigma_{yz}')
% 
% colormap bluewhitered(20)
% set(findobj(gcf,'type','axes'),'FontSize',15,'LineWidth', 1);
% 
% % residuals
% figure(12),clf
% subplot(3,1,1)
% pcolor(x,y,reshape(u1_shz-u1_bem,ny,nx)), shading interp
% axis tight equal
% colorbar
% clim([-1 1]*0.5)
% 
% subplot(3,1,2)
% pcolor(x,y,reshape(s12_shz-s12_bem,ny,nx)), shading interp
% axis tight equal
% clim([-1 1]*0.5)
% colorbar
% 
% subplot(3,1,3)
% pcolor(x,y,reshape(s13_shz-s13_bem,ny,nx)), shading interp
% axis tight equal
% clim([-1 1]*0.5)
% colorbar
% 
% colormap bluewhitered(100)
