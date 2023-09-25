clear
addpath functions/

% Elastic parameters
G = 1;
% construct mesh for triangle source
A =[0,-2];

% dimensions of box
L_x = 2*abs(x0);
L_y = 2*abs(y0);
dx = 0.01;
rcv = construct_box(x0,y0,L_x,L_y,dx);


%% use cuboid kernels to calculate displacement and stress

obs = rcv.xc - rcv.nv.*1e-12;% -ve normal vector means outside the domain

source_strain = 0.5;

% [Stress_12,~] = calc_stressgreensfunctions_antiplaneshz(G,...
%                         obs(:,1) - (x0 + L_x/2),obs(:,2),...
%                         y0,L_x,L_y);
[Stress_12,Stress_13]=computeStressAntiplaneTriangleShearZone(obs(:,1),obs(:,2),A,B,C,1,0,G);

s12_shz = Stress_12(:,1).*source_strain;
s13_shz = Stress_12(:,2).*source_strain;
% s12_shz = Stress_13(:,1).*source_strain;
% s13_shz = Stress_13(:,2).*source_strain;

% calculate displacements
[Disp_shz] = calc_dispgreensfunctions_antiplaneshz(obs(:,1) - (x0 + L_x/2),obs(:,2),...
                                               y0,L_x,L_y);
u1_shz = Disp_shz(:,1).*source_strain;
% u1_shz = Disp_shz(:,2).*source_strain;

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







