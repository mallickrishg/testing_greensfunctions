clear
addpath ~/Dropbox/scripts/utils/

% compute stresses in a full space
G = 1;

% rcv = construct_box(-1,-1,2,2,0.5);

% source location
x2_bot = 0;
x3_bot = -1; % flip signs
dip = 90;
% index = 12;
% x2_bot = rcv.x2(index);
% x3_bot = rcv.y2(index);
% dip = rcv.dip(index);

% dimensions of fault/shearzone
% L_x3 = rcv.W(index);
L_x3 = 2;
L_x2 = 1;

% grid points
nx = 200;
ny = nx/2;

x = linspace(-3,3,nx);
y = linspace(-2,2,ny);
[X,Y] = meshgrid(x,y);

% calculate stresses from a fault
[s12_fault,s13_fault] = calcstress_antiplanefault(G,X,Y,x2_bot,x3_bot,L_x3,dip);
u1_fault = calcdisp_antiplanefault(X,Y,x2_bot,x3_bot,L_x3,dip);

% calculate stresses from a shear zone
[Stress_12,Stress_13] = calc_stressgreensfunctions_antiplaneshz(G,X(:),Y(:),x3_bot,L_x2,L_x3);
s12_shz = Stress_12(:,1);
s13_shz = Stress_12(:,2);
[Disp] = calc_dispgreensfunctions_antiplaneshz(X(:),Y(:),x3_bot,L_x2,L_x3);
u1_shz = Disp(:,1)*1/2/L_x2 - Disp(:,2)*0/2/L_x3;

figure(1),clf
subplot(2,2,1)
pcolor(x,y,s12_fault), shading interp
axis tight equal
clim([-1 1]*1)
cb=colorbar; cb.Label.String = '\sigma_{12} (fault)';

subplot(2,2,2)
pcolor(x,y,s13_fault), shading interp
axis tight equal
clim([-1 1]*1)
cb=colorbar; cb.Label.String = '\sigma_{13} (fault)';

subplot(2,2,3)
pcolor(x,y,reshape(s12_shz,ny,nx)), shading interp
axis tight equal
clim([-1 1]*1)
cb=colorbar; cb.Label.String = '\sigma_{12} (shz)';

subplot(2,2,4)
pcolor(x,y,reshape(s13_shz,ny,nx)), shading interp
axis tight equal
clim([-1 1]*1)
cb=colorbar; cb.Label.String = '\sigma_{13} (shz)';
colormap bluewhitered(20)

figure(2),clf
subplot(2,1,1)
pcolor(x,y,u1_fault), shading interp
% hold on
% plot([rcv.x1(index),rcv.x2(index)],[rcv.y1(index),rcv.y2(index)],'ko-')
% alpha 0.5
axis tight equal
grid on
clim([-1 1]*0.5)
cb=colorbar; cb.Label.String = 'u_1 (fault)';

subplot(2,1,2)
pcolor(x,y,reshape(u1_shz,ny,nx)), shading interp
axis tight equal
clim([-1 1]*0.5)
cb=colorbar; cb.Label.String = 'u_1 (shz)';
colormap bluewhitered(20)


