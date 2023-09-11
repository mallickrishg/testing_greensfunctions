% testing how stress greens functions work in a half-space
% Rishav Mallick, Caltech, 2022

clear

clear
addpath ~/Dropbox/scripts/utils/

% shear modulus
G = 1;
nu = 0.25;

% grid points
nx2 = 100;
nx3 = nx2/2;

x2 = linspace(-3,3,nx2);
x3 = linspace(-2,2,nx3);
[X2,X3] = meshgrid(x2,x3);

% compute half space stress
Y2 = 0e3;
W = 2;
dip = 20;
Y3 = 0;
slip = -1;
open = 0;

%% compute displacements in half-space
ng = 5;
[Disp] = LDdispFS(X2,X3,Y2,Y3,W/2,-deg2rad(dip),slip,open,nu);
ue = Disp(:,1);uz = Disp(:,2);

figure(11),clf
subplot(211)
pcolor(x2,x3,reshape(ue,nx3,nx2)), hold on, shading interp
contour(x2,x3,reshape(ue,nx3,nx2),[-1:.1:1].*0.5,'k-')
% quiver(X2(1:ng:end),X3(1:ng:end),ue(1:ng:end)',uz(1:ng:end)','k-','Linewidth',1)
plot(x2,0.*x2,'k-','LineWidth',1)
axis tight equal, box on
cb = colorbar;cb.Label.String='u_e';
caxis([-1,1].*0.5)
set(gca,'Fontsize',15,'YDir','normal','LineWidth',2)

subplot(212)
pcolor(x2,x3,reshape(uz,nx3,nx2)), hold on, shading interp
contour(x2,x3,reshape(uz,nx3,nx2),[-1:.1:1].*0.5,'k-')
plot(x2,0.*x2,'k-','LineWidth',1)
axis tight equal, box on
cb = colorbar;cb.Label.String='u_z';
caxis([-1,1].*0.5)
colormap(bluewhitered(40))
set(gca,'Fontsize',15,'YDir','normal','LineWidth',2)

%% compute stresses
[Stress] = LDstressFS(X2,X3,Y2,Y3,W/2,-deg2rad(dip),slip,open,nu,2*G*(1+nu));
Sxx = Stress(:,1);
Szz = Stress(:,2);
Sxz = Stress(:,3);

figure(12),clf
subplot(3,1,1)
imagesc(x2,x3,reshape(Sxx,nx3,nx2)), hold on
plot(x2,0.*x2,'k-','LineWidth',1)
axis tight equal, box on
cb = colorbar;cb.Label.String='\sigma_{xx}';
caxis([-1,1].*1)
set(gca,'Fontsize',15,'YDir','normal','LineWidth',2)

subplot(3,1,2)
imagesc(x2,x3,reshape(Szz,nx3,nx2)), hold on
plot(x2,0.*x2,'k-','LineWidth',1)
axis tight equal, box on
cb = colorbar;cb.Label.String='\sigma_{zz}';
caxis([-1,1].*1)
set(gca,'Fontsize',15,'YDir','normal','LineWidth',2)

subplot(3,1,3)
imagesc(x2,x3,reshape(Sxz,nx3,nx2)), hold on
plot(x2,0.*x2,'k-','LineWidth',1)
axis tight equal, box on
cb = colorbar;cb.Label.String='\sigma_{xz}';
caxis([-1,1].*1)
set(gca,'Fontsize',15,'YDir','normal','LineWidth',2)
colormap bluewhitered(100)
