
clear

% shear modulus
mu = 1;
nu = 0.25;

% grid points
nx = 200;
ny = nx/2 + 1;

x_obs = linspace(-3,3,nx);
y_obs = linspace(-1.5,1.5,ny);
[X_obs,Y_obs] = meshgrid(x_obs,y_obs);

% construct a constant traction source and compute [u],[σ] in the medium
L = 1;
xs = 0;
ys = 0;
Fx = 1;
Fy = 0;

[Disp,Stress] = LTkernelFS(X_obs(:),Y_obs(:),xs,ys,L,0,Fx,Fy,nu,mu);

% plot results in 2-d bulk
nskip = 19;
figure(1),clf
toplot = sqrt(Disp(:,1).^2 + Disp(:,2).^2);
pcolor(x_obs,y_obs,reshape(toplot,ny,nx)), shading interp, hold on
quiver(X_obs(1:nskip:end)',Y_obs(1:nskip:end)',Disp(1:nskip:end,1),Disp(1:nskip:end,2),'k')
axis tight equal
colorbar
colormap('hot(20)')
clim([0 0.5])

figure(2),clf
subplot(3,1,1)
toplot = Stress(:,1);
pcolor(x_obs,y_obs,reshape(toplot,ny,nx)), shading interp
axis tight equal
cb=colorbar;cb.Label.String = '\sigma_{xx}';
clim([-1 1]*max(abs(toplot)))

subplot(3,1,2)
toplot = Stress(:,2);
pcolor(x_obs,y_obs,reshape(toplot,ny,nx)), shading interp
axis tight equal
cb=colorbar;cb.Label.String = '\sigma_{xy}';
clim([-1 1]*max(abs(toplot)))

subplot(3,1,3)
toplot = Stress(:,3);
pcolor(x_obs,y_obs,reshape(toplot,ny,nx)), shading interp
axis tight equal
cb=colorbar;cb.Label.String = '\sigma_{yy}';
colormap('bluewhitered(50)')
clim([-1 1]*max(abs(toplot)))

% evaluate results exactly ON the source
x_obs = linspace(min(x_obs),max(x_obs),1000)';
y_obs = zeros(size(x_obs));
[Disp,Stress] = LTkernelFS(x_obs,y_obs,xs,ys,L,0,Fx,Fy,nu,mu);

figure(3),clf
subplot(2,1,1)
plot(x_obs,Disp(:,1),'LineWidth',2), hold on
plot(x_obs,Disp(:,2),'LineWidth',2)
legend('u_x','u_y')
xlabel('x'), ylabel('u')
axis tight

subplot(2,1,2)
plot(x_obs,Stress(:,1),'.','LineWidth',2), hold on
plot(x_obs,Stress(:,2),'.','LineWidth',2)
plot(x_obs,Stress(:,3),'.','LineWidth',2)
legend('\sigma_{xx}','\sigma_{xy}','\sigma_{yy}')
axis tight
xlabel('x'), ylabel('\sigma')

%% try to construct same solution using slip elements
% method - 
% need to figure out a way to use slip elements, which produce
% spatially continuous tractions (apart from singularities) to recreate
% stresses that arise from traction elements, which produce discontinuous
% stresses











