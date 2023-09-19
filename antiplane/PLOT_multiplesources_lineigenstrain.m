
clear
addpath functions/

% Set model parameters 
% Elasticity parameters
mu = 1;

% source strength terms (alpha_0 + alpha_1*x)
Nsources = 2;
alpha_0 = 1;
alpha_1 = 1;
sources = [1,1;0,0];

% dimensions of source
% specify location and dimensions of cources
xloc = [-1,1]*0.5;
yloc = ones(Nsources,1).*0;
Rx = ones(Nsources,1).*0.5;
Ry = ones(Nsources,1).*1.5;

% discretize evaluation points
nx = 100;
ny = nx;

x_vec = linspace(-4, 4, nx);
y_vec = linspace(-4, 4, ny);
[x_mat, y_mat] = meshgrid(x_vec, y_vec);

%% compute displacement and stress
Ku = zeros(nx*ny,2,Nsources);
Ks12 = zeros(nx*ny,2,Nsources);
Ks13 = zeros(nx*ny,2,Nsources);

for i = 1:Nsources
    [u1,s12,s13] = calc_disp_stress_lineigenstrain(x_mat(:),y_mat(:),...
        xloc(i),yloc(i),Rx(i),Ry(i),[alpha_0,alpha_1]);
    Ku(:,:,i) = u1;
    Ks12(:,:,i) = s12;
    Ks13(:,:,i) = s13;
end

%% plot displacement and stresses
u1 = zeros(nx*ny,1);
s12 = zeros(nx*ny,1);
s13 = zeros(nx*ny,1);

alpha0_vec = [1,1];
alpha1_vec = [0,0];
sources = [alpha0_vec;alpha1_vec];

for i = 1:Nsources
    u1 = u1 + Ku(:,:,i)*sources(:,i);
    s12 = s12 + Ks12(:,:,i)*sources(:,i);
    s13 = s13 + Ks13(:,:,i)*sources(:,i);
end


figure(1),clf
pcolor(x_vec,y_vec,reshape(u1,ny,nx)), shading interp, hold on
contour(x_vec,y_vec,reshape(u1,ny,nx),linspace(-1,1,21).*max(abs(u1(:))),'k-')
axis tight equal
colorbar
clim([-1 1].*max(abs(u1(:))))
colormap bluewhitered(20)

figure(2),clf
subplot(2,1,1)
pcolor(x_vec,y_vec,reshape(s12,ny,nx)), shading interp, hold on
contour(x_vec,y_vec,reshape(s12,ny,nx),linspace(-1,1,21).*max(abs(s12(:))),'k-')
axis tight equal
clim([-1 1].*max(abs(s12(:))))
colorbar

subplot(2,1,2)
pcolor(x_vec,y_vec,reshape(s13,ny,nx)), shading interp, hold on
contour(x_vec,y_vec,reshape(s13,ny,nx),linspace(-1,1,21).*max(abs(s13(:))),'k-')
axis tight equal
clim([-1 1].*max(abs(s13(:))))
colorbar

colormap bluewhitered(20)

figure(11),clf
toplot = reshape(u1,ny,nx);
plot(x_vec,toplot(ny/2,:),'.-','LineWidth',2)
axis tight, grid on
xlabel('x'), ylabel('u_1')
