
clear
addpath functions/

% Set model parameters 
% Elasticity parameters
mu = 1;

% source strength terms (alpha_0 + alpha_1*x)
Nsources = 3;
alpha_0 = 1;
alpha_1 = 1;
alpha_2 = 1;

% dimensions of source
% specify location and dimensions of sources
w = ones(Nsources,1).*1;
% xloc = ones(Nsources,1).*0;
% yloc = [2,0,-2];
% dip = 90.*ones(Nsources,1);
dip = [0,90,90];
xloc = [0-cosd(dip(1)),0,0];
yloc = [1+sind(dip(1)),0,-2];


% discretize evaluation points
nx = 500;
ny = nx;

x_vec = linspace(-4, 4, nx);
y_vec = linspace(-6, 6, ny);
[x_mat, y_mat] = meshgrid(x_vec, y_vec);

%% compute displacement and stress
Ku = zeros(nx*ny,3,Nsources);
Ks12 = zeros(nx*ny,3,Nsources);
Ks13 = zeros(nx*ny,3,Nsources);

tic
for i = 1:Nsources
    [u1,s12,s13] = calc_disp_stress_quadfault(x_mat(:),y_mat(:),...
        xloc(i),yloc(i),w(i),[alpha_0,alpha_1,alpha_2],dip(i));
    Ku(:,:,i) = u1;
    Ks12(:,:,i) = s12;
    Ks13(:,:,i) = s13;
end
toc
%% plot displacement and stresses

% % for linear slip elements
% alpha0_vec = [1,1.75,0.75];
% alpha1_vec = [-1,0.25,0.75];
% sources = [alpha0_vec;alpha1_vec];

% for quadratic slip elements
alpha0_vec = [1,6,1];
alpha1_vec = [-2,0,2];
alpha2_vec = [1,-2,1];
% alpha0_vec = [1,0,-1];
% alpha1_vec = [0,0,0];
% alpha2_vec = -[1,0,-1];

sources = [alpha0_vec;alpha1_vec;alpha2_vec].*0.5;

% use tensor products to contract 3-d matrices
u1 = tensorprod(Ku,sources,[2 3],[1 2]);
s12 = tensorprod(Ks12,sources,[2 3],[1 2]);
s13 = tensorprod(Ks13,sources,[2 3],[1 2]);

figure(1),clf
pcolor(x_vec,y_vec,reshape(u1,ny,nx)), shading interp, hold on
contour(x_vec,y_vec,reshape(u1,ny,nx),linspace(-1,1,11).*max(abs(u1(:))),'k-')
axis tight equal
cb=colorbar; cb.Location='northoutside';
clim([-1 1].*max(abs(u1(:))))
colormap(ttscm('vik',1000))

figure(2),clf
subplot(1,2,1)
pcolor(x_vec,y_vec,reshape(s12,ny,nx)), shading interp, hold on
contour(x_vec,y_vec,reshape(s12,ny,nx),linspace(-1,1,21).*max(abs(s12(:))),'k-')
axis tight equal
% clim([-1 1].*max(abs(s12(:))))
clim([-1 1].*1)
cb=colorbar;cb.Label.String = '\sigma_{12}';cb.Location='northoutside';
xlabel('x'), ylabel('y')
set(gca,'FontSize',15,'LineWidth',1.5)

subplot(1,2,2)
pcolor(x_vec,y_vec,reshape(s13,ny,nx)), shading interp, hold on
contour(x_vec,y_vec,reshape(s13,ny,nx),linspace(-1,1,21).*max(abs(s13(:))),'k-')
axis tight equal
% clim([-1 1].*max(abs(s13(:))))
clim([-1 1].*1)
cb=colorbar;cb.Label.String = '\sigma_{13}';cb.Location='northoutside';
xlabel('x'), ylabel('y')
set(gca,'FontSize',15,'LineWidth',1.5)
colormap(ttscm('bam',100))
%% calculate slip and stress change ON the fault

x_mat = zeros(nx*Nsources,1);
y_mat = zeros(nx*Nsources,1);
deps = -1e-8;
for i = 1:Nsources
    x1 = xloc(i)-w(i).*cosd(dip(i)) ;% these are not correct - but work for this example
    x2 = xloc(i)+w(i).*cosd(dip(i)) ;
    y1 = yloc(i)+w(i).*sind(dip(i)) ;
    y2 = yloc(i)-w(i).*sind(dip(i)) ;

    x_mat((i-1)*nx+1:i*nx) = linspace(x1,x2,nx)' + deps*sind(dip(i));
    y_mat((i-1)*nx+1:i*nx) = linspace(y1,y2,nx)' + deps*cosd(dip(i));
end

Ku = zeros(nx*Nsources,3,Nsources);
Ks12 = zeros(nx*Nsources,3,Nsources);
Ks13 = zeros(nx*Nsources,3,Nsources);

tic
for i = 1:Nsources
    [u1,s12,s13] = calc_disp_stress_quadfault(x_mat(:),y_mat(:),...
        xloc(i),yloc(i),w(i),[alpha_0,alpha_1,alpha_2],dip(i));
    Ku(:,:,i) = u1;
    Ks12(:,:,i) = s12;
    Ks13(:,:,i) = s13;
end
toc

% use tensor products to contract 3-d matrices
u1 = tensorprod(Ku,sources,[2 3],[1 2]);
s12 = tensorprod(Ks12,sources,[2 3],[1 2]);
s13 = tensorprod(Ks13,sources,[2 3],[1 2]);

figure(3),clf
subplot(1,2,1)
% toplot = reshape(u1,nx,nx);
plot(u1,'.','LineWidth',3)
grid on, box on
ylabel('displacement')
set(gca,'FontSize',15,'LineWidth',1.5)
subplot(1,2,2)
% toplot = reshape(s12,nx,nx);
plot(s12,'.','LineWidth',3), hold on
% toplot = reshape(s13,nx,nx);
plot(s13,'.','LineWidth',3)
grid on, box on
ylim([-1 1]*2)
% xlim([-1 1].*max(abs(get(gca,'XLim'))))
ylabel('\Delta\tau')
set(gca,'FontSize',15,'LineWidth',1.5)

