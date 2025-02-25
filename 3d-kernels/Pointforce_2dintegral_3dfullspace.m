clear
% Elastic material parameters
mu = 1;
nu = 0.25;

% provide observation points
n_pts = 200;
x_vec = linspace(-3, 3, n_pts);
y_vec = linspace(-3, 3, n_pts);

% create a box grid of points
[x_mat,y_mat] = meshgrid(x_vec, y_vec);
z_mat = zeros(numel(x_mat),1)+1e-4;

% triangle source (the triangle lies in the x-y plane)
% triangle domain: (0,0),(l_x,dl_y),(0,l_y)
l_x = 1;
l_y = 2;
dl_y = 0.0;
% need to provide this triangle in terms of a double integral over x and f(x)
ymax = @(x) l_y - (l_y-dl_y)*(x/l_x);
ymin = @(x) dl_y*x/l_x;

% plot point source solution
uxvec = zeros(numel(x_mat),1);
uyvec = zeros(numel(x_mat),1);
uzvec = zeros(numel(x_mat),1);
Fvec = zeros(numel(x_mat),1);
parfor i = 1:length(uxvec)
    displacement = kelvinGreen([x_mat(i),y_mat(i),z_mat(i)],mu,nu);
    Fvec(i) = pointsource(x_mat(i),y_mat(i),z_mat(i),0,0,0,mu,nu);
    uxvec(i) = displacement(1,1);
    uyvec(i) = displacement(2,1);
    uzvec(i) = displacement(3,1);
end

figure(1),clf
subplot(1,3,1)
pcolor(x_vec,y_vec,reshape(uxvec,n_pts,n_pts)), shading interp
axis tight equal
colorbar

subplot(1,3,2)
pcolor(x_vec,y_vec,reshape(uyvec,n_pts,n_pts)), shading interp
axis tight equal
colorbar

subplot(1,3,3)
pcolor(x_vec,y_vec,reshape(uzvec,n_pts,n_pts)), shading interp
axis tight equal
colorbar

figure(2),clf
pcolor(x_vec,y_vec,reshape(Fvec,n_pts,n_pts)), shading interp
axis tight equal
colorbar
clim([-1 1]*0.1)

%% numerical integration %%%%%
Fintegral = zeros(numel(x_mat),1);
parfor k = 1:numel(x_mat)
    funint = @(x0,y0) pointsource(x_mat(k),y_mat(k),z_mat(k),x0,y0,0,mu,nu);
    Fintegral(k) = integral2(funint,0,l_x,ymin,ymax);
end

%% plot results of integration

figure(3),clf
subplot(3,1,[1 2])
toplot = reshape(Fintegral,n_pts,n_pts);
pcolor(x_vec,y_vec,toplot), 
shading interp, hold on
contourf(x_vec,y_vec,toplot,21)
axis tight equal
clim([-1 1]*0.2)
cb=colorbar;cb.Label.String='F(x,y,0)';cb.Location='northoutside';
colormap(ttscm('vik',1000))
% colormap(turbo)
xlabel('x'), ylabel('y')
set(gca,'Fontsize',15,'TickDir','out','LineWidth',1.5)

subplot(3,1,3)
index = find(y_mat(:)==min(y_vec(y_vec>0)));
plot(x_vec,toplot(index),'LineWidth',3)
ylim([-1 1]*0.2)
xlabel('x'), ylabel('F(x,y,0)')
set(gca,'Fontsize',15,'TickDir','out','LineWidth',1.5)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = pointsource(x,y,z,x0,y0,z0,mu,nu)
% define stress Greens function
rx = x-x0;
ry = y-y0;
rz = z-z0;
rnorm = sqrt(rx.^2 + ry.^2 + rz.^2);

k = 3;
j = 3;
i = 1;
delta_ij = double(i==j);
delta_ik = double(i==k);
delta_jk = double(j==k);
prefactor = 1/(16*pi*mu*(1-nu));

F = prefactor * ( ...
    -(3-4*nu)*delta_ij.*rz./rnorm.^3 ...
    - (delta_ik.*rz + delta_jk.*rx)./rnorm.^3 ...
    + 3.*rx.*rz.*rz./rnorm.^5 ); %s_xz for fz input
% dGdx(i,j,k) = prefactor * ...
%   ( - (3-4*nu)*delta_ij*r(k)/rnorm^3 ...
%     - (delta_ik*r(j)+delta_jk*r(i))/rnorm^3 ...
%     + 3*r(i)*r(j)*r(k)/rnorm^5 );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = kelvinGreen(r, mu, nu)
    % Computes the Kelvin Green's function (3x3 tensor) at vector r.
    rnorm = norm(r);
    if rnorm < 1e-8
        G = zeros(3,3);
        return;
    end
    prefactor = 1/(16*pi*mu*(1-nu));
    G = prefactor * ( (3-4*nu)/rnorm * eye(3) + (r'*r)/(rnorm^3) );
end