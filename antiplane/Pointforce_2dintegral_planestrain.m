% using matlab's adaptive integral to solve the Kelvin integral over a 2-d
% triangular source
% 
% Rishav Mallick, 2023, Caltech Seismolab

clear 
addpath /Users/mallickrishg/Dropbox/scripts/topotoolbox/colormaps/

%% Set model parameters 
% Elasticity parameters
mu_val = 1;
nu_val = 0.25;

% Kelvin force vector (for gravity set fx,fy = 0,-1
fx_val = 0;
fy_val = 1;

% provide plotting type (2-d grid or 1-d line)
% 0 - line, 
% 1 - xy grid
eval_type = 1;

n_pts = 100;
x_vec = linspace(-3, 3, n_pts);
y_vec = linspace(-2, 2, n_pts);

if eval_type == 1
    % create a box grid of points
    [x_mat, y_mat] = meshgrid(x_vec, y_vec);
else
    % only evluate function along a line at y=0
    x_mat = x_vec;
    y_mat = x_vec.*0;
end

%% rectangle or triangle domain
shape = 1; % shape = 0 (rectangle), shape = 1 (triangle)
rectangle_x = 1.5;
rectangle_y = 1;

%% numerical integration (with matlab integral)

ux_vals = zeros(length(x_mat(:,1)),length(x_mat(1,:)));
uy_vals = zeros(length(x_mat(:,1)),length(x_mat(1,:)));
sxx_vals = zeros(length(x_mat(:,1)),length(x_mat(1,:)));
syy_vals = zeros(length(x_mat(:,1)),length(x_mat(1,:)));
sxy_vals = zeros(length(x_mat(:,1)),length(x_mat(1,:)));

ux_numeric_int = zeros(size(x_mat));
uy_numeric_int = zeros(size(x_mat));
sxx_numeric_int = zeros(size(x_mat));
syy_numeric_int = zeros(size(x_mat));
sxy_numeric_int = zeros(size(x_mat));

tic
area_source = rectangle_x*rectangle_y;
Rx = rectangle_x;
Ry = rectangle_y;
% need to provide this triangle in terms of a double integral over x and f(x)
% triangle is defined from 0<=x<=Rx, -Ry<=y<=Ry
ymax = @(x) Ry*(1 - 2*x/Rx);

parfor k=1:numel(x_mat)
    fun_ux = @(x0,y0) gf_ux(x_mat(k),y_mat(k),x0, y0, fx_val, fy_val, mu_val, nu_val);
    fun_uy = @(x0,y0) gf_uy(x_mat(k),y_mat(k),x0, y0, fx_val, fy_val, mu_val, nu_val);
    fun_sxx = @(x0,y0) gf_sxx(x_mat(k),y_mat(k),x0, y0, fx_val, fy_val, mu_val, nu_val);
    fun_syy = @(x0,y0) gf_syy(x_mat(k),y_mat(k),x0, y0, fx_val, fy_val, mu_val, nu_val);
    fun_sxy = @(x0,y0) gf_sxy(x_mat(k),y_mat(k),x0, y0, fx_val, fy_val, mu_val, nu_val);
    
    if shape == 0
        % rectangle source
        ux_numeric_int(k) = integral2(fun_ux,-Rx/2,Rx/2,-Ry/2,Ry/2)./area_source;
        uy_numeric_int(k) = integral2(fun_uy,-Rx/2,Rx/2,-Ry/2,Ry/2)./area_source;
        sxx_numeric_int(k) = integral2(fun_sxx,-Rx/2,Rx/2,-Ry/2,Ry/2)./area_source;
        syy_numeric_int(k) = integral2(fun_syy,-Rx/2,Rx/2,-Ry/2,Ry/2)./area_source;
        sxy_numeric_int(k) = integral2(fun_sxy,-Rx/2,Rx/2,-Ry/2,Ry/2)./area_source;
    elseif shape == 1
        % triangle source
        ux_numeric_int(k) = integral2(fun_ux,0,Rx,-Ry,ymax)./area_source;
        uy_numeric_int(k) = integral2(fun_uy,0,Rx,-Ry,ymax)./area_source;
        sxx_numeric_int(k) = integral2(fun_sxx,0,Rx,-Ry,ymax)./area_source;
        syy_numeric_int(k) = integral2(fun_syy,0,Rx,-Ry,ymax)./area_source;
        sxy_numeric_int(k) = integral2(fun_sxy,0,Rx,-Ry,ymax)./area_source;
    else
        error('Unknwon Shape')
    end
end

toc

%% plot solutions

if eval_type==1
    
    figure(11),clf
    n_skip = 23;

    ux = ux_numeric_int;
    uy = uy_numeric_int;

    toplot_n = sqrt(ux.^2 + uy.^2);
    contourf(x_mat, y_mat, toplot_n,5), hold on
    quiver(x_mat(1:n_skip:end), y_mat(1:n_skip:end),ux(1:n_skip:end),uy(1:n_skip:end),'r','Linewidth',1)
    cb=colorbar;cb.Label.String = 'Displacement |U|';
    clim([0,1].*0.1)
    axis("equal")
    colormap(sky(10))
    xlabel('x'), ylabel('y')
    set(gca,'Fontsize',15)

    % plot stress components
    figure(12),clf
    subplot(3,1,1)
    toplot_n = sxx_numeric_int;
    contourf(x_mat, y_mat, toplot_n,11), hold on
    cb=colorbar;cb.Label.String = '\sigma_{xx}';
    % clim([-1,1].*max(abs(sxx_numeric_int(:))))
    clim([-1,1].*0.5)
    axis("equal")
    xlabel('x'), ylabel('y')
    set(gca,'Fontsize',15)

    subplot(3,1,2)
    toplot_n = syy_numeric_int;
    contourf(x_mat, y_mat, toplot_n,11), hold on
    cb=colorbar;cb.Label.String = '\sigma_{yy}';
    clim([-1,1].*0.5)
    axis("equal")
    xlabel('x'), ylabel('y')
    set(gca,'Fontsize',15)

    subplot(3,1,3)
    toplot_n = sxy_numeric_int;
    contourf(x_mat, y_mat, toplot_n,11), hold on
    cb=colorbar;cb.Label.String = '\sigma_{xy}';
    clim([-1,1].*0.5)
    axis("equal")
    colormap(ttscm('vik',20))
    xlabel('x'), ylabel('y')
    set(gca,'Fontsize',15)
    % print('force2d_planestrain_fy','-djpeg','-r200')

else
    cspec = cool(n_eval);
    for i = 1:n_eval 
        figure(11)
        ux = squeeze(ux_vals(:,:,i));
        uy = squeeze(uy_vals(:,:,i));
        subplot(2,1,1)
        plot(x_mat,ux,'-','LineWidth',2,'Color',cspec(i,:)), hold on
        axis tight, grid on
        xlabel('x'), ylabel('u_x')
        set(gca,'FontSize',15)

        subplot(2,1,2)
        plot(x_mat,uy,'-','LineWidth',2,'Color',cspec(i,:)), hold on
        axis tight, grid on
        xlabel('x'), ylabel('u_y')
        set(gca,'FontSize',15)

        figure(12)
        subplot(3,1,1)
        toplot_n = squeeze(sxx_vals(:,:,i));
        plot(x_mat,toplot_n,'-','LineWidth',2,'Color',cspec(i,:)), hold on
        axis tight, grid on
        xlabel('x'), ylabel('\sigma_{xx}')
        set(gca,'FontSize',15)

        subplot(3,1,2)
        toplot_n = squeeze(syy_vals(:,:,i));
        plot(x_mat,toplot_n,'-','LineWidth',2,'Color',cspec(i,:)), hold on
        axis tight, grid on
        xlabel('x'), ylabel('\sigma_{yy}')
        set(gca,'FontSize',15)

        subplot(3,1,3)
        toplot_n = squeeze(sxy_vals(:,:,i));
        plot(x_mat,toplot_n,'-','LineWidth',2,'Color',cspec(i,:)), hold on
        axis tight, grid on
        xlabel('x'), ylabel('\sigma_{xy}')
        set(gca,'FontSize',15)
        

    end
end


%% define kelvin point source function
function ux = gf_ux(x0, y0, xoffset, yoffset, fx, fy, mu, nu)
[ux, ~, ~, ~, ~] = kelvin_point(x0, y0, xoffset, yoffset, fx, fy, mu, nu);
end
function uy = gf_uy(x0, y0, xoffset, yoffset, fx, fy, mu, nu)
[~, uy, ~, ~, ~] = kelvin_point(x0, y0, xoffset, yoffset, fx, fy, mu, nu);
end
function sxx = gf_sxx(x0, y0, xoffset, yoffset, fx, fy, mu, nu)
[~, ~, sxx, ~, ~] = kelvin_point(x0, y0, xoffset, yoffset, fx, fy, mu, nu);
end
function syy = gf_syy(x0, y0, xoffset, yoffset, fx, fy, mu, nu)
[~, ~, ~, syy, ~] = kelvin_point(x0, y0, xoffset, yoffset, fx, fy, mu, nu);
end
function sxy = gf_sxy(x0, y0, xoffset, yoffset, fx, fy, mu, nu)
[~, ~, ~, ~, sxy] = kelvin_point(x0, y0, xoffset, yoffset, fx, fy, mu, nu);
end

function [ux, uy, sxx, syy, sxy] = kelvin_point(x0, y0, xoffset, yoffset, fx, fy, mu, nu)
    x = x0 - xoffset;
    y = y0 - yoffset;
    C = 1 / (4 * pi * (1 - nu));
    r = sqrt(x.^2 + y.^2);
    g = -C .* log(r);
    gx = -C .* x ./ (x.^2 + y.^2);
    gy = -C .* y ./ (x.^2 + y.^2);
    gxy = C .* 2 .* x .* y ./ (x.^2 + y.^2).^2;
    gxx = C .* (x.^2 - y.^2) ./ (x.^2 + y.^2).^2;
    gyy = -gxx;

    ux = fx / (2 * mu) * ((3 - 4 * nu) .* g - x .* gx) + fy / (2 * mu) .* (-y .* gx);
    uy = fx / (2 * mu) * (-x .* gy) + fy / (2 * mu) * ((3 - 4 * nu) .* g - y .* gy);
    
    sxx = fx .* (2 * (1 - nu) .* gx - x .* gxx) + fy .* (2 * nu .* gy - y .* gxx);
    syy = fx .* (2 * nu .* gx - x .* gyy) + fy .* (2 * (1 - nu) .* gy - y .* gyy);
    sxy = fx .* ((1 - 2 * nu) .* gy - x .* gxy) + fy .* ((1 - 2 * nu) .* gx - y .* gxy);
end

