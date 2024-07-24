% using matlab's adaptive integral to solve the point force GF 
% integrated over a rectangular source
% and compare solution with G-L & T-S quadratures
% Rishav Mallick, 2023, Caltech Seismolab

clear 

%% Set model parameters 
% Elasticity parameters
mu = 1;

% provide plotting type (2-d grid or 1-d line)
% 0 - line, 
% 1 - xy grid
eval_type = 1;

% rectangle domain
rectangle_x = 2;
rectangle_y = 0.5;

nx = 40;
ny = 40;

x_vec = linspace(-2, 2, nx);
y_vec = linspace(-2, 2, ny);

if eval_type == 1
    % create a box grid of points
    [x_mat, y_mat] = meshgrid(x_vec, y_vec);
else
    % only evluate function along a line at y=0
    % x_mat = x_vec;
    % y_mat = x_vec.*0;
    rcv = construct_box(-rectangle_x/2,-rectangle_y/2,rectangle_x,rectangle_y,rectangle_y/nx); 
    x_mat = rcv.xc(:,1);
    y_mat = rcv.xc(:,2);
end

%% numerical integration (with matlab integral)

u1 = zeros(size(x_mat));
s12 = zeros(size(x_mat));
s13 = zeros(size(x_mat));

tic
area_source = 1;%rectangle_x*rectangle_y;
Rx = rectangle_x;
Ry = rectangle_y;

parfor k=1:numel(x_mat)
    fun_u = @(x0,y0) pointforce_disp(x_mat(k),y_mat(k),x0,y0);
    fun_s12 = @(x0,y0) pointforce_s12(x_mat(k),y_mat(k),x0,y0);
    fun_s13 = @(x0,y0) pointforce_s13(x_mat(k),y_mat(k),x0,y0);
    
    % rectangle source
    u1(k) = integral2(fun_u,-Rx/2,Rx/2,-Ry/2,Ry/2)./area_source;
    s12(k) = mu.*integral2(fun_s12,-Rx/2,Rx/2,-Ry/2,Ry/2)./area_source;
    s13(k) = mu.*integral2(fun_s13,-Rx/2,Rx/2,-Ry/2,Ry/2)./area_source;    

end
toc

%% numerical integration (GL quadrature)

u1_GL = zeros(size(x_mat));
s12_GL = zeros(size(x_mat));
s13_GL = zeros(size(x_mat));

% Gauss-Legendre over a 2-d square domain
N_gl = 29;
[xkvec,wxkvec] = calc_gausslegendre_weights(N_gl);
[ykvec,wykvec] = calc_gausslegendre_weights(N_gl);
[xk,yk] = meshgrid(xkvec.*rectangle_x/2,ykvec.*rectangle_y/2);
[wxk,wyk] = meshgrid(wxkvec,wykvec);
wk = wxk.*wyk;

tic
parfor k=1:numel(xk)
    n_u = pointforce_disp(x_mat,y_mat,xk(k),yk(k));
    n_s12 = mu.*pointforce_s12(x_mat,y_mat,xk(k),yk(k));
    n_s13 = mu.*pointforce_s13(x_mat,y_mat,xk(k),yk(k));
    u1_GL = u1_GL + n_u.*wk(k);
    s12_GL = s12_GL + n_s12.*wk(k);
    s13_GL = s13_GL + n_s13.*wk(k);
end
toc

%% compute solution using tanh-sinh quadrature

u1_TS = zeros(size(x_mat));
s12_TS = zeros(size(x_mat));
s13_TS = zeros(size(x_mat));

% numerical solution
h=0.05;% step size for tanh-sinh
n=fix(2/h);

kvec = (-n:n)';
wkvec=(0.5*h*pi*cosh(kvec*h))./(cosh(0.5*pi*sinh(kvec*h))).^2;
xkvec=tanh(0.5*pi*sinh(kvec*h));
[xk,yk] = meshgrid(xkvec.*rectangle_x/2,xkvec.*rectangle_y/2);
[wxk,wyk] = meshgrid(wkvec,wkvec);
wk = wxk.*wyk;

tic
parfor k=1:numel(xk)        
    n_u = pointforce_disp(x_mat,y_mat,xk(k),yk(k));
    n_s12 = mu.*pointforce_s12(x_mat,y_mat,xk(k),yk(k));
    n_s13 = mu.*pointforce_s13(x_mat,y_mat,xk(k),yk(k));
    u1_TS = u1_TS + n_u.*wk(k);
    s12_TS = s12_TS + n_s12.*wk(k);
    s13_TS = s13_TS + n_s13.*wk(k);
end
toc

%% plot results and compare solutions
figure(1),clf
subplot(3,1,1)
pcolor(x_vec,y_vec,u1), shading interp
colorbar
axis tight equal
title('Displacement - integral2')

subplot(3,1,2)
pcolor(x_vec,y_vec,u1_GL), shading interp
axis tight equal
colorbar
title('Gauss-Legendre')

subplot(3,1,3)
pcolor(x_vec,y_vec,u1_TS), shading interp
axis tight equal
colorbar
title('Tanh-sinh')

colormap('hot(20)')

% stress components
figure(2),clf
subplot(3,2,1)
pcolor(x_vec,y_vec,s12), shading interp, hold on
contour(x_vec,y_vec,s12,14,'k-')
cb=colorbar; cb.Label.String = '\sigma_{12}';
axis tight equal
xlabel('x'), ylabel('y')
title('\sigma_{12} integral2')

subplot(3,2,3)
pcolor(x_vec,y_vec,s12_GL), shading interp, hold on
% contour(x_vec,y_vec,s12_GL,14,'k-')
cb=colorbar; cb.Label.String = '\sigma_{12}';
axis tight equal
title('Gauss-Legendre')

subplot(3,2,5)
pcolor(x_vec,y_vec,s12_TS), shading interp, hold on
contour(x_vec,y_vec,s12_TS,14,'k-')
cb=colorbar; cb.Label.String = '\sigma_{12}';
axis tight equal
title('Tanh-sinh')

subplot(3,2,2)
pcolor(x_vec,y_vec,s13), shading interp, hold on
contour(x_vec,y_vec,s13,14,'k-')
cb=colorbar; cb.Label.String = '\sigma_{13}';
axis tight equal
title('\sigma_{13} integral2')

subplot(3,2,4)
pcolor(x_vec,y_vec,s13_GL), shading interp, hold on
contour(x_vec,y_vec,s13_GL,14,'k-')
cb=colorbar; cb.Label.String = '\sigma_{13}';
axis tight equal
title('Gauss-Legendre')

subplot(3,2,6)
pcolor(x_vec,y_vec,s13_TS), shading interp, hold on
contour(x_vec,y_vec,s13_TS,14,'k-')
cb=colorbar; cb.Label.String = '\sigma_{13}';
axis tight equal
title('Tanh-sinh')
colormap('turbo(15)')

set(findobj(gcf,'type','axes'),'Fontsize',15,'Linewidth',1.5,'TickDir','both')%,'clim',[-1 1]*0.5)
%% compare residuals 
% Δ = 100*(1-σ_numeric/σ_true)

figure(10),clf
subplot(2,2,1)
pcolor(x_vec,y_vec,(1-s12_GL./s12).*100), shading interp
cb = colorbar; cb.Label.String = '% error';
axis tight equal
title('\Delta\sigma_{12} Gauss-Legendre')

subplot(2,2,3)
pcolor(x_vec,y_vec,(1-s12_TS./s12).*100), shading interp
cb = colorbar; cb.Label.String = '% error';
axis tight equal
title('\Delta\sigma_{12} Tanh-sinh')

subplot(2,2,2)
pcolor(x_vec,y_vec,(1-s13_GL./s13).*100), shading interp
cb = colorbar; cb.Label.String = '% error';
axis tight equal
title('\Delta\sigma_{13} Gauss-Legendre')

subplot(2,2,4)
pcolor(x_vec,y_vec,(1-s13_TS./s13).*100), shading interp
cb = colorbar; cb.Label.String = '% error';
axis tight equal
title('\Delta\sigma_{13} Tanh-sinh')
colormap(turbo(100))

set(findobj(gcf,'type','axes'),'Clim',[-1 1]*10,'TickDir','both','FontSize',15,'Linewidth',1)

figure(11),clf
subplot(2,2,1)
pcolor(x_vec,y_vec,abs(s12-s12_GL)), shading interp
colorbar
axis tight equal
title('\delta\sigma_{12} Gauss-Legendre')

subplot(2,2,3)
pcolor(x_vec,y_vec,abs(s12-s12_TS)), shading interp
axis tight equal
colorbar
title('\delta\sigma_{12} Tanh-sinh')

subplot(2,2,2)
pcolor(x_vec,y_vec,abs(s13-s13_GL)), shading interp
axis tight equal
colorbar
title('\delta\sigma_{13} Gauss-Legendre')

subplot(2,2,4)
pcolor(x_vec,y_vec,abs(s13-s13_TS)), shading interp
axis tight equal
colorbar
title('\delta\sigma_{13} Tanh-sinh')
colormap(turbo(100))

set(findobj(gcf,'type','axes'),'ColorScale','log','Clim',10.^[-4 -1],'TickDir','both','FontSize',15,'Linewidth',1)

%% define function
function u = pointforce_disp(x,y,xs,ys)
r = sqrt((x-xs).^2 + (y-ys).^2);
u = 1/2/pi.*log(r);
end
function s12 = pointforce_s12(x,y,xs,ys)
r = sqrt((x-xs).^2 + (y-ys).^2);
s12 = 1/2/pi.*(x-xs)./(r.^2);
end
function s13 = pointforce_s13(x,y,xs,ys)
r = sqrt((x-xs).^2 + (y-ys).^2);
s13 = 1/2/pi.*(y-ys)./(r.^2);
end