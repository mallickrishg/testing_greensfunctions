% using matlab's adaptive integral to solve the point force GF 
% integrated over a rectangular source
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
rectangle_x = 1;
rectangle_y = 1;

% triangle domain: (0,0),(l_x,dl_y),(0,l_y)
l_x = 1;
l_y = 1;
dl_y = -0.5;

nx = 100;
ny = 100;

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
% area_source = 0.5*rectangle_x*rectangle_y;
% Rx = rectangle_x;
% Ry = rectangle_y;
% Calculate the lengths of the sides of the triangle
a = sqrt((l_x - 0)^2 + (dl_y - 0)^2);
b = sqrt((0 - l_x)^2 + (l_y - dl_y)^2);
c = sqrt((0 - 0)^2 + (0 - l_y)^2);

% Calculate the semi-perimeter of the triangle
s = (a + b + c) / 2;

% Calculate the area of the triangle using Heron's formula
area_source = sqrt(s * (s - a) * (s - b) * (s - c));

ymax = @(x) l_y - (l_y-dl_y)*(x/l_x);
ymin = @(x) dl_y*x/l_x;
% fun_u = @(x,y,x0,y0) pointforce_disp(x,y,x0,y0);
% fun_s12 = @(x,y,x0,y0) pointforce_s12(x,y,x0,y0);
% fun_s13 = @(x,y,x0,y0) pointforce_s13(x,y,x0,y0);
parfor k=1:numel(x_mat)
    fun_u = @(x0,y0) pointforce_disp(x_mat(k),y_mat(k),x0,y0);
    fun_s12 = @(x0,y0) pointforce_s12(x_mat(k),y_mat(k),x0,y0);
    fun_s13 = @(x0,y0) pointforce_s13(x_mat(k),y_mat(k),x0,y0);
    
    % rectangle source
    % u1(k) = integral2(fun_u,-Rx/2,Rx/2,-Ry/2,Ry/2)./area_source;
    % s12(k) = mu.*integral2(fun_s12,-Rx/2,Rx/2,-Ry/2,Ry/2)./area_source;
    % s13(k) = mu.*integral2(fun_s13,-Rx/2,Rx/2,-Ry/2,Ry/2)./area_source;

    % triangle source
    % u1(k) = integral2(@(x0,y0) fun_u(x_mat(k),y_mat(k),x0,y0),0,Rx,-Ry,ymax)./area_source;
    % s12(k) = mu.*integral2(@(x0,y0) fun_s12(x_mat(k),y_mat(k),x0,y0),0,Rx,-Ry,ymax)./area_source;
    % s13(k) = mu.*integral2(@(x0,y0) fun_s13(x_mat(k),y_mat(k),x0,y0),0,Rx,-Ry,ymax)./area_source;

    u1(k) = integral2(fun_u,0,l_x,ymin,ymax)./area_source;
    s12(k) = mu.*integral2(fun_s12,0,l_x,ymin,ymax)./area_source;
    s13(k) = mu.*integral2(fun_s13,0,l_x,ymin,ymax)./area_source;

end
toc

%% plot results

if eval_type == 1
    figure(21),clf
    pcolor(x_mat,y_mat,reshape(u1(:),ny,nx)),shading interp
    hold on
    contour(x_mat,y_mat,reshape(u1(:),ny,nx),9,'k-')
    % plot(rectangle_x/2.*[-1,1,1,-1,-1],rectangle_y/2.*[1,1,-1,-1,1],'k-','Linewidth',2)
    plot([0,l_x,0,0],[0,dl_y,l_y,0],'k-','Linewidth',2)
    axis tight equal
    cb=colorbar;cb.Label.String = 'u';
    colormap('sky(10)')
    set(gca,'Fontsize',15,'LineWidth',1)

    figure(22),clf
    subplot(211)
    pcolor(x_mat,y_mat,reshape(s12(:),ny,nx)),shading interp
    hold on
    contour(x_mat,y_mat,reshape(s12(:),ny,nx),[-1:0.1:1]*0.5,'k-')
    % plot(rectangle_x/2.*[-1,1,1,-1,-1],rectangle_y/2.*[1,1,-1,-1,1],'k-','Linewidth',2)
    plot([0,l_x,0,0],[0,dl_y,l_y,0],'k-','Linewidth',1)
    axis tight equal
    cb=colorbar;cb.Label.String = '\sigma_{xz}';
    set(gca,'Fontsize',15,'LineWidth',1)
    % clim([-1 1]*max(abs(s12(:))))
    clim([-1 1]*0.5)
    subplot(212)
    pcolor(x_mat,y_mat,reshape(s13(:),ny,nx)),shading interp
    hold on
    contour(x_mat,y_mat,reshape(s13(:),ny,nx),[-1:0.1:1]*0.5,'k-')
    % plot(rectangle_x/2.*[-1,1,1,-1,-1],rectangle_y/2.*[1,1,-1,-1,1],'k-','Linewidth',2)
    plot([0,l_x,0,0],[0,dl_y,l_y,0],'k-','Linewidth',1)
    axis tight equal
    cb=colorbar;cb.Label.String = '\sigma_{yz}';
    % clim([-1 1]*max(abs(s13(:))))
    clim([-1 1]*0.5)
    colormap('turbo(20)')
    set(gca,'Fontsize',15,'LineWidth',1)
    % print('force2d_antiplane','-djpeg','-r200')
else
    figure(31),clf
    subplot(2,1,1)
    plot(u1,'linewidth',2)
    axis tight
    xlabel('node'),ylabel('u_1')
    subplot(2,1,2)
    plot(s12,'LineWidth',2), hold on
    plot(s13,'LineWidth',2)
    axis tight
    legend('\sigma_{12}','\sigma_{13}')
    ylabel('Stress component'),xlabel('node')
end


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