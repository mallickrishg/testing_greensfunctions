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
rectangle_x = 0.5;
rectangle_y = 1.5;

nx = 100;
ny = 50;

x_vec = linspace(-4, 4, nx);
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
area_source = 0.5*rectangle_x*rectangle_y;
Rx = rectangle_x;
Ry = rectangle_y;
ymax = @(x) Ry*(1 - 2*x/Rx);

% fun_u = @(x,y,x0,y0) pointforce_disp(x,y,x0,y0);
% fun_s12 = @(x,y,x0,y0) pointforce_s12(x,y,x0,y0);
% fun_s13 = @(x,y,x0,y0) pointforce_s13(x,y,x0,y0);
parfor k=1:numel(x_mat)
    fun_u = @(x0,y0) pointforce_disp(x_mat(k),y_mat(k),x0,y0);
    fun_s12 = @(x0,y0) pointforce_s12(x_mat(k),y_mat(k),x0,y0);
    fun_s13 = @(x0,y0) pointforce_s13(x_mat(k),y_mat(k),x0,y0);
    
    % rectangle source
    u1(k) = integral2(fun_u,-Rx/2,Rx/2,-Ry/2,Ry/2)./area_source;
    s12(k) = mu.*integral2(fun_s12,-Rx/2,Rx/2,-Ry/2,Ry/2)./area_source;
    s13(k) = mu.*integral2(fun_s13,-Rx/2,Rx/2,-Ry/2,Ry/2)./area_source;

    % triangle source
    % u1(k) = integral2(@(x0,y0) fun_u(x_mat(k),y_mat(k),x0,y0),0,Rx,-Ry,ymax)./area_source;
    % s12(k) = mu.*integral2(@(x0,y0) fun_s12(x_mat(k),y_mat(k),x0,y0),0,Rx,-Ry,ymax)./area_source;
    % s13(k) = mu.*integral2(@(x0,y0) fun_s13(x_mat(k),y_mat(k),x0,y0),0,Rx,-Ry,ymax)./area_source;

    % u1(k) = integral2(fun_u,0,Rx,-Ry,ymax)./area_source;
    % s12(k) = mu.*integral2(fun_s12,0,Rx,-Ry,ymax)./area_source;
    % s13(k) = mu.*integral2(fun_s13,0,Rx,-Ry,ymax)./area_source;

end
toc

%% plot results

if eval_type == 1
    figure(21),clf
    pcolor(x_mat,y_mat,reshape(abs(u1(:)),ny,nx)),shading interp
    axis tight equal
    colorbar
    colormap('sky(10)')

    figure(22),clf
    subplot(211)
    pcolor(x_mat,y_mat,reshape(s12(:),ny,nx)),shading interp
    hold on
    contour(x_mat,y_mat,reshape(s12(:),ny,nx),9,'k-')
    axis tight equal
    cb=colorbar;cb.Label.String = '\sigma_{xz}';
    set(gca,'Fontsize',15,'LineWidth',1)
    % clim([-1 1]*max(abs(s12(:))))
    clim([-1 1]*0.5)
    subplot(212)
    pcolor(x_mat,y_mat,reshape(s13(:),ny,nx)),shading interp
    hold on
    contour(x_mat,y_mat,reshape(s13(:),ny,nx),9,'k-')
    axis tight equal
    cb=colorbar;cb.Label.String = '\sigma_{yz}';
    % clim([-1 1]*max(abs(s13(:))))
    clim([-1 1]*0.5)
    colormap('turbo(1000)')
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