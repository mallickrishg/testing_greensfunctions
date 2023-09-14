%% Set model parameters 
% Elasticity parameters
mu = 1;

% provide plotting type (2-d grid or 1-d line)
% 0 - line, 
% 1 - xy grid
eval_type = 1;

% specify location and dimensions of cources
source_strength = [1,-1]*2;
xloc = [-0.5,0.5].*1;
yloc = [0,0];
rectangle_x = [1,1].*1;
rectangle_y = [1.5,1.5];


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
u1_eval = zeros(size(x_mat));
s12_eval = zeros(size(x_mat));
s13_eval = zeros(size(x_mat));

tic

for i = 1:length(rectangle_x)
    Rx = rectangle_x(i);
    Ry = rectangle_y(i);
    pos0 = [xloc(i),yloc(i)];
    
    parfor k=1:numel(x_mat)
        fun_u = @(x0,y0) pointforce_disp(x_mat(k)-pos0(1),y_mat(k)-pos0(2),x0,y0);
        fun_s12 = @(x0,y0) pointforce_s12(x_mat(k)-pos0(1),y_mat(k)-pos0(2),x0,y0);
        fun_s13 = @(x0,y0) pointforce_s13(x_mat(k)-pos0(1),y_mat(k)-pos0(2),x0,y0);

        % rectangle source
        u1_eval(k) = integral2(fun_u,-Rx/2,Rx/2,-Ry/2,Ry/2);
        s12_eval(k) = mu.*integral2(fun_s12,-Rx/2,Rx/2,-Ry/2,Ry/2);
        s13_eval(k) = mu.*integral2(fun_s13,-Rx/2,Rx/2,-Ry/2,Ry/2);

    end
    u1 = u1 + u1_eval.*source_strength(i);
    s12 = s12 + s12_eval.*source_strength(i);
    s13 = s13 + s13_eval.*source_strength(i);
end
toc

%% plot results

if eval_type == 1
    figure(21),clf
    pcolor(x_mat,y_mat,reshape(abs(u1(:)),ny,nx)),shading interp
    axis tight equal
    colorbar
    clim([0 0.5])
    colormap('sky(10)')

    figure(22),clf
    subplot(211)
    pcolor(x_mat,y_mat,reshape(s12(:),ny,nx)),shading interp
    hold on
    contour(x_mat,y_mat,reshape(s12(:),ny,nx),[-0.5:0.1:0.5],'k-')
    axis tight equal
    cb=colorbar;cb.Label.String = '\sigma_{xz}';
    set(gca,'Fontsize',15,'LineWidth',1)
    clim([-1 1]*0.5)
    subplot(212)
    pcolor(x_mat,y_mat,reshape(s13(:),ny,nx)),shading interp
    hold on
    contour(x_mat,y_mat,reshape(s13(:),ny,nx),[-0.5:0.1:0.5],'k-')
    axis tight equal
    cb=colorbar;cb.Label.String = '\sigma_{yz}';
    clim([-1 1]*0.5)
    colormap('turbo(1000)')
    set(gca,'Fontsize',15,'LineWidth',1)
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