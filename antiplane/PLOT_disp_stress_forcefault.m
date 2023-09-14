clear
addpath functions/

u1 = @(x,y,a) (-4*a + ...
                2.*x.*(atan2(2*a.*x,x.^2+y.^2-a^2)) +...
                (a - y).*log(x.^2 + (a - y).^2) + ...
                (a + y).*log(x.^2 + (a + y).^2))./(4*pi);

s12 = @(x,y,a) (atan2(2*a.*x,x.^2+y.^2-a^2))./(2*pi);

s13 = @(x,y,a) (log(x.^2 + (a + y).^2) - log(x.^2 + (a - y).^2))./(4*pi);

% calculate displacements and stresses in the medium
nx = 100;
ny = round(nx/2);
x = linspace(-6,6,nx);
y = linspace(-3,3,ny);
[X,Y] = meshgrid(x,y);

% half-length of source
a = 1;

% uplot = u1(X(:),Y(:),a);
% s12plot = s12(X(:),Y(:),a);
% s13plot = s13(X(:),Y(:),a);

% for a dipping source
dip = 90;
[uplot,s12plot,s13plot] = calc_disp_stress_forcefault(X(:),Y(:),a,dip);

figure(1),clf
pcolor(x,y,reshape(uplot,ny,nx)), shading interp
axis tight equal
colorbar
xlabel('x'),ylabel('y')
title('Displacement')
colormap bluewhitered

figure(2),clf
subplot(211)
pcolor(x,y,reshape(s12plot,ny,nx)), shading interp
axis tight equal
clim([-1 1]*.5)
colorbar
xlabel('x'),ylabel('y')
title('\sigma_{12}')

subplot(212)
pcolor(x,y,reshape(s13plot,ny,nx)), shading interp
axis tight equal
colorbar
xlabel('x'),ylabel('y')
title('\sigma_{13}')
colormap bluewhitered(20)