% script to convert Mathematica expressions to MATLAB
% for antiplane point force integrals with 
% linearly varying body forces along a line segment
% 
% AUTHORS:
% Rishav Mallick, JPL, 2024

clear

mu = 1;

Nobs = 1000;
x = linspace(-2,2,Nobs)';
y = zeros(Nobs,1);

% Displacement kernels copied over from Mathematica for a horizontal line
% element present at y = 0, -w <= x <= w
% 
% ux,uy kernels have two linear basis functions
% f1: goes from (-w,0) to (w,1)
% f2: goes from (-w,1) to (w,0)

u_1 = @(w) 0.5.*((1/8).*pi.^(-1).*w.^(-1).*((-4).*w.*(2.*w+x)+4.*(w+x).*y.*atan2(w-x,y)+...
    4.*(w+x).*y.*atan2(w+x,y)+((w+(-1).*x) ...
  .*(3.*w+x)+y.^2).*log((w+(-1).*x).^2+y.^2)+(w+x+(-1).*y).*(w+x+y) ...
  .*log((w+x).^2+y.^2)));

u_2 = @(w) 0.5.*((1/8).*pi.^(-1).*w.^(-1).*(4.*w.*((-2).*w+x)+4.*(w+(-1).*x).*y.* ...
  atan2(w-x,y)+4.*(w+(-1).*x).*y.*atan2(w+x,y) ...
  +(w+(-1).*x+(-1).*y).*(w+(-1).*x+y).*log((w+(-1).*x).^2+y.^2)+(( ...
  3.*w+(-1).*x).*(w+x)+y.^2).*log((w+x).^2+y.^2)));

ex_1 = @(w) (1/8).*pi.^(-1).*w.^(-1).*((-4).*w+2.*y.*(atan2(w-x,y)+atan2(w+x,y))+...
    (-1).*(w+x).*log((w+(-1).*x).^2+y.^2)+(w+x).*log((w+x).^2+y.^2));

ex_2 = @(w) (-1/8).*pi.^(-1).*w.^(-1).*((-4).*w+2.*y.*(atan2(w-x,y)+atan2(w+x,y))+(w+(-1).*x).*(log((w+(-1).*x).^2+y.^2)+ ...
  (-1).*log((w+x).^2+y.^2)));

ey_1 = @(w) (1/8).*pi.^(-1).*w.^(-1).*(2.*(w+x).*(atan2(w-x,y)+ ...
  atan2(w+x,y))+y.*(log((w+(-1).*x).^2+y.^2)+(-1).*log((w+x) ...
  .^2+y.^2)));

ey_2 = @(w) (1/8).*pi.^(-1).*w.^(-1).*(2.*(w+(-1).*x).*(atan2(w-x,y)+atan2(w+x,y))+...
    y.*((-1).*log((w+(-1).*x).^2+y.^2)+log((w+x).^2+y.^2)));

% set element length
w = 1;

% Displacements
figure(1),clf
plot(x,u_1(w),'Linewidth',2), hold on
plot(x,u_2(w),'Linewidth',2)
plot([1 1]*w,get(gca,'YLim'),'k--')
plot([1 1]*-w,get(gca,'YLim'),'k--')
xlabel('x')
title('u (f)','FontWeight','normal')
set(findobj(gcf,'type','axes'),'FontSize',20,'LineWidth',1,'TickDir','both','YLim',[-1 1]*0.23);

% Stress
figure(2),clf
subplot(2,1,1)
plot(x,mu.*ex_1(w),'Linewidth',2), hold on
plot(x,mu.*ex_2(w),'Linewidth',2)
plot([1 1]*w,get(gca,'YLim'),'k--')
plot([1 1]*-w,get(gca,'YLim'),'k--')
xlabel('x')
title('e_x (f)','FontWeight','normal')

subplot(2,1,2)
plot(x,mu.*ey_1(w),'Linewidth',2), hold on
plot(x,mu.*ey_2(w),'Linewidth',2)
plot([1 1]*w,get(gca,'YLim'),'k--')
plot([1 1]*-w,get(gca,'YLim'),'k--')
xlabel('x')
title('e_y (f)','FontWeight','normal')
set(findobj(gcf,'type','axes'),'FontSize',20,'LineWidth',1,'TickDir','both','YLim',[-1 1]*0.5);

%% use kernel calculation function

[Disp,Stress] = LinForceKernelFS(x,y,0,0,w,mu);
% Disp_kernels - [Nobs x 2 basis functions]
% Stress_kernels - [Nobs x (sx or sy) x 2 basis functions]

% Displacements
figure(1)
plot(x,Disp(:,1),'k--','Linewidth',1)
plot(x,Disp(:,2),'k--','Linewidth',1)
axis tight

% Stress
figure(2)
subplot(2,1,1)
for i = 1:2
    plot(x,Stress(:,1,i),'k--','Linewidth',1)
end

subplot(2,1,2)
for i = 1:2
    plot(x,Stress(:,2,i),'k--','Linewidth',1)
end
