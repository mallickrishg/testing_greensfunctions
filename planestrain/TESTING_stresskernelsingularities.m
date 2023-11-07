
clear

% Elastic parameters
mu = 1;
nu = 0.25;
a = 1.0;
Beta = 0;
xe = 0;
ye = 0;

n_pts = 1000;
x_obs = linspace(-2, 2, n_pts);
y_obs = zeros(size(x_obs)) + 0;

Fx = 1;
Fy = 0;

[Disp,Stress] = LTkernelFS(x_obs,y_obs,xe,ye,a,Beta,Fx,Fy,nu,mu);

figure(1),clf
plot(x_obs,Disp(:,1),'LineWidth',2), hold on
plot(x_obs,Disp(:,2),'LineWidth',2)
