% script to run viscoelastic simulations (semi-analytical solver) and
% demonstrate mesh size invariance in results
% 
% Authors:
% Rishav Mallick, JPL/Caltech

clear

addpath ~/Documents/GitHub/utils/
addpath ~/Documents/GitHub/topotoolbox/colormaps/

Nvec = 20;

% time steps to evaluate solution
tvec = logspace(0,6,100).*86400;% in seconds

% construct mesh
x2extent = 10e3;
x3extent = 20e3;
x3shift = 2e3;

%% compute eigen values and plot positive values
Nx2 = Nvec;
Nx3 = x3extent*(Nx2)/2/x2extent;

shz = create_shzmesh(x2extent,Nx2,x3extent,Nx3,x3shift);

e22dot = zeros(shz.N,length(tvec));
e23dot = zeros(shz.N,length(tvec));
e22dot_reg = zeros(shz.N,length(tvec));
e23dot_reg = zeros(shz.N,length(tvec));

%% construct deviatoric stress kernels
tic
evl_kernel = compute_shzstresskernels_planestrain(shz,1);
toc

% prescribe viscosity
eta_matrix = 1e18*1e-6; % in MPa-s

% store original deviatoric stress kernels
L2222o = (evl_kernel.LL2222 - evl_kernel.LL3322 - evl_kernel.LL2233 + evl_kernel.LL3333)./(2.*eta_matrix);
L2322o = (evl_kernel.LL2322 - evl_kernel.LL2333)./(2.*eta_matrix);
L2223o = (evl_kernel.LL2223 - evl_kernel.LL3323)./eta_matrix;
L2323o = evl_kernel.LL2323./eta_matrix;

%  combine all individual deviatorickernels
rheoparam = [L2222o L2322o;...
             L2223o L2323o];

% eigen-value decomposition of the rheological parameter
[Evector,Evals] = eig(rheoparam);
lambda = diag(Evals);
lambda_reg = lambda;
lambda_positive = (real(lambda) >= 0);
lambda_reg(lambda_positive) = -Inf;

%% provide stress perturbation as an IC (Gaussian source)  
x20 = 0;
x30 = -(x3shift + x3extent/2);
r = sqrt((shz.x3-x30).^2 + (shz.x2-x20).^2);
r0 = 5e3;
strainmax = 1;
deltastrainrate = repmat(strainmax.*exp(-r.^2/(2*r0^2)),2,1);
% compute solution
parfor ti = 1:length(tvec)
    tval = tvec(ti);
    % time-integrate dynammics
    sol = real(Evector*diag(exp(lambda.*tval))/Evector*deltastrainrate);    
    e22dot(:,ti) = sol(1:end/2);
    e23dot(:,ti) = sol(1+end/2:end);
    sol = real(Evector*diag(exp(lambda_reg.*tval))/Evector*deltastrainrate);
    e22dot_reg(:,ti) = sol(1:end/2);
    e23dot_reg(:,ti) = sol(1+end/2:end);
end

% circle source
theta = linspace(0,360,100);
xcircle = r0.*cosd(theta)./1e3;
ycircle = (r0.*sind(theta)+x30)./1e3;

figure(1),clf
% set(gcf,'Position',[0 0 3 2]*500)
Nx2 = Nvec;
Nx3 = x3extent*(Nx2)/2/x2extent;
shz = create_shzmesh(x2extent,Nx2,x3extent,Nx3,x3shift);
figure(1)
toplot = exp(-r.^2/(2*r0^2));
plotshz(shz,toplot,1), hold on
plot(xcircle,ycircle,'k-','LineWidth',2)
axis tight equal, box on
xlabel('x'), ylabel('depth')
cb = colorbar; cb.Label.String = '\Delta\epsilon_v';cb.LineWidth=1.5;
clim([0 1])
set(gca,'FontSize',20,'Linewidth',1)

%% plot solution at time steps
tnorm = eta_matrix/30e3;

tsnapshots = [0.1,0.5,1,2,5,10,15,20,22].*tnorm;

figure(10),clf
set(gcf,'Position',[0 0 1.5 1.5]*500)
subplot(2,1,1)
sourceindex = r<=r0;
plot(tvec./tnorm,mean(sqrt(e22dot(sourceindex,:).^2 + e23dot(sourceindex,:).^2)),'-','LineWidth',3), hold on
plot(tvec./tnorm,mean(sqrt(e22dot_reg(sourceindex,:).^2 + e23dot_reg(sourceindex,:).^2)),'r--','LineWidth',2)
xlim([min(tvec) max(tvec)]./tnorm)
% ylim([1e-2 1e2]), grid on
ylim([0 1.2])
xlabel('t*'), ylabel('\epsilon_v (inside)')
set(gca,'Xscale','log','YScale','lin','LineWidth',1.5,'FontSize',20)

subplot(2,1,2)
sourceindex = r>r0 & r<1.5*r0;
plot(tvec./tnorm,mean(sqrt(e22dot(sourceindex,:).^2 + e23dot(sourceindex,:).^2)),'-','LineWidth',3), hold on
plot(tvec./tnorm,mean(sqrt(e22dot_reg(sourceindex,:).^2 + e23dot_reg(sourceindex,:).^2)),'r--','LineWidth',2)
xlim([min(tvec) max(tvec)]./tnorm)
% ylim([1e-2 1e2]), grid on
ylim([0 1.2])
xlabel('t*'), ylabel('\epsilon_v (outside)')
set(gca,'Xscale','log','YScale','lin','LineWidth',1.5,'FontSize',20)

print('relaxation_plots_timeseries','-djpeg','-r300')

figure(11),clf
set(gcf,'Position',[0 0 3 2]*500)
for i = 1:length(tsnapshots)  
    tindex = abs(tsnapshots(i)-tvec) == min(abs(tsnapshots(i)-tvec));
    subplot(3,3,i)
    toplot = sqrt(e22dot(:,tindex).^2 + e23dot(:,tindex).^2);
    plotshz(shz,toplot,1)
    axis tight equal, box on
    cb=colorbar;cb.LineWidth=1.5;cb.Label.String = '\epsilon_v';
    clim([0 1])
    colormap(ttscm('bilbao',10))
    title(['t* = ' num2str(tsnapshots(i)/tnorm)],'FontWeight','normal')
    set(gca,'FontSize',25,'LineWidth',1.5,'TickDir','out','TickLength',[0.03 0.1])
end
print('relaxation_plots_positiveeigen','-djpeg','-r300')
figure(12),clf
set(gcf,'Position',[0 0 3 2]*500)
for i = 1:length(tsnapshots)  
    tindex = abs(tsnapshots(i)-tvec) == min(abs(tsnapshots(i)-tvec));
    subplot(3,3,i)
    toplot = sqrt(e22dot_reg(:,tindex).^2 + e23dot_reg(:,tindex).^2);
    plotshz(shz,toplot,1)
    axis tight equal, box on
    cb=colorbar;cb.LineWidth=1.5;cb.Label.String = '\epsilon_v';
    clim([0 1])
    colormap(ttscm('bilbao',10))
    title(['t* = ' num2str(tsnapshots(i)/tnorm)],'FontWeight','normal')
    set(gca,'FontSize',25,'LineWidth',1.5,'TickDir','out','TickLength',[0.03 0.1])
end
print('relaxation_plots_regularizedeigen','-djpeg','-r300')
