% script to run viscoelastic simulations (semi-analytical solver) and
% demonstrate mesh size invariance in results
% 
% Authors:
% Rishav Mallick, JPL/Caltech

clear

addpath ~/Documents/GitHub/utils/
addpath ~/Documents/GitHub/topotoolbox/colormaps/

Nvec = [10,14,18,22,26,30,34,38,42];
% Nvec = [10,20];

% time steps to evaluate solution
tvec = logspace(0,6,40).*86400;% in seconds

% solution matrix (evaluated inside the source domain)
sol22_in_matrix = zeros(length(Nvec),length(tvec));
sol23_in_matrix = zeros(length(Nvec),length(tvec));
sol22_out_matrix = zeros(length(Nvec),length(tvec));
sol23_out_matrix = zeros(length(Nvec),length(tvec));

% construct mesh
x2extent = 10e3;
x3extent = 20e3;
x3shift = 2e3;
% provide source properties
x20 = 0;
x30 = -(x3shift + x3extent/2);
r0 = 5e3;
strainmax = 1;
% circle source
theta = linspace(0,360,100);
xcircle = r0.*cosd(theta)./1e3;
ycircle = (r0.*sind(theta)+x30)./1e3;

figure(1),clf
set(gcf,'Position',[0 0 3 2]*500)
for i = 1:length(Nvec)
    Nx2 = Nvec(i);
    Nx3 = x3extent*(Nx2)/2/x2extent;
    shz = create_shzmesh(x2extent,Nx2,x3extent,Nx3,x3shift);
    figure(1)
    subplot(3,3,i)
    r = sqrt((shz.x3-x30).^2 + (shz.x2-x20).^2);
    toplot = exp(-r.^2/(2*r0^2));
    plotshz(shz,toplot,1), hold on
    plot(xcircle,ycircle,'k-','LineWidth',2)
    axis tight equal, box on
    xlabel('x'), ylabel('depth')
    cb = colorbar; cb.Label.String = '\Delta\epsilon_v';cb.LineWidth=1.5;
    clim([0 1])
    colormap(ttscm('bilbao',10))
    title(['N = ' num2str(Nx2)])
    set(gca,'FontSize',20,'Linewidth',1)
end
print('meshtest_meshdomain','-djpeg','-r300')
%% loop through various meshes
for i = 1:length(Nvec)
    Nx2 = Nvec(i);
    Nx3 = x3extent*(Nx2)/2/x2extent;

    shz = create_shzmesh(x2extent,Nx2,x3extent,Nx3,x3shift);
    
    e22dot = zeros(shz.N,length(tvec));
    e23dot = zeros(shz.N,length(tvec));

    % load kernels for displacement and stress
    tic
    evl_kernel = compute_shzstresskernels_planestrain(shz,1);
    toc
    % construct deviatoric stress kernels and remove positive eigen values
    
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

    % remove eigen values that cause instabilities   
    lambda = diag(Evals);
    lambda_positive = (real(lambda) >= 0);
    lambda(lambda_positive) = -Inf;

    % provide stress perturbation as an IC (Gaussian source)    
    r = sqrt((shz.x3-x30).^2 + (shz.x2-x20).^2);
    % Gaussian source
    % deltastrainrate = [r.*0;strainmax.*exp(-r.^2/(2*r0^2))];   
    deltastrainrate = repmat(strainmax.*exp(-r.^2/(2*r0^2)),2,1);

    % compute solution
    parfor ti = 1:length(tvec)
        tval = tvec(ti);
        % time-integrate dynammics
        sol = real(Evector*diag(exp(lambda.*tval))/Evector*deltastrainrate);
        e22dot(:,ti) = sol(1:end/2);
        e23dot(:,ti) = sol(1+end/2:end);
    end    
    
    % compute strain inside domain
    sourceindex = r<=r0;
    sol22_in_matrix(i,:) = mean(e22dot(sourceindex,:),1);
    sol23_in_matrix(i,:) = mean(e23dot(sourceindex,:),1);
    sourceindex = r>r0 & r<1.5*r0;
    sol22_out_matrix(i,:) = mean(e22dot(sourceindex,:),1);
    sol23_out_matrix(i,:) = mean(e23dot(sourceindex,:),1);
end

%% %%%%%% plot time series of strain rate evolution averaged over center
cspec = jet(length(Nvec));
tnorm = eta_matrix/30e3;

figure(2),clf
set(gcf,'Position',[0 0 1.5 1.5]*500)
subplot(2,1,1)
for i = 1:length(Nvec)
    semilogx(tvec./tnorm,sqrt(sol22_in_matrix(i,:).^2 + sol23_in_matrix(i,:).^2),'-','LineWidth',2,'Color',cspec(i,:)), hold on
end
xlabel('t*'), ylabel('\epsilon_v (inside)')
axis tight, ylim([0 1.1])
legend(num2str(Nvec'),'box','off')
set(gca,'FontSize',20,'Linewidth',1.5)
subplot(2,1,2)
for i = 1:length(Nvec)
    semilogx(tvec./tnorm,sqrt(sol22_out_matrix(i,:).^2 + sol23_out_matrix(i,:).^2),'-','LineWidth',2,'Color',cspec(i,:)), hold on
end
xlabel('t*'), ylabel('\epsilon_v (outside)')
axis tight, ylim([0 1.1])
set(gca,'FontSize',20,'Linewidth',1.5)
print('meshtest_timeseries','-djpeg','-r300')