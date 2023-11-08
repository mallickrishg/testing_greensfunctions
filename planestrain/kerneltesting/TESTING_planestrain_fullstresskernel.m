% script to test how stress singularities manifest in the deviatoric stress
% kernels and the scale at which they are most important 
% for eq-cycle calculations
% 
% Authors:
% Rishav Mallick, JPL/Caltech

clear

addpath ~/Documents/GitHub/utils/

% construct mesh
x2extent = 10e3;
x3extent = 20e3;

Nx2 = 10;
Nx3 = x3extent*(Nx2)/2/x2extent;

x3shift = 100e3;

shz = create_shzmesh(x2extent,Nx2,x3extent,Nx3,x3shift);

nobs = 20;
xobs = linspace(-x2extent,x2extent,nobs);
zobs = linspace(-(x3shift+x3extent),-x3shift,nobs);
[xobs,zobs] = meshgrid(xobs,zobs);
obs = [xobs(:), zobs(:)];

%% load kernels for displacement and stress
tic
stresskernel = compute_shzstresskernels_planestrain(shz,1);
dispkernel = compute_shzdispkernels_planestrain(shz,obs);
toc
%% construct full stress kernels and remove positive eigen values

bigL =  [stresskernel.LL2222, stresskernel.LL2322, stresskernel.LL3322;...
         stresskernel.LL2223, stresskernel.LL2323, stresskernel.LL3323;...
         stresskernel.LL2233, stresskernel.LL2333, stresskernel.LL3333];

[Evec,eig_ps] = eig(bigL);

% remove eigen values that cause instabilities 
eig_ps_c = diag(eig_ps);
eig_ps_c(real(eig_ps_c)>0) = 0;
bigLmod = real(Evec*diag(eig_ps_c)/Evec);

% plot eigen values
figure(1),clf
subplot(211)
stem(sort(real(diag(eig_ps))),'k.'), hold on
stem(sort(real(eig_ps_c)),'r.')
axis tight, grid on

subplot(212)
plot(real(eig_ps),imag(eig_ps),'ko'), hold on
plot(real(eig_ps_c),imag(eig_ps_c),'r.')
axis tight, grid on
xlim([-1 1].*max(abs(get(gca,'XLim'))))
set(findobj(gcf,'type','axes'),'Fontsize',15,'Linewidth',1)
