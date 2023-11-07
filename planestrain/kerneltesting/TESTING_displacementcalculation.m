mu = 30e3;
nu = 0.25;

index = 45;
W = shz.W(index);
T = shz.T(index);
q2 = shz.x2(index);
q3 = -shz.x3(index) - shz.W(index)/2;
phi = 90;

epsv22p = 1/W;
epsv23p = 0;
epsv33p = -1/W;

nobs = 60;
x = linspace(-x2extent,x2extent,nobs);
z = linspace(-(x3shift+x3extent),-x3shift,nobs);
[x2,x3] = meshgrid(x,z);

[u2,u3]=computeDisplacementPlaneStrainShearZone(x2(:),-x3(:),q2,q3,T,W,phi,epsv22p,epsv23p,epsv33p,mu,nu);

figure(2),clf
dummy = nan(shz.N,1);
dummy(index) = 1;
plotshz(shz,dummy,1), hold on
axis tight equal
scatter(x2(:)./1e3,x3(:)./1e3,100,sqrt(u2.^2 + u3.^2))
quiver(x2(:)./1e3,x3(:)./1e3,u2,-u3,'k','LineWidth',2)
colorbar;
colormap turbo(20)
clim([0 1]*0.5)
