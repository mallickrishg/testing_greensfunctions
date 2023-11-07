function evl = compute_shzstresskernels_planestrain(shz,scf)
% compute stress interaction kernel for shear zones in plane strain
% compute displacement kernels at the surface

G = 30e3;
nu = 0.25;

% always use vertical shear zones
phi = 90;

evl = [];

LL2222 = zeros(shz.N,shz.N);
LL2223 = zeros(shz.N,shz.N);
LL2233 = zeros(shz.N,shz.N);

if ~exist('scf','var')
    scf = 1;
end

parfor i = 1:shz.N
    % rescaling factor for far-field stresses    
    
    % evluate stresses at these points (shear zone centers)
    X2 = shz.x2;
    X3 = shz.x3;
    q2 = shz.x2(i);
    q3 = -shz.x3(i) - shz.W(i)/2;
        
    % for phi=90, T - x2; W - x3
    T = shz.T(i);
    W = shz.W(i);
    
    % eps22
    epsv22p = 1*scf^2;
    epsv23p = 0;
    epsv33p = 0;
    [s22,s23,s33]=computeStressPlaneStrainShearZone( ...
        X2(:),-X3(:),q2,q3,T/scf,W/scf,phi,epsv22p,epsv23p,epsv33p,G,nu);    
    LL2222(:,i) = s22;
    LL2223(:,i) = s23;
    LL2233(:,i) = s33;      

    % eps23
    epsv22p = 0;
    epsv23p = 1*scf^2;
    epsv33p = 0;
    [s22,s23,s33]=computeStressPlaneStrainShearZone( ...
        X2(:),-X3(:),q2,q3,T/scf,W/scf,phi,epsv22p,epsv23p,epsv33p,G,nu);
    LL2322(:,i) = s22;
    LL2323(:,i) = s23;
    LL2333(:,i) = s33;

    % eps33
    epsv22p = 0;
    epsv23p = 0;
    epsv33p = 1*scf^2;
    [s22,s23,s33]=computeStressPlaneStrainShearZone( ...
        X2(:),-X3(:),q2,q3,T/scf,W/scf,phi,epsv22p,epsv23p,epsv33p,G,nu);    
    LL3322(:,i) = s22;
    LL3323(:,i) = s23;
    LL3333(:,i) = s33;

end

evl.LL2222 = LL2222;
evl.LL2223 = LL2223;
evl.LL2233 = LL2233;
evl.LL2322 = LL2322;
evl.LL2323 = LL2323;
evl.LL2333 = LL2333;
evl.LL3322 = LL3322;
evl.LL3323 = LL3323;
evl.LL3333 = LL3333;

end