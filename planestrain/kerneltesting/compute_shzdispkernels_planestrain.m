function evl = compute_shzdispkernels_planestrain(shz,obs)
% compute displacement kernel for shear zones in plane strain

mu = 30e3;
nu = 0.25;

% always use vertical shear zones
phi = 90;

Nobs = length(obs(:,1));

evl = [];

Gx22 = zeros(Nobs,shz.N);
Gx23 = zeros(Nobs,shz.N);
Gx33 = zeros(Nobs,shz.N);
Gz22 = zeros(Nobs,shz.N);
Gz23 = zeros(Nobs,shz.N);
Gz33 = zeros(Nobs,shz.N);

parfor i = 1:shz.N
    % rescaling factor for far-field stresses    
    
    % evluate stresses at these points (shear zone centers)
    X2 = obs(:,1);
    X3 = obs(:,2);
    q2 = shz.x2(i);
    q3 = -shz.x3(i) - shz.W(i)/2;
        
    % for phi=90, T - x2; W - x3
    T = shz.T(i);
    W = shz.W(i);

    [u2,u3]=computeDisplacementPlaneStrainShearZone(X2(:),-X3(:),...
        q2,q3,T,W,phi,1,0,0,mu,nu);
    Gx22(:,i) = u2;
    Gz22(:,i) = -u3;

    [u2,u3]=computeDisplacementPlaneStrainShearZone(X2(:),-X3(:),...
        q2,q3,T,W,phi,0,1,0,mu,nu);
    Gx23(:,i) = u2;
    Gz23(:,i) = -u3;

    [u2,u3]=computeDisplacementPlaneStrainShearZone(X2(:),-X3(:),...
        q2,q3,T,W,phi,0,0,1,mu,nu);
    Gx33(:,i) = u2;
    Gz33(:,i) = -u3;

end

evl.Gx22 = Gx22;
evl.Gx23 = Gx23;
evl.Gx33 = Gx33;
evl.Gz22 = Gz22;
evl.Gz23 = Gz23;
evl.Gz33 = Gz33;

end