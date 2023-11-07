function evl = compute_shzdispkernels_planestrain(shz)
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