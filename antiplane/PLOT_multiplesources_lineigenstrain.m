
clear

% Set model parameters 
% Elasticity parameters
mu = 1;

% source strength terms (alpha_0 + alpha_1*x)
alpha_0 = 1;
alpha_1 = 1;

% dimensions of source
% specify location and dimensions of cources
xloc = 0;
yloc = 0;
Rx = 1;
Ry = 1.5;

% discretize evaluation points
nx = 100;
ny = nx;

x_vec = linspace(-4, 4, nx);
y_vec = linspace(-4, 4, ny);
[x_mat, y_mat] = meshgrid(x_vec, y_vec);

%% compute displacement and stress

for i = 1:Nsources
    [u1,s12,s13] = calc_disp_stress_lineigenstrain(x_mat(:),y_mat(:),xloc,yloc,Rx,Ry,[alpha_0,alpha_1]);
end
