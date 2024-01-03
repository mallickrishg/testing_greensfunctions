function bimat = create_horizontalinterfaces(x2min,x2max,nx2,D)
% function to create a data structure of interface mesh
% x2min -> x2max - horizontal extent of layer
% nx2 - number of BE to mesh interface
% D - vector of depth levels
% Rishav Mallick, Caltech SeismoLab, 2022

Nint = length(D);

bimat = [];
x2h = linspace(x2min,x2max,nx2);

nxmesh = (length(x2h)-1);
bimat.Nxmesh = nxmesh;
bimat.N = Nint*nxmesh;
bimat.x2c = zeros(bimat.N,1);
bimat.x3c = zeros(bimat.N,1);
bimat.x2 = zeros(bimat.N,1);
bimat.x3 = zeros(bimat.N,1);

bimat.W = zeros(bimat.N,1);
bimat.dip = zeros(bimat.N,1);
bimat.nvec = zeros(bimat.N,2);

for i = 1:Nint
    
    index = [(i-1)*nxmesh + 1:i*nxmesh];
    
    x3h = D(i).*ones(1,nx2);
            
    bimat.x2c(index) = mean([x2h(1:end-1)',x2h(2:end)'],2);
    bimat.x3c(index) = D(i).*ones(1,nxmesh);

    bimat.x2(index) = x2h(1:end-1)';
    bimat.x3(index) = D(i).*ones(1,nxmesh);
    
    bimat.W(index) = sqrt((x2h(1:end-1)-x2h(2:end)).^2 + (x3h(1:end-1)-x3h(2:end)).^2);
    
    bimat.dip(index) = -atan2d(-x3h(1:end-1)+x3h(2:end),-x2h(1:end-1)+x2h(2:end));
    
    bimat.nvec(index,[1,2]) = [sind(bimat.dip(index)),cosd(bimat.dip(index))];

end