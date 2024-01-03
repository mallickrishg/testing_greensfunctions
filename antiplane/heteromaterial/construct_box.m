function rcv = construct_box(Lx, Lz, Nx, Nz)
% Input:
%   width, height: dimensions of the box
%   resolution: number of elements along each dimension

% Generate node coordinates
x = linspace(-Lx/2, Lx/2, Nx+1);
y = linspace(-Lz/2, Lz/2, Nz+1);

% Create 2D grid
[X, Y] = meshgrid(x, y);

% Reshape coordinates
vertices = [X(:), Y(:)];

% Create faces
faces = [];
for i = 1:Nx
    for j = 1:Nz
        % Define the indices of the four vertices of each rectangular element
        v1 = (Nz+1)*(i-1) + j;
        v2 = v1 + 1;
        v3 = v1 + Nz + 2;
        v4 = v1 + Nz + 1;

        % Create two triangular faces for each rectangular element
        faces = [faces; v1, v2, v3, v4];
    end
end

rcv = [];
rcv.N = length(faces(:,1));
rcv.xc = zeros(rcv.N,2);
rcv.Lx = zeros(rcv.N,1);
rcv.Lz = zeros(rcv.N,1);
for i = 1:rcv.N
    rcv.xc(i,:) = mean(vertices(faces(i,:),:));
    rcv.Lz(i) = abs(diff(vertices(faces(i,1:2),2)));
    rcv.Lx(i) = abs(diff(vertices(faces(i,2:3),1)));
end
rcv.A = rcv.Lx.*rcv.Lz;

end
