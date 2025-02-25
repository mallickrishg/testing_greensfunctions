function KelvinTriangleIntegration
    %% Material and source parameters
    mu = 30e9;            % shear modulus [Pa]
    nu = 0.25;            % Poisson's ratio
    stress_drop = 1e6;    % stress drop [Pa]
    
    %% Define a triangular patch (each row is a vertex in 3D)
    vertices = [0, 0, 0;
                5000, 0, 0;
                0, 5000, 0];
    % For a patch in the x-y plane, choose an upward unit normal:
    n = [0, 0, 1];
    
    %% Define observation grid in the x-y plane at fixed depth z_obs
    x_min = -10000; x_max = 10000;
    y_min = -10000; y_max = 10000;
    nx = 30; ny = 30;
    [X, Y] = meshgrid(linspace(x_min, x_max, nx), linspace(y_min, y_max, ny));
    z_obs = 0;
    
    tol = 1e-6;
    
    %% Compute displacement field on the grid using Kelvin kernel integration
    u_field = zeros(nx, ny, 3);
    for i = 1:nx
        for j = 1:ny
            x_obs = [X(i,j), Y(i,j), z_obs];
            u_field(i,j,:) = integrateTriangleKelvin(x_obs, vertices, n, stress_drop, mu, nu, tol);
        end
    end
    u_magnitude = sqrt(sum(u_field.^2,3));
    
    %% Compute stress field using the displacement gradient
    stress_field = zeros(nx, ny, 3, 3);
    for i = 1:nx
        for j = 1:ny
            x_obs = [X(i,j), Y(i,j), z_obs];
            grad_u = integrateTriangleKelvinGradient(x_obs, vertices, n, stress_drop, mu, nu, tol);
            % Compute strain as the symmetric gradient
            eps = 0.5*(grad_u + grad_u');
            lam = 2*mu*nu/(1-2*nu);
            I3 = eye(3);
            sigma = lam*trace(eps)*I3 + 2*mu*eps;
            stress_field(i,j,:,:) = sigma;
        end
    end
    vonMises = zeros(nx, ny);
    for i = 1:nx
        for j = 1:ny
            sigma = squeeze(stress_field(i,j,:,:));
            sxx = sigma(1,1); syy = sigma(2,2); szz = sigma(3,3);
            sxy = sigma(1,2); sxz = sigma(1,3); syz = sigma(2,3);
            vonMises(i,j) = sqrt(0.5*((sxx-syy)^2 + (syy-szz)^2 + (szz-sxx)^2 + 6*(sxy^2+sxz^2+syz^2)));
        end
    end
    
    %% Plot displacement magnitude and von Mises stress
    figure;
    subplot(1,2,1);
    contourf(X, Y, u_magnitude, 20, 'LineColor', 'none');
    colorbar; title('Displacement Magnitude');
    xlabel('X (m)'); ylabel('Y (m)'); hold on;
    tri = [vertices(:,1:2); vertices(1,1:2)]; % close the triangle
    plot(tri(:,1), tri(:,2), 'r-', 'LineWidth', 2);
    
    subplot(1,2,2);
    contourf(X, Y, vonMises, 20, 'LineColor', 'none');
    colorbar; title('Von Mises Stress');
    xlabel('X (m)'); ylabel('Y (m)'); hold on;
    plot(tri(:,1), tri(:,2), 'r-', 'LineWidth', 2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = integrateTriangleKelvin(x, vertices, n, stress_drop, mu, nu, tol)
    % Computes the displacement at observation point x due to a uniform stress drop
    % over the triangle defined by vertices.
    v0 = vertices(1,:);
    v1 = vertices(2,:);
    v2 = vertices(3,:);
    J = norm(cross(v1-v0, v2-v0));
    
    % Integrate component-by-component using integral3.
    u = zeros(1,3);
    for i = 1:3
        % Wrap helperKelvinComponent with arrayfun so that scalar evaluation is forced.
        fi = @(u,w,z) arrayfun(@(u,w,z) helperKelvinComponent(x, v0, v1, v2, mu, nu, n, J, u, w, i), u, w, z);
        u(i) = stress_drop * integral3(fi, 0, 1, 0, 1, 0, 1, 'AbsTol', tol, 'RelTol', tol);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grad_u = integrateTriangleKelvinGradient(x, vertices, n, stress_drop, mu, nu, tol)
    % Computes the displacement gradient (3x3 matrix) at x due to the uniform stress drop
    % over the triangle defined by vertices.
    v0 = vertices(1,:);
    v1 = vertices(2,:);
    v2 = vertices(3,:);
    J = norm(cross(v1-v0, v2-v0));
    
    grad_u = zeros(3,3);
    for i = 1:3
        for k = 1:3
            fk = @(u,w,z) arrayfun(@(u,w,z) helper_dG_contract(x, v0, v1, v2, mu, nu, n, J, u, w, i, k), u, w, z);
            grad_u(i,k) = stress_drop * integral3(fk, 0, 1, 0, 1, 0, 1, 'AbsTol', tol, 'RelTol', tol);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = helperKelvinComponent(x, v0, v1, v2, mu, nu, n, J, u, w, i)
    % Helper function for integrateTriangleKelvin.
    % Computes the i-th component of G(x-y)*n, where:
    %   y = v0 + u.*(v1-v0) + (1-u).*w.*(v2-v0)
    % and multiplies by the Jacobian factor (J*(1-u)).
    y = v0 + u .* (v1 - v0) + (1 - u) .* w .* (v2 - v0);
    G = kelvinGreen(x - y, mu, nu);
    temp = G * n.';  
    val = temp(i) * (J * (1 - u));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = helper_dG_contract(x, v0, v1, v2, mu, nu, n, J, u, w, i, k)
    % Helper function for integrateTriangleKelvinGradient.
    % Computes the contraction: sum_j dG(i,j,k)*n(j),
    % where y = v0 + u.*(v1-v0) + (1-u).*w.*(v2-v0) and dG is the derivative of the Kelvin Green's function.
    y = v0 + u .* (v1 - v0) + (1 - u) .* w .* (v2 - v0);
    dG = dkelvinGreen(x - y, mu, nu);
    temp = squeeze(dG(i,:,k)) * n.'; 
    val = temp * (J * (1 - u));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = kelvinGreen(r, mu, nu)
    % Computes the Kelvin Green's function (3x3 tensor) at vector r.
    rnorm = norm(r);
    if rnorm < 1e-8
        G = zeros(3,3);
        return;
    end
    prefactor = 1/(16*pi*mu*(1-nu));
    G = prefactor * ( (3-4*nu)/rnorm * eye(3) + (r'*r)/(rnorm^3) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dG = dkelvinGreen(r, mu, nu)
    % Computes the spatial derivative of the Kelvin Green's function.
    % dG is a 3x3x3 tensor with dG(i,j,k) = dG_{ij}/dx_k.
    rnorm = norm(r);
    if rnorm < 1e-8
        dG = zeros(3,3,3);
        return;
    end
    prefactor = 1/(16*pi*mu*(1-nu));
    dG = zeros(3,3,3);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                delta_ij = double(i==j);
                delta_ik = double(i==k);
                delta_jk = double(j==k);
                dG(i,j,k) = prefactor * ( - (3-4*nu)*delta_ij*r(k)/rnorm^3 ...
                                 - (delta_ik*r(j)+delta_jk*r(i))/rnorm^3 ...
                                 + 3*r(i)*r(j)*r(k)/rnorm^5 );
            end
        end
    end
end