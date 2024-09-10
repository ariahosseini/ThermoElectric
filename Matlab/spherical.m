function [mag_kpoint, E, tau, N] = spherical(ro, nk, uo, m_frac, v_frac, ko, del_k, n)
%% Electron scattering rate from spherical symmetry potential wall in ellipsoid conduction band;
% q = |k'-k|;
% Matrix element: M = 4*pi*u0*(1./q.*sin(r_inc*q)-r_inc*cos(r_inc*q))./(q.^2);
% SR matrix ignoring delta(E-E'): SR = 2*pi/hbar*M.*conj(M);
% E = Ec + hbar^2/2*((kx-k0x)^2/ml+(ky^2+kz^2)/mt);
% Ec = 0;

% Constants
hbar = 6.582119514e-16;    % Reduced Planck constant (eV.s)
eV2J = 1.60218e-19;        % Unit conversion from eV to Joules
me = 9.10938356e-31;       % Electron rest mass (Kg)
m = me * m_frac;           % Effective mass of electrons (Kg)
N = 3 * v_frac / (4 * pi * ro^3);  % Number of particles per unit volume

%% Define kpoints
kx = linspace(ko(1), ko(1) + del_k(1), nk(1));  
ky = linspace(ko(2), ko(2) + del_k(2), nk(2));  
kz = linspace(ko(3), ko(3) + del_k(3), nk(3));  

[xk, yk, zk] = meshgrid(kx, ky, kz);             
kpoint = [xk(:), yk(:), zk(:)];  % Matrix of the meshing kpoints [kx, ky, kz]
mag_kpoint = vecnorm(kpoint, 2, 2);  % Wavevector magnitude

%% Energy calculation (eV)
E = hbar^2 / 2 * (...
    (kpoint(:, 1) - ko(1)).^2 / m(1) + ...
    (kpoint(:, 2) - ko(2)).^2 / m(2) + ...
    (kpoint(:, 3) - ko(3)).^2 / m(3)) * eV2J;

%% Preallocate Scattering Rate
scattering_rate = zeros(1, size(E, 1));

%% Ellipsoid and Numerical Integration (Centroid Method)
% Number of triangles : 2*n*(n-1);
% [x,y,z] = ellipsoid(xc,yc,zc,xr,yr,zr,n) generates a surface mesh described by three n+1-by-n+1 matrices;
% Enable surf(x,y,z) to plot an ellipsoid centered at (xc,yc,zc) with semi-axis lengths of (xr,yr,zr);
% Numerical integration on ellipsoid surface (iso energy surface) using cetroid method with triangle meshes;
% Here, a, b and c are the vertices of the triangles;
% s = (a+b+c)/2
% A = sqrt(s*(s-a)*(s-b)*(s-c)); % triangle surface are in 3D    

for u = 1:size(E, 1)
    % Initialize centroids and surface areas
    Q = zeros(2 * n * (n - 1), 3);  % Centroid of the triangles
    A = zeros(2 * n * (n - 1), 1);  % Surface area of the triangles
    k = 1;
    
    % Generate ellipsoid surface mesh    
    [x, y, z] = ellipsoid(ko(1), ko(2), ko(3), ...
        sqrt(2 / (hbar^2 * eV2J) * m(1) * E(u, 1)), ...
        sqrt(2 / (hbar^2 * eV2J) * m(2) * E(u, 1)), ...
        sqrt(2 / (hbar^2 * eV2J) * m(3) * E(u, 1)), n);

    % Compute centroids and surface areas for the triangles
    [Q, A, k] = computeTriangleCentroidsAndAreas(x, y, z, n, Q, A, k);

    % Find q = |k' - k| for each centroid
    qx = kpoint(u, 1) - Q(:, 1);
    qy = kpoint(u, 2) - Q(:, 2);
    qz = kpoint(u, 3) - Q(:, 3);
    q = sqrt(qx.^2 + qy.^2 + qz.^2);

    % Calculate cos(theta)
    cosTheta = dot(repmat(kpoint(u, :), size(Q, 1), 1), Q, 2) ./ ...
        (norm(kpoint(u, :)) * sqrt(sum(Q.^2, 2)));

    % Matrix element M
    M = 4 * pi * uo * (1 ./ q .* sin(ro * q) - ro * cos(ro * q)) ./ (q.^2);

    % Scattering rate ignoring delta(E(k') - E(k))
    SR = 2 * pi / hbar * M .* conj(M);

    % delta E(k') - E(k)
    delE = abs(hbar^2 * ((Q(:, 1) - ko(1)) / m(1) + (Q(:, 2) - ko(2)) / m(2) + (Q(:, 3) - ko(3)) / m(3)));

    % Integrand for scattering rate
    f = SR ./ delE .* (1 - cosTheta);

    % Scattering rate
    scattering_rate(1, u) = N / (2 * pi)^3 * sum(f .* A);
end

% Electron inclusion relaxation time (s)
tau = 1 ./ scattering_rate * eV2J;
end

%% Helper Function to Compute Centroids and Areas of Triangles

function [Q, A, k] = computeTriangleCentroidsAndAreas(x, y, z, n, Q, A, k)
    % Loop through the ellipsoid to compute centroids and areas for each triangle
    for j = 2:n
        for i = 3:n+1
            S1 = [x(i,j), y(i,j), z(i,j)] + [x(i-1,j), y(i-1,j), z(i-1,j)] + [x(i-1,j-1), y(i-1,j-1), z(i-1,j-1)];
            % Centroid
            Q(k, :) = S1 / 3;
            % Surface area
            A(k, 1) = calculateTriangleArea(x, y, z, i, j);
            k = k + 1;
        end
    end

    for j = 2:n
        for i = 2:n
            S2 = [x(i,j-1), y(i,j-1), z(i,j-1)] + [x(i,j), y(i,j), z(i,j)] + [x(i-1,j-1), y(i-1,j-1), z(i-1,j-1)];
            % Centroid
            Q(k, :) = S2 / 3;
            % Surface area
            A(k, 1) = calculateTriangleArea(x, y, z, i, j);
            k = k + 1;
        end
    end

    % Handle boundary triangles
    for i = 3:n+1
        S1 = [x(i,1), y(i,1), z(i,1)] + [x(i-1,1), y(i-1,1), z(i-1,1)] + [x(i-1,end-1), y(i-1,end-1), z(i-1,end-1)];
        Q(k, :) = S1 / 3;
        A(k, 1) = calculateTriangleArea(x, y, z, i, 1);
        k = k + 1;
    end

    for i = 2:n
        S2 = [x(i,end-1), y(i,end-1), z(i,end-1)] + [x(i,1), y(i,1), z(i,1)] + [x(i-1,end-1), y(i-1,end-1), z(i-1,end-1)];
        Q(k, :) = S2 / 3;
        A(k, 1) = calculateTriangleArea(x, y, z, i, end-1);
        k = k + 1;
    end
end

% Function to calculate the surface area of a triangle
function area = calculateTriangleArea(x, y, z, i, j)
    a = norm([x(i,j), y(i,j), z(i,j)] - [x(i-1,j), y(i-1,j), z(i-1,j)]);
    b = norm([x(i-1,j), y(i-1,j), z(i-1,j)] - [x(i-1,j-1), y(i-1,j-1), z(i-1,j-1)]);
    c = norm([x(i-1,j-1), y(i-1,j-1), z(i-1,j-1)] - [x(i,j), y(i,j), z(i,j)]);
    s = (a + b + c) / 2;
    area = sqrt(s * (s - a) * (s - b) * (s - c));
end



