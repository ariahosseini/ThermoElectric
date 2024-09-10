function [mag_kpoint, E, tau, N] = cylinder(l, ro, nk, uo, m_frac, v_frac, ko, del_k, n)
%% Electron scattering rate from cylindrical potential wall in ellipsoid conduction band

%% Constants
hbar = 6.582119514e-16;        % Reduced Planck constant (eV.s)
eV2J = 1.60218e-19;            % Conversion factor from eV to Joules
me = 9.10938356e-31;           % Electron rest mass (Kg)
m = me * m_frac;               % Effective mass of electrons (Kg)
N = v_frac / l / (pi * ro^2);  % Number of particles in unit volume

%% Define k-points
kx = linspace(ko(1), ko(1) + del_k(1), nk(1));  % kx points
ky = linspace(ko(2), ko(2) + del_k(2), nk(2));  % ky points
kz = linspace(ko(3), ko(3) + del_k(3), nk(3));  % kz points
[xk, yk, zk] = meshgrid(kx, ky, kz);            % 3D grid of k-points

kpoint = [xk(:), yk(:), zk(:)];                 % Flatten into a list of k-points
mag_kpoint = vecnorm(kpoint, 2, 2);             % Magnitude of k-points

% Calculate the energy at each k-point
E = hbar^2 / 2 * ((kpoint(:, 1) - ko(1)).^2 / m(1) + ...
                  (kpoint(:, 2) - ko(2)).^2 / m(2) + ...
                  (kpoint(:, 3) - ko(3)).^2 / m(3)) * eV2J;

%% Numerical integration on the ellipsoid surface
% Preallocate arrays for scattering rate and relaxation time
scattering_rate = zeros(1, size(E, 1));

% Loop over each k-point and compute the scattering rate
for u = 1:size(E, 1)
    % Generate ellipsoid surface for iso-energy
    [x, y, z] = generateEllipsoid(ko, m, E(u), n, hbar, eV2J);
    
    % Calculate centroids and areas of the ellipsoid triangles
    [Q, A] = calculateCentroidsAndAreas(x, y, z, n);
    
    % Compute the scattering rate for this energy level
    scattering_rate(u) = computeScatteringRate(Q, kpoint(u, :), A, ro, uo, l, ko, hbar, m, N);
end

%% Relaxation time
tau = 1 ./ scattering_rate * eV2J; % Electron inclusion relaxation time
end

%% Helper Functions

function [x, y, z] = generateEllipsoid(ko, m, E, n, hbar, eV2J)
    % Generates ellipsoid surface points
    xr = sqrt(2 / (hbar^2 * eV2J) * m(1) * E);
    yr = sqrt(2 / (hbar^2 * eV2J) * m(2) * E);
    zr = sqrt(2 / (hbar^2 * eV2J) * m(3) * E);
    [x, y, z] = ellipsoid(ko(1), ko(2), ko(3), xr, yr, zr, n);
end

function [Q, A] = calculateCentroidsAndAreas(x, y, z, n)
    % Calculates centroids and areas of ellipsoid triangles
    numTriangles = 2 * n * (n - 1);
    Q = zeros(numTriangles, 3); % Centroids of the triangles
    A = zeros(numTriangles, 1); % Surface area of the triangles
    
    k = 1;
    
    % First loop: from (2:n, 3:n+1) and (2:n, 2:n)
    for j = 2:n
        for i = 3:n+1
            [Q(k, :), A(k)] = computeTriangleProperties(x, y, z, i, j, i-1, j, i-1, j-1);
            k = k + 1;
        end
    end
    for j = 2:n
        for i = 2:n
            [Q(k, :), A(k)] = computeTriangleProperties(x, y, z, i, j-1, i, j, i-1, j-1);
            k = k + 1;
        end
    end
    
    % Second loop: from (3:n+1, 1) and (2:n, end-1)
    for i = 3:n+1
        [Q(k, :), A(k)] = computeTriangleProperties(x, y, z, i, 1, i-1, 1, i-1, end-1);
        k = k + 1;
    end
    for i = 2:n
        [Q(k, :), A(k)] = computeTriangleProperties(x, y, z, i, end-1, i, 1, i-1, end-1);
        k = k + 1;
    end
end

function [centroid, area] = computeTriangleProperties(x, y, z, i1, j1, i2, j2, i3, j3)
    % Computes centroid and area of a triangle on the ellipsoid surface
    V1 = [x(i1, j1), y(i1, j1), z(i1, j1)];
    V2 = [x(i2, j2), y(i2, j2), z(i2, j2)];
    V3 = [x(i3, j3), y(i3, j3), z(i3, j3)];
    centroid = (V1 + V2 + V3) / 3;
    
    % Compute the area of the triangle
    a = norm(V1 - V2);
    b = norm(V2 - V3);
    c = norm(V3 - V1);
    s = (a + b + c) / 2;
    area = sqrt(s * (s - a) * (s - b) * (s - c));
end

function scattering_rate = computeScatteringRate(Q, kpoint, A, ro, uo, l, ko, hbar, m, N)
    % Computes the scattering rate using matrix elements and q vectors
    q = kpoint - Q;
    rq = sqrt(q(:, 1).^2 + q(:, 2).^2); % Magnitude of q (radial components)

    % Matrix element (using cylindrical Bessel function)
    J = besselj(1, ro * rq);
    M = 4 * pi * uo * ro * J ./ rq .* sin(l * q(:, 3) / 2) ./ q(:, 3);
    
    % Scattering rate functional (ignoring delta(E-E'))
    SR = 2 * pi / hbar * M .* conj(M);
    
    % Energy difference factor
    delE = abs(hbar^2 * (Q(:, 1) - ko(1)) / m(1) + ...
               hbar^2 * (Q(:, 2) - ko(2)) / m(2) + ...
               hbar^2 * (Q(:, 3) - ko(3)) / m(3));
    
    % Compute the integrand and scattering rate
    f = SR ./ delE .* (1 - dot(kpoint, Q, 2) ./ (norm(kpoint) * sqrt(sum(Q.^2, 2))));
    scattering_rate = N / (2 * pi)^3 * sum(f .* A);
end
