function [mag_kpoint, E, tau, N] = cubic(l, nk, uo, m_frac, v_frac, ko, del_k, n)
%% Electron scattering rate from spherical potential wall in ellipsoid conduction band

%% Constants
hbar = 6.582119514e-16;        % Reduced Planck constant (eV.s)
eV2J = 1.60218e-19;            % Conversion factor from eV to Joules
me = 9.10938356e-31;           % Electron rest mass (Kg)
m = me * m_frac;               % Effective mass of electrons (Kg)
N = v_frac / prod(l);          % Number of particles in unit volume

%% Define k-points
kx = linspace(ko(1), ko(1) + del_k(1), nk(1));  % kx points
ky = linspace(ko(2), ko(2) + del_k(2), nk(2));  % ky points
kz = linspace(ko(3), ko(3) + del_k(3), nk(3));  % kz points
[xk, yk, zk] = meshgrid(kx, ky, kz);            % 3D grid of k-points

kpoint = [xk(:), yk(:), zk(:)];                 % Flatten into a list of k-points
mag_kpoint = vecnorm(kpoint, 2, 2);             % Magnitude of k-points

% Energy at each k-point
E = hbar^2 / 2 * (...
    (kpoint(:,1) - ko(1)).^2 / m(1) + ...
    (kpoint(:,2) - ko(2)).^2 / m(2) + ...
    (kpoint(:,3) - ko(3)).^2 / m(3)) * eV2J;

%% Precompute ellipsoid surface and scatter rate for each energy level
scattering_rate = zeros(1, size(E, 1));  % Preallocate for scattering rates

for u = 1:size(E, 1)
    % Generate ellipsoid surface for iso-energy
    [x, y, z] = generateEllipsoid(ko, m, E(u), n, hbar, eV2J);
    
    % Calculate centroids and areas of triangles on the ellipsoid surface
    [Q, A] = calculateCentroidsAndAreas(x, y, z, n);
    
    % Compute the scattering rate for this energy level
    scattering_rate(u) = computeScatteringRate(Q, kpoint(u, :), A, uo, l, ko, hbar, m, N);
end

%% Compute relaxation time (tau)
tau = 1 ./ scattering_rate * eV2J;  % Relaxation time (s)
end

%% Helper Functions

% Function to generate ellipsoid surface points based on energy level
function [x, y, z] = generateEllipsoid(ko, m, E, n, hbar, eV2J)
    xr = sqrt(2 / (hbar^2 * eV2J) * m(1) * E);
    yr = sqrt(2 / (hbar^2 * eV2J) * m(2) * E);
    zr = sqrt(2 / (hbar^2 * eV2J) * m(3) * E);
    [x, y, z] = ellipsoid(ko(1), ko(2), ko(3), xr, yr, zr, n);
end

% Function to calculate centroids and areas of the ellipsoid triangles
function [Q, A] = calculateCentroidsAndAreas(x, y, z, n)
    numTriangles = 2 * n * (n - 1);
    Q = zeros(numTriangles, 3);  % Centroids of the triangles
    A = zeros(numTriangles, 1);  % Surface area of the triangles
    
    k = 1;
    % Loop through ellipsoid surface to compute centroids and areas
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
end

% Function to compute the scattering rate
function scattering_rate = computeScatteringRate(Q, kpoint, A, uo, l, ko, hbar, m, N)
    % Calculate q = k' - k for each centroid
    qx = kpoint(1) - Q(:, 1);
    qy = kpoint(2) - Q(:, 2);
    qz = kpoint(3) - Q(:, 3);
    
    % Calculate cos(theta) for each triangle centroid
    cosTheta = dot(repmat(kpoint, size(Q, 1), 1), Q, 2) ./ (norm(kpoint) * sqrt(sum(Q.^2, 2)));
    
    % Compute matrix element M
    M = 8 * uo * (sin(l(1) * qx / 2) ./ qx) .* (sin(l(2) * qy / 2) ./ qy) .* (sin(l(3) * qz / 2) ./ qz);
    
    % Compute scattering rate ignoring delta(E - E')
    SR = 2 * pi / hbar * M .* conj(M);
    
    % Energy difference factor, del(E) = E(k') - E(k)
    delE = abs(hbar^2 * ((Q(:, 1) - ko(1)) / m(1) + (Q(:, 2) - ko(2)) / m(2) + (Q(:, 3) - ko(3)) / m(3)));
    
    % Compute the scattering rate integral
    f = SR ./ delE .* (1 - cosTheta);
    scattering_rate = N / (2 * pi)^3 * sum(f .* A);
end

% Function to calculate triangle properties: centroids and areas
function [centroid, area] = computeTriangleProperties(x, y, z, i1, j1, i2, j2, i3, j3)
    % Vertices of the triangle
    V1 = [x(i1, j1), y(i1, j1), z(i1, j1)];
    V2 = [x(i2, j2), y(i2, j2), z(i2, j2)];
    V3 = [x(i3, j3), y(i3, j3), z(i3, j3)];
    
    % Centroid of the triangle
    centroid = (V1 + V2 + V3) / 3;
    
    % Side lengths of the triangle
    a = norm(V1 - V2);
    b = norm(V2 - V3);
    c = norm(V3 - V1);
    
    % Area of the triangle using Heron's formula
    s = (a + b + c) / 2;
    area = sqrt(s * (s - a) * (s - b) * (s - c));
end
