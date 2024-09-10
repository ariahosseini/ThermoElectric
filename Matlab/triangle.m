function [mag_kpoint, E, tau, N] = triangle(l, nk, uo, m_frac, v_frac, ko, del_k, n, kind)
%% electron scattering rate from spherical symmetry potential wall in ellipsoid conduction band
% q = |k'-k|;
% Matrix element: M;
% SR matrix ignoring delta(E-E'): 2*pi/hbar*M.*conj(M);
% E = Ec + hbar^2/2*((kx-k0x)^2/ml+(ky^2+kz^2)/mt);
% Ec = 0;

% Constants
hbar = 6.582119514e-16;    % Reduced Planck constant (eV.s)
eV2J = 1.60218e-19;        % Unit conversion from eV to Joules
me = 9.10938356e-31;       % Electron rest mass (Kg)
m = me * m_frac;           % Effective mass of electrons (Kg)
N = v_frac / l(3) / (l(1) * l(2) / 2);  % Number of particles per unit volume

%% Define kpoints
kx = linspace(ko(1), ko(1) + del_k(1), nk(1));  
ky = linspace(ko(2), ko(2) + del_k(2), nk(2));  
kz = linspace(ko(3), ko(3) + del_k(3), nk(3));  

[xk, yk, zk] = meshgrid(kx, ky, kz);             
kpoint = [xk(:), yk(:), zk(:)];  % Matrix of kpoints
mag_kpoint = vecnorm(kpoint, 2, 2);  % Wavevector magnitude

%% Energy calculation (eV)
E = hbar^2 / 2 * (...
    (kpoint(:, 1) - ko(1)).^2 / m(1) + ...
    (kpoint(:, 2) - ko(2)).^2 / m(2) + ...
    (kpoint(:, 3) - ko(3)).^2 / m(3)) * eV2J; 

%% Preallocate Scattering Rate
scattering_rate = zeros(1, size(E, 1));

%% Ellipsoid and Numerical Integration (Centroid Method)
% number of triangles : 2*n*(n-1)
% [x,y,z] = ellipsoid(xc,yc,zc,xr,yr,zr,n) generates a surface mesh described by three n+1-by-n+1 matrices, 
% enabling surf(x,y,z) to plot an ellipsoid with center (xc,yc,zc) and semi-axis lengths (xr,yr,zr).
% Numerical integration on ellipsoid surface (iso energy surface) using cetroid method with triangle meshes
% a, b and c are the vertices of the triangles
% s = (a+b+c)/2
% A = sqrt(s*(s-a)*(s-b)*(s-c)); % triangle surface are in 3D  
for u = 1:size(E, 1)
    [x, y, z] = ellipsoid(ko(1), ko(2), ko(3), ...
        sqrt(2 / (hbar^2 * eV2J) * m(1) * E(u, 1)), ...
        sqrt(2 / (hbar^2 * eV2J) * m(2) * E(u, 1)), ...
        sqrt(2 / (hbar^2 * eV2J) * m(3) * E(u, 1)), n);
    
    % Triangle centroids and surface areas
    [Q, A] = computeCentroidsAndAreas(x, y, z, n);

    % Compute q = k' - k for each centroid
    qx = kpoint(u, 1) - Q(:, 1);
    qy = kpoint(u, 2) - Q(:, 2);
    qz = kpoint(u, 3) - Q(:, 3);

    % Cosine of scattering angle
    cosTheta = dot(repmat(kpoint(u, :), size(Q, 1), 1), Q, 2) ./ ...
        (norm(kpoint(u, :)) * sqrt(sum(Q.^2, 2)));

    % Calculate Matrix Element based on scattering type (kind)
    M = calculateMatrixElement(kind, l, qx, qy, qz, uo);

    % Scattering rate ignoring delta(E - E')
    SR = 2 * pi / hbar * M .* conj(M);

    % delta E(k')-E(k)
    delE = abs(hbar^2 * ((Q(:, 1) - ko(1)) / m(1) + (Q(:, 2) - ko(2)) / m(2) + (Q(:, 3) - ko(3)) / m(3)));

    % Integrand for scattering rate
    f = SR ./ delE .* (1 - cosTheta);

    % Scattering rate calculation
    scattering_rate(1, u) = N / (2 * pi)^3 * sum(f .* A);
end

% Electron inclusion relaxation time (s)
tau = 1 ./ scattering_rate * eV2J;
tau = real(tau);  % Ensure the output is real
end

%% Helper Functions

% Function to compute centroids and surface areas
function [Q, A] = computeCentroidsAndAreas(x, y, z, n)
    numTriangles = 2 * n * (n - 1);
    Q = zeros(numTriangles, 3);
    A = zeros(numTriangles, 1);
    k = 1;

    for j = 2:n
        for i = 3:n + 1
            [Q(k, :), A(k)] = calculateTriangle(x(i, j), y(i, j), z(i, j), x(i-1, j), y(i-1, j), z(i-1, j), x(i-1, j-1), y(i-1, j-1), z(i-1, j-1));
            k = k + 1;
        end
    end
end

% Function to calculate the centroid and area of a triangle
function [centroid, area] = calculateTriangle(x1, y1, z1, x2, y2, z2, x3, y3, z3)
    centroid = ( [x1, y1, z1] + [x2, y2, z2] + [x3, y3, z3] ) / 3;
    a = norm([x1, y1, z1] - [x2, y2, z2]);
    b = norm([x2, y2, z2] - [x3, y3, z3]);
    c = norm([x3, y3, z3] - [x1, y1, z1]);
    s = (a + b + c) / 2;
    area = sqrt(s * (s - a) * (s - b) * (s - c));
end

% Function to calculate the matrix element based on scattering kind
function M = calculateMatrixElement(kind, l, qx, qy, qz, uo)
    switch kind
        case 3
            M = 2 * uo * (...
                -1 * ((exp(l(1) * qx * 1i) * 1i - 1i) * 1i) ./ (qx .* qy) ...
                + (l(1) * (exp(l(1) * qx * 1i) .* exp(l(2) * qy * 1i) * 1i - 1i) * 1i) ./ ...
                (qy .* (l(1) * qx + l(2) * qy))) .* (sin(l(3) * qz / 2) ./ qz);
        case 2
            M = 2 * uo * (...
                (-((exp((l(1) * qx * 1i) / 2) * 1i - 1i) * 1i) ./ (qx .* qy) + ...
                (l(1) * (exp((l(1) * qx * 1i) / 2) .* exp(l(2) * qy * 1i) * 1i - 1i) * 1i) ./ ...
                (qy .* (l(1) * qx + 2 * l(2) * qy))) + ...
                ((exp((l(1) * qx * 3i) / 4) .* sin((l(1) * qx) / 4) * 2i) ./ (qx .* qy) + ...
                (l(1) * exp((l(1) * qx * 1i) / 2) .* (exp((l(1) * qx * 1i) / 2) * 1i - ...
                exp(l(2) * qy * 1i) * 1i) * 1i) ./ (qy .* (l(1) * qx - 2 * l(2) * qy)))) .* ...
                (sin(l(3) * qz / 2) ./ qz);
        case 1
            M = 2 * uo * (...
                -(qx .* (l(1) - l(1) .* exp(l(2) * qy * 1i)) - ...
                qy .* (l(2) - l(2) * exp(l(1) .* qx * 1i))) ./ ...
                (qx .* qy .* (l(1) * qx - l(2) * qy))) .* (sin(l(3) * qz / 2) ./ qz);
    end
end
