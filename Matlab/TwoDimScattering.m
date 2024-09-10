function [mag_kpoint, E, tau, N] = TwoDimScattering(l, nk, uo, m_frac, v_frac, ko, del_k, n, geom_type, kind, ro)
% Electron scattering rate for different geometries (cubic, cylinder, triangle)

% Constants
hbar = 6.582119514e-16;  % Reduced Planck constant (eV.s)
eV2J = 1.60218e-19;      % Conversion factor from eV to Joules
me = 9.10938356e-31;     % Electron rest mass (kg)
m = me * m_frac;         % Effective electron mass (kg)

% Particle number density based on geometry
N = compute_particle_density(l, v_frac, geom_type, ro);

% k-points mesh and magnitude calculation
[kpoint, mag_kpoint] = define_k_points(ko, del_k, nk);

% Energy calculation
E = compute_energy(kpoint, ko, m, hbar, eV2J);

% Compute scattering time based on geometry type
tau = compute_scattering_time(geom_type, kpoint, mag_kpoint, uo, N, l, m, E, hbar, eV2J, ko, n, kind, ro);

end

function N = compute_particle_density(l, v_frac, geom_type, ro)
% Calculate particle number density based on geometry
switch geom_type
    case 'cubic'
        N = v_frac / prod(l);  % For cubic
    case 'cylinder'
        N = v_frac / (pi * ro^2);  % For cylinder
    case 'triangle'
        N = v_frac / (l(1) * l(2) / 2);  % For triangular
    otherwise
        error('Unknown geometry type');
end
end

function [kpoint, mag_kpoint] = define_k_points(ko, del_k, nk)
% Define k-points mesh and compute magnitude of k-points
kx = linspace(ko(1), ko(1) + del_k(1), nk(1));  % kx mesh
ky = linspace(ko(2), ko(2) + del_k(2), nk(2));  % ky mesh
kz = linspace(ko(3), ko(3) + del_k(3), nk(3));  % kz mesh

[xk, yk, zk] = meshgrid(kx, ky, kz);            % Create 3D grid of k-points
kpoint = [xk(:), yk(:), zk(:)];                 % Flatten k-point matrix
mag_kpoint = vecnorm(kpoint, 2, 2);             % Magnitude of wavevector
end

function E = compute_energy(kpoint, ko, m, hbar, eV2J)
% Compute energy at k-points based on effective mass and constants
E = (hbar^2 / 2) * (...
    (kpoint(:, 1) - ko(1)).^2 / m(1) + ...
    (kpoint(:, 2) - ko(2)).^2 / m(2) + ...
    (kpoint(:, 3) - ko(3)).^2 / m(3)) * eV2J;   % Energy in Joules
end

function tau = compute_scattering_time(geom_type, kpoint, mag_kpoint, uo, N, l, m, E, hbar, eV2J, ko, n, kind, ro)
% Compute scattering time based on geometry

% Parametric angles and ellipse parameters
t = linspace(0, 2 * pi, n);                     
a = sqrt(2 * m(1) / hbar^2 .* E / eV2J);        
b = sqrt(2 * m(2) / hbar^2 .* E / eV2J);        
ds = sqrt((a .* sin(t)).^2 + (b .* cos(t)).^2); % Ellipse area element

% Calculate cos(theta) and energy difference
cos_theta = compute_cos_theta(a, b, kpoint, ko, mag_kpoint, t);
delE = compute_energy_difference(a, b, kpoint, ko, m, hbar);

% Geometry-specific scattering rates
switch geom_type
    case 'cubic'
        qx = kpoint(:, 1) - a .* cos(t);
        qy = kpoint(:, 2) - b .* sin(t);
        qx(qx == 0) = eps;
        qy(qy == 0) = eps;
        SR = 2 * pi / hbar * uo^2 .* (sin(l(1) * qx / 2) ./ qx).^2 .* ...
                                      (sin(l(2) * qy / 2) ./ qy).^2;
    
    case 'cylinder'
        qx = kpoint(:, 1) - a .* cos(t);
        qy = kpoint(:, 2) - b .* sin(t);
        qr = sqrt(qx.^2 + qy.^2);
        J = besselj(1, ro * qr);  % Bessel function of the first kind
        SR = 2 * pi / hbar * uo^2 * (2 * pi)^3 * (ro * J ./ qr).^2;
    
    case 'triangle'
        qx = kpoint(:,1) - a .* cos(t);
        qy = kpoint(:,2) - b .* sin(t);
        SR = scattering_rate_triangle(qx, qy, l, uo, kind);
        
    otherwise
        error('Unknown geometry type');
end

% Scattering function and numerical integration
func = SR .* (1 - cos_theta) ./ delE .* ds;
Int = trapz(t, func, 2);  % Numerical integration
tau = (N / (2 * pi)^3 .* Int).^-1 * eV2J;  % Scattering time
end


function cos_theta = compute_cos_theta(a, b, kpoint, ko, mag_kpoint, t)
% Calculate cos(theta)
cos_theta = (a .* kpoint(:, 1) .* cos(t) + b .* kpoint(:, 2) .* sin(t) + kpoint(:, 3).^2) ./ ...
    sqrt(a.^2 .* cos(t).^2 + b.^2 .* sin(t).^2 + kpoint(:, 3).^2) ./ mag_kpoint;
end


function delE = compute_energy_difference(a, b, kpoint, ko, m, hbar)
% Calculate energy difference between initial and scattered states
delE = hbar^2 * abs((a .* cos(t) - ko(1)) / m(1) + ...
    (b .* sin(t) - ko(2)) / m(2) + (kpoint(:, 3) - ko(3)) / m(3));
end


function SR = scattering_rate_triangle(qx, qy, l, uo, kind)
% Calculate scattering rate for triangular geometry
switch kind
    case 1
        tmp = uo * (-(qx .* (l(1) - l(1) .* exp(l(2) * qy * 1i)) - ...
            qy .* (l(2) - l(2) * exp(l(1) .* qx * 1i))) ./ ...
            (qx .* qy .* (l(1) * qx - l(2) * qy)));

    case 2
        tmp = uo * ((-((exp((l(1) * qx * 1i) / 2) * 1i - 1i) * 1i) ./ (qx .* qy)) + ...
            (l(1) * (exp((l(1) * qx * 1i) / 2) .* exp(l(2) * qy * 1i) * 1i - 1i) * 1i) ./ ...
            (qy .* (l(1) * qx + 2 * l(2) * qy))) + ...
            ((exp((l(1) * qx * 3i) / 4) .* sin((l(1) * qx) / 4) * 2i) ./ (qx .* qy)) + ...
            (l(1) * exp((l(1) * qx * 1i) / 2) .* (exp((l(1) * qx * 1i) / 2) * 1i - ...
            exp(l(2) * qy * 1i) * 1i) * 1i) ./ (qy .* (l(1) * qx - 2 * l(2) * qy)));
        
    case 3
        tmp = uo * (-1 * ((exp(l(1) * qx * 1i) * 1i - 1i) * 1i) ./ (qx .* qy) + ...
            (l(1) * (exp(l(1) * qx * 1i) .* exp(l(2) * qy * 1i) * 1i - 1i) * 1i) ./ ...
            (qy .* (l(1) * qx + l(2) * qy)));
        
    otherwise
        error('Unknown kind for triangular geometry');
end

% Scattering rate
SR = 2 * pi / hbar * tmp .* conj(tmp);
end

