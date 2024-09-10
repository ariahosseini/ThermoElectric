function [mag_kpoint, E, tau, N] = inf_cubic(l, nk, uo, m_frac, v_frac, ko, del_k, n)
% Electron scattering rate from spherical symmetry potential in ellipsoid conduction band

%% Constants
hbar = 6.582119514e-16;  % Reduced Planck constant (eV.s)
eV2J = 1.60218e-19;      % Conversion factor from eV to Joules
me = 9.10938356e-31;     % Electron rest mass (kg)
m = me * m_frac;         % Effective electron mass (kg)
N = v_frac / prod(l);    % Particle number density in unit volume

%% Define k-points mesh
kx = linspace(ko(1), ko(1) + del_k(1), nk(1));  % kx mesh
ky = linspace(ko(2), ko(2) + del_k(2), nk(2));  % ky mesh
kz = linspace(ko(3), ko(3) + del_k(3), nk(3));  % kz mesh

[xk, yk, zk] = meshgrid(kx, ky, kz);            % Create 3D grid of k-points
kpoint = [xk(:), yk(:), zk(:)];                 % Flatten k-point matrix
mag_kpoint = vecnorm(kpoint, 2, 2);             % Magnitude of wavevector

%% Energy calculation
E = (hbar^2 / 2) * (...
    (kpoint(:, 1) - ko(1)).^2 / m(1) + ...
    (kpoint(:, 2) - ko(2)).^2 / m(2) + ...
    (kpoint(:, 3) - ko(3)).^2 / m(3)) * eV2J;   % Energy in Joules

%% Parametrize scattering ellipse
t = linspace(0, 2 * pi, n);                     % Parametric angles for ellipse
a = sqrt(2 * m(1) / hbar^2 .* E / eV2J);        % Semi-major axis
b = sqrt(2 * m(2) / hbar^2 .* E / eV2J);        % Semi-minor axis
ds = sqrt((a .* sin(t)).^2 + (b .* cos(t)).^2); % Ellipse area element

%% Calculate cos(theta) and energy difference
cos_theta = (a .* kpoint(:, 1) .* cos(t) + b .* kpoint(:, 2) .* sin(t) + kpoint(:, 3).^2) ./ ...
    sqrt(a.^2 .* cos(t).^2 + b.^2 .* sin(t).^2 + kpoint(:, 3).^2) ./ mag_kpoint;

delE = hbar^2 * abs((a .* cos(t) - ko(1)) / m(1) + ...
    (b .* sin(t) - ko(2)) / m(2) + (kpoint(:, 3) - ko(3)) / m(3));

%% Scattering rate calculation
qx = kpoint(:, 1) - a .* cos(t);
qy = kpoint(:, 2) - b .* sin(t);

% Avoid division by zero using a small tolerance (eps)
qx(qx == 0) = eps;
qy(qy == 0) = eps;

SR = 2 * pi / hbar * uo^2 .* (sin(l(1) * qx / 2) ./ qx).^2 .* ...
                             (sin(l(2) * qy / 2) ./ qy).^2; % Scattering rate

%% Calculate the scattering function and integrate
func = SR .* (1 - cos_theta) ./ delE .* ds;
Int = trapz(t, func, 2);  % Numerical integration

tau = (N / (2 * pi)^3 .* Int).^-1 * eV2J;  % Scattering time

end
