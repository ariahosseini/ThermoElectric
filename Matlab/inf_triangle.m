function [mag_kpoint, E, tau, N] = inf_triangle(l, nk, uo, m_frac, v_frac, ko, del_k, n, kind)
% Calculate electron scattering rate in an ellipsoidal conduction band

%% Constants
hbar    = 6.582119514e-16;       % Reduced Planck constant (eV.s)
eV2J    = 1.60218e-19;           % Conversion from eV to Joules
me      = 9.10938356e-31;        % Electron rest mass (kg)
m       = me * m_frac;           % Effective electron mass (kg)
N       = v_frac / (l(1) * l(2) / 2); % Number of particles per unit volume

%% Define k-points mesh
kx = linspace(ko(1), ko(1) + del_k(1), nk(1)); % kx mesh
ky = linspace(ko(2), ko(2) + del_k(2), nk(2)); % ky mesh
kz = linspace(ko(3), ko(3) + del_k(3), nk(3)); % kz mesh

[xk, yk, zk] = meshgrid(kx, ky, kz); % 3D grid of k-points
kpoint = [xk(:), yk(:), zk(:)];      % Flattened k-point matrix
mag_kpoint = vecnorm(kpoint, 2, 2);  % Magnitude of wavevector

%% Energy calculation
E = hbar^2 / 2 * ( ...
    (kpoint(:,1) - ko(1)).^2 / m(1) + ...
    (kpoint(:,2) - ko(2)).^2 / m(2) + ...
    (kpoint(:,3) - ko(3)).^2 / m(3)) * eV2J;  % Energy (Joules)

%% Parametrize the ellipse
t = linspace(0, 2 * pi, n);        % Parametric angles
a = sqrt(2 * m(1) / hbar^2 .* E / eV2J); % Semi-major axis
b = sqrt(2 * m(2) / hbar^2 .* E / eV2J); % Semi-minor axis
ds = sqrt((a .* sin(t)).^2 + (b .* cos(t)).^2); % Ellipse area element

%% Calculate cos(theta) and energy difference
cos_theta = (a .* kpoint(:,1) .* cos(t) + b .* kpoint(:,2) .* sin(t) + kpoint(:,3).^2) ./ ...
    sqrt(a.^2 .* cos(t).^2 + b.^2 .* sin(t).^2 + kpoint(:,3).^2) ./ mag_kpoint;

delE = hbar^2 * abs((a .* cos(t) - ko(1)) / m(1) + ...
    (b .* sin(t) - ko(2)) / m(2) + (kpoint(:,3) - ko(3)) / m(3));

%% Calculate qx and qy
qx = kpoint(:,1) - a .* cos(t);
qy = kpoint(:,2) - b .* sin(t);

%% Scattering rate based on potential type (kind)
switch kind
    case 3
        tmp = uo * (-1 * ((exp(l(1) * qx * 1i) * 1i - 1i) * 1i) ./ (qx .* qy) + ...
            (l(1) * (exp(l(1) * qx * 1i) .* exp(l(2) * qy * 1i) * 1i - 1i) * 1i) ./ ...
            (qy .* (l(1) * qx + l(2) * qy)));
        
    case 2
        tmp = uo * ((-((exp((l(1) * qx * 1i) / 2) * 1i - 1i) * 1i) ./ (qx .* qy)) + ...
            (l(1) * (exp((l(1) * qx * 1i) / 2) .* exp(l(2) * qy * 1i) * 1i - 1i) * 1i) ./ ...
            (qy .* (l(1) * qx + 2 * l(2) * qy))) + ...
            ((exp((l(1) * qx * 3i) / 4) .* sin((l(1) * qx) / 4) * 2i) ./ (qx .* qy)) + ...
            (l(1) * exp((l(1) * qx * 1i) / 2) .* (exp((l(1) * qx * 1i) / 2) * 1i - ...
            exp(l(2) * qy * 1i) * 1i) * 1i) ./ (qy .* (l(1) * qx - 2 * l(2) * qy))));

    case 1
        tmp = uo * (-(qx .* (l(1) - l(1) .* exp(l(2) * qy * 1i)) - ...
            qy .* (l(2) - l(2) * exp(l(1) .* qx * 1i))) ./ ...
            (qx .* qy .* (l(1) * qx - l(2) * qy)));
end

%% Scattering rate calculation
SR = 2 * pi / hbar * tmp .* conj(tmp);
func = SR .* (1 - cos_theta) ./ delE .* ds;

%% Numerical integration and scattering time (tau)
Int = trapz(t, func, 2);
tau = (N / (2 * pi)^3 .* Int).^-1 * eV2J;

end
