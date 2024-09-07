function [mag_kpoint, E, tau, N] = inf_cylinder(ro,nk,uo,m_frac,v_frac,ko,del_k,n)
%% electron scattering rate from spherical symmetry potential wall in ellipsoid conduction band
% S. Aria Hosseini
% 03/29/2020
% q = |k'-k|;
% Matrix element
% M = 4*pi*u0*(1./q.*sin(r_inc*q)-r_inc*cos(r_inc*q))./(q.^2)
% SR matrix ignoring delta(E-E')
% SR = 2*pi/hbar*M.*conj(M): Scattering rate functional
% E = Ec + hbar^2/2*((kx-k0x)^2/ml+(ky^2+kz^2)/mt)
% Ec = 0;
%% define variables
hbar    =        6.582119514e-16;              % reduced Planck constant (eV.s)
eV2J    =            1.60218e-19;              % unit change from eV to Joul
me      =         9.10938356e-31;              % electron rest mass (Kg)
m       =              me*m_frac;              % effective mass of electrons (Kg)
N       =         v_frac/pi/ro^2;              % number of particles in unit volume
%% define kpoints
kx           =          linspace(ko(1),ko(1)+del_k(1),nk(1));       % kpoints mesh
ky           =          linspace(ko(2),ko(2)+del_k(2),nk(2));       % kpoints mesh
kz           =          linspace(ko(3),ko(3)+del_k(3),nk(3));       % kpoints mesh
[xk,yk,zk]   =                            meshgrid(kx,ky,kz);             
kpoint       =                           [xk(:),yk(:),zk(:)];       % matrix of the meshing kpoints [kx, ky, kz]
mag_kpoint   =                           vecnorm(kpoint,2,2);       % wavevector magnitude

E            =         hbar^2/2*(...
                                 (kpoint(:,1)-ko(1)).^2/m(1)...
                                +(kpoint(:,2)-ko(2)).^2/m(2)...
                                +(kpoint(:,3)-ko(3)).^2/m(3)...
                                )*eV2J;                            % energy (eV)

t            =                            linspace(0,2*pi,n);      % meshing the ellipse
a            =                   sqrt(2*m(1)/hbar^2.*E/eV2J);      % ellispe semi major axis
b            =                   sqrt(2*m(2)/hbar^2.*E/eV2J);      % ellispe semi minor axis
ds           =           sqrt((a.*sin(t)).^2+(b.*cos(t)).^2);      % parametrize ellipse area

% find cos(theta)
cos_theta    = (a.*kpoint(:,1).*cos(t)+b.*kpoint(:,2).*sin(t)+kpoint(:,3).^2)./...
               sqrt(a.^2.*cos(t).^2+b.^2.*sin(t).^2+kpoint(:,3).^2)./...
                                                  mag_kpoint;
                                              
                                              
delE = hbar^2*abs((a.*cos(t)-ko(1))/m(1)+...
           (b.*sin(t)-ko(2))/m(2)+...
           (kpoint(:,3)-ko(3))/m(3));
       
qx = kpoint(:,1)-a.*cos(t);
qy = kpoint(:,2)-b.*sin(t);
qr = sqrt(qx.^2+qy.^2);

J = besselj(1,ro*qr);

SR              =              2*pi/hbar*uo^2*(2*pi)^3*(ro*J./qr).^2;
func = SR.*(1-cos_theta)./delE.*ds;

Int = trapz(t,func,2);
tau = (N/(2*pi)^3.*Int).^-1*eV2J;

end
