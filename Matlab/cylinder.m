function [mag_kpoint, E, tau, N] = cylinder(l,ro,nk,uo,m_frac,v_frac,ko,del_k,n)
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
N       =     v_frac/l/(pi*ro^2);              % number of particles in unit volume
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
%% ellipsoid
% number of triangles : 2*n*(n-1)
% [x,y,z] = ellipsoid(xc,yc,zc,xr,yr,zr,n) generates a surface mesh described by three n+1-by-n+1 matrices, 
% enabling surf(x,y,z) to plot an ellipsoid with center (xc,yc,zc) and semi-axis lengths (xr,yr,zr).
% Numerical integration on ellipsoid surface (iso energy surface) using
% cetroid method with triangle meshes
% a, b and c are the vertices of the triangles
% s = (a+b+c)/2
% A = sqrt(s*(s-a)*(s-b)*(s-c)); % triangle surface are in 3D    
scattering_rate = zeros(1,size(E,1));
for u = 1:size(E,1)

    Q           =                    zeros(2*n*(n-1),3);          % centroid of the triangles
    A           =                    zeros(2*n*(n-1),1);          % surface area of the triangles
    k           =                                     1;
[x, y, z]       =           ellipsoid(ko(1),ko(2),ko(3),...
                      sqrt(2/(hbar^2*eV2J)*m(1)*E(u,1)),...
                      sqrt(2/(hbar^2*eV2J)*m(2)*E(u,1)),...
                      sqrt(2/(hbar^2*eV2J)*m(3)*E(u,1)),...
                                 n);
                                  

for j = 2:n
    for i = 3:n+1
        S1 = [x(i,j),y(i,j),z(i,j)]+...
             [x(i-1,j),y(i-1,j),z(i-1,j)]+...
             [x(i-1,j-1),y(i-1,j-1),z(i-1,j-1)];
        % centroid
        Q(k,:) = S1/3;
        % surface area
        a = norm([x(i,j),y(i,j),z(i,j)]-[x(i-1,j),y(i-1,j),z(i-1,j)]);
        b = norm([x(i-1,j),y(i-1,j),z(i-1,j)]-[x(i-1,j-1),y(i-1,j-1),z(i-1,j-1)]);
        c = norm([x(i-1,j-1),y(i-1,j-1),z(i-1,j-1)]-[x(i,j),y(i,j),z(i,j)]);
        s = a+b+c;
        s = s/2;
        A(k,1) = sqrt(s*(s-a)*(s-b)*(s-c));
        % k
        k = k+1;
    end
end
for j = 2:n
    for i = 2:n

        S2 = [x(i,j-1),y(i,j-1),z(i,j-1)]+[x(i,j),y(i,j),z(i,j)]+...
            [x(i-1,j-1),y(i-1,j-1),z(i-1,j-1)];
        % centroid
        Q(k,:) = S2/3;
        % surface area
        a = norm([x(i,j-1),y(i,j-1),z(i,j-1)]-[x(i,j),y(i,j),z(i,j)]);
        b = norm([x(i,j),y(i,j),z(i,j)]-[x(i-1,j-1),y(i-1,j-1),z(i-1,j-1)]);
        c = norm([x(i-1,j-1),y(i-1,j-1),z(i-1,j-1)]-[x(i,j-1),y(i,j-1),z(i,j-1)]);
        s = a+b+c;
        s = s/2;
        A(k,1) = sqrt(s*(s-a)*(s-b)*(s-c));
        % k
        k = k+1;
    end
end
for i = 3:n+1
    
        S1 = [x(i,1),y(i,1),z(i,1)]+...
             [x(i-1,1),y(i-1,1),z(i-1,1)]+...
             [x(i-1,end-1),y(i-1,end-1),z(i-1,end-1)];
        % centroid
        Q(k,:) = S1/3;  
        % surface area
        a = norm([x(i,1),y(i,1),z(i,1)]-[x(i-1,1),y(i-1,1),z(i-1,1)]);
        b = norm([x(i-1,1),y(i-1,1),z(i-1,1)]-[x(i-1,end-1),y(i-1,end-1),z(i-1,end-1)]);
        c = norm([x(i-1,end-1),y(i-1,end-1),z(i-1,end-1)]-[x(i,1),y(i,1),z(i,1)]);
        s = a+b+c;
        s = s/2;
        A(k,1) = sqrt(s*(s-a)*(s-b)*(s-c));
        % k
        k = k+1;       
end
for i = 2:n

        S2 = [x(i,end-1),y(i,end-1),z(i,end-1)]+...
             [x(i,1),y(i,1),z(i,1)]+...
             [x(i-1,end-1),y(i-1,end-1),z(i-1,end-1)];   
        % centroid     
        Q(k,:) = S2/3;
        % surface area
        a = norm([x(i,end-1),y(i,end-1),z(i,end-1)]-[x(i,1),y(i,1),z(i,1)]);
        b = norm([x(i,1),y(i,1),z(i,1)]-[x(i-1,end-1),y(i-1,end-1),z(i-1,end-1)]);
        c = norm([x(i-1,end-1),y(i-1,end-1),z(i-1,end-1)]-[x(i,end-1),y(i,end-1),z(i,end-1)]);
        s = a+b+c;
        s = s/2;
        A(k,1) = sqrt(s*(s-a)*(s-b)*(s-c));
        % k
        k = k+1;
end

% find q = k'-k for each centroid
qx             =                 kpoint(u,1) - Q(:,1);
qy             =                 kpoint(u,2) - Q(:,2);
qz             =                 kpoint(u,3) - Q(:,3);
rq             =                  sqrt(qx.^2 + qy.^2);

% find cos(theta)
cosTheta       =                dot(repmat(kpoint(u,:)...
                                    ,size(Q,1),1),Q,2)...
                                  ./(norm(kpoint(u,:))...
                                   *sqrt(sum(Q.^2,2)));
% find matrix element
J              =                     besselj(1,ro*rq);
M              =  4*pi*uo*ro*J./rq.*(sin(l*qz/2)./qz);

% find scattering rate ingnoring del(E(k')-E(k))
SR             =                2*pi/hbar*M.*conj(M);
% G = E(k')-E(k), delG use to map volume integral to surface integral
delE           =           abs(hbar^2*(...
                                      (Q(:,1)-ko(1))/m(1)+...
                                      (Q(:,2)-ko(2))/m(2)+...
                                      (Q(:,3)-ko(3))/m(3) ...
                                      ));
% integrand
f                       =      SR./delE.*(1-cosTheta);
% scattering rate
scattering_rate(1,u)    =        N/(2*pi)^3*sum(f.*A);

end
tau            =       1./scattering_rate*eV2J; % Electron inclusion relaxation time 
end




