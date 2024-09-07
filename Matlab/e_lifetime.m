%% clear workspace
clear;
clc;
close all;
%% define variables
a           =              5.43e-10;        % lattice parameter (m)
nk          =            [12,12,12];        % number of kpoints along x y z
uo          =                     2;        % electron affinity (eV)
m_frac      =      [0.89,0.19,0.19];        % electron effective mass ratio (Kg)
v_frac      =                  0.05;        % volume fraction
ko          =     2*pi/a*[0.85,0,0];        % conduction band valley along (1/m)
del_k       =   2*pi/a*0.15*[1,1,1];        % BZ edge from valley 
n           =                   123;

%% spherical pore life time
% ro          =                  3e-9;        % pore radius (m)
% [~, E, tau_sphirical, N_spherical] = spherical(ro,nk,uo,m_frac,v_frac,ko,del_k,n);
%% cubic pore life time
% l           =          [1,2,3]*1e-9;        % pore length (m)
% [~, E, tau_cubic, N_cubic]     =     cubic(l,nk,uo,m_frac,v_frac,ko,del_k,n);
%% cylinder pore life time
% l           =                1*1e-9;        % pore length (m)
% ro          =                  3e-9;        % pore radius (m)
% [~, E, tau_cylinder, N_cylinder] = cylinder(l,ro,nk,uo,m_frac,v_frac,ko,del_k,n);
%% triangle pore life time
% l           =          [1,2,3]*1e-9;        % pore length (m)
% kind        =                     1;
% 
% [mag_kpoint, E, tau_triangle, N] = triangle(l,nk,uo,m_frac,v_frac,ko,del_k,n,kind);
%% inf cubic pore life time
% l           =              [1,1]*1e-9;        % pore length (m)
% nt          =                     2e4;
% [~, E, tau_inf_cubic, N_inf_cubic]     =     inf_cubic(l,nk,uo,m_frac,v_frac,ko,del_k,nt);
%% inf cylinder pore life time
% ro           =                 1*1e-9;        % pore length (m)
% nt          =                     2e4;
% [~, E, tau_inf_cylinder, N_inf_cylinder]     =    inf_cylinder(ro,nk,uo,m_frac,v_frac,ko,del_k,nt);
%% inf triangle pore life time
l           =          [1,2]*1e-9;        % pore length (m)
nt          =                 1e4;
kind        =                   2;
[mag_kpoint, E, tau_triangle, N] = inf_triangle(l,nk,uo,m_frac,v_frac,ko,del_k,nt,kind);
%% phonon life time
T           =                   300;
A0          =             220.5e-15;
tau_p       =         A0./sqrt(E*T);
semilogy(E,tau_triangle,'o')
hold on
semilogy(E,tau_p,'o')