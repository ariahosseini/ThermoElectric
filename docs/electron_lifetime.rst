Electron Lifetime
=========================================================

Intrinsic scattering terms
------------------------------------------------------------------

The electron-phonon and electron-impurity scattering rates are computed as

::

    bulk_module = 98  # Bulk module (GPA)
    rho = 2329  # Mass density (Kg/m3)
    speed_sound = np.sqrt(bulk_module/rho)  # Speed of sound
    num_vally = 6

    tau_ph = TE.tau_p(energy = engr, alpha_term = nonparabolic_term , D_v = 2.94, D_a = 9.5, temp = tmpr,
                      vel_sound = speed_sound, DoS = e_density, rho = rho)
    tau_imp = TE.tau_strongly_screened_coulomb(DoS = e_density, screen_len = screen_len, n_imp = carrier_con,
                                                    dielectric = dielectric)
    tau = TE.matthiessen(engr, num_vally * tau_ph['nonparabolic_ph_lifetime'], num_vally * tau_imp)


The following equations are used

.. math::

    \mathrm{\tau_p(E)=\frac{\rho \nu_s^2 \hbar}{\pi \Phi_A^2 k_B T D(E)} \left ( \left[1-\frac{\alpha E}{1+2\alpha E}
            \left(1-\frac{\Phi_v}{\Phi_A} \right)\right]^2-\frac{8}{3} \frac{\alpha E(1+ \alpha E)}{(1+2 \alpha E)^2}
            \frac{D_v}{D_A} \right)^{-1}}

.. math::

    \mathrm{\tau_{ion}(E)=\frac{\hbar}{\pi N_i \left(\frac{e^2 L_D^2}{4\pi \epsilon \epsilon_o}\right)^2 D(E)}}

The Matthiessen's rule is used to add the scattering rates.

Scattering terms raised from nanoengineering
------------------------------------------------------------------

To compute the electron lifetime from the scattering processes raised from nano-engineering we used Fermi’s golden rule
to relate the transmission probability from the initial energy state to the distribution of final energy states for a
given time-invariant potential. In the case where energy conservation is imposed (elastic scattering)
the scattering rate in Born approximation can be written as

.. math::
    \tau^{-1}(s) = \frac{N}{(2\pi)^2\hbar}\int_{E(k')=0}\frac{M_{kk'}\overline{M}_{kk'}}{\nabla E(k')}(1-\cos\theta)dS(k').

Here, :math:`M_{kk'}` is the matrix element operator shows the coupling strength between initial and final
wave-functions and the number of ways the transmission may happen, :math:`N` is the number density of scattering source
and :math:`\theta` is the angle through which the electron's momentum is turned between the initial and scattered
states. For Bloch waves, :math:`M_{kk'}` is
defined as :math:`M_{kk'}= \int e^{i(k'-k).r} U(r)dr`. :math:`S(k')` represents iso-energic surface of electronic
states in k-space. For semiconductors with Si-like  band structures with indirect degenerate band gap, the contour
surfaces of k-states around the conduction band valley with energy :math:`E(k)` above the conduction band edge is
approximated well by an ellipsoid :math:`E(k)=\hbar^2 \left[\frac{(k_l-k_{ol} )^2}{2m_l^*} +\frac{k_t^2}{m_t^*}\right]`,
where :math:`k_l` and :math:`k_t` are the components of the wavevector that are parallel ant transverse to the long axis
of the conduction band valley. The term :math:`k_{ol}` describes the location of the conduction band minimum, while
:math:`m_l^*` and :math:`m_t^*` are the effective masses of electrons traveling along and transverse to the conduction
band valley, respectively. For silicon, :math:`m_l^*=0.98m_o` and :math:`m_t^*=0.19m_o` where :math:`m_o` is free
electron rest mass, and :math:`k_{ol}=0.85 2\pi/a` where :math:`a` is silicon's lattice parameter of 0.543 nm.

Model of Electron Lifetime for Scattering by Nanoparticles
------------------------------------------------------------------

The band alignment at the interface of nano-particles presents a barrier to electron transport equal to the conduction
band offset, :math:`\Delta\!E_c` between bulk silicon and the inclusions. For spherical nano-particles, the scattering
potential term, given as, :math:`U(r)=\Delta\!E_c \Pi(r_o-r)`, where :math:`r_o` is the nano-particle’s radius and
:math:`\Pi(r)` is a dimensionless boxcar function equal to unit inside and zero outside of the particles. For the
spherical symmetric potential, :math:`M_{kk'}` only depends on :math:`q=k'-k` and is defined as

.. math::
    M_{kk'}=\frac{4\pi \Delta\!E_c}{|q|^2}\left( \frac{1}{|q|}\sin\left(r_o|q|\right)-r_o\cos\left(r_o|q|\right)\right).

At equilibrium, the Fermi energy level of nano-particles and parent material aligned leaving the band offset between
SiC nano-particles and silicon, :math:`\Delta\!E_c`, equal to the difference between Fermi energy level and conduction
band edge of the SiC. For intrinsic semiconductors Fermi energy level is located at the middle of band gap so that
:math:`\Delta\!E_c=\frac{1}{2}E_g`. The SiC band gap varies from 2.36 eV at 300 K down to 2.036 eV at 1200 K following
(:math:`E_g  = 2.39-6.0\times10^{-4}\times \frac{T^2}{T+1200}`). Such a variation has negligible effect on scattering
rate so that we used temperature independent value of :math:`E_g =2.19\ eV` (and therefore :math:`\Delta\!E_c = 1.095\ eV`)
to model electron-nanoparticle scattering rate. Note that :math:`N` is the number density of nano-particles and is equal
to :math:`N=\frac{3\phi}{4\pi r_o^3}`, with :math:`\phi` the volume fraction of nano-particle. We have computed the
rates of electron scattering from SiC nano-particles

Model of Electron Lifetime for Scattering by Grain Boundaries
------------------------------------------------------------------
Along with the change in dopant concentration, the addition of 1% and 5% of SiC nano-particles results in a 22% and 40%
reduction in the grain size, respectively. It is known that grain boundaries can cause an electron filtering effect,
particularly if the boundaries include segregated species such as oxygen that provide centers for trapping charge
carriers. However, this effect only becomes significant in much smaller grain sizes. For our Si/SiC nanocomposites,
even with a 40% size reduction, the grains are still an order of magnitude larger than the average electron mean free
path in P-doped Si (which is only a few nanometers only at room temperature for carrier concentrations in excess of
:math:`10^{20}\ 1/cm^3`). Furthermore, we have computed the rate of electron scattering from grains (this is of special
importance in the next section where we evaluate the scope of enhancement in power factor in Si nanocomposites) using
the approach by Minnich et al. in which they have modeled grain boundaries as
decomposition of many local regions, each interacts independently with charge carriers and coherently scatters electron
waves. The model potential for grain boundaries in their work described as

.. math::
    U_{GB} =\left\{\begin{matrix}
    U_g e^{\frac{-|z|}{z_o}}& r<r_{GB} \\
    0& r>r_{GB}
    \end{matrix}\right.

In this equation, :math:`z` is the direction normal to the grain with :math:`z=0` at the center of the grain boundary,
:math:`r_{GB}` is a constant on the order of the screening length, and :math:`z_o` is a constant related to the size of
the depletion region. :math:`U_g` in this model is proposed as, :math:`U_g=\frac{e^2 N_t^2}{8 \epsilon \epsilon_o N_c}`,
where :math:`\epsilon` and :math:`\epsilon_o` are relative and vacuum permittivity, respectively, :math:`N_c` is the
doping concentration, and :math:`N_t` is the area density of traps. The matrix element of this potential is

.. math::
    M_{kk'}=4\pi U_g \left[ \frac{z_o}{1+(q_zz_o)^2} \right]r_o^2\left[ \frac{J_1(q_rr_o)}{q_rr_o} \right]

where :math:`J_1 (q_r r_o )` is the first-order Bessel function of the first kind, :math:`q=k-k'`, :math:`q_r` and
:math:`q_z` are the :math:`r` and :math:`z` component of :math:`q`, respectively. Unfortunately, there is a limit
information about the trap area density (:math:`N_t`) and the exact value of :math:`z_o` and :math:`r_o`. Nevertheless,
we know that depletion regime and the screening length are on the order of few nm.
We used :math:`N_t  = 10^{13}  \mathrm{\frac{1}{cm^2}}` for trap density of doped silicon,
:math:`z_o=1\ nm` and :math:`r_o=1\ nm`.