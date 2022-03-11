def tau_p(self, energyRange, alpha, Dv, DA, T, vs, D, rho):
    # Electron-phonon scattering rate using Ravich model
    # See the manual for the reference

    nonparabolic_term = (1 - ((alpha.T * energyRange) / (1 + 2 * alpha.T * energyRange) * (1 - Dv / DA))) ** 2 \
                        - 8 / 3 * (alpha.T * energyRange) * (1 + alpha.T * energyRange) / (
                                    1 + 2 * alpha.T * energyRange) ** 2 * (Dv / DA)  # Nonparabolic term in Ravich model

    tau = rho * vs ** 2 * thermoelectricProperties.hBar \
          / np.pi / thermoelectricProperties.kB / T.T / DA / DA * 1e9 / thermoelectricProperties.e2C / D  # Lifetime for parabolic band

    tau_p = tau / nonparabolic_term  # Lifetime in nonparabolic band

    return [tau, tau_p]  # The first row does not count for nonparabolicity, the second row does

    """
    In the following lines three models to predict electron-ion scattering rate is defined : 
    "tau_Screened_Coulomb", "tau_Unscreened_Coulomb", "tau_Strongly_Screened_Coulomb", ...
    the first one is the Brook-Herring model, the second one is for shallow dopants concentrationn up to  ~10^18 1/cm^3 ...
    (no screening effect is considered), and the last one is for strongly doped dielectrics.

    Note that for highly doped semiconductors, screen length plays a significant role, ...
    therefor should be computed carefully. Highly suggest to use following matlab file "Fermi.m" from:
    https://www.mathworks.com/matlabcentral/fileexchange/13616-fermi

    If committed to use python, the package "dfint" works with python2
    pip install fdint

    See the manual for details
    A good reference book on this topic: Fundamentals of Carrier Transport by Mark Lundstrom
    """


def tau_Screened_Coulomb(self, energyRange, m_c, LD, N):
    # Electron-ion scattering rate following Brook-Herring model

    g = 8 * m_c.T * LD.T ** 2 * energyRange / thermoelectricProperties.hBar ** 2 / thermoelectricProperties.e2C  # Gamma term

    var_tmp = np.log(1 + g) - g / (1 + g)  # tmp var.

    tau = 16 * np.pi * np.sqrt(2 * m_c.T) * (4 * np.pi * self.dielectric * thermoelectricProperties.e0) ** 2 \
          / N.T / var_tmp * energyRange ** (3 / 2) / thermoelectricProperties.e2C ** (
                      5 / 2)  # Brook-Herring model for electron-impurity scattering

    where_are_NaNs = np.isnan(tau)
    tau[where_are_NaNs] = 0

    return tau  # The array size is [1, numEnergySampling]


def tau_Unscreened_Coulomb(self, energyRange, m_c, N):
    # Electron-ion scattering rate for shallow dopants ~10^18 1/cm^3 (no screening effect is considered)

    g = 4 * np.pi * (4 * np.pi * self.dielectric * thermoelectricProperties.e0) * energyRange / N.T ** (
                1 / 3) / thermoelectricProperties.e2C  # Gamma term

    var_tmp = np.log(1 + g ** 2)  # tmp var.

    tau = 16 * np.pi * np.sqrt(2 * m_c.T) * (4 * np.pi * self.dielectric * thermoelectricProperties.e0) ** 2 \
          / N.T / var_tmp * energyRange ** (3 / 2) / thermoelectricProperties.e2C ** (
                      5 / 2)  # Electron-impurity scattering model fpr shallow doping

    where_are_NaNs = np.isnan(tau)
    tau[where_are_NaNs] = 0

    return tau  # The array size is [1, numEnergySampling]


def tau_Strongly_Screened_Coulomb(self, D, LD, N):
    tau = thermoelectricProperties.hBar / N.T / np.pi / D / \
          (LD.T ** 2 / (4 * np.pi * self.dielectric * thermoelectricProperties.e0)) ** 2 \
          * 1 / thermoelectricProperties.e2C ** 2  # Electron-impurity scattering model in highly doped dielectrics

    return tau  # The array size is [1, numEnergySampling]


def tau2D_cylinder(self, energyRange, nk, Uo, m, vfrac, valley, dk_len, ro, n=2000):
    """
    This is a fast algorithm that uses Fermi’s golden rule to compute the energy dependent electron scattering rate
    due cylindrical nanoparticles or pores extended perpendicular to the electrical current

    See manual for the detail
    """

    meff = np.array(m) * thermoelectricProperties.men  # Electron conduction nband effective mass
    ko = 2 * np.pi / self.latticeParameter * np.array(valley)
    del_k = 2 * np.pi / self.latticeParameter * dk_len * np.array([1, 1, 1])
    N = vfrac / np.pi / ro ** 2  # volume fraction/ porosity

    kx = np.linspace(ko[0], ko[0] + del_k[0], nk[0], endpoint=True)  # kpoints mesh
    ky = np.linspace(ko[1], ko[1] + del_k[1], nk[1], endpoint=True)  # kpoints mesh
    kz = np.linspace(ko[2], ko[2] + del_k[2], nk[2], endpoint=True)  # kpoints mesh
    [xk, yk, zk] = np.meshgrid(kx, ky, kz)
    xk_ = np.reshape(xk, -1)
    yk_ = np.reshape(yk, -1)
    zk_ = np.reshape(zk, -1)
    kpoint = np.array([xk_, yk_, zk_])  # Kpoint mesh sampling
    mag_kpoint = norm(kpoint, axis=0)

    E = thermoelectricProperties.hBar ** 2 / 2 * \
        ((kpoint[0, :] - ko[0]) ** 2 / meff[0] + (kpoint[1, :] - ko[1]) ** 2 / meff[1] +
         (kpoint[2, :] - ko[2]) ** 2 / meff[
             2]) * thermoelectricProperties.e2C  # Energy levels in ellipsoidal band structure

    # Write the ellips shape in parametric form

    t = np.linspace(0, 2 * np.pi, n)
    a = np.expand_dims(np.sqrt(2 * meff[1] / thermoelectricProperties.hBar ** 2 *
                               E / thermoelectricProperties.e2C), axis=0)
    b = np.expand_dims(np.sqrt(2 * meff[2] / thermoelectricProperties.hBar ** 2 *
                               E / thermoelectricProperties.e2C), axis=0)

    ds = np.sqrt((a.T * np.sin(t)) ** 2 + (b.T * np.cos(t)) ** 2)

    cos_theta = ((a * kpoint[0]).T * np.cos(t) + (b * kpoint[1]).T * np.sin(t) +
                 np.expand_dims(kpoint[2] ** 2, axis=1)) / \
                np.sqrt(a.T ** 2 * np.cos(t) ** 2 + b.T ** 2 * np.sin(t) ** 2 +
                        np.expand_dims(kpoint[2] ** 2, axis=1)) / np.expand_dims(mag_kpoint, axis=1)

    delE = thermoelectricProperties.hBar ** 2 * \
           np.abs((a.T * np.cos(t) - ko[0]) / meff[0] +
                  (b.T * np.sin(t) - ko[1]) / meff[1] + (
                              np.expand_dims(kpoint[2] ** 2, axis=1) - ko[2] / meff[2]))  # Energy increment

    # qpints
    qx = np.expand_dims(kpoint[0], axis=1) - a.T * np.cos(t)
    qy = np.expand_dims(kpoint[1], axis=1) - b.T * np.sin(t)
    qr = np.sqrt(qx ** 2 + qy ** 2)

    tau = np.empty((len(ro), len(E)))

    for r_idx in np.arange(len(ro)):
        J = jv(1, ro[r_idx] * qr)  # Bessel func.
        SR = 2 * np.pi / thermoelectricProperties.hBar * Uo ** 2 * (2 * np.pi) ** 3 * (
                    ro[r_idx] * J / qr) ** 2  # Scattering rate
        f = SR * (1 - cos_theta) / delE * ds
        int_ = np.trapz(f, t, axis=1)
        tau[r_idx] = 1 / (N[r_idx] / (2 * np.pi) ** 3 * int_) * thermoelectricProperties.e2C

    Ec, indices, return_indices = np.unique(E, return_index=True, return_inverse=True)

    tau_c = np.empty((len(ro), len(indices)))

    tauFunctionEnergy = np.empty((len(ro), len(energyRange[0])))

    for r_idx in np.arange(len(ro)):
        tau_c[r_idx] = accum(return_indices, tau[r_idx], func=np.mean, dtype=float)

    # Map lifetime to desired energy range
    for tau_idx in np.arange(len(tau_c)):
        ESpline = PchipInterpolator(Ec[30:], tau_c[tau_idx, 30:])
        tauFunctionEnergy[tau_idx] = ESpline(energyRange)

    return tauFunctionEnergy


def tau3D_spherical(self, energyRange, nk, Uo, m, vfrac, valley, dk_len, ro, n=32):
    """
    This is a fast algorithm that uses Fermi’s golden rule to compute the energy dependent electron scattering rate
    due spherical nanoparticles or pores.

    See manual for the detail
    """

    meff = np.array(m) * thermoelectricProperties.me  # Electron conduction nband effective mass
    ko = 2 * np.pi / self.latticeParameter * np.array(valley)
    del_k = 2 * np.pi / self.latticeParameter * dk_len * np.array([1, 1, 1])

    N = 3 * vfrac / 4 / np.pi / ro ** 3  # volume fraction/ porosity

    kx = np.linspace(ko[0], ko[0] + del_k[0], nk[0], endpoint=True)  # kpoints mesh
    ky = np.linspace(ko[1], ko[1] + del_k[1], nk[1], endpoint=True)  # kpoints mesh
    kz = np.linspace(ko[2], ko[2] + del_k[2], nk[2], endpoint=True)  # kpoints mesh
    [xk, yk, zk] = np.meshgrid(kx, ky, kz)
    xk_ = np.reshape(xk, -1)
    yk_ = np.reshape(yk, -1)
    zk_ = np.reshape(zk, -1)

    kpoint = np.array([xk_, yk_, zk_])  # Kpoint mesh sampling
    mag_kpoint = norm(kpoint, axis=0)

    # Energy levels in ellipsoidal band structure
    E = thermoelectricProperties.hBar ** 2 / 2 * \
        ((kpoint[0, :] - ko[0]) ** 2 / meff[0] +
         (kpoint[1, :] - ko[1]) ** 2 / meff[1] +
         (kpoint[2, :] - ko[2]) ** 2 / meff[2]) * thermoelectricProperties.e2C

    scattering_rate = np.zeros((len(ro), len(E)))

    nu = np.linspace(0, np.pi, n)
    z_ = -1 * np.cos(nu)

    r = np.sqrt(1.0 - z_ ** 2)[:, None]
    theta = np.linspace(0, 2 * np.pi, n)[None, :]

    x_ = r * np.cos(theta)
    y_ = r * np.sin(theta)

    # Mesh energy ellipsiod in triangular elements

    for u in np.arange(len(E)):

        Q = np.zeros((2 * (n - 2) * (n - 1), 3))
        A = np.zeros((2 * (n - 2) * (n - 1), 1))
        k = 0
        a_axis = np.sqrt(2 / (thermoelectricProperties.hBar ** 2 * thermoelectricProperties.e2C) * meff[0] * E[u])
        b_axis = np.sqrt(2 / (thermoelectricProperties.hBar ** 2 * thermoelectricProperties.e2C) * meff[1] * E[u])
        c_axis = np.sqrt(2 / (thermoelectricProperties.hBar ** 2 * thermoelectricProperties.e2C) * meff[2] * E[u])

        y = -1 * b_axis * y_ + ko[1]
        x = -1 * a_axis * x_ + ko[0]
        Z_ = c_axis * z_ + ko[2]
        z = np.tile(Z_[:, None], (1, n))
        for j in np.arange(1, n - 1):
            for i in np.arange(2, n):
                S = np.array(np.array([x[i, j], y[i, j], z[i, j]]) +
                             np.array([x[i - 1, j], y[i - 1, j], z[i - 1, j]]) +
                             np.array([x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]]))
                Q[k] = S / 3
                a = norm(np.array([x[i, j], y[i, j], z[i, j]]) - np.array([x[i - 1, j], y[i - 1, j], z[i - 1, j]]))
                b = norm(np.array([x[i - 1, j], y[i - 1, j], z[i - 1, j]]) - np.array(
                    [x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]]))
                c = norm(np.array([x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]]) - np.array(
                    [x[i, j], y[i, j], z[i, j]]))
                s = a + b + c
                s = s / 2
                A[k] = np.sqrt(s * (s - a) * (s - b) * (s - c))  # surface area of the triangular mesh elements
                k += 1
        for j in np.arange(1, n - 1):
            for i in np.arange(1, n - 1):
                S = np.array([x[i, j - 1], y[i, j - 1], z[i, j - 1]]) + \
                    np.array([x[i, j], y[i, j], z[i, j]]) + \
                    np.array([x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]])

                Q[k] = S / 3

                a = norm(np.array([x[i, j - 1], y[i, j - 1], z[i, j - 1]]) - np.array([x[i, j], y[i, j], z[i, j]]))
                b = norm(np.array([x[i, j], y[i, j], z[i, j]]) - np.array(
                    [x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]]))
                c = norm(np.array([x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]]) - np.array(
                    [x[i, j - 1], y[i, j - 1], z[i, j - 1]]))
                s = a + b + c
                s = s / 2

                A[k] = np.sqrt(s * (s - a) * (s - b) * (s - c))
                k += 1

        for i in np.arange(2, n):
            S = np.array([x[i, 0], y[i, 0], z[i, 0]]) + np.array([x[i - 1, 0], y[i - 1, 0], z[i - 1, 0]]) + np.array(
                [x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]])
            Q[k] = S / 3

            a = norm(np.array([x[i, 0], y[i, 0], z[i, 0]]) - np.array([x[i - 1, 0], y[i - 1, 0], z[i - 1, 0]]))
            b = norm(np.array([x[i - 1, 0], y[i - 1, 0], z[i - 1, 0]]) - np.array(
                [x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]]))
            c = norm(np.array([x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]]) - np.array([x[i, 0], y[i, 0], z[i, 0]]))
            s = a + b + c
            s = s / 2

            A[k] = np.sqrt(s * (s - a) * (s - b) * (s - c))
            k += 1

        for i in np.arange(1, n - 1):
            S = np.array([x[i, -2], y[i, -2], z[i, -2]]) + np.array([x[i, 0], y[i, 0], z[i, 0]]) + np.array(
                [x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]])
            Q[k] = S / 3

            a = norm(np.array([x[i, -2], y[i, -2], z[i, -2]]) - np.array([x[i, 0], y[i, 0], z[i, 0]]))
            b = norm(np.array([x[i, 0], y[i, 0], z[i, 0]]) - np.array([x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]]))
            c = norm(np.array([x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]]) - np.array([x[i, -2], y[i, -2], z[i, -2]]))
            s = a + b + c
            s = s / 2

            A[k] = np.sqrt(s * (s - a) * (s - b) * (s - c))
            k += 1

        qx = kpoint[0, u] - Q[:, 0]
        qy = kpoint[1, u] - Q[:, 1]
        qz = kpoint[2, u] - Q[:, 2]
        q = np.sqrt(qx ** 2 + qy ** 2 + qz ** 2)

        cosTheta = np.matmul(kpoint[:, u][None, :], Q.T) / norm(kpoint[:, u]) / np.sqrt(np.sum(Q ** 2, axis=1))

        delE = np.abs(thermoelectricProperties.hBar ** 2 * (
                    (Q[:, 0] - ko[0]) / meff[0] + (Q[:, 1] - ko[1]) / meff[1] + (Q[:, 2] - ko[2]) / meff[2]))

        for ro_idx in np.arange(len(ro)):
            M = 4 * np.pi * Uo * (1 / q * np.sin(ro[ro_idx] * q) - ro[ro_idx] * np.cos(ro[ro_idx] * q)) / (
                        q ** 2)  # Matrix element
            SR = 2 * np.pi / thermoelectricProperties.hBar * M * np.conj(M)  # Scattering rate
            f = SR / delE * (1 - cosTheta)
            scattering_rate[ro_idx, u] = N[ro_idx] / (2 * np.pi) ** 3 * np.sum(f * A.T)

    return scattering_rate  # Electorn scattering rate from the spherical pores/ nanoparticles