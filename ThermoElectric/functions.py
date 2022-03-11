"""Provide the primary functions."""

from .util import *
import numpy as np


def band_gap(Eg_o: float, Ao: float, Bo: float, temp: np.ndarray = None):

    """
    This method uses Eg(T)=Eg(T=0)-Ao*T**2/(T+Bo) to approximate the temperature dependency of the dielectrics band gap.
    A good reference is "Properties of Advanced Semiconductor Materials" by Michael E. Levinshtein.

    Parameters
    ----------
    Eg_o: float
        The band gap at zero Kelvin
    Ao: float
        Experimentally fitted parameter
    Bo: float
        Experimentally fitted parameter
    temp: np.ndarray
        Temperature range

    Returns
    -------
    Eg: np.ndarray
        Temperature-dependent band gap
    """

    if temp is None:
        '''Use default temperature range (300K up to 1300K with step of 100K) '''
        T = temperature()
    else:
        T = temp

    Eg = Eg_o - Ao * np.divide(T ** 2, T + Bo)

    return Eg


def analyticalDoS(range_energy: np.ndarray, electron_eff_mass: float, nonparabolic_term: np.ndarray):

    """
       This function approximate the electron density of state for parabolic and non-parabolic bands
       in case DFT calculation is not available.

    Parameters
    ----------
    range_energy: np.ndarray
        Electron energy range
    electron_eff_mass: float
        Electron effective mass
    nonparabolic_term: np.ndarray
        Non-parabolic term (shows the mixture of S and P orbitals)

    Returns
    -------
    DoS: np.ndarray
        First row is phonon density of states for non-parabolic band and
        the second row is the phonon density of states for non-parabolic.
        The array size is [2, numEnergySampling].
    """

    hBar = 6.582119e-16  # Reduced Planck constant in eV.s
    e2C = 1.6021765e-19  # e to Coulomb unit change

    DoS_nonparabolic = 1 / np.pi ** 2 * np.sqrt(2 * range_energy * (1 + range_energy * nonparabolic_term.T)) \
                       * np.sqrt(electron_eff_mass / hBar ** 2) ** 3 \
                       * (1 + (2 * range_energy * nonparabolic_term.T)) / e2C ** (3.0 / 2)

    DoS_parabolic = np.sqrt(range_energy) / np.pi ** 2 * np.sqrt(2) / hBar ** 3 \
                    * electron_eff_mass ** (3.0 / 2) / e2C ** (3.0 / 2)

    DoS = np.array(DoS_nonparabolic, DoS_parabolic)

    return DoS


def fermi_level(carrier, energy, density, Nc=None, Ao=None, temp=None):

    """
    This function uses Joice Dixon approximation to predict Ef and thereby the carrier concentration at each temperature
    A good reference book is "Principles of Semiconductor Devices" by Sima Dimitrijev.

    Parameters
    ----------
    carrier: np.ndarray
        Total carrier concentration
    energy: np.ndarray
        The electron energy level
    density: np.ndarray
        The electron density of states
    Ao: float
        Experimentally fitted parameter (Nc ~ Ao*T^(3/2))
    Nc: float
        The effective densities of states in the conduction band
    temp: np.ndarray
        Temperature range

    Returns
    -------
    output: np.ndarray
        The first roa is the Fermi level and the second one is the carrier concentration
    """

    k_bolt = 8.617330350e-5  # Boltzmann constant in eV/K

    if temp is None:
        T = temperature()
    else:
        T = temp
    if Ao is None and Nc is None:
        raise Exception("Either Ao or Nc should be defined")
    if Nc is None:
        Nc = Ao * temp ** (3.0 / 2)

    JD_CC = np.log(carrier/ Nc) + 1 / np.sqrt(8) * carrier / Nc - (3. / 16 - np.sqrt(3) / 9) * (carrier/Nc)**2
    fermi_energy = k_bolt * (T * JD_CC)  # Joice Dixon approximation of Ef
    f, _ = fermi_distribution(energy, fermi_energy, temp=T)  # Fermi distribution
    n = np.trapz(np.multiply(density, f), energy, axis=1)  # Carrier concentration

    output = np.array([fermi_energy, n])

    return output


def fermi_level_self_consistent(carrier, temp, energy, density, fermi_level):

    """
    A tool for self-consistent calculation of the Fermi level from a given carrier concentration
    using Joyce Dixon approximation as the initial guess for degenerate semiconductors.
    As a default value of 4000 sampling points in energy range from Ef(JD)-0.4 eV up to Ef(JD)+0.2 is considered.
    This looks reasonable in most cases. The index is printed out if it reaches the extreme index of (0) or (4000),
    increase energy range. Then increase sampling point number to finner results.

    Parameters
    ----------



        :arg
                carrierConcentration                           : Function object, total carrier concentration
                energyRange                                    : Function object, the electron energy level
                DoS                                            : Function object, the electron density of state
                fermilevel                                     : Function object, Joyce Dixon approximation as the initial guess
                Temp                                           : Function object, temperature range
        :returns
                [Ef,n]                                         : A 1 by 2 list, The first element is a NumPy array of Fermi level for each temperature ...
                                                                 while the second element is a Numpy array of the corresponding carrier concentration
        """


    fermi = np.linspace(fermi_level[0] - 0.4, fermi_level[0] + 0.2, 4000, endpoint=True).T  # Range of energy arounf Ef(JD )to consider

    result_array = np.empty((np.shape(temp)[1], np.shape(fermi)[1]))
    idx_j = 0
    for j in temp[0]:
        idx_i = 0
        for i in fermi[idx_j]:
            f, _ = fermi_distribution(energy, np.expand_dims(np.array([i]), axis=0),
                                      np.expand_dims(np.array([j]), axis=0))
            tmp = np.trapz(np.multiply(density, f), energy, axis=1)
            result_array[idx_j, idx_i] = tmp
            idx_i += 1
        idx_j += 1

    diff = np.tile(np.transpose(carrier), (1, np.shape(fermi)[1])) - abs(result_array)

    min_idx = np.argmin(np.abs(diff), axis=1)
    print("Fermi Level Self Consistent Index ", min_idx)

    Ef = np.empty((1, np.shape(temp)[1]))

    for Ef_idx in np.arange(len(min_idx)):
        Ef[0, Ef_idx] = fermi[Ef_idx, min_idx[Ef_idx]]
    elm = 0
    n = np.empty((1, np.shape(temp)[1]))
    for idx in min_idx:
        n[0, elm] = result_array[elm, idx]
        elm += 1

    return [Ef,n]  # The array size is [2, size(temp)], The first row is the Fermi and the second row is the carrier concentration


    def electronGroupVelocity(self, kp, energy_kp, energyRange):

        # This is the derivation of band structure from DFT.
        # BTE needs single band data. Reciprocal lattice vector is needed, ...
        # See the example (si.py) or the manual for the details.

        dE = np.roll(energy_kp, -1, axis=0) - np.roll(energy_kp, 1, axis=0)
        dk = np.roll(kp, -1, axis=0) - np.roll(kp, 1, axis=0)
        dEdk = np.divide(dE, dk)
        dEdk[0] = (energy_kp[1] - energy_kp[0]) / (kp[1] - kp[0])
        dEdk[-1] = (energy_kp[-1] - energy_kp[-2]) / (kp[-1] - kp[-2])

        dEdkSpline = InterpolatedUnivariateSpline(energy_kp, np.array(dEdk))
        dEdkFunctionEnergy = dEdkSpline(energyRange)

        groupVel = dEdkFunctionEnergy / thermoelectricProperties.hBar

        return groupVel

    def analyticalGroupVelocity(self, energyRange, nk, m, valley, dk_len, alpha):

        """
        If no DFT calculation is availble this function approximate the group velocity near the conduction band edge.
        This works well up to few hundreds of mev.
        """

        meff = np.array(m)
        ko = 2 * np.pi / self.latticeParameter * np.array(valley)
        del_k = 2 * np.pi / self.latticeParameter * dk_len * np.array([1, 1, 1])
        kx = np.linspace(ko[0], ko[0] + del_k[0], nk[0], endpoint=True)  # kpoints mesh
        ky = np.linspace(ko[1], ko[1] + del_k[1], nk[1], endpoint=True)  # kpoints mesh
        kz = np.linspace(ko[2], ko[2] + del_k[2], nk[2], endpoint=True)  # kpoints mesh
        [xk, yk, zk] = np.meshgrid(kx, ky, kz)
        xk_ = np.reshape(xk, -1)
        yk_ = np.reshape(yk, -1)
        zk_ = np.reshape(zk, -1)

        kpoint = np.array([xk_, yk_, zk_])
        mag_kpoint = norm(kpoint, axis=0)

        mc = 3 / (1 / meff[0] + 1 / meff[1] + 1 / meff[2])  # Conduction band effective mass

        E = thermoelectricProperties.hBar ** 2 / 2 * \
            ((kpoint[0] - ko[0]) ** 2 / meff[0] + (kpoint[1] - ko[1]) ** 2 / meff[1] + (kpoint[2] - ko[2]) ** 2 / meff[
                2]) \
            * thermoelectricProperties.e2C  # Ellipsoidal energy band shape

        vel = thermoelectricProperties.hBar * np.sqrt((kpoint[0] - ko[0]) ** 2 + (kpoint[1] - ko[1]) ** 2
                                                      + (kpoint[2] - ko[2]) ** 2) / mc / (
                          1 + 2 * alpha * E) * thermoelectricProperties.e2C

        Ec, indices, return_indices = np.unique(E, return_index=True, return_inverse=True)  # Smooth data

        vg = accum(return_indices, vel, func=np.mean, dtype=float)

        ESpline = PchipInterpolator(Ec, vg)
        velFunctionEnergy = ESpline(energyRange)

        return velFunctionEnergy  # The array size is [1, numEnergySampling]




if __name__ == "__main__":

    print('ThermoElectric')

