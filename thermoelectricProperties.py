import numpy as np
import os
import numpy.matlib
import math
from scipy.interpolate import InterpolatedUnivariateSpline
from functools import wraps
import time
import cmath
import matplotlib
import matplotlib.pyplot as plt
from scipy import interpolate, optimize
from scipy.interpolate import UnivariateSpline, InterpolatedUnivariateSpline, interp1d
from scipy.optimize import anderson, broyden1, broyden2, excitingmixing, diagbroyden, newton_krylov, linearmixing
import yaml


class thermoelectricProperties:

    hBar = 6.582119e-16     # Reduced Planck constant in eV.s/rad
    kB = 8.617330350e-5     # Boltzmann constant in eV/K
    e2C = 1.6021765e-19     # e to Coulomb unit change
    e0 = 8.854187817e-12    # Permittivity in vacuum F/m
    Ang2meter = 1e-10       # Unit conversion from Angestrom to meter

    def __init__(self, latticeParameter, dopantElectricCharge, electronEffectiveMass, energyMin, energyMax, kS, numKpoints, numBands, numQpoints, electronDispersian=None, kpoints=None, numEnergySampling=1e3):

        self.latticeParameter = latticeParameter            # Lattice parameter in A
        self.dopantElectricCharge = dopantElectricCharge
        self.electronEffectiveMass = electronEffectiveMass
        self.energyMax = energyMax                          # Maximum energy in eV
        self.energyMin = energyMin                          # Minimum energy in eV
        self.kS = kS                                        # Relative permittivity
        self.numEnergySampling = numEnergySampling          # Number of energy space samples to generate in eV, the defaulf value is 1000
        self.numKpoints = numKpoints
        self.numBands = numBands
        self.electronDispersian = electronDispersian
        self.numQpoints = numQpoints

    # def loggerFanction(orig_func):
    #     import logging
    #     logging.basicconf(filename='{}.log'.format(orig_func.__name__), level=logging.info)

    #     @wraps(orig_func)
    #     def wrapper(*args, **kwargs):
    #         tStart = time.time()
    #         result = orig_func(*args, **kwargs)
    #         elapsedTime = time.time() - tStart
    #         logging.info(
    #             '{} ran in {} with args: {}, and kwargs: {} '.format(org_func.__name__, elapsedTime, args, kwargs))
    #         return result
    #     return wrapper

    def energyRange(self):                                  # Create an array of energy space sampling
        energyRange = np.linspace(self.energyMin, self.energyMax, self.numEnergySampling)
        return np.expand_dims(energyRange, axis=0)

    def kpoints(self, path2kpoints, delimiter=None, skiprows=0):
        kpoints = np.loadtxt(expanduser(path2kpoints), delimiter=None, skiprows=0)
        return kpoints

    def temp(self, TempMin=300, TempMax=1301, dT=100):
        temperature = np.arange(TempMin, TempMax, dT)
        return np.expand_dims(temperature, axis=0)

    def bandGap(self, Eg_o, Ao, Bo, Temp=None):
        if Temp is None:
            T = self.temp()
        else:
            T = Temp
        Eg = Eg_o - Ao * np.divide(T**2, T + Bo)
        return Eg

    def analyticalDoS(self, energyRange, alpha, unitcell_volume):
        print(np.shape(alpha))
        DoS = 1 / np.pi**2 * np.sqrt(2 * energyRange * (1 + energyRange / np.transpose(alpha))) * np.sqrt(self.electronEffectiveMass / thermoelectricProperties.hBar**2)**3 * (1 + (2 * energyRange / np.transpose(alpha))) / unitcell_volume / thermoelectricProperties.e2C**(3. / 2)
        return DoS

    def carrierConcentration(self, path2extrinsicCarrierConcentration, bandGap, Ao=None, Bo=None, Nc=None, Nv=None, Temp=None):
        if Temp is None:
            T = self.temp()
        else:
            T = Temp

        if Ao is None and Nc is None:
            raise Exception("Either Ao or Nc should be defined")
        if Bo is None and Nv is None:
            raise Exception("Either Bo or Nv should be defined")

        if Nc is None:
            Nc = Ao * Temp**(3. / 2)
        if Nv is None:
            Nv = Bo * Temp**(3. / 2)

        exCarrierFile = np.loadtxt(expanduser(path2extrinsicCarrierConcentration), delimiter=None, skiprows=0)
        extrinsicCarrierConcentration_tmp = InterpolatedUnivariateSpline(exCarrierFile[0, :], exCarrierFile[1, :] * 1e6)
        extrinsicCarrierConcentration = extrinsicCarrierConcentration_tmp(T)
        intrinsicCarrierConcentration = np.multiply(np.sqrt(np.multiply(Nc, Nv)), np.exp(-(np.divide(bandGap, (2 * thermoelectricProperties.kB * T)))))
        totalCarrierConcentration = intrinsicCarrierConcentration + abs(extrinsicCarrierConcentration)
        # return totalCarrierConcentration
        return totalCarrierConcentration

    def fermiLevel(self, carrierConcentration, Nc=None, Ao=None, Temp=None):

        if Temp is None:
            T = self.temp()
        else:
            T = Temp

        if Ao is None and Nc is None:
            raise Exception("Either Ao or Nc should be defined")

        if Nc is None:
            Nc = Ao * Temp**(3. / 2)

        JD_CC = np.log(np.divide(carrierConcentration, Nc)) + 1 / np.sqrt(8) * np.divide(carrierConcentration, Nc) - (3. / 16 - np.sqrt(3) / 9) * np.power(np.divide(carrierConcentration, Nc), 2)
        fermiLevelEnergy = thermoelectricProperties.kB * np.multiply(T, JD_CC)
        return fermiLevelEnergy

    def fermiDistribution(self, energyRange, fermiLevel, Temp=None):

        if Temp is None:
            T = self.temp()
        else:
            T = Temp

        fermiDirac = np.divide(1, np.exp(np.divide((repmat(energyRange, np.shape(T)[1], 1) - repmat(np.transpose(fermiLevel), 1, np.shape(energyRange)[1])), thermoelectricProperties.kB * repmat(np.transpose(T), 1, np.shape(energyRange)[1]))) + 1)  # Fermi Dirac distribution
        dE = energyRange[0, 2] - energyRange[0, 1]
        df = np.roll(fermiDirac, -1, axis=1) - np.roll(fermiDirac, 1, axis=1)
        dfdE = df / dE / 2
        dfdE[:, 0] = (fermiDirac[:, 1] - fermiDirac[:, 0]) / dE
        dfdE[:, -1] = (fermiDirac[:, -1] - fermiDirac[:, -2]) / dE
        fermi = np.array([fermiDirac, dfdE])
        return fermi

    def electronBandStructure(self, path2eigenval, skipLines):

        with open(expanduser(path2eigenval)) as eigenvalFile:
            for _ in range(skipLines):
                next(eigenvalFile)
            block = [[float(_) for _ in line.split()] for line in eigenvalFile]
        eigenvalFile.close()
        electronDispersian = [range(1, self.numBands + 1)]  # First line is atoms id
        kpoints = np.asarray(block[1::self.numBands + 2])[:, 0:3]
        # print(type(kpoints))
        for _ in range(self.numKpoints):
            binary2Darray = []
            for __ in range(self.numBands):
                binary2Darray = np.append(binary2Darray, block[__ + 2 + (self.numBands + 2) * _][1])
            electronDispersian = np.vstack([electronDispersian, binary2Darray])  # Next lines are eigenvalues in eV
        dispersian = [kpoints, electronDispersian]
        return dispersian

    def electronDoS(self, path2DoS, headerLines, numDoSpoints, unitcell_volume, valleyPoint, energyRange):

        DoS = np.loadtxt(expanduser(path2DoS), delimiter=None, skiprows=headerLines, max_rows=numDoSpoints)
        # return DoS
        valleyPointEnergy = DoS[valleyPoint, 0]  # Energy of vally in conduction band
        # print(valleyPointEnergy)
        DoSSpline = InterpolatedUnivariateSpline(DoS[valleyPoint:, 0] - valleyPointEnergy, DoS[valleyPoint:, 1] / unitcell_volume)
        DoSFunctionEnergy = DoSSpline(energyRange)  # Density of state
        return DoSFunctionEnergy

    def fermiLevelSelfConsistent(self, carrierConcentration, Temp, energyRange, DoS):

        fermi = np.linspace(-1.5, 1.5, 1000, endpoint=True)
        result_array = np.empty((np.shape(Temp)[1], np.shape(fermi)[0]))
        idx_j = 0
        for j in Temp[0]:
            idx_i = 0
            for i in fermi:
                f, _ = self.fermiDistribution(energyRange=energyRange, fermiLevel=np.expand_dims(np.array([i]), axis=0), Temp=np.expand_dims(np.array([j]), axis=0))
                tmp = np.trapz(np.multiply(DoS, f), energyRange, axis=1)
                result_array[idx_j, idx_i] = tmp
                idx_i += 1
            idx_j += 1

        diff = np.tile(np.transpose(carrierConcentration), (1, np.shape(fermi)[0])) - abs(result_array)
        # print('diff', diff)
        # print('re', result_array)
        # print(np.argmin(np.abs(diff), axis=1))
        min_idx = np.argmin(np.abs(diff), axis=1)
        # print(min_idx)
        Ef = fermi[min_idx]
        elm = 0
        n = np.empty((1, np.shape(Temp)[1]))
        # print(n, np.shape(n))
        # print([result_array[0, 0]], np.shape([result_array[0, 0]]))
        for idx in min_idx:
            n = np.append(n, [np.array(np.expand_dims(result_array[elm, idx], axis=0))], axis=1)
            elm += 1
        return np.expand_dims(Ef, axis=0)

    def electronGroupVelocity(self, kp, energy_kp, energyRange):

        dE = np.roll(energy_kp, -1, axis=0) - np.roll(energy_kp, 1, axis=0)
        dk = np.roll(kp, -1, axis=0) - np.roll(kp, 1, axis=0)
        dEdk = np.divide(dE, dk)
        dEdk[0] = (energy_kp[1] - energy_kp[0]) / (kp[1] - kp[0])
        dEdk[-1] = (energy_kp[-1] - energy_kp[-2]) / (kp[-1] - kp[-2])
        dEdkSpline = InterpolatedUnivariateSpline(energy_kp, np.array(dEdk))
        dEdkFunctionEnergy = dEdkSpline(energyRange)
        groupVel = dEdkFunctionEnergy / thermoelectricProperties.hBar
        # groupVel = dEdk / thermoelectricProperties.hBar
        return groupVel

    def matthiessen(self, energyRange, *args):
        tau = 1. / sum([1. / arg for arg in args])
        tau[np.isinf(tau)] = 0
        return tau

    def electricalProperties(self, E, DoS, vg, Ef, dfdE, Temp, tau):

        X = DoS * vg**2 * dfdE
        print(np.shape(DoS), np.shape(vg), np.shape(dfdE), np.shape(X * tau))
        Y = (E - np.transpose(Ef)) * X
        Sigma = -1 * np.trapz(X * tau, E, axis=1) / 3 * thermoelectricProperties.e2C
        print(Sigma)
        # X = np.multiply(numpy.matlib.repmat(np.multiply(np.power(electronGroupVelocity[0], 2), electronDoS), len(temp), 1), fermiDistribution[1])
        # Y = np.multiply(np.multiply(numpy.matlib.repmat(np.multiply(np.power(electronGroupVelocity[0], 2), electronDoS), len(temp), 1), (numpy.matlib.repmat(Erange, len(temp), 1) - numpy.matlib.repmat((np.array([fermiLevel]).T), 1, len(Erange)))), fermiDistribution[1])
        # conductivity = factor * (-1.0 / 3) * ((q**2) * thermoelectricProperties.e2C) * np.trapz(np.multiply(X, numpy.matlib.repmat(matthiessen, len(temp), 1)), energyRange, dx=dE, axis=1)
        # resistivity = 1. / conductivity
        # seebeck = factor * q * (-1.0 / 3) * thermoelectricProperties.e2C * np.divide(np.divide(np.trapz(np.multiply(Y, numpy.matlib.repmat(matthiessen, len(temp), 1)), energyRange, dx=dE, axis=1), conductivity), temp)
        # powerFactor = np.multiply(conductivity, np.power(seebeck, 2))
        return np.expand_dims(Sigma, axis=0)

    # def electricalProperties(self, energyRange, electronBandStructure, temp, electronDoS, electronGroupVelocity, matthiessen, bandGap, carrierConcentration, fermiLevel, fermiDistribution, factor, q):
    #     X = np.multiply(numpy.matlib.repmat(np.multiply(np.power(electronGroupVelocity[0], 2), electronDoS), len(temp), 1), fermiDistribution[1])
    #     Y = np.multiply(np.multiply(numpy.matlib.repmat(np.multiply(np.power(electronGroupVelocity[0], 2), electronDoS), len(temp), 1), (numpy.matlib.repmat(Erange, len(temp), 1) - numpy.matlib.repmat((np.array([fermiLevel]).T), 1, len(Erange)))), fermiDistribution[1])
    #     conductivity = factor * (-1.0 / 3) * ((q**2) * thermoelectricProperties.e2C) * np.trapz(np.multiply(X, numpy.matlib.repmat(matthiessen, len(temp), 1)), energyRange, dx=dE, axis=1)
    #     resistivity = 1. / conductivity
    #     seebeck = factor * q * (-1.0 / 3) * thermoelectricProperties.e2C * np.divide(np.divide(np.trapz(np.multiply(Y, numpy.matlib.repmat(matthiessen, len(temp), 1)), energyRange, dx=dE, axis=1), conductivity), temp)
    #     powerFactor = np.multiply(conductivity, np.power(seebeck, 2))
    #     return seebeck, resistivity, conductivity, powerFactor, X, Y

    def filteringEffect(self, U0, tau0, tauOff, energyRange, electronBandStructure, temp, electronDoS, electronGroupVelocity, bandGap, carrierConcentration, fermiLevel, fermiDistribution, factor, q, uIncrement=0.05, tauIncrement=1e-15, tempIndex=0):
        n = 0
        m = 0
        sMatrix = np.array([])
        # rMatrix = np.array([])
        # pfMatrix = np.array([])

        for _ in np.arange(0.1, U0, uIncrement):
            m += 1
            for __ in np.arange(tauIncrement, tau0, tauIncrement):
                tauInc = np.ones(len(Erange))
                tauInc[np.where(Erange < _)] = __
                tau_on = self.matthiessen(energyRange, tauOff, tauInc)
                s, r, c, pf, X, Y = self.electricalProperties(energyRange, electronBandStructure, temp, electronDoS, electronGroupVelocity, tau_on, bandGap, carrierConcentration, fermiLevel, fermiDistribution, factor, q)
                sMatrix = np.append(sMatrix, s[tempIndex])
                n += 1
        sMatrix = np.reshape(sMatrix, (n, m))
        return sMatrix

    # def qpoints(self):
    #     qpoints = np.array([np.zeros(self.numQpoints), np.zeros(self.numQpoints), np.linspace(-math.pi / self.latticeParameter, math.pi / self.latticeParameter, num=self.numQpoints)])
    #     return qpoints

    def dynamicalMatrix(self, path2massWeightedHessian, path2atomsPositions, skipLines, numAtoms, baseLatticePoint, numAtomsInUnitCell, qpoints):
        with open(os.path.expanduser(path2massWeightedHessian)) as hessianFile:
            hessianMatrix = hessianFile.readlines()
        hessianMatrix = [line.split() for line in hessianMatrix]
        hessianMatrix = np.array([[float(_) for _ in __] for __ in hessianMatrix])
        hessianFile.close()
        hessianSymmetry = (np.triu(hessianMatrix) + np.tril(hessianMatrix).transpose()) / 2
        hessianMatrix = hessianSymmetry + np.triu(hessianSymmetry, 1).transpose()
        #
        with open(os.path.expanduser(path2atomsPositions)) as atomsPositionsFile:
            atomsPositions = atomsPositionsFile.readlines()
        atomsPositions = [line.split() for line in atomsPositions]
        [atomsPositions.pop(0) for _ in range(skipLines)]
        atomsPositions = np.array([[float(_) for _ in __] for __ in atomsPositions[0:numAtoms]])
        atomsPositions = atomsPositions[atomsPositions[:, 0].argsort()]
        # atomsPositions = np.sort(atomsPositions.view('i8,i8,f8,f8,f8,f8,f8,f8'), order=['f0'], axis=0).view(np.float)
        latticePoints = np.array([_[2:5] for _ in atomsPositions[::numAtomsInUnitCell]])
        latticePointsVectors = latticePoints - numpy.matlib.repmat(latticePoints[baseLatticePoint], len(latticePoints), 1)
        dynamicalMatrix = np.zeros((numAtomsInUnitCell * 3, numAtomsInUnitCell * 3))
        for _ in range(self.numQpoints):
            dynamMatPerQpoint = np.zeros((numAtomsInUnitCell * 3, numAtomsInUnitCell * 3))
            for __ in range(len(latticePointsVectors)):
                sumMatrix = hessianMatrix[__ * numAtomsInUnitCell * 3: (__ + 1) * numAtomsInUnitCell * 3, baseLatticePoint * numAtomsInUnitCell * 3: (baseLatticePoint + 1) * numAtomsInUnitCell * 3] * cmath.exp(-1j * np.dot(latticePointsVectors[__], qpoints[:, _]))
                dynamMatPerQpoint = dynamMatPerQpoint + sumMatrix
            dynamicalMatrix = np.append(dynamicalMatrix, dynamMatPerQpoint, axis=0)
        dynamicalMatrix = dynamicalMatrix[numAtomsInUnitCell * 3:]
        eigVal = np.array([])
        eigVec = np.zeros((numAtomsInUnitCell * 3, numAtomsInUnitCell * 3))
        for _ in range(self.numQpoints):
            dynmat = dynamicalMatrix[_ * numAtomsInUnitCell * 3:(_ + 1) * numAtomsInUnitCell * 3]
            eigvals, eigvecs, = np.linalg.eigh(dynmat)
            eigVal = np.append(eigVal, eigvals).reshape(-1, numAtomsInUnitCell * 3)
            eigVec = np.append(eigVec, eigvecs, axis=0)
        eigVec = eigVec[numAtomsInUnitCell * 3:]
        frequencies = np.sqrt(np.abs(eigVal.real)) * np.sign(eigVal.real)
        # conversion_factor_to_THz = 15.633302
        # frequencies = frequencies * conversion_factor_to_THz
        return eigVec

    def phonopyQpointYamlInterface(self, path2QpointYaml):
        qpointsData = yaml.load(open("qpoints.yaml"))
        nqpoint = qpointsData['nqpoint']
        natom = qpointsData['natom']
        qpoints = []
        qpoints = np.append(qpoints, [qpointsData['phonon'][_]['q-position'] for _ in range(nqpoint)]).reshape(-1, 3)
        frequency = []
        frequency = np.append(frequency, [[qpointsData['phonon'][_]['band'][__]['frequency'] for __ in range(3 * natom)] for _ in range(nqpoint)]). reshape(-1, 3 * natom)
        eigVal = np.array([])
        eigVec = np.zeros((natom * 3, natom * 3))
        for _ in range(nqpoint):
            dynmat = []
            dynmat_data = qpointsData['phonon'][_]['dynamical_matrix']
            for row in dynmat_data:
                vals = np.reshape(row, (-1, 2))
                dynmat.append(vals[:, 0] + vals[:, 1] * 1j)
            dynmat = np.array(dynmat)
            eigvals, eigvecs, = np.linalg.eigh(dynmat)
            eigVal = np.append(eigVal, eigvals).reshape(-1, natom * 3)
            eigVec = np.append(eigVec, eigvecs, axis=0)
        eigVec = eigVec[natom * 3:]
        frequencies = np.sqrt(np.abs(eigVal.real)) * np.sign(eigVal.real)
        conversion_factor_to_THz = 15.633302
        frequencies = frequencies * conversion_factor_to_THz
        return eigVec

    def gaussianDestribution(self, sigma, expectedValue, qpoints):
        gauss = (1.0 / np.sqrt(2 * pi) / sigma) * np.exp((-1.0 / 2) * np.power(((qpoints - expectedValue) / sigma), 2))
        return gauss

    # def singleWave(self, path2atomsPosition, numberOfAtomsInChain, skipLines)
    # def __str(self):

    # def __repr(self):


silicon = thermoelectricProperties(5.42, 1, 0.2, 0.001, 2, 4.1, 800, 8, 201)
# Erange = thermoelectricProperties.energyRange(silicon)
# Disper = thermoelectricProperties.electronBandStructure(silicon, '~/Desktop/Si-Vasp-5_1_17/Si-Band-Structure-05_01_17/EIGENVAL', 6)
# print('Good so far')
# kpoints = thermoelectricProperties.kpoints(silicon, '~/Desktop/Si-Vasp-5_1_17/KvE-in-CB-valley.txt')
# print(kpoints[1])
# plt.plot(kpoints[1], Disper[401:601, 3], kpoints[1], Disper[401:601, 4])
# plt.show()
# T = thermoelectricProperties.temp(silicon, 400, 1201)
# print(len(T))
# bandgap = thermoelectricProperties.bandGap(silicon, T, 1.170 - np.divide(4.73e-4 * np.power(T, 2), (T + 636)))
# carrierCon = thermoelectricProperties.carrierConcentration(silicon, 2.81e19 * np.power((T / 300.00), 3.0 / 2), 1.831e19 * np.power((T / 300.00), 3.0 / 2), '~/Desktop/data', bandgap, T)
# carrierCon_5 = thermoelectricProperties.carrierConcentration(silicon, 2.81e19 * np.power((T / 300.00), 3.0 / 2), 1.831e19 * np.power((T / 300.00), 3.0 / 2), '~/Desktop/data_5', bandgap, T)

# plt.figure(1)
# plt.plot(T, carrierCon, T, carrierCon_5)
# plt.show()

# fermilevelenergy = thermoelectricProperties.fermiLevel(silicon, 2.81e19 * np.power((T / 300.00), 3.0 / 2), T, carrierCon)
# fermilevelenergy_5 = thermoelectricProperties.fermiLevel(silicon, 2.81e19 * np.power((T / 300.00), 3.0 / 2), T, carrierCon_5)

# print(fermilevelenergy_5)
# plt.figure(2)
# plt.plot(T, fermilevelenergy, T, fermilevelenergy_5)
# plt.show()
# fermiDist = thermoelectricProperties.fermiDistribution(silicon, Erange, T, fermilevelenergy)
# fermiDist_5 = thermoelectricProperties.fermiDistribution(silicon, Erange, T, fermilevelenergy_5)

# plt.plot(Erange, fermiDist[1][1], Erange, fermiDist[1][50], Erange, fermiDist_5[1][1], Erange, fermiDist_5[1][50])
# plt.show()
# plt.plot(Erange, fermiDist[1][1], Erange, fermiDist[1][50], Erange, fermiDist_5[1][1], Erange, fermiDist_5[1][50])
# plt.show()
# print(len(fermiWindow[1]))
# DoS = thermoelectricProperties.electronDoS(silicon, '~/Desktop/Si-Vasp-5_1_17/Si-DoS-05_04_17/DOSCAR', 6, 2000, 1118, Erange)
# plt.plot(Erange, DoS)
# plt.show()
# Vel = thermoelectricProperties.electronGroupVelocity(silicon, Disper, 4, kpoints, 401, 601, Erange)
# print len(Vel[0]), type(Vel[0])
# plt.plot(kpoints[1], Vel[1])
# plt.show()
# plt.plot(Erange, abs(Vel[0]))
# plt.show()
# # test = np.trapz(np.divide(DoS, (1 + np.exp((-0.1) / (thermoelectricProperties.kB * 400)))), Erange, dx=dE, axis=0)
# # print(test)
# # Ef = thermoelectricProperties.fermiLevelSelfConsistent(silicon, carrierCon, DoS, fermiDist, T, Erange)
# # print(Ef)
# tauPhonon = 2e-15 * Erange**-0.5
# tauInc = np.ones(len(Erange))
# tauInc[np.where(Erange < 0.2)] = 1e-14
# print tauInc
# tauOff = thermoelectricProperties.matthiessen(silicon, Erange, tauPhonon)
# tauOn = thermoelectricProperties.matthiessen(silicon, Erange, tauPhonon, tauInc)
# plt.plot(Erange, tauOff, Erange, tauOn)
# plt.show()
# s, r, c, pf, X, Y = thermoelectricProperties.electricalProperties(silicon, Erange, Disper, T, DoS, Vel, tauOff, bandgap, carrierCon, fermilevelenergy, fermiDist, 3, -1)
# s_5, r_5, c_5, pf_5, X_5, Y_5 = thermoelectricProperties.electricalProperties(silicon, Erange, Disper, T, DoS, Vel, tauOn, bandgap, carrierCon_5, fermilevelenergy_5, fermiDist_5, 3, -1)
# ss, rr, cc, pff, XX, YY = thermoelectricProperties.electricalProperties(silicon, Erange, Disper, T, DoS, Vel, tauOff, bandgap, carrierCon_5, fermilevelenergy_5, fermiDist_5, 3, -1)

# T_exp = np.array([374.1, 472.9, 575.6, 675.7, 775.6, 874, 973.7, 1073.6, 1173.6, 1273.5])
# R_exp = np.array([7.086e-06, 7.6208e-06, 8.295e-06, 8.9915e-06, 1.0114e-05, 1.1529e-05, 1.3554e-05, 1.3404e-05, 1.0918e-05, 1.1011e-05])
# S_exp = [-0.00008005, -0.00009907, -0.00011806, -0.00013041, -0.000147142, -0.00015356, -0.00016366, -0.00015978, -0.000149326, -0.00014794]  # Seebeck 5% inclusion
# PF_exp = np.divide(np.power(S_exp, 2), R_exp)
# T5_exp = np.array([373, 473.2, 573.3, 673.5, 773.3, 873.3, 973.3, 1073.3, 1173.3, 1273.4])
# R5_exp = np.array([1.2949e-05, 1.3688e-05, 1.458e-05, 1.5659e-05, 1.6999e-05, 1.8321e-05, 2.049e-05, 2.0107e-05, 1.761e-05, 1.6853e-05])
# S5_exp = np.array([-0.00010974, -0.00014027, -0.00016151, -0.00017394, -0.000186419, -0.00019944, -0.00021742, -0.000218135, -0.00021185, -0.0002107])  # Seebeck 5% inclusion
# PF5_exp = np.divide(np.power(S5_exp, 2), R5_exp)
# print s[1]
# print len(T), len(Erange), len(x), len(x[1])
# print T[50]
# plt.plot(Erange, X[0], Erange, X[50])
# plt.show()
# plt.plot(Erange, Y[0], Erange, Y[50])
# plt.show()
# plt.plot(Erange, X[0], Erange, X[50], Erange, X_5[0], Erange, X_5[50])
# plt.show()
# plt.plot(Erange, Y[0], Erange, Y[50], Erange, Y_5[0], Erange, Y_5[50])
# plt.show()
# plt.figure(1)
# plt.subplot(311)
# plt.plot(T, s, T, s_5, T, ss, T_exp, S_exp, T5_exp, S5_exp)
# plt.ylabel("Seebeck")
# plt.xlabel('Temperature ($^o$K)')
# plt.subplot(312)
# plt.plot(T, r, T, r_5, T, rr, T_exp, R_exp, T5_exp, R5_exp)
# plt.ylabel("Resistivity")
# plt.xlabel('Temperature ($^o$K)')
# plt.subplot(313)
# plt.plot(T, pf, T, pf_5, T, pff, T_exp, PF_exp, T5_exp, PF5_exp)
# plt.ylabel("Power factor")
# plt.xlabel('Temperature ($^o$K)')
# plt.legend(['No inclusion', 'With 5% of inclusion', 'cc'])
# plt.show()

# matrixPro = thermoelectricProperties.filteringEffect(silicon, 0.2, 5e-15, tauOff, Erange, Disper, T, DoS, Vel, bandgap, carrierCon, fermilevelenergy, fermiDist, 3, -1, uIncrement=0.05, tauIncrement=1e-15)
# print matrixPro
# Qpoint = np.array([np.zeros(silicon.numQpoints), np.zeros(silicon.numQpoints), np.linspace(-math.pi / silicon.latticeParameter, math.pi / silicon.latticeParameter, num=silicon.numQpoints)])
# print len(Qpoint[1])
# dynamicalMatrix = thermoelectricProperties.dynamicalMatrix(silicon, '~/Desktop/Notes/Box_120a_Lambda_10a/Si-hessian-mass-weighted-hessian.d', '~/Desktop/Notes/Box_120a_Lambda_10a/data.Si-3x3x3', 15, 216, 14, 8, Qpoint)
# print dynamicalMatrix
# print(dynamicalMatrix)
# Eig = thermoelectricProperties.phonopyQpointYamlInterface(silicon, '~/Desktop/qpoints.yaml')
# np.savetxt('EigVec', Eig.real, fmt='%10.5f', delimiter=' ')
print('done')
