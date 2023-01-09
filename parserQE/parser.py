import re
from munch import Munch
import numpy as np
import sys

def _unitVec(vecA):
    return vecA / np.linalg.norm(vecA)

def angle2Vec(vecA, vecB):
    vAu = _unitVec(vecA)
    vBu = _unitVec(vecB)
    angle = np.arccos(np.dot(vAu, vBu))
    angle *= (180 / np.pi)
    return angle

'''
    This function 'checkandextract' helping to check a list of data and extracting time.
    If the list is none, the function return 0.
'''
def checkandextract(data_list):
    data = 0
    if data_list:
        data = float (data_list[0])
    return data
'''
    This function return the total time in seconds
'''
def myTime(hours, minutes, seconds):
    h = checkandextract(hours)
    m = checkandextract(minutes)
    s = checkandextract(seconds)

    return h * 3600 + m * 60 + s

'''
    This functions create a list for checking the number of iterations and the time that cpu spent in each scf cycle
'''
def getTotalIterations(text):
    total_iter = 0
    for line in text:
        x = re.search(r"iteration #",line)
        if x is not None:
            total_iter += 1
    return total_iter

def getListTimesSCF(text):
    listTime=[]
    for line in text:
        x = re.search(r"total cpu time spent up to now is", line)
        if x is not None:
            y = x.string.split()
            listTime.append(float(y[8]))

    listscf = []
    for i in range(len(listTime)-1):
        listscf.append(listTime[i+1]-listTime[i])

    return listscf

'''
    Functions to get the time in PWSCF
    The PWSCF is the total time. It is reported by CPU time and Wall time.
'''

def get_pwscf_cpu_time_line(line):
    x = re.search(r"PWSCF  ", line)
    if x is not None:
        no_spaces = []
        y = re.split('CPU',x.string)
        y = re.split(' ',y[0])

        for string in y:
            if string != '':
                no_spaces.append(string)
        string = no_spaces[0]
        for i in range(len(no_spaces)-1):
            string = string + no_spaces[i+1]

        ho = re.findall(r"(\d+)h",string)
        mi = re.findall(r"(\d+)m",string)
        se = re.findall(r"([0-9.0-9]*\d+)s",string)

        time = myTime(ho, mi, se)

        return time
    return None

def get_pwscf_cpu_time(text):
    for line in text:
        time = get_pwscf_cpu_time_line(line)
        if time is not None:
            return time

    return None

def get_pwscf_wall_time_line(line):
    pass
'''
    x = re.search(r"PWSCF  ", line)
    if x is not None:
        no_spaces = []
        y = re.split('CPU',x.string)
        y = re.split(' ',y[0])

        for string in y:
            if string != '':
                no_spaces.append(string)
        string = no_spaces[0]
        for i in range(len(no_spaces)-1):
            string = string + no_spaces[i+1]

        ho = re.findall(r"(\d+)h",string)
        mi = re.findall(r"(\d+)m",string)
        se = re.findall(r"([0-9.0-9]*\d+)s",string)

        time = myTime(ho, mi, se)

        return time
    return None
'''

def get_pwscf_wall_time(text):
    for line in text:
        time = get_pwscf_wall_time_line(line)
        if time is not None:
            return time

    return None

'''
    Functions to get the time of "init_run"
    The time is reported by CPU time and Wall time.
'''

def get_init_wall_time_line(line):
    x = re.search(r"init_run   ", line)
    if x is not None:
        no_spaces = []
        y = re.split('CPU',x.string)
        y = re.split(' ',y[1])

        for string in y:
            if string != '':
                no_spaces.append(string)
        string = no_spaces[0]
        for i in range(len(no_spaces)-1):
            string = string + no_spaces[i+1]

        ho = re.findall(r"(\d+)h",string)
        mi = re.findall(r"(\d+)m",string)
        se = re.findall(r"([0-9.0-9]*\d+)s",string)

        time = myTime(ho, mi, se)

        return time
    return None

def get_pwscf_wall_time(text):
    for line in text:
        time = get_init_wall_time_line(line)
        if time is not None:
            return time

    return None


def get_init_cpu_time_line(line):
    x = re.search(r"init_run   ", line)
    if x is not None:
        y = x.string.split()
        se = re.findall(r"([0-9.0-9]*\d+)s",y[2])

        time = myTime([],[],se)

        return time
    return None

def get_init_cpu_time(text):
    for line in text:
        time = get_init_cpu_time_line(line)
        if time is not None:
            return time

    return None
'''
    Functions to get the time of "electrons".
    "electrons" calculate the self-consistent solution.
    The time is reported by CPU time and Wall time.
'''

def get_electrons_cpu_time_line(line):
    x = re.search(r"electrons    :", line)
    if x is not None:
        y = x.string.split()
        se = re.findall(r"([0-9.0-9]*\d+)s",y[2])

        time = myTime([],[],se)

        return time
    return None

def get_electrons_cpu_time(text):
    for line in text:
        time = get_electrons_cpu_time_line(line)
        if time is not None:
            return time
    return None

def get_electrons_wall_time_line(line):
    x = re.search(r"electrons    :", line)
    if x is not None:
        y = x.string.split()
        se = re.findall(r"([0-9.0-9]*\d+)s",y[4])

        time = myTime([],[],se)

        return time
    return None

def get_electrons_wall_time(text):
    for line in text:
        time = get_electrons_wall_time_line(line)
        if time is not None:
            return time
    return None

'''
    Functions to get the time of "c_bands".
    "c_bands" subroutine is call by "electrons", it calculates Kohn-Sham states.
    The time is reported by CPU time and Wall time.
'''

def get_cbands_cpu_time_line(line):
    x = re.search(r"c_bands      :", line)
    if x is not None:
        y = x.string.split()
        se = re.findall(r"([0-9.0-9]*\d+)s",y[2])

        time = myTime([],[],se)

        return time
    return None

def get_cbands_cpu_time(text):
    for line in text:
        time = get_cbands_cpu_time_line(line)
        if time is not None:
            return time
    return None

def get_cbands_wall_time_line(line):
    x = re.search(r"c_bands      :", line)
    if x is not None:
        y = x.string.split()
        se = re.findall(r"([0-9.0-9]*\d+)s",y[4])

        time = myTime([],[],se)

        return time
    return None

def get_cbands_wall_time(text):
    for line in text:
        time = get_cbands_wall_time_line(line)
        if time is not None:
            return time
    return None

'''
    Functions to get the time of "sum_band".
    "sum_band" subroutine is call by "electrons", it calculates charge density.
    The time is reported by CPU time and Wall time.
'''

def get_sumband_cpu_time_line(line):
    x = re.search(r"sum_band     :", line)
    if x is not None:
        y = x.string.split()
        se = re.findall(r"([0-9.0-9]*\d+)s",y[2])

        time = myTime([],[],se)

        return time
    return None

def get_sumband_cpu_time(text):
    for line in text:
        time = get_sumband_cpu_time_line(line)
        if time is not None:
            return time
    return None

def get_sumband_wall_time_line(line):
    x = re.search(r"sum_band     :", line)
    if x is not None:
        y = x.string.split()
        se = re.findall(r"([0-9.0-9]*\d+)s",y[4])

        time = myTime([],[],se)

        return time
    return None

def get_sumband_wall_time(text):
    for line in text:
        time = get_sumband_wall_time_line(line)
        if time is not None:
            return time
    return None

def isJobDone(text):
    done = False
    for line in text:
        x = re.search(r"JOB DONE", line)
        if x is not None:
            done = True
            break

    return done

def getFunctional(text):
    for line in text:
        x = re.search(r"Exchange-correlation= ", line)
        if x is not None:
            y = re.split(r'=',line)
            functional = y[1].strip()
            return functional
    return None

def getNumAtomsPerCell(text):
    for line in text:
        x = re.search(r"number of atoms/cell      =", line)
        if x is not None:
            y = re.split(r'=',line)
            natom = int(y[1].strip())
            return natom
    return None

def getLatticeParam(text):
    for line in text:
        x = re.search(r"lattice parameter \(alat\)  =", line)
        if x is not None:
            y = x.string.split()
            alat = float(y[4])
            return alat
    return None

def getNumElectrons(text):
    for line in text:
        x = re.search(r"number of electrons       =", line)
        if x is not None:
            y = re.split(r'=',line)
            nelectron = float(y[1].strip())
            return nelectron
    return None

def getNumCores(text):
    for line in text:
        x = re.search(r"Parallel version", line)
        if x is not None:
            y = x.string.split()
            if isinstance(int(y[-1]),int):
                return int(y[-1])
            else:
                return None
    return None

def getMPIprocesses(text):
    for line in text:
        x = re.search(r"Number of MPI processes:", line)
        if x is not None:
            y = x.string.split()
            return int(y[4])
    return None

def getThreadsPerMPI(text):
    for line in text:
        x = re.search(r"Threads/MPI process:", line)
        if x is not None:
            y = x.string.split()
            return int(y[2])
    return None

def getNumNodes(text):
    for line in text:
        x = re.search(r"MPI processes distributed on", line)
        if x is not None:
            y = x.string.split()
            return int(y[4])
    return None

def getVersion(text):
    for line in text:
        x = re.search(r"Program PWSCF",line)
        if x is not None:
            y = x.string.split()
            return y[2]
    return None

def getAccuracy(text):
    listAccuracy = []
    for line in text:
        x = re.search(r"estimated scf accuracy", line)
        if x is not None:
            y = x.string.split()
            listAccuracy.append(float(y[4]))

    return listAccuracy;
'''
def getListTimesSCF(text):
    listTime=[]
    for line in text:
        x = re.search(r"total cpu time spent up to now is", line)
        if x is not None:
            y = x.string.split()
            listTime.append(float(y[8]))

    listscf = []
    for i in range(len(listTime)-1):
        listscf.append(listTime[i+1]-listTime[i])

    return listscf
'''

def getNameInput(text):
    for line in text:
        x = re.search(r"Reading input from ", line)
        if x is not None:
            y = x.string.split()
            return y[3]

def getGPUAcceleration(text):
    for line in text:
        x = re.search(r"GPU acceleration",line)
        if x is not None:
            y = x.string.split()
            return y[3].rstrip(y[3][-1])
    return None

def getCutOffEnergy(text):
    for line in text:
        x = re.search(r"kinetic-energy cutoff     =",line)
        if x is not None:
            y = x.string.split()
            return float(y[3])
    return None

def getCutOffDensity(text):
    for line in text:
        x = re.search(r"charge density cutoff     =",line)
        if x is not None:
            y = x.string.split()
            return float(y[4])
    return None

def getGeomChem(text):
    natom = getNumAtomsPerCell(text)
    alat = getLatticeParam(text)
    idx = 0;
    symb = []
    coords = []
    for line in text:
        x = re.search(r'positions \(alat units\)', line)
        idx += 1
        if x is not None:
            for j in range(natom):
                y = text[idx + j]
                y = y.split()
                symb.append(y[1])
                xyz = [float(y[6]) * alat, float(y[7]) * alat, float(y[8]) * alat]
                coords.append(xyz)

            break
    return symb, coords

def getCell(text):
    alat = getLatticeParam(text)
    idx = 0
    cell = []
    for line in text:
        x = re.search('crystal axes:', line)
        idx += 1
        if x is not None:
            for j in range(3):
                y = text[idx + j]
                y = y.split()
                vec = [float (y[3]) * alat, float (y[4]) * alat, float (y[5]) * alat]
                cell.append(vec)
            break

    alpha = angle2Vec(cell[1], cell[2])
    beta = angle2Vec(cell[0], cell[2])
    gamma = angle2Vec(cell[0], cell[1])
    a = np.sqrt(np.dot(cell[0], cell[0]))
    b = np.sqrt(np.dot(cell[1], cell[1]))
    c = np.sqrt(np.dot(cell[2], cell[2]))

    axes = [a, b, c]
    angles = [alpha, beta, gamma]

    return cell, axes, angles

'''
    Functions to get energies:
    Fermi Energy in eV
    Total energy F = E - TS  in Ry
    Internal Energy :  E = F + TS
    Smearing Contribution:  TS
    Energy E is the sum of: One-electron contribution + XC Contribution + Ewald
    Dispersion Correction
'''

def getFermiEnergy(text):
    for line in text:
        x = re.search(r"the Fermi energy is", line)
        if x is not None:
            y = x.string.split()
            return float(y[4])
    return None

def getTotalEnergy(text):
    for line in text:
        x = re.search(r"!    total energy", line)
        if x is not None:
            y = x.string.split()
            return float(y[4])
    return None

def getInternalEnergy(text):
    for line in text:
        x = re.search(r"internal energy E\=F\+TS", line)
        if x is not None:
            y = x.string.split()
            return float(y[4])
    return None

def getSmearingContrib(text):
    for line in text:
        x = re.search(r"smearing contrib. \(-TS\)", line)
        if x is not None:
            y = x.string.split()
            return float(y[4])
    return None

def getOneElectronEnergy(text):
    for line in text:
        x = re.search(r"one-electron contribution =", line)
        if x is not None:
            y = x.string.split()
            return float(y[3])
    return None

def getHartreeEnergy(text):
    for line in text:
        x = re.search(r"hartree contribution      =", line)
        if x is not None:
            y = x.string.split()
            return float(y[3])
    return None

def getXCEnergy(text):
    for line in text:
        x = re.search(r"xc contribution           =", line)
        if x is not None:
            y = x.string.split()
            return float(y[3])
    return None

def getEwaldContrib(text):
    for line in text:
        x = re.search(r"ewald contribution        =", line)
        if x is not None:
            y = x.string.split()
            return float(y[3])
    return None

def getDispersion(text):
    for line in text:
        x = re.search(r"Dispersion Correction     =", line)
        if x is not None:
            y = x.string.split()
            return float(y[3])
    return None

def dumpXYZ(data, fname):
    bohr2angs = 0.529177249
    try:
        with open(fname,"w+") as fout:
            print(" {0}".format(data.chem.natoms), file=fout)
            print(" File generate by ParserQE in Angstroms:{0: .6f} {1: .6f} {2: .6f} {3: .6f} {4: .6f} {5: .6f}".format(
                      data.chem.axes[0] * bohr2angs,
                      data.chem.axes[1] * bohr2angs,
                      data.chem.axes[2] * bohr2angs,
                      data.chem.angles[0],
                      data.chem.angles[1],
                      data.chem.angles[2]
                  ), file=fout)
            for i in range(data.chem.natoms):
                print(" {0:3} {1: 10.6f} {2: 10.6f} {3: 10.6f}".format(data.chem.symb[i],
                                                data.chem.coors[i][0] * bohr2angs,
                                                data.chem.coors[i][1] * bohr2angs,
                                                data.chem.coors[i][2] * bohr2angs),
                                                file=fout)

    except IOError as error:
        print('Error to open file: {0}'.format(fname))

def loadQE(fname):

    qe = Munch()
    qe.init = Munch()
    qe.status = Munch()
    qe.chem = Munch()
    qe.profile = Munch()
    qe.energy = Munch()

    ncores = None
    try:
        with open(fname,"r") as file:
            text = file.readlines()

            try:
                qe.status.done = isJobDone(text)
                if qe.status.done is False:
                    print('The file {0} did not finish correctly!'.format(fname))
            except IOError as error:
                print('The file {0} did not finish correctly!'.format(fname))
                qe = None
            else:
                qe.profile.tot_cpu_time = get_pwscf_cpu_time(text)
                qe.profile.tot_wall_time = get_pwscf_wall_time(text)
                qe.status.iter = getTotalIterations(text)
                qe.status.scftimes = getListTimesSCF(text)
                qe.status.nameinp = getNameInput(text)
                qe.accuracy = getAccuracy(text)
                #qe.status.accuracy = accuracy[-1]
                qe.chem.funct = getFunctional(text)
                qe.chem.natoms = getNumAtomsPerCell(text)
                qe.chem.nelect = getNumElectrons(text)
                qe.chem.symb, qe.chem.coors = getGeomChem(text)
                qe.chem.alat = getLatticeParam(text)
                qe.chem.cell, qe.chem.axes, qe.chem.angles = getCell(text)
                #ncores = getNumCores(text)
                qe.init.nmpi = getMPIprocesses(text)
                qe.init.nthr = getThreadsPerMPI(text)
                qe.init.nnodes = getNumNodes(text)
                qe.energy.fermi = getFermiEnergy(text)
                qe.energy.total = getTotalEnergy(text)
                qe.energy.internal = getInternalEnergy(text)
                qe.energy.smearing = getSmearingContrib(text)
                qe.energy.hcore = getOneElectronEnergy(text)
                qe.energy.hartree = getHartreeEnergy(text)
                qe.energy.xc = getXCEnergy(text)
                qe.energy.ewald = getEwaldContrib(text)
                qe.energy.dispersion = getDispersion(text)
                qe.profile.init_cpu_time = get_init_cpu_time(text)
                qe.profile.electrons_cpu_time = get_electrons_cpu_time(text)
                qe.profile.electrons_wall_time = get_electrons_wall_time(text)
                qe.profile.cbands_cpu_time = get_cbands_cpu_time(text)
                qe.profile.cbands_wall_time = get_cbands_wall_time(text)
                qe.profile.sumband_cpu_time = get_sumband_cpu_time(text)
                qe.profile.sumband_wall_time = get_sumband_wall_time(text)
                qe.chem.ecut = getCutOffEnergy(text)
                qe.chem.dencut = getCutOffDensity(text)
                qe.init.version = getVersion(text)
                qe.status.gpuacc = getGPUAcceleration(text)
            #print(qe)
            #print(qe.keys())
    except IOError as error:
        print('Error to open file: {0}'.format(fname))
    else:
       return qe


def printQE_info(qe):
    if qe is None:
        print('Error no information')
        return
            
    if qe.status.done is True:
        print(" STATUS")
        print(" The job is done : ", qe.status.done)
        print(" Iterations      : ", qe.status.iter)
        print(" Size in list scf: ", len(qe.status.scftimes))
        print(" Times in scf    : ", qe.status.scftimes)
        print(" Name's inp file : ", qe.status.nameinp)
        print(" Acceleration GPU: ", qe.status.gpuacc)
        print("")

        print(" INITIALIZATION ")
        print(" PWSCF Version   : ", qe.init.version)
        print(" MPI processes   : ", qe.init.nmpi)
        print(" Threads / MPI   : ", qe.init.nthr)
        print(" Num nodes       : ", qe.init.nnodes)
        #print(" Num cores      : ", ncores)
        print("")

        print(" CHEMISTRY ")
        print(" XC Functional   : ", qe.chem.funct)
        print(" Num of Atoms    : ", qe.chem.natoms)
        print(" Num electrons   : ", qe.chem.nelect)
        print(" Cutoff Energy   : ", qe.chem.ecut)
        print(" Cutoff Density  : ", qe.chem.dencut)
        print(" List of Symbols : ", qe.chem.symb)
        print(" List of Coords  : ", qe.chem.coors)
        print(" Lattice param   : ", qe.chem.alat)
        print(" Cell param      : ", qe.chem.cell)
        print(" Cell axes       : ", qe.chem.axes)
        print(" Cell angles     : ", qe.chem.angles)
        print("")

        print(" ENERGY PROFILE ")
        print(" Fermi energy    : ", qe.energy.fermi)
        print(" Total energy    : ", qe.energy.total)
        print(" Internal energy : ", qe.energy.internal)
        print(" Smearing contrib: ", qe.energy.smearing)
        print(" One-electron    : ", qe.energy.hcore)
        print(" Hartree contrib : ", qe.energy.hartree)
        print(" XC contrib      : ", qe.energy.xc)
        print(" Ewald contrib   : ", qe.energy.ewald)
        print(" Dispersion Corre: ", qe.energy.dispersion)
        print("")

        print(" PROFILE TIMING  :")
        print("   CPU     time  : ", qe.profile.tot_cpu_time)
        print("   WALL    time  : ", qe.profile.tot_wall_time)
        print("   Init    time  : ", qe.profile.init_cpu_time)
        print(" Elecs CPU  time : ", qe.profile.electrons_cpu_time)
        print(" Elecs Wall time : ", qe.profile.electrons_wall_time)
        print(" cbands CPU time : ", qe.profile.cbands_cpu_time)
        print(" cbands Wall time: ", qe.profile.cbands_wall_time)
        print(" sum ba CPU time : ", qe.profile.sumband_cpu_time)
        print(" sum ba Wall time: ", qe.profile.sumband_wall_time)

