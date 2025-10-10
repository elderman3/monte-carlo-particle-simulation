import numpy as np
import matplotlib.pyplot as plt

Na = 6.02214076e23
barn = 1e-24

class Material():
    def __init__(self):
        self.sym = ""
        self.z = 0
        self.a = 0
        self.aw = 0
        self.T = 0
        self.neutrons = []
        self.mt = {}

def read_file(filename):
    mat = Material()
    with open(filename, "r") as f:
        l = f.readline().rstrip().split(" ")
        mat.sym = l[0]
        mat.z = l[1]
        mat.a = l[2]
        mat.aw = float(l[3])
        mat.T = l[4]
        neu_num = int(f.readline())
        for i in range(neu_num):
            l = f.readline().split(" ")
            mat.neutrons.append((float(l[0]), float(l[1])))
        while l := f.readline():
            l = l.rstrip().split(" ")
            mat.mt[l[0]] = []
            for j in range(int(l[2])):
                ln = f.readline().split(" ")
                mat.mt[l[0]].append((float(ln[0]), float(ln[1])))
    return mat

def value_interp(data, target):
    for (E1, s1), (E2, s2) in zip(data, data[1:]):
        if E1 <= target <= E2:
            return (s2 - s1) / (E2 - E1) * (target - E1) + s1
    return 0 # if a value is not found in this increment

def sigma(sigs, weights, molar_masses, rho):
    sigs = [np.asarray(s, float) for s in sigs]
    n = len(sigs)
    total = np.zeros_like(sigs[0])
    for i in range(n):
        Ni = rho * weights[i] / molar_masses[i] * Na
        total += Ni * sigs[i] * barn
    return total

def drawH1O16():
    matH1 = read_file("data/H1.dat")
    matO16 = read_file("data/O16.dat")
    x = np.logspace(-11, np.log10(20.0), 500)
    mt = matH1.mt.keys()
    yH1 = []
    for i in x:
        sum = 0
        for m in mt:
            sum += value_interp(matH1.mt[m], i)
        yH1.append(sum)
    mt = matO16.mt.keys()
    yO16 = []
    for i in x:
        sum = 0
        for m in mt:
            sum += value_interp(matO16.mt[m], i)
        yO16.append(sum)
    plt.figure()
    plt.title("Cross section of H1 and O16")
    plt.xlabel("Energy (MeV)")
    plt.ylabel("Cross section (barn)")
    plt.loglog(x, yH1, label="H1")
    plt.loglog(x, yO16, label="O16")
    plt.xlim(10**(-11), 20)
    plt.legend()
    plt.grid(True)
    plt.show()

#drawH1O16()

def drawU235U238():
    matU235 = read_file("data/U235.dat")
    matU238 = read_file("data/U238.dat")
    x = np.logspace(-11, np.log10(20.0), 500)
    yU235R = []
    yU238R = []
    yU235F = []
    yU238F = []  
    for i in x:
        yU235R.append(value_interp(matU235.mt["102"], i))
        yU238R.append(value_interp(matU238.mt["102"], i))
        yU235F.append(value_interp(matU235.mt["18"], i))
        yU238F.append(value_interp(matU238.mt["18"], i))

    plt.figure()
    plt.title("Radiative and Total cross section of U238 and U235")
    plt.xlabel("Energy (MeV)")
    plt.ylabel("Cross section (barn)")
    plt.loglog(x, yU235F, label="Fission U235")
    plt.loglog(x, yU235R, label="Radiative U235")
    plt.loglog(x, yU238F, label="Fission U238")
    plt.loglog(x, yU238R, label="Radiative U238")
    plt.xlim(10**(-11), 20)
    plt.legend()
    plt.grid(True)
    plt.show()

#drawU235U238()

def drawU238Inelastic():
    matU238 = read_file("data/U238.dat")
    x = np.logspace(-11, np.log10(20.0), 500)
    yU238 = []
    for i in x:
        sum = 0
        for mt in range(51, 77):
            sum += value_interp(matU238.mt[str(mt)], i)
        yU238.append(sum)
    plt.figure()
    plt.title("Inelastic Scattering Cross section of U238")
    plt.xlabel("Energy (MeV)")
    plt.ylabel("Cross section (barn)")
    plt.loglog(x, yU238, label="U238 Inelastic")
    plt.xlim(10**(-11), 20)
    plt.legend()
    plt.grid(True)
    plt.show()

#drawU238Inelastic()

def drawH2ONatU():
    matH1 = read_file("data/H1.dat")
    matO16 = read_file("data/O16.dat")
    matU235 = read_file("data/U235.dat")
    matU238 = read_file("data/U238.dat")
    x = np.logspace(-11, np.log10(20.0), 500)

    MH2O = 2*matH1.aw + matO16.aw
    wH = 2*matH1.aw / MH2O
    wO = matO16.aw / MH2O
    rhoH2O = 1

    MU = 0.0072 * matU235.aw + 0.9928 * matU238.aw
    w235 = 0.0072 * matU235.aw / MU
    w238 = 0.9928 * matU238.aw / MU
    rhoU = 19.1

    yH1 = []
    yO16 = []
    yU235 = []
    yU238 = []
    mt = matH1.mt.keys()
    for i in x:
        sum = 0
        for m in mt:
            sum += value_interp(matH1.mt[m], i)
        yH1.append(sum)
    mt = matO16.mt.keys()
    for i in x:
        sum = 0
        for m in mt:
            sum += value_interp(matO16.mt[m], i)
        yO16.append(sum)
    mt = matU235.mt.keys()
    for i in x:
        sum = 0
        for m in mt:
            sum += value_interp(matU235.mt[m], i)
        yU235.append(sum)
    mt = matU238.mt.keys()
    for i in x:
        sum = 0
        for m in mt:
            sum += value_interp(matU238.mt[m], i)
        yU238.append(sum)

    sigmaH20 = sigma((yH1, yO16), (wH, wO), (matH1.aw, matO16.aw), rhoH2O)
    sigmaU = sigma((yU235, yU238), (w235, w238), (matU235.aw, matU238.aw), rhoU)

    plt.figure()
    plt.title("Macroscopic total of H20 and U")
    plt.xlabel("Energy (MeV)")
    plt.ylabel("Macroscopic Σ (cm$^{-1}$)")
    plt.loglog(x, sigmaH20, label="H2O Σ (cm$^{-1}$)")
    plt.loglog(x, sigmaU, label="Nat U Σ (cm$^{-1}$)")
    plt.xlim(10**(-11), 20)
    plt.legend()
    plt.grid(True)
    plt.show()

#drawH2ONatU()

def printCollisionEnergy():
    x = []
    y = []
    with open("output/col_1758851991.csv", "r") as f:
        for l in f:
            l = l.rstrip().split(",")
            x.append(int(l[0]))
            y.append(float(l[1]))
    x = x[1:]
    y = y[1:]
    plt.figure()
    plt.title("Average Energy per Collision H2")
    plt.xlabel("Collision number")
    plt.ylabel("Energy (MeV)")
    plt.loglog(x, y)
    plt.legend()
    plt.grid(True)
    plt.show()

#printCollisionEnergy()

def printkeff():
    x = []
    y = []
    with open("output/keff_1758847988.csv", "r") as f:
        for l in f:
            l = l.rstrip().split(",")
            x.append(float(l[0]))
            y.append(float(l[1]))
    plt.figure()
    plt.title("k_eff per Enrichment")
    plt.xlabel("% of U235")
    plt.ylabel("k_eff")
    plt.plot(x, y)
    plt.legend()
    plt.grid(True)
    plt.show()

#printkeff()

def printHistTime():
    x = []
    y = []
    with open("output/time_1758854746.csv", "r") as f:
        for l in f:
            l = l.rstrip().split(",")
            x.append(float(l[0]))
            y.append(int(l[1]))
    plt.figure()
    plt.title("Timing of Fission (in neutron time)")
    plt.xlabel("Time form birth (s)")
    plt.ylabel("Collision #")
    plt.plot(x, y)
    plt.yscale("log")
    plt.legend()
    plt.grid(True)
    plt.show()

#printHistTime()
