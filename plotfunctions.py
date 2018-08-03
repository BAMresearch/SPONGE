import matplotlib
import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')

def simPlot(dataset, RayleighR = None, title = None):
    plt.errorbar(dataset.Q, dataset.I * 10. / 1.1, dataset.IError * 10. / 1.1, label = title)
    if RayleighR is not None:
        q = dataset.Q
        f = 3. * (np.sin(q*RayleighR) - q*RayleighR * np.cos(q*RayleighR)) / ((q*RayleighR)**3.) # Rayleigh form factor
        plt.loglog(q, f**2, label = "Rayleigh function, R = {}".format(RayleighR))
    plt.xscale("log")
    plt.yscale("log")
    plt.grid("on")
    plt.legend(loc = 0)
    plt.xlabel("q (1/nm)")
    plt.ylabel("I (A.U.)")

