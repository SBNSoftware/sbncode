import numpy as np
from scipy.interpolate import interp2d
import landau
from matplotlib import gridspec
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple

# Calculate the MPV dEdx v. R.R. for muons

# From: https://pdg.lbl.gov/2012/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.pdf
LAr_density_gmL = 1.3973

# Reference data
CSDA_RR_REF = np.array([
    9.833e-1,
    1.786e0,
    3.321e0,
    6.598e0,
    1.058e1,
    3.084e1,
    4.250e1,
    6.732e1,
    1.063e2,
    1.725e2,
    2.385e2,
    4.934e2,
    6.163e2
]) / LAr_density_gmL

KE_REF = np.array([
    10.,
    14.,
    20.,
    30.,
    40.,
    80.,
    100.,
    140.,
    200.,
    300.,
    400.,
    800.,
    1000.
])

# Constants
mass_electron = 0.5109989461 # MeV https://pdg.lbl.gov/2020/listings/rpp2020-list-K-plus-minus.pdf
mass = 105.6583745 # MeV https://pdg.lbl.gov/2020/listings/rpp2020-list-muon.pdf
Ival = 188.0e-6
Zval = 18.0
Aval = 39.948
Kfactor = 0.307075

def Calc_MPV_DEDX(pitch, T):
    gamma = (mass+T)/mass
    beta = np.power(1.0-np.power(gamma,-2.0),0.5)
    Wmax = (2.0*mass_electron*np.power(beta,2.0)*np.power(gamma,2.0))/(1.0+2.0*gamma*(mass_electron/mass)+np.power(mass_electron/mass,2.0))

    # Medium energy 
    dens_factor = 2.0*np.log(10)*np.log10(beta*gamma)-5.2146+0.19559*np.power(3.0-np.log10(beta*gamma),3.0)
    # low energy
    dens_factor[np.log10(beta*gamma) < 0.2] = 0.
    # high energy
    dens_factor[np.log10(beta*gamma) > 3.0] = (2.0*np.log(10)*np.log10(beta*gamma)-5.2146)[np.log10(beta*gamma) > 3.0]
    xi = (Kfactor/2.0)*(Zval/Aval)*np.power(beta,-2.0)*LAr_density_gmL*pitch
    kappa = xi/Wmax
    dEdx_MPV = xi*(np.log((2.0*mass_electron*np.power(beta*gamma,2.0))/Ival)+np.log(xi/Ival)+0.200-np.power(beta,2.0)-dens_factor)/pitch
  
    return dEdx_MPV

def Calc_MEAN_DEDX(T):
    gamma = (mass+T)/mass
    beta = np.power(1.0-np.power(gamma,-2.0),0.5)
    Wmax = (2.0*mass_electron*np.power(beta,2.0)*np.power(gamma,2.0))/(1.0+2.0*gamma*(mass_electron/mass)+np.power(mass_electron/mass,2.0))

    # Medium energy 
    dens_factor = 2.0*np.log(10)*np.log10(beta*gamma)-5.2146+0.19559*np.power(3.0-np.log10(beta*gamma),3.0)
    # low energy
    dens_factor[np.log10(beta*gamma) < 0.2] = 0.
    dens_factor[beta < 1e-6] = 0.
    # high energy
    dens_factor[np.log10(beta*gamma) > 3.0] = (2.0*np.log(10)*np.log10(beta*gamma)-5.2146)[np.log10(beta*gamma) > 3.0]
    dEdx_mean = LAr_density_gmL*Kfactor*(Zval/Aval)*np.power(beta,-2.0)*(0.5*np.log(2.0*mass_electron*np.power(beta,2.0)*np.power(gamma,2.0)*Wmax*np.power(Ival,-2.0))-np.power(beta,2.0)-dens_factor/2.0)

    return dEdx_mean

# Map R.R. to KE
def make_mpv_map():
    KE_points_max = 1000.
    dRR = 0.01
    thisKE = KE_points_max
    
    KE_points = [thisKE]
    RR_points = [0.]
    
    while thisKE > 0.0:
        deltaKE = Calc_MEAN_DEDX(np.array([thisKE])) * dRR
        RR_points.append(RR_points[-1] + dRR)
        thisKE -= deltaKE[0]
        KE_points.append(thisKE)
    
    KE_points = np.array(list(reversed(KE_points[:-1])))
    RR_points = np.array(RR_points[:-1])
    
    
    # Map KE to MPV dE/dx
    
    # pitches
    PITCH_points = np.linspace(0.2, 3, 28*100+1)
    
    KE_points_2d = np.tile(KE_points, (PITCH_points.size, 1))
    RR_points_2d = np.tile(RR_points, (PITCH_points.size, 1))
    PITCH_points_2d = np.tile(PITCH_points, (KE_points.size, 1)).T
    
    MPV_dEdx_points_2d = Calc_MPV_DEDX(PITCH_points_2d, KE_points_2d)
    
    RRpitch2dEdx = interp2d(RR_points, PITCH_points, MPV_dEdx_points_2d, kind="cubic")

    return RRpitch2dEdx

RRpitch2dEdx = make_mpv_map()

# ArgoNeuT params
MODA = 0.930
MODB = 0.212
Wion = 1e3 / 4.237e7
Efield = 0.5

def recombination(dEdx, A=MODA, B=MODB, E=Efield):
    alpha = A
    beta = B / (LAr_density_gmL * E)
    
    dQdx = np.log(alpha + dEdx*beta) / (Wion * beta)
    return dQdx

# ICARUS params
k = 0.0486
A = 0.8

def Birks_recombination(dEdx):
    R =  A / (1 + k*dEdx / (Efield*LAr_density_gmL))
    return R * dEdx / Wion

def landau_gaus(X, *p):
    mpv, eta, sigma, A = p
    sigma = np.minimum(sigma, 100*eta)
    return landau.landau.gauss_landau(X, mpv, eta, sigma, A)

def langau_chi2(x, y, yerr, popt):
    return np.sum(((landau_gaus(x, *popt) - y) / yerr)**2)

def gain_predicted_MPV(RRs, CAL, pitch, A=MODA, B=MODB, E=Efield):
    dEdxs = RRpitch2dEdx(RRs, pitch)
    dQdxs = recombination(dEdxs, A, B, E)
    return dQdxs / CAL

def gain_predicted_MPV_Birks(RRs, CAL, pitch):
    dEdxs = RRpitch2dEdx(RRs, pitch)
    dQdxs = Birks_recombination(dEdxs)
    return dQdxs / CAL

def gain_chi2(RRs, CAL, MPV, err, pitch, when, A=MODA, B=MODB):
    dEdxs = RRpitch2dEdx(RRs, pitch)
    dQdxs = recombination(dEdxs, A, B)
    dQdxs_ADC = np.outer(1. / CAL, dQdxs[when])
    chi2s = (MPV[when] - dQdxs_ADC)**2 / err[when]**2
    return np.sum(chi2s, axis=-1)

def gain_chi2_Birks(RRs, CAL, MPV, err, pitch, when):
    dEdxs = RRpitch2dEdx(RRs, pitch)
    dQdxs = Birks_recombination(dEdxs)
    dQdxs_ADC = np.outer(1. / CAL, dQdxs[when])
    chi2s = (MPV[when] - dQdxs_ADC)**2 / err[when]**2
    return np.sum(chi2s, axis=-1)

def valid_mpv(RRs, MPV, err, minRR=2.):
    return np.isfinite(err) & (RRs > minRR) & (err/MPV < 0.1)

def calibrate_plot(fig, Xlist, preds, datas, errs, text=None, title=None, labels=None):
    if not isinstance(preds, list):
        preds = [preds]

    if not isinstance(datas, list):
        datas = [datas]

    if not isinstance(Xlist, list):
        Xlist = [Xlist]

    assert(len(preds) == len(datas) == len(Xlist))

    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1]) 
    ax1 = fig.subplot(gs[0])    
    ax2 = fig.subplot(gs[1], sharex = ax1)

    colors = fig.rcParams['axes.prop_cycle'].by_key()['color']

    ps = []
    ds = []
    
    for i, (Xs,pred) in enumerate(zip(Xlist,preds)):
        color = None if len(preds) == 1 else colors[i]
        label = "Cal. Fit"
        if labels is not None: 
            label = label + " " + labels[i]
        p = ax1.plot(Xs, pred, label=label, color=color)
        ps.append(p[0])

    for i, (Xs, data, err) in enumerate(zip(Xlist, datas, errs)):
        color = None if len(preds) == 1 else colors[i]
        label = "Data M.P.V."
        if labels is not None: 
            label = label + " " + labels[i]
        d = ax1.errorbar(Xs, data, yerr=err, ls="none", 
                     color=color, label=label, marker=".", markersize=5.)
        ds.append(d)

    leg_labels = ["Fit/Data" for _ in ps]
    if labels is not None:
        leg_labels = [ll+l for ll,l in zip(leg_labels, labels)]

    ax1.legend(list(zip(ps, ds)), leg_labels, numpoints=1, 
        handler_map={tuple: HandlerTuple(ndivide=None)})

    ax2.set_label("Residual Range [cm]")
    ax1.set_ylabel("dQ/dx [ADDC/cm]")
    
    if text:
        ax1.text(0.5, 2.25, text, transform=fig.gca().transAxes, verticalalignment="top")
    if title:
        ax1.set_title(title)
        
    ax1.get_shared_x_axes().join(ax1, ax2)
    fig.subplots_adjust(hspace=.0)
    
    for i, (Xs, data, pred, err) in enumerate(zip(Xlist, datas, preds, errs)):
        color = None if len(preds) == 1 else colors[i]
        ax2.plot(Xs, (data - pred) / err, color=color) #, linestyle="None", marker=".", markersize=5)

    ax2.set_ylim([-4, 4])
    ax2.axhline(0, color="gray")
    ax2.axhline(1, color="gray", linestyle="--")
    ax2.axhline(-1, color="gray", linestyle="--")
    ax2.axhline(2, color="gray", linestyle=":")
    ax2.axhline(-2, color="gray", linestyle=":")
    ax2.set_ylabel("Data - Pred. [$\\sigma$]")
    ax2.set_xlabel("Residual Range [cm]")
    
    yticks = ax2.yaxis.get_major_ticks()
    yticks[-1].label1.set_visible(False)
    return ax1, ax2
