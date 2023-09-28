#! /usr/bin/env python
######################################################################
#
# Name: generate_simple_weighted_template.py
#
# Purpose: Create templates
#
# Created: February-2020  Iker Loïc de Icaza Astiz (icaza@fnal.gov)
#
# Usage:
#
# generate_simple_weighted_template.py (--sbnd | --icarus) file
#
# Options:
#
# [-h|--help] - Print help message.
# (--sbnd or --icarus) to select for which experiment generate metrics
# Arguments:
#
# file  ... - Input file.
#
######################################################################

import sys
import os
import string
import argparse
import numpy as np
from time import sleep
from array import array

from ROOT import TStyle, TCanvas, TColor, TGraph, TGraphErrors
from ROOT import TH1D, TH2D, TProfile, TProfile3D, TFile, TF1
from ROOT import gROOT, TList, TTree, TDirectoryFile
import ROOT

try:
    # import project_utilities
    import larbatch_posix
    import fhicl
except ImportError:
    print("Failed to import 'larbatch_posix' or 'fhicl' modules")
    print("Setup 'fhiclpy' first:\n\tsetup fhiclpy vV_VV_VV -q QQQ:QQQQ")
    exit(1)


def pretty_print(string_to_prettify):
    print("\n")
    print("*" * 50)
    print(string_to_prettify)
    print("*" * 50)
    print("\n")
    return

class dotDict(dict):
    def __getattr__(self,val):
        return self[val]


class metrics_stuff:
    def __init__(self, name, pset):
        profile_option = 's'  # errors are the standard deviation
        self.name    = name
        self.flash_type = pset.FlashType
        self.x_bins  = pset.XBins
        self.x_bins_ = pset.x_bins_
        self.x_low   = pset.x_low
        self.x_up    = pset.x_up
        self.y_bins  = pset.y_bins
        self.y_low   = pset.y_low
        self.y_up    = pset.y_up
        self.z_bins  = pset.z_bins
        self.z_low   = pset.z_low
        self.z_up    = pset.z_up
        self.xvals   = np.arange(xbin_width/2., drift_distance, xbin_width)
        self.xerrs   = np.array([xbin_width/2.] * len(self.xvals))
        self.bins    = int(pset[name]['bins'])
        self.low     = pset[name]['low']
        self.up      = pset[name]['up']
        self.means   = [None] * self.x_bins
        self.spreads = [None] * self.x_bins
        self.h2      = TH2D(name+"_h2", "", self.x_bins, self.x_low, self.x_up,
                            self.bins, self.low, self.up)
        self.prof    = TProfile(name+"_prof", "",
                                self.x_bins, self.x_low, self.x_up,
                                self.low, self.up, profile_option)
        self.prof3   = TProfile3D(name+"_prof3", "",
                                  self.x_bins_, self.x_low, self.x_up,
                                  self.y_bins, self.y_low, self.y_up,
                                  self.z_bins, self.z_low, self.z_up,
                                  # self.low, self.up,
                                  profile_option)
        self.h1      = TH1D(name+"_h1", "", self.x_bins, self.x_low, self.x_up)

    def update_metrics(self):
        # fill histograms for match score calculation from profile histograms
        for ib in list(range(0, self.x_bins)):
            ibp = ib + 1
            self.h1.SetBinContent(ibp, self.prof.GetBinContent(ibp))
            self.h1.SetBinError(ibp, self.prof.GetBinError(ibp))
            self.means[ib] = self.prof.GetBinContent(ibp)
            self.spreads[ib] = self.prof.GetBinError(ibp)
            if self.spreads[ib] <= 0.001:
                print(f"Warning {self.name} spread close to zero.\n",
                      f"index: {ib}. spread: {self.spreads[ib]} \n",
                      "Setting to 0.001")
                self.spreads[ib] = 0.001

    def write_metrics(self):
        self.h2.Write()
        self.prof.Write()
        self.prof3.Write()
        self.h1.Write()

    def draw_metrics(self, canv, directory):
        self.h2.Draw()
        crosses = TGraphErrors(x_bins,
                             array('f', self.xvals), array('f', self.means),
                             array('f', self.xerrs), array('f', self.spreads))
        crosses.SetLineColor(ROOT.kAzure+9)
        crosses.SetLineWidth(3)
        crosses.Draw("Psame")
        canv.Print(directory.removesuffix("yzmaps") + f"{self.name}.pdf")
        canv.Update()

    def draw_3D(self, canv, directory):
        ROOT.gStyle.SetPalette(104)
        ROOT.gStyle.SetOptStat(0)
        m_min =  1.e8
        m_max = -1.e8
        s_min =  1.e8
        s_max = -1.e8
        for xb in list(range(1, self.x_bins_+1)):
            # get mins and maxes
            for yb in list(range(1, self.y_bins+1)):
                for zb in list(range(1, self.z_bins+1)):
                    mean = self.prof3.GetBinContent(xb,yb,zb)
                    sprd = self.prof3.GetBinError(xb,yb,zb)
                    if mean < m_min: m_min = mean
                    if mean > m_max: m_max = mean
                    if sprd < s_min: s_min = sprd
                    if sprd > s_max: s_max = sprd

        for xb in list(range(1, self.x_bins_+1)):
            yzmap_m = TH2D(f"{self.name}_yzmap_m_xb_{xb}",
                           f"{self.name}_yzmap_mean_xb_{xb}",
                           self.z_bins, self.z_low, self.z_up,
                           self.y_bins, self.y_low, self.y_up)
            if self.name == "dy" or self.name == "dz":
                zrange = max(abs(m_min), abs(m_max))
                yzmap_m.GetZaxis().SetRangeUser(-1.*zrange, zrange)
            else:
                yzmap_m.GetZaxis().SetRangeUser(m_min, m_max)
            yzmap_s = TH2D(f"{self.name}_yzmap_s_xb_{xb}",
                           f"{self.name}_yzmap_sprd_xb_{xb}",
                           self.z_bins, self.z_low, self.z_up,
                           self.y_bins, self.y_low, self.y_up)
            yzmap_s.GetZaxis().SetRangeUser(0., s_max)

            for yb in list(range(1, self.y_bins+1)):
                y_row_means = ''
                y_row_spdrs = ''
                for zb in list(range(1, self.z_bins+1)):
                    mean = self.prof3.GetBinContent(xb,yb,zb)
                    sprd = self.prof3.GetBinError(xb,yb,zb)
                    # y_row_means += str(f"{self.prof3.GetBinContent(xb,yb,zb):.1E}  ")
                    y_row_means += str(f"{mean:.1E}+-{sprd:.1E}  ")
                    yzmap_m.SetBinContent(zb, yb, mean)
                    yzmap_s.SetBinContent(zb, yb, sprd)
                # print(y_row_means)
            yzmap_m.Draw("colz")
            canv.Print(directory + f"/{self.name}_yzmap_mean_xb_{xb}.pdf")
            canv.Update()
            yzmap_s.Draw("colz")
            canv.Print(directory + f"/{self.name}_yzmap_sprd_xb_{xb}.pdf")
            canv.Update()
        print(f"Finish printing 3D maps of {self.name}")


def sig_fig_round(number, digits=3):
    power = "{:e}".format(number).split('e')[1]
    return round(number, -(int(power) - digits))

# Globally turn off root warnings.
# Don't let root see our command line options.
myargv = sys.argv
sys.argv = myargv[0:1]
if 'TERM' in os.environ:
    del os.environ['TERM']
ROOT.gErrorIgnoreLevel = ROOT.kError
sys.argv = myargv

# TODO: work out the best values for these
ROOT.gStyle.SetStatX(0.4)
ROOT.gStyle.SetStatY(0.88)
ROOT.gStyle.SetStatW(0.15)
# ROOT.gStyle.SetStatH(0.2)
ROOT.gStyle.SetOptStat(1111)
ROOT.gStyle.SetOptFit(1111)

gROOT.SetBatch(True) # to not show plots

detector = "experiment"
drift_distance = 0.
xbins = 0
xbin_width = 0.
time_delay = 0.
tolerable_time_diff = 0.1

# # Print help
# def help():

#     filename = sys.argv[0]
#     file = open(filename)

#     doprint = False

#     for line in file.readlines():
#         if line[2:9] == 'stat.py':
#             doprint = True
#         elif line[0:6] == '######' and doprint:
#             doprint = False
#         if doprint:
#             if len(line) > 2:
#                 print(line[2:], end=' ')
#             else:
#                 print()



def hypo_flashx_from_H2(flash_rr, rr_h2, flash_ratio, ratio_h2):
    rr_hypoX, rr_hypoXRMS = x_estimate_and_rms(flash_rr, rr_h2);
    ratio_hypoX, ratio_hypoXRMS = x_estimate_and_rms(flash_ratio, ratio_h2);

    drr2 = rr_hypoXRMS * rr_hypoXRMS;
    dratio2 = ratio_hypoXRMS * ratio_hypoXRMS;
    rr_hypoXWgt = 1./drr2;
    ratio_hypoXWgt = 1./dratio2;

    sum_weights = rr_hypoXWgt + ratio_hypoXWgt
    if sum_weights < 0.0002: return (-10., -10., rr_hypoX, rr_hypoXRMS, ratio_hypoX, ratio_hypoXRMS)
    hypo_x = (rr_hypoX*rr_hypoXWgt + ratio_hypoX*ratio_hypoXWgt) / sum_weights
    # consistent estimates, resulting error is smaller
    hypo_x_err = np.sqrt( 1. / sum_weights)
    if np.abs(rr_hypoX - ratio_hypoX) > 2.*np.sqrt(drr2+dratio2):
        # inconsistent estimates, resulting error is larger
        hypo_x_err = np.sqrt(drr2 + dratio2)
    return (hypo_x, hypo_x_err, rr_hypoX, rr_hypoXRMS,
            ratio_hypoX, ratio_hypoXRMS)


def x_estimate_and_rms(metric_value, metric_h2):
    kMinEntriesInProjection = 100
    bin_ = metric_h2.GetYaxis().FindBin(metric_value);
    bins = metric_h2.GetNbinsY();
    metric_hypoX = -1.;
    bin_buff = 0;
    while 0 < bin_-bin_buff or bin_+bin_buff <= bins :
        low_bin = bin_-bin_buff if 0 < bin_-bin_buff else 0
        high_bin = bin_+bin_buff if bin_+bin_buff <= bins else -1
        metric_px = metric_h2.ProjectionX("metric_px", low_bin, high_bin);
        if metric_px.GetEntries() > kMinEntriesInProjection :
            metric_hypoX = metric_px.GetRandom();
            metric_rmsX = metric_px.GetRMS();
            if metric_rmsX < xbin_width/2.: # something went wrong
                print(f"{metric_h2.GetName()} projected on metric_value: {metric_value}, "
                      f"bin: {bin_}, bin_buff: {bin_buff}; has {metric_px.GetEntries()} entries.")
                print(f"  metric_hypoX: {metric_hypoX}, metric_rmsX: {metric_rmsX}")
                return (-10., drift_distance); # no estimate
            return (metric_hypoX, metric_rmsX);
        bin_buff += 1;
    return (-10., drift_distance); # no estimate


def quality_checks(e, beam_spill_time_end):
    if e.slices != 1: return False
    if e.is_nu != 1: return False
    if e.mcT0 < 0. or beam_spill_time_end < e.mcT0 : return False
    if (e.flash_time-time_delay - e.mcT0) < 0. or (e.flash_time-time_delay - e.mcT0) > tolerable_time_diff : False
    if e.charge_x < 0. or e.charge_x > drift_distance: return False
    return True


def polynomial_correction(skew, hypo_x, pol_coeffs, skew_high_limit=10.):
    # TODO: maybe some other condition to prevent corrections?
    if np.abs(skew) > skew_high_limit or np.isnan(skew) or np.isnan(hypo_x):
        return 0.
    correction = 0.
    exponent = 1.
    for coeff in pol_coeffs:
        correction += coeff * exponent
        exponent *= hypo_x
    return correction * skew


def parameters_correction_fitter(nuslice_tree, var, flash_type, profile_bins,
                                 x_low, x_up, fit_func, beam_spill_time_end,
                                 skew_high_limit=10., skew_low_limit=0.05):
    fit_prof = TProfile(f"fit_prof_{var}", "", profile_bins,
                        x_low, x_up)
    draw_expression = (f"((flash_{var}b-charge_{var})/{var}_skew):new_hypo_x"
                       f">>fit_prof_{var}")
    # Another option for the fit: using charge_x instead of the hypoX estimate,
    # it might be justified to use the "truth" to make the fit,
    # however since the flashX estimate is biased towards shorter
    # distances near the cathode it's better to use the same as in runtime
    #
    # draw_expression = (f"((flash_{var}b-charge_{var})/{var}_skew):charge_x"
    #                    f">>fit_prof_{var}")
    #
    # HACK: Filter out from the fit the regions where discrepancy  is
    # not as bad, to give more weight to the edges where there large discrepancy
    filter_tolerable = "true"
    if detector == "sbnd":
        filter_tolerable = "abs(charge_y) > 60." if var=="y" \
            else "(charge_z<120. || 380.<charge_z)"
    elif detector == "icarus":
        filter_tolerable = "(charge_y<-65. || 19.<charge_y)" if var=="y" \
            else "(charge_z<-800. || 800.<charge_z)"

    # the filters are the ones in quality_checks(), plus the one above
    # to give more weight to the edges
    draw_filters = (f"abs({var}_skew)>{skew_low_limit} && "
                    f"abs({var}_skew)<{skew_high_limit} && "
                    f"is_nu==1 && slices==1 && "
                    f"0.<=mcT0 && mcT0<={beam_spill_time_end} && "
                    f"(flash_time-{time_delay} - mcT0) >= 0. && (flash_time-{time_delay} - mcT0) <= {tolerable_time_diff} && "
                    f"charge_x >= {x_low} && charge_x <= {x_up} && new_hypo_x >= 0. &&"
                    f"{filter_tolerable}"
                    )
    draw_option = "prof"
    print("Draw expression: ", draw_expression,
          "\nfilters: ", draw_filters,
          "\noptions: ", draw_option)
    can = TCanvas("can")
    nuslice_tree.Draw(draw_expression, draw_filters, draw_option)
    fit_result = fit_prof.Fit(fit_func, "S")
    fit_prof.Write()
    # fit_result.Print("V")
    can.Print(f"plots/{flash_type}/{var}_correction_fit.pdf")
    params = []
    for p in fit_result.Parameters():
        params.append(sig_fig_round(p, 3))
    print("The fitted and rounded correction parameters for ", var, " are: ", params)
    print("These are now stored in the metrics file.\n")
    return params


def generator(nuslice_tree, rootfile, pset, suffix, flash_type):
    # BIG TODO: Metrics should depend on X,Y,Z.
    # Many changes needed everywhere
    half_bin_width = xbin_width/2.

    directory = "plots/" + flash_type + "/yzmaps"
    if not os.path.exists(directory):
        os.makedirs(directory)

    if detector == "sbnd":
        x_gl_low = -215
        x_gl_up = 215
    elif detector == "icarus":
        if pset.Cryostat == 0:
            x_gl_low = -380
            x_gl_up = -50
        if pset.Cryostat == 1:
            x_gl_low = 380
            x_gl_up = 50
    beam_spill_time_end = pset.BeamSpillTimeEnd - pset.BeamSpillTimeStart

    md = {
        'dy':    metrics_stuff('dy', pset),
        'dz':    metrics_stuff('dz', pset),
        'rr':    metrics_stuff('rr', pset),
        'ratio': metrics_stuff('ratio', pset),
        'slope': metrics_stuff('slope', pset),
        'petoq': metrics_stuff('petoq', pset)
    }

    unfolded_score_scatter = TH2D("unfolded_score_scatter", "Scatter plot of match scores",
                                  2*pset.XBins, x_gl_low, x_gl_up,
                                  pset['score']['bins'], pset['score']['low'], pset['score']['up'])
    unfolded_score_scatter.GetXaxis().SetTitle("X global (cm)")
    unfolded_score_scatter.GetYaxis().SetTitle("match score (arbitrary)")
    oldunfolded_score_scatter = TH2D("oldunfolded_score_scatter", "Scatter plot of match scores, old metrics",
                                     2*pset.XBins, x_gl_low, x_gl_up,
                                     pset['score']['bins'], pset['score']['low'], pset['score']['up'])
    oldunfolded_score_scatter.GetXaxis().SetTitle("X global (cm)")
    oldunfolded_score_scatter.GetYaxis().SetTitle("match score (arbitrary)")
    unfolded_score_scatter_3D = TH2D("unfolded_score_scatter_3D", "Scatter plot of match scores 3D",
                                     2*pset.XBins, x_gl_low, x_gl_up,
                                     pset['score']['bins'], pset['score']['low'], pset['score']['up'])
    unfolded_score_scatter_3D.GetXaxis().SetTitle("X global (cm)")
    unfolded_score_scatter_3D.GetYaxis().SetTitle("match score (arbitrary)")

    match_score_scatter = TH2D("match_score_scatter", "Scatter plot of match scores",
                               pset.XBins, pset.x_low, pset.x_up,
                               pset['score']['bins'], pset['score']['low'], pset['score']['up'])
    match_score_scatter.GetXaxis().SetTitle("distance from anode (cm)")
    match_score_scatter.GetYaxis().SetTitle("match score (arbitrary)")
    match_score_h1 = TH1D("match_score", "Match Score",
                          pset['score']['bins'], pset['score']['low'], pset['score']['up'])
    metrics_filename = 'fm_metrics_' + detector + '.root'
    hfile_top = TFile(metrics_filename, 'UPDATE',
                  'Simple flash matching metrics for ' + detector.upper())
    keys = ROOT.gDirectory.GetListOfKeys()
    print(keys)
    print(ftype_long)
    dir_exists = False
    for key in keys:
        if ftype_long == key.GetName():
            dir_exists = True
    hfile = hfile_top.Get(metrics_filename+":/" + ftype_long) if(dir_exists==True) else \
            TDirectoryFile(ftype_long, f"Metrics directory for {ftype_long}")
    hfile.cd()

    # fill rr_h2 and ratio_h2 first
    for e in nuslice_tree:
        oldunfolded_score_scatter.SetFillColor(ROOT.kBlue)
        unfolded_score_scatter.SetFillColor(ROOT.kBlue)
        unfolded_score_scatter_3D.SetFillColor(ROOT.kBlue)
        match_score_scatter.SetFillColor(ROOT.kBlue)

        if not quality_checks(e, beam_spill_time_end):
            oldunfolded_score_scatter.SetFillColor(ROOT.kRed)
            unfolded_score_scatter.SetFillColor(ROOT.kRed)
            unfolded_score_scatter_3D.SetFillColor(ROOT.kRed)
            match_score_scatter.SetFillColor(ROOT.kRed)

        qX = e.charge_x
        md['rr'].h2.Fill(qX, e.flash_rr)
        md['rr'].prof.Fill(qX, e.flash_rr)
        md['rr'].prof3.Fill(e.charge_x, e.charge_y, e.charge_z, e.flash_rr)
        md['ratio'].h2.Fill(qX, e.flash_ratio)
        md['ratio'].prof.Fill(qX, e.flash_ratio)
        md['ratio'].prof3.Fill(e.charge_x, e.charge_y, e.charge_z, e.flash_ratio)

    # Use rr_h2 and ratio_h2 to compute flash drift distance
    # estimates, and store them as 'new_'...
    new_hypo_x = array('d',[0])
    new_hypo_x_branch = nuslice_tree.Branch("new_hypo_x", new_hypo_x, "new_hypo_x/D");
    new_hypo_x_err = array('d',[0])
    new_hypo_x_err_branch = nuslice_tree.Branch("new_hypo_x_err", new_hypo_x_err, "new_hypo_x_err/D");
    new_hypo_x_rr = array('d',[0])
    new_hypo_x_rr_branch = nuslice_tree.Branch("new_hypo_x_rr", new_hypo_x_rr, "new_hypo_x_rr/D");
    new_hypo_x_rr_err = array('d',[0])
    new_hypo_x_rr_err_branch = nuslice_tree.Branch("new_hypo_x_rr_err", new_hypo_x_rr_err, "new_hypo_x_rr_err/D");
    new_hypo_x_ratio = array('d',[0])
    new_hypo_x_ratio_branch = nuslice_tree.Branch("new_hypo_x_ratio", new_hypo_x_ratio, "new_hypo_x_ratio/D");
    new_hypo_x_ratio_err = array('d',[0])
    new_hypo_x_ratio_err_branch = nuslice_tree.Branch("new_hypo_x_ratio_err", new_hypo_x_ratio_err, "new_hypo_x_ratio_err/D");
    for e in nuslice_tree:
        # No need to check quality in this loop
        hypo_x, hypo_x_err, rr_hypoX, rr_hypoXRMS, ratio_hypoX, ratio_hypoXRMS = \
            hypo_flashx_from_H2(e.flash_rr, md['rr'].h2,
                                e.flash_ratio, md['ratio'].h2)
        new_hypo_x[0] = hypo_x
        new_hypo_x_branch.Fill()
        new_hypo_x_err[0] = hypo_x_err
        new_hypo_x_err_branch.Fill()
        new_hypo_x_rr[0] = rr_hypoX
        new_hypo_x_rr_branch.Fill()
        new_hypo_x_rr_err[0] = rr_hypoXRMS
        new_hypo_x_rr_err_branch.Fill()
        new_hypo_x_ratio[0] = ratio_hypoX
        new_hypo_x_ratio_branch.Fill()
        new_hypo_x_ratio_err[0] = ratio_hypoXRMS
        new_hypo_x_ratio_err_branch.Fill()

    # Fit and create std::vector objects for polynomial correction
    # coefficients
    y_pol_coeffs = parameters_correction_fitter(nuslice_tree, "y", flash_type, pset.XBins,
                                                0., pset.DriftDistance,
                                                pset.fit_func_y,
                                                beam_spill_time_end,
                                                pset.SkewLimitY)
    y_pol_coeffs_vec = ROOT.std.vector['double']()
    for yp in y_pol_coeffs: y_pol_coeffs_vec.push_back(yp)

    z_pol_coeffs = parameters_correction_fitter(nuslice_tree, "z", flash_type, pset.XBins,
                                                0., pset.DriftDistance,
                                                pset.fit_func_z,
                                                beam_spill_time_end,
                                                pset.SkewLimitZ)
    if detector == "icarus":
        z_pol_coeffs = [0.] # No Z corrections for ICARUS
    z_pol_coeffs_vec = ROOT.std.vector['double']()
    for zp in z_pol_coeffs: z_pol_coeffs_vec.push_back(zp)

    # Using the new estimation new_hypo_x, and the just fitted
    # polynomial coefficients; get the corrected new_flash_y and new_flash_z
    new_flash_y = array('d',[0])
    new_flash_y_branch = nuslice_tree.Branch("new_flash_y", new_flash_y, "new_flash_y/D");
    new_flash_z = array('d',[0])
    new_flash_z_branch = nuslice_tree.Branch("new_flash_z", new_flash_z, "new_flash_z/D");
    for e in nuslice_tree:
        # No need to check quality in this loop
        new_flash_y[0] = e.flash_yb - polynomial_correction(
            e.y_skew, e.new_hypo_x, y_pol_coeffs, pset.SkewLimitY)
        new_flash_y_branch.Fill()
        new_flash_z[0] = e.flash_zb - polynomial_correction(
            e.z_skew, e.new_hypo_x, z_pol_coeffs, pset.SkewLimitZ)
        new_flash_z_branch.Fill()

    # Update the file
    rootfile.Write()
    hfile_top.ReOpen("UPDATE")
    hfile.cd()
    hfile.Write()

    # Use the new corrected terms to fill the rest of H2s and Profs
    for e in nuslice_tree:
        if not quality_checks(e, beam_spill_time_end): continue
        qX = e.charge_x
        md['dy'].h2.Fill(qX, e.new_flash_y - e.charge_y)
        md['dy'].prof.Fill(qX, e.new_flash_y - e.charge_y)
        md['dy'].prof3.Fill(e.charge_x, e.charge_y, e.charge_z, e.new_flash_y - e.charge_y)
        md['dz'].h2.Fill(qX, e.new_flash_z - e.charge_z)
        md['dz'].prof.Fill(qX, e.new_flash_z - e.charge_z)
        md['dz'].prof3.Fill(e.charge_x, e.charge_y, e.charge_z, e.new_flash_z - e.charge_z)
        # these H2s below use no new corrections
        # md['slope'].h2.Fill(qX, e.flash_slope - e.charge_slope)
        # md['slope'].prof.Fill(qX, e.flash_slope - e.charge_slope)
        # md['slope'].prof3.Fill(e.charge_x, e.charge_y, e.charge_z, e.flash_slope - e.charge_slope)
        md['slope'].h2.Fill(qX, e.flash_xw)
        md['slope'].prof.Fill(qX, e.flash_xw)
        md['slope'].prof3.Fill(e.charge_x, e.charge_y, e.charge_z, e.flash_xw)
        md['petoq'].h2.Fill(qX, e.petoq)
        md['petoq'].prof.Fill(qX, e.petoq)
        md['petoq'].prof3.Fill(e.charge_x, e.charge_y, e.charge_z, e.petoq)

    # fill histograms for match score calculation from profile histograms
    for m in md.values():
        m.update_metrics()

    for e in nuslice_tree:
        # calculate match score
        # if not quality_checks(e, beam_spill_time_end): continue
        qX = e.charge_x
        qXGl = e.charge_x_gl
        isl = int(qX/xbin_width)
        score = 0.
        score += abs((e.new_flash_y-e.charge_y) - md['dy'].means[isl])/md['dy'].spreads[isl]
        score += abs((e.new_flash_z-e.charge_z) - md['dz'].means[isl])/md['dz'].spreads[isl]
        score += abs(e.flash_ratio-md['ratio'].means[isl])/md['ratio'].spreads[isl]
        # score += abs((e.flash_slope - e.charge_slope) - md['slope'].prof3.GetBinContent(xb,yb,zb))/md['slope'].prof3.GetBinError(xb,yb,zb) # TODO: if useful add it to the total score
        # score += abs(e.flash_xw - md['slope'].means[isl])/md['slope'].spreads[isl] # TODO: if useful add it to the total score
        score += abs(e.petoq-md['petoq'].means[isl])/md['petoq'].spreads[isl]

        xb  = md['dy'].prof3.GetXaxis().FindBin(e.charge_x)
        if xb < 1: xb = 1
        elif xb > pset.x_bins_: xb = pset.x_bins_
        yb  = md['dy'].prof3.GetYaxis().FindBin(e.charge_y)
        if yb < 1: yb = 1
        elif yb > pset.y_bins: yb = pset.y_bins
        zb  = md['dy'].prof3.GetZaxis().FindBin(e.charge_z)
        if zb < 1: zb = 1
        elif zb > pset.z_bins: zb = pset.z_bins
        try:
            score_3D = 0.
            score_3D += abs((e.new_flash_y-e.charge_y) - md['dy'].prof3.GetBinContent(xb,yb,zb))/md['dy'].prof3.GetBinError(xb,yb,zb)
            score_3D += abs((e.new_flash_z-e.charge_z) - md['dz'].prof3.GetBinContent(xb,yb,zb))/md['dz'].prof3.GetBinError(xb,yb,zb)
            score_3D += abs(e.flash_rr-md['rr'].prof3.GetBinContent(xb,yb,zb))/md['rr'].prof3.GetBinError(xb,yb,zb)
            score_3D += abs(e.flash_ratio-md['ratio'].prof3.GetBinContent(xb,yb,zb))/md['ratio'].prof3.GetBinError(xb,yb,zb)
            # score_3D += abs((e.flash_slope - e.charge_slope) - md['slope'].prof3.GetBinContent(xb,yb,zb))/md['slope'].prof3.GetBinError(xb,yb,zb) # TODO: if useful add it to the total score
            # score_3D += abs(e.flash_xw - md['slope'].prof3.GetBinContent(xb,yb,zb))/md['slope'].prof3.GetBinError(xb,yb,zb) # TODO: if useful add it to the total score
            score_3D += abs(e.petoq-md['petoq'].prof3.GetBinContent(xb,yb,zb))/md['petoq'].prof3.GetBinError(xb,yb,zb)
        except:
            print(xb, f"{e.charge_x:.3E}", yb, f"{e.charge_y:.3E}", zb, f"{e.charge_z:.3E}", md['dy'].prof3.GetBinContent(xb,yb,zb), md['dy'].prof3.GetBinError(xb,yb,zb))
            pass


        oldunfolded_score_scatter.Fill(qXGl, e.score)
        unfolded_score_scatter.Fill(qXGl, score)
        unfolded_score_scatter_3D.Fill(qXGl, score_3D)
        match_score_scatter.Fill(qX, score)
        match_score_h1.Fill(score)

    for m in md.values():
        m.write_metrics()

    hfile.WriteObject(y_pol_coeffs_vec, "pol_coeffs_y")
    hfile.WriteObject(z_pol_coeffs_vec, "pol_coeffs_z")
    match_score_scatter.Write()
    oldunfolded_score_scatter.Write()
    unfolded_score_scatter.Write()
    unfolded_score_scatter_3D.Write()
    match_score_h1.Write()
    hfile.Close()

    canv = TCanvas("canv")
    for m in md.values():
        m.draw_metrics(canv, directory)
        m.draw_3D(canv, directory)

    oldunfolded_score_scatter.Draw()
    canv.Print(f"plots/{flash_type}/oldunfolded_score_scatter.pdf")
    canv.Update()

    unfolded_score_scatter.Draw()
    canv.Print(f"plots/{flash_type}/unfolded_score_scatter.pdf")
    canv.Update()

    unfolded_score_scatter_3D.Draw()
    canv.Print(f"plots/{flash_type}/unfolded_score_scatter_3D.pdf")
    canv.Update()

    match_score_scatter.Draw()
    canv.Print(f"plots/{flash_type}/match_score_scatter.pdf")
    canv.Update()

    match_score_h1.Draw()
    canv.Print(f"plots/{flash_type}/match_score.pdf")
    canv.Update()
    sleep(20)

# Main program.
def main():

    # Parse arguments.
    parser = argparse.ArgumentParser(prog='generate_simple_weighted_template.py')
    parser.add_argument('file')
    # parser.add_argument('--help', '-h',
    #                     action='store_true',
    #                     help='help flag' )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--sbnd', action='store_true', help='Generate metrics for SBND')
    group.add_argument('--icarus', action='store_true', help='Generate metrics for ICARUS')
    args = parser.parse_args()

    # if args.help:
    #     print("To run do:\n"/
    #           "generate_simple_weighted_template.py file.root\n"/
    #           "where file.root has a fmatch/nuslicetree")
    #     return(0)
    if args.sbnd :
        print("Generate metrics for SBND")
    elif args.icarus :
        print("Generate metrics for ICARUS")

    if not larbatch_posix.exists(args.file):
        print('Input file %s does not exist.' % args.file)
        return 1

    print('\nOpening ', args.file)
    rootfile_orig = TFile.Open(args.file, 'READ')
    if not rootfile_orig.IsOpen() or rootfile_orig.IsZombie():
        print('Failed to open ', args.file)
        return 1

    file_updated = "updated_" + args.file
    print('\nCopying  ', args.file, ' to ', file_updated)
    if not rootfile_orig.Cp(file_updated):
        print('Failed to copy ', args.file)
        return 1
    rootfile = TFile.Open(file_updated, 'UPDATE')
    if not rootfile.IsOpen() or rootfile.IsZombie():
        print('Failed to open ', file_updated)
        return 1

    global detector
    global drift_distance
    global x_bins
    global xbin_width
    global time_delay # TODO: Improve timing res to make this obsolete
    global tolerable_time_diff
    global suffix_list
    global ftype_long
    if args.sbnd:
        detector = "sbnd"
        metrics_filename = 'fm_metrics_' + detector + '.root'
        hfile_top = TFile(metrics_filename, 'RECREATE',
                      'Simple flash matching metrics for ' + detector.upper())
        hfile_top.Close()
        suffix_list = ["", "_op", "_ara", "_opara"]
        time_delays = [0.15, 0, 0.04, 0]
        long_names = ["SimpleFlash_PMT", "OpFlash_PMT", "SimpleFlash_ARA", "OpFlash_ARA"]
        for (dir_, t_delay, l_name) in zip(suffix_list, time_delays, long_names):
            fcl_params = fhicl.make_pset('flashmatch_sbnd.fcl')
            fhicl_table = "sbnd_simple_flashmatch" + dir_
            pset = dotDict(fcl_params[fhicl_table])
            time_delay = t_delay
            ftype_long = pset.FlashType
            dir = rootfile.Get(file_updated+":/fmatch"+dir_.replace("_", ""))
            nuslice_tree = dir.Get("nuslicetree")
            drift_distance = pset.DriftDistance
            x_bins = pset.XBins
            xbin_width = drift_distance/x_bins
            pretty_print(l_name)
            generator(nuslice_tree, rootfile, pset, dir_, l_name)
    elif args.icarus:
        fcl_params = fhicl.make_pset('flashmatch_simple_icarus.fcl')
        time_delay = 0.
        # by default merge the trees from both Cryos
        metrics_filename = 'fm_metrics_' + detector + '.root'
        hfile_top = TFile(metrics_filename, 'RECREATE',
                      'Simple flash matching metrics for ' + detector.upper())
        hfile_top.Close()
        suffix_list = ["", "_op"]
        long_names = ["SimpleFlash_PMT", "OpFlash_PMT"]
        for (dir_, l_name) in zip(suffix_list, long_names):
            detector = "icarus"
            pset = dotDict(fcl_params['icarus_simple_flashmatch_E' + dir_])
            fhicl_table_E = "icarus_simple_flashmatch_E" + dir_
            fhicl_table_W = "icarus_simple_flashmatch_W" + dir_
            dir0 = rootfile.Get(file_updated+":/fmatch"+dir_.replace("_", "")+"CryoE")
            nuslice_tree0 = dir0.Get("nuslicetree")
            dir1 = rootfile.Get(file_updated+":/fmatch"+dir_.replace("_", "")+"CryoW")
            nuslice_tree1 = dir1.Get("nuslicetree")
            treelist = TList()
            treelist.Add(nuslice_tree0)
            treelist.Add(nuslice_tree1)
            nuslice_tree = TTree.MergeTrees(treelist);
            nuslice_tree.SetName("nuslice_tree");
            # nuslice_tree.Write();
            ftype_long = pset.FlashType
            drift_distance = pset.DriftDistance
            x_bins = pset.XBins
            xbin_width = drift_distance/x_bins
            pretty_print(l_name)
            generator(nuslice_tree, rootfile, pset, dir_, l_name)


#    generator(nuslice_tree, rootfile, pset)


if __name__ == '__main__':
    sys.exit(main())
