## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Plotting Outputs
#
#  (from covariance and sensitivity
#  calculations)
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import ROOT
from ROOT import TFile, TCanvas, TH2D, TH1D, TGraph, TGraph2D, TStyle, TLegend, THStack, TPad
import argparse
import os
import math


## Main function
## ~~~~~~~~~~~~~

def make_count_plot(args):
    
    # Count plots
    
    countfile = TFile(args.cntfile)
    
    dets = []
    for key in countfile.GetListOfKeys():
        key = key.GetName()
        if key[:key.find("_")] not in dets:
            dets.append(key[:key.find("_")])
    
    samples = []
    for key in countfile.GetListOfKeys():
        key = key.GetName()
        if key[key.find("_")+1:key.find("_")+1+key[key.find("_")+1:].find("_")] not in samples:
            samples.append(key[key.find("_")+1:key.find("_")+1+key[key.find("_")+1:].find("_")])
    
    samplename = {"numu": "#nu_{#mu}", "nue": "#nu_{e}"}
    canvases = []; stacks = []; ct_hists = []; bkg_hists = []; legends = []
    for s, sample in enumerate(samples):
        
        canvases.append(TCanvas(sample+"_canvas"))
        canvases[s].Divide(len(dets), 1)
        
        ct_hists.append([])
        bkg_hists.append([])
        stacks.append([])
        
        legends.append([])
        
        for d, det in enumerate(dets):
            
            ct_hists[s].append(countfile.Get(det+"_"+sample+"_cts"))
            ct_hists[s][d].SetLineColor(38)
            ct_hists[s][d].SetFillColor(38)
            
            bkg_hists[s].append(countfile.Get(det+"_"+sample+"_bkg"))
            bkg_hists[s][d].SetLineColor(30)
            bkg_hists[s][d].SetFillColor(30)
            
            for b in range(ct_hists[s][d].GetNbinsX()):
                ct_hists[s][d].SetBinContent(1+b, ct_hists[s][d].GetBinContent(1+b) / ct_hists[s][d].GetBinWidth(1+b))
                bkg_hists[s][d].SetBinContent(1+b, bkg_hists[s][d].GetBinContent(1+b) / bkg_hists[s][d].GetBinWidth(1+b))
             
            stacks[s].append(THStack(det+"_"+sample+"_stack", det+"; Counts; Energy (GeV)"))
            stacks[s][d].Add(bkg_hists[s][d])
            stacks[s][d].Add(ct_hists[s][d])
            
            legends[s].append(TLegend())
            legends[s][d].AddEntry(ct_hists[s][d], 'CC -> #mu + X', 'f')
            legends[s][d].AddEntry(bkg_hists[s][d], 'NC -> #pi^{#pm} + X', 'f')
            
            canvases[s].cd(d+1)
            stacks[s][d].Draw("hist")
            legends[s][d].Draw()
            
        # NOTE: This will have to be changed when nue samples are added; this assumed a simple
        #       count + background histogram stack, whereas nue plots in the proposal have the 
        #       background broken up into many different 'types'
        
        canvases[s].SaveAs(args.outdir + sample + "_selection.pdf")
    
    
def plot_cov_output(args):
    
    # Covariance, fractional covariance and correlation

    covfile = TFile(args.covfile)

    covcanvas = TCanvas()

    gStyle = TStyle()
    gStyle.SetPadLeftMargin(0.25); gStyle.SetPadRightMargin(0.15)
    gStyle.SetPalette(56)

    for matname in ('cov', 'fcov', 'corr'):

        mat = covfile.Get(matname)

        mat.GetYaxis().LabelsOption('v')    # doesn't work...
        mat.GetYaxis().SetLabelSize(0.07)
        mat.GetXaxis().SetLabelSize(0.07)

        mat.SetTitleSize(0.3, 't')  # doesn't work...

        if matname == 'corr': 
            mat.GetZaxis().SetRangeUser(-0.4, 1)
            mat.SetTitle("Flux Correlation Matrix")

        mat.Draw("colz")
        mat.SetStats(False)
        covcanvas.SetLeftMargin(0.25)
        covcanvas.Update()
        covcanvas.SaveAs(args.outdir + matname + "_plot.pdf")


def plot_chi2_output(args):

    # Chi squareds
    
    chi2file = TFile(args.chifile)
    
    gStyle = TStyle()
    
    chi2 = chi2file.Get('chisq')
    
    chi2canvas = TCanvas()
    
    chi2.SetTitle('#chi^{2}; log_{10}(sin^{2}(2#theta)); log_{10}(#Delta m^{2}); #chi^{2}');
    gStyle.SetPalette(1)
    chi2.Draw('surf1')
    chi2canvas.SaveAs(args.outdir + "chisq.pdf")
    
    
    # Contours
    contcanvas = TCanvas('cont_canvas', '', 1020, 990)
    gStyle.SetPadLeftMargin(0.15); gStyle.SetPadRightMargin(0.15)
    
    colours = [30, 38, 46]
    contours = [chi2file.Get('90pct'), 
                chi2file.Get('3sigma'), 
                chi2file.Get('5sigma')]
    
    for g in range(len(contours)):
        
        contours[g].SetMarkerStyle(20)
        contours[g].SetMarkerSize(0.25)
        contours[g].SetMarkerColor(colours[g])
        contours[g].SetLineColor(colours[g])
    
    gr_range = TGraph()
    gr_range.SetPoint(0, 0.001, 0.01)
    gr_range.SetPoint(1, 1, 100)
    gr_range.SetMarkerColor(0)
    
    bestfit = TGraph()
    bestfit.SetPoint(0, 0.062, 1.7)
    bestfit.SetMarkerStyle(29)
    bestfit.SetMarkerSize(1.6)
    bestfit.SetMarkerColor(40)
    
    gr_range.SetTitle('SBN Sensitivity; sin^{2}(2#theta); #Delta m^{2} (eV^{2})')
    
    legend = TLegend()
    legend.AddEntry(contours[0], '90% CL', 'l')
    legend.AddEntry(contours[1], '3#sigma CL', 'l')
    legend.AddEntry(contours[2], '5#sigma CL', 'l')
    legend.AddEntry(bestfit, 'Best Fit Point', 'p')
    
    contcanvas.SetLogy()
    contcanvas.SetLogx()
    
    gr_range.Draw('AP')
    gr_range.GetXaxis().SetRangeUser(0.001, 1)
    gr_range.GetYaxis().SetRangeUser(0.01, 100)
    
    contours[0].Draw('P same')
    contours[1].Draw('P same')
    contours[2].Draw('P same')
    
    legend.Draw()
    bestfit.Draw('P same')
    
    contcanvas.SetLeftMargin(0.15)
    contcanvas.Update()
    contcanvas.SaveAs(args.outdir+'Sensitivity.pdf')
    
    
def compare_w_proposal(args):
    
    chi2file = TFile(args.chifile)
    
    gStyle = TStyle()
    gStyle.SetPadLeftMargin(0.15); gStyle.SetPadRightMargin(0.15)
    
    colours = [30, 38, 46]
    contours = [chi2file.Get('90pct'), 
                chi2file.Get('3sigma'), 
                chi2file.Get('5sigma')]
    
    propcontours = []
    contournames = ['90pct', '3s', '5s']
    contourtitles = ['90% Confidence Level', '3#sigma Confidence Level', '5#sigma Confidence Level']
    
    gr_range = TGraph()
    gr_range.SetPoint(0, 0.001, 0.01)
    gr_range.SetPoint(1, 1, 100)
    gr_range.SetMarkerColor(0)
    
    bestfit = TGraph()
    bestfit.SetPoint(0, 0.062, 1.7)
    bestfit.SetMarkerStyle(29)
    bestfit.SetMarkerSize(1.6)
    bestfit.SetMarkerColor(40)
    
    print("contours has length " + str(len(contours)))
    
    for i in range(len(contours)):
        
        x = []
        y = []
        
        with open(args.compdir+'numu'+contournames[i]+'.txt') as f:
            for l, line in enumerate(f):
                if l == 0: continue
                x.append(float(line.split(', ')[0]))
                y.append(float(line.split(', ')[1].replace("\n", "")))
        
        propcontours.append(TGraph())
        for j in range(len(x)):
            propcontours[i].SetPoint(j, x[j], y[j])

        tempcanvas = TCanvas('temp_canvas', '', 1020, 990)

        templegend = TLegend()
        templegend.AddEntry(contours[i], 'Our contour', 'l')
        templegend.AddEntry(propcontours[i], 'From proposal', 'l')
        templegend.AddEntry(bestfit, 'Best Fit Point', 'p')

        tempcanvas.SetLogy()
        tempcanvas.SetLogx()

        gr_range.SetTitle(contourtitles[i]+' Comparison; sin^{2}(2#theta); #Delta m^{2} (eV^{2})')

        gr_range.Draw('AP')
        gr_range.GetXaxis().SetRangeUser(0.001, 1)
        gr_range.GetYaxis().SetRangeUser(0.01, 100)
        
        for lst in (contours, propcontours):
            lst[i].SetMarkerStyle(20)
            lst[i].SetMarkerSize(0.25)
            lst[i].SetMarkerColor(colours[i] if lst == contours else 1)
            lst[i].SetLineColor(colours[i] if lst == contours else 1)
        
        contours[i].Draw('P same')
        propcontours[i].Draw('P same')

        templegend.Draw()
        bestfit.Draw('P same')

        tempcanvas.SaveAs(args.outdir+contournames[i]+'_comparison.pdf')
    
    
def dm2_chi2_slice(args):
    
    dm2s = [float(dm2) for dm2 in args.dm2list.split(",")]
    
    chi2file = TFile(args.chifile)
    chi2 = chi2file.Get('chisq')
    
    min_sin = chi2.GetXmin(); max_sin = chi2.GetXmax()
    min_dm2 = chi2.GetYmin(); max_dm2 = chi2.GetYmax()
    
    chi2_vals = [float(chi) for chi in chi2.GetZ()]
    NP = int(math.sqrt(len(chi2_vals)))
    
    reusable_canvas = TCanvas()
    print(dm2s)
    print("min_sin = " + str(min_sin) + " and max_sin = " + str(max_sin))
    print("min_dm2 = " + str(min_dm2) + " and max_dm2 + " + str(max_dm2))
    print("NP + " + str(NP))
    slices = []
    for m, dm2 in enumerate(dm2s):
        
        jval = int((math.log10(dm2) - min_dm2)*(NP-1)/(max_dm2-min_dm2))
        tempslice = [chi2_vals[i*NP + jval] for i in range(NP)]
        
        slices.append(TGraph())
        for c, chi in enumerate(tempslice):
            slices[m].SetPoint(c, (min_sin + c*(max_sin - min_sin)/(NP-1)), chi)
        
        if args.chilims:
            chilims = [float(lim) for lim in args.chilims.split(",")]
            if len(chilims) == 1:
                slices[m].GetYaxis().SetRangeUser(0, chilims[0])
            else:
                slices[m].GetYaxis().SetRangeUser(chilims[0], chilims[1])
        
        slices[m].SetTitle("#chi^{2} @ #Delta m^{2} = " + str(dm2) + "; log_{10}(sin^{2}(2#theta)); #chi^{2}")
        slices[m].Draw("AP")
        
        reusable_canvas.SaveAs(args.outdir + "dm2_"+str(dm2).replace(".", "-")+"_slice.pdf")
    
def sin_chi2_slice(args):
    
    sins = [float(sin) for sin in args.sinlist.split(",")]
    
    chi2file = TFile(args.chifile)
    chi2 = chi2file.Get('chisq')
    
    min_sin = chi2.GetXmin(); max_sin = chi2.GetXmax()
    min_dm2 = chi2.GetYmin(); max_dm2 = chi2.GetYmax()
    
    chi2_vals = [float(chi) for chi in chi2.GetZ()]
    NP = int(math.sqrt(len(chi2_vals)))
    
    reusable_canvas = TCanvas()
    
    slices = []
    for s, sin in enumerate(sins):
        
        ival = int((math.log10(sin) - min_sin)*(NP-1)/(max_sin-min_sin))
        tempslice = [chi2_vals[ival*NP + j] for j in range(NP)]
        
        slices.append(TGraph())
        for c, chi in enumerate(tempslice):
            slices[s].SetPoint(c, 10**(min_dm2 + c*(max_dm2 - min_dm2)/(NP-1)), chi)
        
        if args.chilims:
            chilims = [float(lim) for lim in args.chilims.split(",")]
            if len(chilims) == 1:
                slices[s].GetYaxis().SetRangeUser(0, chilims[0])
            else:
                slices[s].GetYaxis().SetRangeUser(chilims[0], chilims[1])
        
        slices[s].SetTitle("#chi^{2} @ #sin^{2}(2#theta) = " + str(sin) + "; log_{10}(#Delta m^{2}); #chi^{2}")
        slices[s].Draw("AP")
        reusable_canvas.SaveAs(args.outdir + "sin_"+str(sin).replace(".", "-")+"_slice.pdf")
    

if __name__ == "__main__":
    
    buildpath = os.environ['SBN_LIB_DIR']
    if not buildpath:
        print('ERROR: SBNDDAQ_ANALYSIS_BUILD_PATH not set')
        sys.exit()
    
    ROOT.gROOT.ProcessLine('.L ' + buildpath + '/libsbnanalysis_Event.so')
    ROOT.gROOT.ProcessLine('.L ' + buildpath + '/libsbnanalysis_SBNOsc_classes.so')
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-chi", "--chifile", default = False)
    parser.add_argument("-cov", "--covfile", default = False)
    parser.add_argument("-cts", "--cntfile", default = False)
    parser.add_argument("-o", "--outdir", required = True)
    parser.add_argument("-comp", "--compdir", default = False)
    parser.add_argument("-dm2slice", "--dm2list", default = False)
    parser.add_argument("-sinslice", "--sinlist", default = False)
    parser.add_argument("-chilims", "--chilims", default = False)
    
    args = parser.parse_args()
    
    if args.cntfile: make_count_plot(args)
    if args.covfile: plot_cov_output(args)
    if args.chifile: plot_chi2_output(args)
    if args.compdir: compare_w_proposal(args)
    if args.dm2list: dm2_chi2_slice(args)
    if args.sinlist: sin_chi2_slice(args)
    
    #if parser.parse_args().compare: 
    #    compare_w_proposal(parser_parse_args())
    #    with open('filename') as f:
    #        for line in f:
    #            data = [float(x) for x in line.split(",")]




