from ROOT import *
import math
import sys
from array import array

gROOT .SetBatch()
gStyle.SetLineWidth(3)
TGaxis.SetMaxDigits(3)

can = TCanvas("can", "can", 800, 800)
can .cd()


def pretty_hierachy_name(hie):
    if hie == "ih": return "IH"
    elif hie == "nh": return "NH"
    return hie


def pretty_dcp_name(dcp):
    new_name = dcp.replace("over", "/")
    new_name = new_name.replace("pi", "#pi")
    return new_name


def get_round_number(hist):
    nevents = hist.Integral()

    if nevents>1e6:
        return str(round(nevents/1e6))+"M"
    if nevents>1e4:
        return str(round(nevents/1e3))+"k"
    else:
        return str(int(round(nevents)))


def get_hist_list(inputFileName, inHistNameList, project=""):

    histList  = []
    inputFile = TFile(inputFileName)
    
    ## First get the plots
    for histName in inHistNameList:
        hist = inputFile.Get(histName)
        hist .SetDirectory(0)
        histList.append(hist)
    inputFile.Close()

    ## Need to deal with projections here
    if project == "x":
        for hist in range(len(histList)):
            histList[hist] = histList[hist].ProjectionX()
    if project == "y":
        for hist in range(len(histList)):
            histList[hist] = histList[hist].ProjectionY()

    return histList


def comparison_plots(outFileName, inputFileName, inHistNameList, nameList, colzList, styleList, project="", show_rate=False):

    histList  = get_hist_list(inputFileName, inHistNameList, project)
    hist0 = histList[0]

    ## Sort out the maximum
    maxVal = hist0.GetMaximum()
    for hist in histList:
        if hist.GetMaximum() > maxVal:
            maxVal = hist.GetMaximum()
    maxVal *= 1.05

    ## Prettify the histogram
    hist0 .GetXaxis().SetTitleOffset(0.95)
    hist0 .GetYaxis().SetTitleOffset(1)
    hist0 .GetYaxis().SetNdivisions(504)
    hist0 .GetXaxis().SetNdivisions(504)
    
    hist0  .Draw("HIST")
    hist0  .SetMaximum(maxVal)
    hist0  .SetMinimum(0)
    ## This needs to be smarter, and state the bin width
    hist0  .GetYaxis().SetTitle("Events")
        
    ## Now the constraints
    for hist in range(len(histList)):
        histList[hist] .SetLineWidth(4)
        histList[hist] .SetMarkerColor(0)
        histList[hist] .SetLineColor(colzList[hist])
        histList[hist] .SetLineStyle(styleList[hist])
        histList[hist] .Draw("HIST SAME")
    hist0  .Draw("HIST SAME")
    hist0  .Draw("AXIS SAME")

    ## Now a legend
    dim = [0.55, 0.55, 0.92, 0.90]
    leg = TLegend(dim[0], dim[1], dim[2], dim[3], "", "NDC")
    leg .SetShadowColor(0)
    leg .SetFillColor(0)
    leg .SetLineWidth(0)
    leg .SetTextSize(0.042)
    leg .SetLineColor(kWhite)
    if show_rate:
        for hist in range(len(histList)): leg .AddEntry(histList[hist], nameList[hist]+" ("+get_round_number(histList[hist])+")", "l")
    else:
        for hist in range(len(histList)): leg .AddEntry(histList[hist], nameList[hist], "l")
    leg .Draw("SAME")

    gPad .SetRightMargin(0.03)
    gPad .SetLeftMargin(0.15)
    gPad .SetTopMargin(0.06)
    gPad .SetBottomMargin(0.15)
    gPad .Update()
    can  .Update()
    can.SaveAs("plots/"+outFileName)


def comparison_plots_ratio(outFileName, inputFileName, inHistNameList, nameList, colzList, styleList, project=""):

    ratioList  = get_hist_list(inputFileName, inHistNameList, project)

    # for now, assume that the central value is the nominal
    #ratioList = []
    nomHist   = ratioList[3].Clone()
    for hist in ratioList:
        hist.Divide(nomHist)
        
    hist0 = ratioList[0]

    nomHist .SetLineStyle(7)

    ## Sort out the maximum
    maxVal = hist0.GetBinContent(hist0.GetMaximumBin())
    for hist in ratioList:
        if hist.GetMaximum() > maxVal:
            maxVal = hist.GetBinContent(hist.GetMaximumBin())
    maxVal += 0.01
    if maxVal > 2: maxVal = 2
    
    ## Prettify the histogram
    hist0 .GetXaxis().SetTitleOffset(0.95)
    hist0 .GetYaxis().SetTitleOffset(1)
    hist0 .GetYaxis().SetNdivisions(504)
    hist0 .GetXaxis().SetNdivisions(504)
    
    hist0  .Draw("HIST")
    hist0  .SetMaximum(maxVal)
    hist0  .SetMinimum(1 - (maxVal - 1))
    ## This needs to be smarter, and state the bin width
    hist0  .GetYaxis().SetTitle("Ratio w.r.t. nominal")
        
    ## Now the constraints
    for hist in range(len(ratioList)):
        ratioList[hist] .SetLineWidth(4)
        ratioList[hist] .SetMarkerColor(0)
        ratioList[hist] .SetLineColor(colzList[hist])
        ratioList[hist] .SetLineStyle(styleList[hist])
        ratioList[hist] .Draw("HIST SAME")
    hist0  .Draw("HIST SAME")
    hist0  .Draw("AXIS SAME")

    ## Now a legend
    dim = [0.86, 0.2, 0.99, 0.8]
    leg = TLegend(dim[0], dim[1], dim[2], dim[3], "", "NDC")
    leg .SetShadowColor(0)
    leg .SetFillColor(0)
    leg .SetLineWidth(0)
    leg .SetTextSize(0.042)
    leg .SetLineColor(kWhite)
    for hist in range(len(ratioList)): leg .AddEntry(ratioList[hist], nameList[hist], "l")
    leg .Draw("SAME")

    gPad .SetRightMargin(0.15)
    gPad .SetLeftMargin(0.15)
    gPad .SetTopMargin(0.06)
    gPad .SetBottomMargin(0.15)
    gPad .Update()
    can  .Update()
    can.SaveAs("plots/"+outFileName)

    
def sample_plot_2D(outFileName, inputFileName, inHistName):

    inputFile = TFile(inputFileName)
    hist = inputFile.Get(inHistName)
    hist .SetDirectory(0)

    ## Prettify the histogram
    hist .GetXaxis().SetTitleOffset(0.95)
    hist .GetYaxis().SetTitleOffset(1)
    hist .GetYaxis().SetNdivisions(504)
    hist .GetXaxis().SetNdivisions(504)
    hist .GetZaxis().SetNdivisions(505)
    hist .Draw("COLZ")

    ## This needs to be smarter, and state the bin width
    hist .GetZaxis().SetTitle("Events")
    hist .GetZaxis().RotateTitle(1)
    hist .Draw("AXIS SAME")

    gPad .SetRightMargin(0.20)
    gPad .SetLeftMargin(0.15)
    gPad .SetTopMargin(0.06)
    gPad .SetBottomMargin(0.15)
    gPad .Update()
    can  .Update()
    can.SaveAs("plots/"+outFileName)


def make_basic_sample_plots(inputFileName, hie="nh", dcp="piover2"):
    
    ## FD NUE
    for mode in ["FHC", "RHC"]:
        inHistNameList = ["FD_"+mode+"_Nue_"+hie+"_"+dcp+"_total",
                          "FD_"+mode+"_Nue_"+hie+"_"+dcp+"_AllSignalNue",                      
                          "FD_"+mode+"_Nue_"+hie+"_"+dcp+"_AllBeamNue",
                          "FD_"+mode+"_Nue_"+hie+"_"+dcp+"_AllNumu",
                          "FD_"+mode+"_Nue_"+hie+"_"+dcp+"_NC"]
        nameList       = ["Total",
                          "Sig. (#nu_{e}+#bar{#nu}_{e}) CC",
                          "Beam (#nu_{e}+#bar{#nu}_{e}) CC",
                          "(#nu_{#mu}+#bar{#nu}_{#mu}) CC",
                          "NC"]
        colzList       = [kBlack, kGray+1, kBlue, kGreen+1, kRed]
        styleList      = [1, 1, 1, 1, 1]
        comparison_plots("FD_"+mode+"_nue_"+hie+"_"+dcp+"_app.png", inputFileName, inHistNameList, nameList, colzList, styleList)
        
        ## FD numu
        inHistNameList = ["FD_"+mode+"_Numu_"+hie+"_"+dcp+"_total",
                          "FD_"+mode+"_Numu_"+hie+"_"+dcp+"_Numu",
                          "FD_"+mode+"_Numu_"+hie+"_"+dcp+"_Numubar",
                          "FD_"+mode+"_Numu_"+hie+"_"+dcp+"_AllBeamNue",
                          "FD_"+mode+"_Numu_"+hie+"_"+dcp+"_NC"]
        nameList       = ["Total",
                          "#nu_{#mu} CC",
                          "#bar{#nu}_{#mu} CC",
                          "(#nu_{e}+#bar{#nu}_{e}) CC",
                          "NC"]
        colzList       = [kBlack, kGray+1, kBlue, kGreen+1, kRed]
        styleList      = [1, 1, 1, 1, 1]
        comparison_plots("FD_"+mode+"_numu_"+hie+"_"+dcp+"_dis.png", inputFileName, inHistNameList, nameList, colzList, styleList)
        
        ## ND
        inHistNameList = ["ND_"+mode+"_total",
                          "ND_"+mode+"_Numu",
                          "ND_"+mode+"_Numubar",
                          "ND_"+mode+"_AllBeamNue",                     
                          "ND_"+mode+"_NC"]
        nameList       = ["Total",
                          "#nu_{#mu} CC",
                          "#bar{#nu}_{#mu} CC",
                          "(#nu_{e}+#bar{#nu}_{e}) CC",
                          "NC"]
        colzList       = [kBlack, kGray+1, kBlue, kGreen+1, kRed]
        styleList      = [1, 1, 1, 1, 1]
        comparison_plots("ND_"+mode+"_Ereco.png", inputFileName, inHistNameList, nameList, colzList, styleList, "x")
        comparison_plots("ND_"+mode+"_yreco.png", inputFileName, inHistNameList, nameList, colzList, styleList, "y")

        ## Need to add 2D ND plots
        sample_plot_2D("ND_"+mode+"_2D.png", inputFileName, "ND_"+mode+"_total")
        
def make_osc_comp_plots(baseName, inputFileName):


    inHistNameList = []
    nameList       = []
    
    for hie in ["nh", "ih"]:
        for dcp in ["0pi", "piover2", "3piover2"]:
            inHistNameList .append(baseName+"_"+hie+"_"+dcp+"_total")
            nameList       .append(pretty_hierachy_name(hie)+" "+pretty_dcp_name(dcp))
    colzList       = [kBlack, kRed, kBlue, kBlack, kRed, kBlue]
    styleList      = [1, 1, 1, 7, 7, 7]
    comparison_plots(baseName+"_oscvar.png", inputFileName, inHistNameList, nameList, colzList, styleList, "", True)

def get_total_rates(inputFileName, hie="nh", dcp="piover2"):

    inHistNameList = ["FD_FHC_Nue_"+hie+"_"+dcp+"_total",
                      "FD_FHC_Numu_"+hie+"_"+dcp+"_total",
                      "FD_RHC_Nue_"+hie+"_"+dcp+"_total",
                      "FD_RHC_Numu_"+hie+"_"+dcp+"_total",
                      "ND_FHC_total",
                      "ND_RHC_total"]
    
    print "Getting rates for hierarchy =", hie, "and dCP =", dcp
    ## First get the plots
    inputFile = TFile(inputFileName)
    
    for histName in inHistNameList:
        hist = inputFile.Get(histName)

        ## Format the rate
        rate = "{:.0f}".format(hist.Integral())

        ## Special case for ND
        if "ND_" in histName:
            rate = "{:.2f}M".format(hist.Integral()/1e6)
        print " -->", histName, "has", rate, "events"

def make_validation_plots(inFileName, systList):

    colzList  = [kRed, kBlue, kGray+1, kBlack, kGray+1, kBlue, kRed]
    styleList = [1, 1, 1, 1, 7, 7, 7]
    nameList  = ["3#sigma", "2#sigma", "1#sigma", "0", "-1#sigma", "-2#sigma", "-3#sigma"]

    for syst in systList:

        for sample in ["FD_FHC_Nue_total", "FD_FHC_Numu_total",
                       "FD_RHC_Nue_total", "FD_RHC_Numu_total",
                       "ND_FHC_1D_total", "ND_RHC_1D_total"]:

            inHistNameList = [sample+"_"+syst+"_"+str(x) for x in range(3, -4, -1)]
            comparison_plots(sample+"_"+syst+"_validation.png", inFileName, inHistNameList, nameList, colzList, styleList)
            comparison_plots_ratio(sample+"_"+syst+"_validation_ratio.png", inFileName, inHistNameList, nameList, colzList, styleList)

        for sample in ["ND_FHC_total", "ND_RHC_total"]:
            
            inHistNameList = [sample+"_"+syst+"_"+str(x) for x in range(3, -4, -1)]
            comparison_plots(sample+"_x_"+syst+"_validation.png", inFileName, inHistNameList, nameList, colzList, styleList, "x")
            comparison_plots(sample+"_y_"+syst+"_validation.png", inFileName, inHistNameList, nameList, colzList, styleList, "y")
            comparison_plots_ratio(sample+"_x_"+syst+"_validation_ratio.png", inFileName, inHistNameList, nameList, colzList, styleList, "x")
            comparison_plots_ratio(sample+"_y_"+syst+"_validation_ratio.png", inFileName, inHistNameList, nameList, colzList, styleList, "y")
 
    
if __name__ == '__main__':

    ## Note that this script expects to be run in a directory where there's a plots folder to dump everything into
    systList  = ["flux0", "flux1", "flux2", "flux3", "flux4", "flux5", "flux6", "flux7",
                 "flux8", "flux9", "eScaleMuLAr", "eScaleMuND", "eScaleE", "ChargedHadCorr",
                 "ChargedHadUncorrFD", "ChargedHadUncorrND", "eScaleN_ND", "eScaleN_FD",
                 "eScalePi0Corr", "Pi0UncorrFD", "Pi0UncorrND", "MaCCQE", "MaNCEL", "EtaNCEL", "MaCCRES",
                 "MvCCRES", "RDecBR1gamma", "RDecBR1eta", "E2p2h_A_nu",
                 "E2p2h_B_nu", "E2p2h_A_nubar", "E2p2h_B_nubar", "NR_nu_n_CC_2Pi",
                 "NR_nu_n_CC_3Pi", "NR_nu_p_CC_2Pi", "NR_nu_p_CC_3Pi", "NR_nu_np_CC_1Pi",
	         "NR_nu_n_NC_1Pi", "NR_nu_n_NC_2Pi", "NR_nu_n_NC_3Pi", "NR_nu_p_NC_1Pi",
                 "NR_nu_p_NC_2Pi", "NR_nu_p_NC_3Pi", "NR_nubar_n_CC_1Pi", "NR_nubar_n_CC_2Pi",
                 "NR_nubar_n_CC_3Pi", "NR_nubar_p_CC_1Pi", "NR_nubar_p_CC_2Pi", "NR_nubar_p_CC_3Pi",
                 "NR_nubar_n_NC_1Pi", "NR_nubar_n_NC_2Pi", "NR_nubar_n_NC_3Pi", "NR_nubar_p_NC_1Pi",
                 "NR_nubar_p_NC_2Pi", "NR_nubar_p_NC_3Pi", "BeRPA_A", "BeRPA_B", "BeRPA_D", "BeRPA_E",
                 "SPPLowQ2Suppression"]
    
    ## make_validation_plots("../output/validation_hists_mcc11.root", systList)
    ## make_basic_sample_plots("../output/spec_hist_mcc11.root")
    ## get_total_rates("../output/spec_hist_mcc11.root")
    ## make_osc_comp_plots("FD_FHC_Nue", "../output/spec_hist_mcc11.root")
    ## make_osc_comp_plots("FD_FHC_Numu", "../output/spec_hist_mcc11.root")
    ## make_osc_comp_plots("FD_RHC_Nue", "../output/spec_hist_mcc11.root")
    ## make_osc_comp_plots("FD_RHC_Numu", "../output/spec_hist_mcc11.root")

    
