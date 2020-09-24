import ROOT
import util
import argparse
from array import array
import math

def main(args):
    util.optstat(args)
    hist = util.get_tobject(args, args.hist)
    util.validate_hists([args.hist], [hist])
    hist = util.resize_histo(args, hist)

    if args.make_percentages:
        if isinstance(hist, ROOT.TH2D):
            for i in range(1, hist.GetNbinsX()+1):
                this_norm = hist.Integral(i, i, 1, hist.GetNbinsY())
                for j in range(1, hist.GetNbinsY()+1):
                    if this_norm < 1e-3:
                        hist.SetBinContent(i, j, 0)
                    else:
                        hist.SetBinContent(i, j, hist.GetBinContent(i, j) / this_norm)
        else:
            for k in range(1, hist.GetNbinsZ()+1):
                for i in range(1, hist.GetNbinsX()+1):
                    this_norm = hist.Integral(i, i, 1, hist.GetNbinsY(), k, k)
                    for j in range(1, hist.GetNbinsY()+1):
                            if this_norm < 1e-3:
                                hist.SetBinContent(i, j, k, 0)
                            else:
                                hist.SetBinContent(i, j, k, hist.GetBinContent(i, j, k) / this_norm)

    hist2 = None
    if args.meanY:
        hist2 = mean_and_err(hist)
    
    canvas = ROOT.TCanvas("canvas", "Canvas", 250,100,700,500)
    drawstr = ""
    if args.drawstr is not None:
        drawstr = args.drawstr
    elif isinstance(hist, ROOT.TH2D):
        drawstr = "COLZ"
        if args.draw_text:
            drawstr += " TEXT"
            hist.SetMarkerSize(3)
            ROOT.gStyle.SetPaintTextFormat("1.3f")
    elif isinstance(hist, ROOT.TH1D):
        drawstr = "HIST"
    elif isinstance(hist, ROOT.TGraph):
        drawstr = "AL"
    hist.Draw(drawstr)
    if hist2:
        hist2.SetLineColor(ROOT.kRed)
        hist2.SetLineWidth(3)
        hist2.SetMarkerSize(1)
        hist2.SetMarkerStyle(21)
        hist2.SetMarkerColor(ROOT.kRed)
        hist2.Draw("P")
    if args.title is not None:
        hist.SetTitle(args.title)
    util.style(args, canvas, hist)

    if args.logy:
        canvas.SetLogy()
    if args.logz:
        canvas.SetLogz()
    box = util.draw_text(args)

    canvas.Update()

    util.wait(args)
    util.write(args, canvas)

def project_variance(hist): 
    # calculate the variance in each y-axis -- this is the y-value 
    var = []
    xs = []
    for i in range(1, hist.GetNbinsX()+1):
        this_sum = 0
        n_entries = 0
        for j in range(1, hist.GetNbinsY()+1):
            this_sum += hist.GetYaxis().GetBinCenter(j) * hist.GetBinContent(i, j)
            n_entries += hist.GetBinContent(i, j)
        if n_entries < 1e-3: continue
        this_mean = this_sum / n_entries if n_entries > 1e-3 else this_sum
        this_var_sum = 0
        for j in range(1, hist.GetNbinsY()+1):
            this_var_sum += hist.GetBinContent(i, j) * (hist.GetYaxis().GetBinCenter(j) - this_mean)**2
        xs.append(hist.GetXaxis().GetBinCenter(i))
        this_var = this_var_sum / n_entries if n_entries > 1e-3 else this_var_sum
        var.append(math.sqrt(this_var))
    return ROOT.TGraph(len(xs), array('d', xs), array('d', var))

def project_mean(hist):
    # calculate the variance in each y-axis -- this is the y-value 
    mean = []
    xs = []
    for i in range(1, hist.GetNbinsX()+1):
        this_sum = 0
        n_entries = 0
        for j in range(1, hist.GetNbinsY()+1):
            this_sum += hist.GetYaxis().GetBinCenter(j) * hist.GetBinContent(i, j)
            n_entries += hist.GetBinContent(i, j)
        if n_entries < 1e-3: continue
        this_mean = this_sum / n_entries if n_entries > 1e-3 else this_sum
        mean.append(this_mean)
        xs.append(hist.GetXaxis().GetBinCenter(i))
    return ROOT.TGraph(len(xs), array('d', xs), array('d', mean))

def mean_and_err(hist):
    hist = hist.Clone()
    hist.RebinX(4)
    g_mean = project_mean(hist)
    g_err = project_variance(hist)
    x_err = array('d', [0 for i in range(len(g_mean.GetX()))])
    return ROOT.TGraphErrors(len(g_mean.GetX()), g_mean.GetX(), g_mean.GetY(), x_err, g_err.GetY())

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = util.with_io_args(parser)
    parser = util.with_histosize_args(parser)
    parser = util.with_histostyle_args(parser)
    parser = util.with_text_args(parser)
    parser.add_argument("-t" ,"--title", default=None)
    parser.add_argument("-hn", "--hist", required=True)
    parser.add_argument("-ly", "--logy", action="store_true")
    parser.add_argument("-lz", "--logz", action="store_true")
    parser.add_argument("-meanY", "--meanY", action="store_true")
    parser.add_argument("-mp", "--make_percentages", action="store_true")
    parser.add_argument("-dt", "--draw_text", action="store_true")
    parser.add_argument("-ds", "--drawstr", default=None)
    main(parser.parse_args())

