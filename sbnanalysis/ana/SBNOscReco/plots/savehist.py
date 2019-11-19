import ROOT
import util
import argparse
from array import array

def main(args):
    util.optstat(args)
    f = ROOT.TFile(args.input)
    hist = f.Get(args.hist)
    util.validate_hists([args.hist], [hist])
    hist = util.resize_histo(args, hist)
    if args.varianceY:
        hist = project_variance(hist)
    elif args.meanY:
        hist = project_mean(hist)
    canvas = ROOT.TCanvas("canvas", "Canvas", 250,100,700,500)
    drawstr = ""
    if args.drawstr is not None:
        drawstr = args.drawstr
    elif isinstance(hist, ROOT.TH2D):
        drawstr = "COLZ"
    elif isinstance(hist, ROOT.TH1D):
        drawstr = "HIST"
    elif isinstance(hist, ROOT.TGraph):
        drawstr = "AL"
    hist.Draw(drawstr)
    if args.title is not None:
        hist.SetTitle(args.title)
    util.style(args, hist)
    if args.logy:
        canvas.SetLogy()
    if args.logz:
        canvas.SetLogz()
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
            this_var_sum += ((hist.GetYaxis().GetBinCenter(j) - this_mean) * hist.GetBinContent(i, j))**2
        xs.append(hist.GetXaxis().GetBinCenter(i))
        this_var = this_var_sum / n_entries if n_entries > 1e-3 else this_var_sum
        var.append(this_var)
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = util.with_io_args(parser)
    parser = util.with_histosize_args(parser)
    parser = util.with_histostyle_args(parser)
    parser.add_argument("-t" ,"--title", default=None)
    parser.add_argument("-hn", "--hist", required=True)
    parser.add_argument("-ly", "--logy", action="store_true")
    parser.add_argument("-lz", "--logz", action="store_true")
    parser.add_argument("-varY", "--varianceY", action="store_true")
    parser.add_argument("-meanY", "--meanY", action="store_true")
    parser.add_argument("-ds", "--drawstr", default=None)
    main(parser.parse_args())

