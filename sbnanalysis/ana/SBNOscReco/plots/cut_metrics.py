import ROOT
import util
import argparse
from array import array

def get_histos(f, variable, cut, names):
    return [f.Get(variable + "_" + n + cut) for n in names]

def get_integral(f, variable, cut, names):
    return sum([h.Integral() for h in get_histos(f, variable, cut, names)])

def main(args):
    f = ROOT.TFile(args.input)
    effs = []
    puritys = []

    signal_base = get_integral(f, args.variable, args.base_cut, args.signal_names)
    background_base = get_integral(f, args.variable, args.base_cut, args.background_names)

    effs.append(1.)
    puritys.append(signal_base / (signal_base + background_base))
    
    for cut in args.cuts:
        this_signal = get_integral(f, args.variable, cut, args.signal_names) 
        this_background = get_integral(f, args.variable, cut, args.background_names) 
        this_eff = this_signal / signal_base
        this_purity = this_signal / (this_signal + this_background) 
        effs.append(this_eff)
        puritys.append(this_purity)

    effs = array('d', effs)
    puritys = array('d', puritys)
    xvals = array('d', range(len(effs)))

    g_eff = ROOT.TGraph(len(effs), xvals, effs)
    g_eff.SetName("Efficiency")
    g_eff.SetTitle("Efficiency")

    g_purity = ROOT.TGraph(len(puritys), xvals, puritys)
    g_purity.SetName("Purity")
    g_purity.SetTitle("Purity")

    canvas = ROOT.TCanvas("canvas", "Canvas", 250,100,700,500)

    g_eff.Draw("AL")
    g_eff.SetLineColor(ROOT.kGreen)
    util.style(args, g_eff)
    for i in range(len(effs)):
        g_eff.GetXaxis().SetBinLabel(g_eff.GetXaxis().FindBin(i), args.names[i])
    g_purity.Draw("SAME")
    g_purity.SetLineColor(ROOT.kRed)

    legend = ROOT.gPad.BuildLegend(0.75,0.75,0.95,0.95,"")
    canvas.Update()

    util.wait(args)
    util.write(args, canvas)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = util.with_io_args(parser)
    parser = util.with_histostyle_args(parser)
    parser.add_argument("-s", "--signal_names", type=util.comma_separated, required=True)
    parser.add_argument("-b", "--background_names", type=util.comma_separated, required=True)
    parser.add_argument("-bc", "--base_cut", required=True)
    parser.add_argument("-c", "--cuts", type=util.comma_separated, required=True)
    parser.add_argument("-n", "--names", type=util.comma_separated, required=True)
    parser.add_argument("-v", "--variable", required=True)
    main(parser.parse_args())

