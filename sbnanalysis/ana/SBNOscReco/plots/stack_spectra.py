import ROOT
import util
import argparse


def main(args):
    f = ROOT.TFile(args.input)
    hists = [util.resize_histo(args, f.Get(h)) for h in args.hstack]
    canvas = ROOT.TCanvas("canvas", "Canvas", 250,100,700,500)
    hstack = ROOT.THStack()
    for i, h in enumerate(hists):
        h.SetFillColor(util.fillcolors(i))
        if args.names: h.SetTitle(args.names[i])
        hstack.Add(h)
        
    hstack.Draw("HIST")
    util.style(args, hstack)
    legend = ROOT.gPad.BuildLegend(0.75,0.75,0.95,0.95,"")
    if args.logy:
        canvas.SetLogy()
    canvas.Update()

    util.wait(args)
    util.write(args, canvas)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = util.with_io_args(parser)
    parser = util.with_histosize_args(parser)
    parser = util.with_histostyle_args(parser)
    parser.add_argument("-ly", "--logy", action="store_true")
    parser.add_argument("-hs", "--hstack", required=True, type=util.comma_separated)
    parser.add_argument("-n", "--names", default=None, type=util.comma_separated)
    main(parser.parse_args())

