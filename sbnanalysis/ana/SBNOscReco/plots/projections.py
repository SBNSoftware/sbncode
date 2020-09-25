import ROOT
import util
import argparse
from array import array
import math

# ROOT.gStyle.SetPalette(21)

def main(args):
    util.optstat(args)
    f = ROOT.TFile(args.input)
    hist = f.Get(args.hist)
    util.validate_hists([args.hist], [hist])
    hist = util.resize_histo(args, hist, "")
    hstack = ROOT.THStack()
    canvas = ROOT.TCanvas("canvas", "Canvas", 250,100,700,500)

    class ArgDict(dict):
        __getattr__ = dict.get
    
    hists = []
    for i,(low,high) in enumerate(zip(args.lows, args.highs)):
        this_args = ArgDict({"y_min": low, "y_max": high, "projectionX": True})
        h = util.resize_histo(this_args, hist.Clone(), str(i))
        h.SetLineColor(util.fillcolors(i)) 
        h.SetLineWidth(3)
        if args.names:
            name = args.names[i]
            h.SetTitle(name)
        hists.append(h)
        hstack.Add(h)

    hstack.Draw("NOSTACK HIST")
    legend = ROOT.gPad.BuildLegend(*(args.legend_position + [""]))
    if args.title is not None:
        hstack.SetTitle(args.title)
    util.style(args, canvas, hstack)

    box = util.draw_text(args)

    canvas.Update()

    util.wait(args)
    util.write(args, canvas)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = util.with_io_args(parser)
    parser = util.with_histostyle_args(parser)
    parser = util.with_text_args(parser)
    parser = util.with_histosize_args(parser)
    parser.add_argument("-t" ,"--title", default=None)
    parser.add_argument("-hn", "--hist", required=True)
    parser.add_argument("-ls", "--lows", nargs="+", type=float)
    parser.add_argument("-hs", "--highs", nargs="+", type=float)
    parser.add_argument("-n", "--names", type=util.comma_separated)
    parser.add_argument("-lp", "--legend_position", default=[0.75,0.75,0.95,0.95], type=util.legend_position)
    main(parser.parse_args())

