import ROOT
import util
import argparse


def main(args):
    f = ROOT.TFile(args.input)
    hist = util.resize_histo(args, f.Get(args.hist))
    canvas = ROOT.TCanvas("canvas", "Canvas", 250,100,700,500)
    hist.Draw("HIST")
    util.style(args, hist)
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
    parser.add_argument("-hn", "--hist", required=True)
    parser.add_argument("-ly", "--logy", action="store_true")
    main(parser.parse_args())

