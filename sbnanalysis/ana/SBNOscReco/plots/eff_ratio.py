import ROOT
import util
import argparse

def main(args):
    f = ROOT.TFile(args.input)
    base_hist = util.resize_histo(args, f.Get(args.base_histo))
    cut_histos = [util.resize_histo(args, f.Get(r)) for r in args.cut_histos]
    canvas = ROOT.TCanvas("canvas", "Canvas", 250,100,700,500)

    all_histos = []
    for i, cut_histo in enumerate(cut_histos):
        ratio = cut_histo.Clone()
        ratio.Divide(base_hist)
        ratio.SetLineColor(util.colors(i))
        ratio.SetLineWidth(2)
        if i == 0:
            util.style(args, ratio)
            ratio.Draw("L")
        else:
            ratio.Draw("L SAME")
        all_histos.append(ratio)

    canvas.Update()

    util.wait(args)
    util.write(args, canvas)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = util.with_io_args(parser)
    parser = util.with_histosize_args(parser)
    parser = util.with_histostyle_args(parser)
    parser.add_argument("-b", "--base_histo", required=True)
    parser.add_argument("-c", "--cut_histos", type=util.comma_separated, required=True)
    main(parser.parse_args())

