import ROOT
import util
import argparse

def main(args):
    f = ROOT.TFile(args.input)

    base_histos = [f.Get(r) for r in args.base_histos]
    util.validate_hists(args.base_histos, base_histos)
    base_histos = [util.resize_histo(args, f.Get(r)) for r in args.base_histos]

    cut_histos = [f.Get(r) for r in args.cut_histos]
    util.validate_hists(args.cut_histos, cut_histos)
    cut_histos = [util.resize_histo(args, f.Get(r)) for r in args.cut_histos]

    canvas = ROOT.TCanvas("canvas", "Canvas", 250,100,700,500)

    all_histos = []
    for i, (base_histo,cut_histo) in enumerate(zip(base_histos,cut_histos)):
        ratio = cut_histo.Clone()
        ratio.Divide(base_histo)
        ratio.SetLineColor(util.fillcolors(i))
        ratio.SetLineWidth(2)
        if i == 0:
            util.style(args, canvas, ratio)
            ratio.Draw("L")
        else:
            ratio.Draw("L SAME")
        all_histos.append(ratio)

        if args.names: 
            name = args.names[i]
            ratio.SetTitle(name)

    legend = ROOT.gPad.BuildLegend(*(args.legend_position + [""]))
    if args.title is not None:
        all_histos[0].SetTitle(args.title)
    canvas.Update()

    util.wait(args)
    util.write(args, canvas)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = util.with_io_args(parser)
    parser = util.with_histosize_args(parser)
    parser = util.with_histostyle_args(parser)
    parser.add_argument("-b", "--base_histos", type=util.comma_separated, required=True)
    parser.add_argument("-c", "--cut_histos", type=util.comma_separated, required=True)
    parser.add_argument("-t" ,"--title", default=None)
    parser.add_argument("-n", "--names", default=None, type=util.comma_separated)
    parser.add_argument("-lp", "--legend_position", default=[0.75,0.75,0.95,0.95], type=util.legend_position)
    main(parser.parse_args())

