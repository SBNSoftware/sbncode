import ROOT
import util
import argparse


def main(args):
    f = ROOT.TFile(args.input)
    hists = [[h for h in hlist] for hlist in args.hstack]
    util.validate_hists(args.hstack, hists)
    hists = [[util.resize_histo(args, f.Get(h)) for h in hlist] for hlist in hists]

    canvas = ROOT.TCanvas("canvas", "Canvas", 250,100,700,500)
    hstack = ROOT.THStack()

    for i, hlist in enumerate(hists):
        h = hlist[0]
        if len(hlist) > 1:
            for hadd in hlist[1:]:
                h.Add(hadd)
  
        if args.area_normalize and h.Integral() > 1e-4: h.Scale(1. / h.Integral())
        if args.stack:
            h.SetFillColor(util.fillcolors(i))
        else:
            h.SetLineColor(util.fillcolors(i))
            h.SetLineWidth(3)
        if args.names: 
            name = args.names[i]
            if args.nevent_in_legend:
                name += " (%i)" % int(h.Integral()) 
            h.SetTitle(name)
        hstack.Add(h)
        
    drawstr = "HIST" if args.stack else "NOSTACK HIST"
    hstack.Draw(drawstr)
    util.style(args, hstack)
    legend = ROOT.gPad.BuildLegend(*(args.legend_position + [""]))
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
    parser.add_argument("-hs", "--hstack", required=True, type=util.histo_list)
    parser.add_argument("-s", "--stack", action="store_true")
    parser.add_argument("-n", "--names", default=None, type=util.comma_separated)
    parser.add_argument("-a", "--area_normalize", action="store_true")
    parser.add_argument("--nevent_in_legend", action="store_true")
    parser.add_argument("-lp", "--legend_position", default=[0.75,0.75,0.95,0.95], type=util.legend_position)
    main(parser.parse_args())

