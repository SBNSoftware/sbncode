import ROOT
import util
import argparse


def main(args):
    f = ROOT.TFile(args.input)

    sigs = [f.Get(h.lstrip("-")) for h in args.histo_sig]
    util.validate_hists(args.histo_sig, sigs)
    sigs = [util.resize_histo(args, h) for h in sigs]

    bkgs = [f.Get(h.lstrip("-")) for h in args.histo_bkg]
    util.validate_hists(args.histo_bkg, bkgs)
    bkgs = [util.resize_histo(args, h) for h in bkgs]


    n_sig = 0.
    for h in sigs:
      n_sig += h.Integral()
   
    n_bkg = 0.
    for h in bkgs:
      n_bkg += h.Integral()

    print n_sig / n_bkg

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = util.with_io_args(parser)
    parser = util.with_histosize_args(parser)
    parser.add_argument("-hs", "--histo_sig", required=True, type=util.comma_separated)
    parser.add_argument("-hb", "--histo_bkg", required=True, type=util.comma_separated)
    main(parser.parse_args())

