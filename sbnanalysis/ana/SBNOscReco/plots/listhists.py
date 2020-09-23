import ROOT
import util
import argparse


def main(args):
    f = ROOT.TFile(args.input)
    f.ls()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = util.with_input_args(parser)
    main(parser.parse_args())

