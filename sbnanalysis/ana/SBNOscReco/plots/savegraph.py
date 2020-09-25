import ROOT
import util
import argparse
from array import array

def main(args):
    f = ROOT.TFile(args.input)
    graph = f.Get(args.graph)
    canvas = ROOT.TCanvas("canvas", "Canvas", 250,100,700,500)
    util.style(args, graph)
    graph.Draw("AC")
    util.resize_graph(args, graph)
    graph.Draw("AC")

    line = ROOT.TGraph(2, array('d', [0., 1.]), array('d', [1., 0.]))
    line.SetLineColor(ROOT.kRed)
    line.Draw("SAME")

    canvas.Update()

    util.wait(args)
    util.write(args, canvas)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = util.with_io_args(parser)
    parser = util.with_histostyle_args(parser)
    parser = util.with_graphsize_args(parser)
    parser.add_argument("-g", "--graph", required=True)
    main(parser.parse_args())

