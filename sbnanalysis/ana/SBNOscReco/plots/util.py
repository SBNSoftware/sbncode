import ROOT
import os
from array import array
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gStyle.SetOptStat(0)

def root_env():
    buildpath = os.environ["SBN_LIB_DIR"]
    if not buildpath:
        print "ERROR: SBNDDAQ_ANALYSIS_BUILD_PATH not set"
        sys.exit()
    ROOT.gROOT.ProcessLine(".L " + buildpath + "/libsbnanalysis_Event.so")
    ROOT.gROOT.ProcessLine(".L " + buildpath + "/libsbnanalysis_PandoraTesting_classes.so")
    ROOT.gROOT.ProcessLine(".L " + buildpath + "/libsbnanalysis_SBNOsc_classes.so")
    ROOT.gROOT.ProcessLine(".L " + buildpath + "/libsbnanalysis_SBNOscReco_classes.so")

def wait(args):
    if args.wait:
        raw_input("Press Enter to continue...")

def write(args, canvas):
    if args.output:
        canvas.SaveAs(args.output)

def with_input_args(parser):
    if "INPUTFILE" in os.environ:
        parser.add_argument("-i", "--input", default=os.environ["INPUTFILE"])
    else:
        parser.add_argument("-i", "--input", required=True)
    return parser

def with_display_args(parser):
    parser.add_argument("-w", "--wait", action="store_true")
    parser.add_argument("-o", "--output", default=None)
    return parser

def with_io_args(parser):
    parser = with_input_args(parser)
    parser = with_display_args(parser)
    return parser

def with_graphsize_args(parser):
    parser.add_argument("-xm", "--x_min", type=float, default=None)
    parser.add_argument("-xh", "--x_max", type=float, default=None)

    parser.add_argument("-ym", "--y_min", type=float, default=None)
    parser.add_argument("-yh", "--y_max", type=float, default=None)

    parser.add_argument("-zm", "--z_min", type=float, default=None)
    parser.add_argument("-zh", "--z_max", type=float, default=None)
    return parser

def with_histosize_args(parser):
    parser.add_argument("-rX", "--rebinX", type=int, default=1)
    parser.add_argument("-rY", "--rebinY", type=int, default=1)
    parser.add_argument("-rZ", "--rebinZ", type=int, default=1)
    parser.add_argument("-pX", "--projectionX", action="store_true")
    parser.add_argument("-pY", "--projectionY", action="store_true") 
    parser.add_argument("-pZ", "--projectionZ", action="store_true") 
    parser.add_argument("-pXY", "--projectionXY", action="store_true")
    parser.add_argument("-pYZ", "--projectionYZ", action="store_true") 
    parser.add_argument("-pXZ", "--projectionXZ", action="store_true") 
    return with_graphsize_args(parser)

def resize_graph(args, hist):
    if args.range_lo is not None and args.range_hi is not None: 
        hist.GetXaxis().SetLimits(args.range_lo, args.range_hi)
    if args.y_min is not None and args.y_max is not None: 
        hist.GetYaxis().SetRangeUser(args.y_min, args.y_max)

def resize_histo(args, hist):
    histdim = 0
    if isinstance(hist, ROOT.TH1D): histdim = 1
    elif isinstance(hist, ROOT.TH2D): histdim = 2
    elif isinstance(hist, ROOT.TH3D): histdim = 3

    axes_ranges = [(args.x_min, args.x_max), (args.y_min, args.y_max), (args.z_min, args.z_max)]
    for i_axis, (lo, hi) in enumerate(axes_ranges):
        if lo is not None and hi is not None:
            if i_axis >= histdim:
                raise Exception("Error: setting axis (%i) for histogram of dimmension (%i)" % (i_axis, histdim))
                
            if i_axis == 0: axis = hist.GetXaxis()
            elif i_axis == 1: axis = hist.GetYaxis()
            elif i_axis == 2: axis = hist.GetZaxis()

            range_lo_ind = None
            for i in range(1, axis.GetNbins()+1):
                if range_lo_ind is None and axis.GetBinLowEdge(i) > lo:
                    range_lo_ind = max(0, i-2)
                if axis.GetBinLowEdge(i) >= hi:
                    range_hi_ind = i
                    break
            else: 
                range_hi_ind = axis.GetNbins()+1
            assert(range_lo_ind < range_hi_ind)

            axis_range = [axis.GetBinLowEdge(i) for i in range(1,axis.GetNbins()+1)] + [axis.GetBinUpEdge(axis.GetNbins())]
            #new_axis_range = array('d', [axis_range[i] for i in range(range_lo_ind, range_hi_ind)])

            new_axis_args = [range_hi_ind - range_lo_ind - 1, axis_range[range_lo_ind], axis_range[range_hi_ind-1]]
            if histdim == 1:
                new_hist = ROOT.TH1D(hist.GetName() + " resized", hist.GetTitle(), *new_axis_args)
            elif histdim == 2:
               x_axis_args =  [hist.GetXaxis().GetNbins(), hist.GetXaxis().GetBinLowEdge(1), hist.GetXaxis().GetBinUpEdge(hist.GetXaxis().GetNbins())]
               y_axis_args = [hist.GetYaxis().GetNbins(), hist.GetYaxis().GetBinLowEdge(1), hist.GetYaxis().GetBinUpEdge(hist.GetYaxis().GetNbins())]
               histo_args = new_axis_args + y_axis_args if i_axis == 0 else x_axis_args + new_axis_args
               new_hist = ROOT.TH2D(hist.GetName() + " resized", hist.GetTitle(), *histo_args)
            elif histdim == 3:
               x_axis_args = [hist.GetXaxis().GetNbins(), hist.GetXaxis().GetBinLowEdge(1), hist.GetXaxis().GetBinUpEdge(hist.GetXaxis().GetNbins())]
               y_axis_args = [hist.GetYaxis().GetNbins(), hist.GetYaxis().GetBinLowEdge(1), hist.GetYaxis().GetBinUpEdge(hist.GetYaxis().GetNbins())]
               z_axis_args = [hist.GetZaxis().GetNbins(), hist.GetZaxis().GetBinLowEdge(1), hist.GetZaxis().GetBinUpEdge(hist.GetZaxis().GetNbins())]
               if i_axis == 0:
                   histo_args = new_axis_args + y_axis_args + z_axis_args
               elif i_axis == 1:
                   histo_args = x_axis_args + new_axis_args + z_axis_args
               elif i_axis == 2:
                   histo_args = x_axis_args + y_axis_args + new_axis_args
               new_hist = ROOT.TH3D(hist.GetName() + " resized", hist.GetTitle(), *histo_args)

            for i in range(range_lo_ind+1, range_hi_ind+1):
                i_set = i - range_lo_ind
                if histdim == 1:
                    new_hist.SetBinContent(i - range_lo_ind, hist.GetBinContent(i))
                elif histdim == 2:
                    other_axis = hist.GetYaxis() if i_axis == 0 else hist.GetXaxis()
                    for j in range(1, other_axis.GetNbins()+1):
                        i_bin_old = hist.GetBin(i, j) if i_axis == 0 else hist.GetBin(j, i)
                        i_bin_new = new_hist.GetBin(i_set, j) if i_axis == 0 else new_hist.GetBin(j, i_set)
                        new_hist.SetBinContent(i_bin_new, hist.GetBinContent(i_bin_old))
                elif histdim == 3:
                    other_axes = [hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis()]
                    other_axes.pop(i_axis)
                    for j in range(1,other_axes[0].GetNbins()+1):
                        for k in range(1, other_axes[1].GetNbins()+1):
                            old_axes_indices = [j,k]
                            old_axes_indices.insert(i_axis, i)

                            new_axes_indices = [j,k]
                            new_axes_indices.insert(i_axis, i_set)

                            i_bin_old = hist.GetBin(*old_axes_indices)
                            i_bin_new = new_hist.GetBin(*new_axes_indices)
                            new_hist.SetBinContent(i_bin_new, hist.GetBinContent(i_bin_old))
            hist = new_hist
    if args.projectionX:
        hist = hist.ProjectionX()
    if args.projectionY:
        hist = hist.ProjectionY()
    if args.projectionZ:
        hist = hist.ProjectionZ()
    if args.projectionXY:
        hist = hist.Project3D("yx")
    if args.projectionYZ:
        hist = hist.Project3D("zy")
    if args.projectionXZ:
        hist = hist.Project3D("zx")

    if args.rebinX != 1:
        hist.RebinX(args.rebinX) 
    if args.rebinY != 1:
        hist.RebinY(args.rebinY) 
    if args.rebinZ != 1:
        hist.RebinZ(args.rebinZ) 
    return hist

def with_histostyle_args(parser):
    parser.add_argument("-xl", "--xlabel", default=None)
    parser.add_argument("-yl", "--ylabel", default=None)
    parser.add_argument("-os", "--optstat", default=None)
    return parser

def optstat(args):
    if args.optstat:
        ROOT.gStyle.SetOptStat(args.optstat)

def style(args, hist):
    hist.GetYaxis().SetTitleSize(20)
    hist.GetYaxis().SetTitleFont(43)
    hist.GetYaxis().SetLabelFont(43)
    hist.GetYaxis().SetLabelSize(20)
    hist.GetYaxis().CenterTitle()

    hist.GetXaxis().SetTitleSize(20)
    hist.GetXaxis().SetTitleFont(43)
    hist.GetXaxis().SetLabelFont(43)
    hist.GetXaxis().SetLabelSize(20)
    hist.GetXaxis().CenterTitle()

    if args.xlabel: hist.GetXaxis().SetTitle(args.xlabel)
    if args.ylabel: hist.GetYaxis().SetTitle(args.ylabel)

def colors(index):
    colors = [ROOT.kRed, ROOT.kGreen]
    return colors[index]

def fillcolors(index):
    colors = [ROOT.kGreen + 1, ROOT.kCyan -9 , ROOT.kBlue+1,  ROOT.kViolet, ROOT.kRed +2, ROOT.kYellow+3, ROOT.kBlack]
    return colors[index]

def comma_separated(inp):
    return inp.split(",")

def legend_position(inp):
    if inp == "ur":
        return [0.75,0.75,0.95,0.95]
    elif inp == "ul":
        return [0.35,0.75,0.15,0.95]
    elif inp == "um":
        return [0.4, 0.69, 0.6, 0.89]
    else:
        return [float(x) for x in inp.split(",")][:4]

def histo_list(inp):
    inp = inp.split(",")
    in_parens = False
    out = []
    for name in inp:
        if name.startswith("("):
            assert(not in_parens)
            in_parens = True
            name = name.lstrip("(")
            out.append([name])
        elif name.endswith(")"):
            assert(in_parens)
            in_parens = False
            name = name.rstrip(")")
            out[-1].append(name)
        elif in_parens:
            out[-1].append(name)
        else:
            out.append([name])
    return out

def validate_hists(names, hists):
    for nn,hh in zip(names, hists):
        if isinstance(hh, list):
            for n,h in zip(nn,hh):
                if not h:
                     raise Exception("Error: invalid histogram with name (%s)" % n)
                    # raise Exception("Error: invalid histogram (%s) with name (%s)" % (h, n))
        else:
            if not hh:
                raise Exception("Error: invalid histogram with name (%s)" % (nn))
                #raise Exception("Error: invalid histogram (%s) with name (%s)" % (hh, nn))


