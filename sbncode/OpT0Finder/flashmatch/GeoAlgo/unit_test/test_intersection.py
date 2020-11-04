from test_msg import debug, info, error, warning
import traceback,sys
from random import *
import numpy as np
from time import *
from test_import import test_import
test_import()
from ROOT import geoalgo

_epsilon = 1E-6

# colors
OK = '\033[92m'
NO = '\033[91m'
BLUE = '\033[94m'
ENDC = '\033[0m'

def test_dAlgo():

    debug()
    debug(BLUE + "Precision Being Required to Consider Two numbers Equal: {0:.2e}".format(_epsilon) + ENDC)
    debug()

    # number of times to test each function
    tests = 10000

    # import Distance Algo
    iAlgo = geoalgo.GeoAlgo()
    
    try:

        # test distance to and back from wall
        info('Testing Intersection Between Half-Line & AABox')
        totSuccess_f = 0
        intersectT_f = 0
        totSuccess_b = 0
        intersectT_b = 0
        for y in range(tests):
            success_f = 1
            success_b = 1
            # generate a unit cubic box from (0,0,0) to (1,1,1)
            box = geoalgo.AABox(0,0,0,1,1,1)
            # generate random point inside box
            # half-line will start from this point
            p = geoalgo.Vector(random(),random(),random())
            # now find random intersection point
            # pick a side of the box that the point is on
            pick = random()
            side = 0
            if ( (pick > 0.33) and (pick < 0.67) ) : side = 1
            if ( pick > 0.67 ) : side = 2
            # pick a direction (Left/Right), (Top/Bottom), (Front/Back) that the point should be on
            # given the direction & side, find the intersection point on the relevant face of the cube. That will be our intersection point. Make the half-line pass through there
            direction = 1
            if ( random() < 0.5 ) : direction = -1
            # pl is a parameter to know numerical value of coordinate for the face we are on
            # if direction = 1  -> pl = +1 (box.Max())
            # if direction = -1 -> pl = -1 (box.Min())
            pl = 1
            if ( direction == -1 ) :
                pl = 0
            if ( side == 0 ) :
                i = geoalgo.Vector(pl,random(),random())
            elif ( side == 1 ) :
                i = geoalgo.Vector(random(),pl,random())
            else :
                i = geoalgo.Vector(random(),random(),pl)
            # now generate half-line passing thorugh p and i & starting from p
            d = i-p
            lin = geoalgo.HalfLine(p,d)
            # answer should be distance between p & i
            # point should be i
            answer = i.SqDist(p)
            pt1 = geoalgo.Vector(3)
            pt2 = geoalgo.Vector(3)
            tim = time()
            pt1_v = iAlgo.Intersection(box,lin)
            intersectT_f += (time()-tim)
            pt2_v = iAlgo.Intersection(lin,box)
            if pt1_v.size(): pt1 = pt1_v[0]
            if pt2_v.size(): pt2 = pt2_v[0]
            a1 = pt1.SqDist(p)
            a2 = pt2.SqDist(p)
            if not ( np.abs(answer-a1) < _epsilon ) : success_f = 0
            if not ( np.abs(answer-a2) < _epsilon ) : success_f = 0
            for x in xrange(3):
                if not ( np.abs(pt1[x]-i[x]) < _epsilon) : success_f = 0
                if not ( np.abs(pt1[x]-i[x]) < _epsilon) : success_f = 0
            totSuccess_f += success_f
            # now backwards:
            # make line go in opposite direction -> intersection point should be the same
            d = p-i
            lin = geoalgo.HalfLine(p,d)
            pt1 = geoalgo.Vector(3)
            pt2 = geoalgo.Vector(3)
            tim = time()
            pt1_v = iAlgo.Intersection(box,lin,1)
            intersectT_b += (time()-tim)
            pt2_v = iAlgo.Intersection(lin,box,1)
            if pt1_v.size(): pt1 = pt1_v[0]
            if pt2_v.size(): pt2 = pt2_v[0]
            a1 = pt1.SqDist(p)
            a2 = pt2.SqDist(p)
            if not ( np.abs(answer-a1) < _epsilon ) : success_b = 0
            if not ( np.abs(answer-a2) < _epsilon ) : success_b = 0
            for x in xrange(3):
                if not ( np.abs(pt1[x]-i[x]) < _epsilon) : success_b = 0
                if not ( np.abs(pt1[x]-i[x]) < _epsilon) : success_b = 0
            if (success_b == 1) : totSuccess_b += 1

        if ( float(totSuccess_f)/tests < 1):
            info(NO + "Success: {0}%".format(100*float(totSuccess_f)/tests) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess_f)/tests) + ENDC)
        info("Time for Intersection                    : {0:.3f} us".format(1E6*intersectT_f/tests))
        if ( float(totSuccess_b)/tests < 1):
            info(NO + "Success: {0}%".format(100*float(totSuccess_b)/tests) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess_b)/tests) + ENDC)
        info("Time for Intersection                    : {0:.3f} us".format(1E6*intersectT_b/tests))


    except Exception:
        error('geoalgo::IntersectAlgo unit test failed.')
        print(traceback.format_exception(*sys.exc_info())[2])
        return 1
    
    info('geoalgo::IntersectAlgo unit test complete.')
    return 0

if __name__ == '__main__':
    import test_msg
    test_msg.test_msg.level = 0
    test_dAlgo()

