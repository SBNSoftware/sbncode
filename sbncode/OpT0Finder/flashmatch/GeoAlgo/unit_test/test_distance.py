from test_msg import debug, info, error, warning
import traceback,sys
from random import *
import numpy as np
from time import *
from test_import import test_import
test_import()
from ROOT import geoalgo

_epsilon = 1E-7

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
    dAlgo = geoalgo.GeoAlgo()
    
    try:

        # test distance from point to infinite line
        info('Testing Point & Infinite Line Distance')
        totSuccess = 0
        sqdistT = 0
        closestT = 0
        for y in range(tests):
            success = 1
            # generate a random point (will be point on line closest to external point
            p1 = geoalgo.Vector(random(),random(),random())
            # second point on line
            p2 = geoalgo.Vector(random(),random(),random())
            # back-project pt1 so that we crate a line starting from an earlier position
            pDir = p2-p1
            p0 = p1 - pDir.Dir() * 3;
            # generate random line
            l = geoalgo.Line(p0,p2)
            # generate a point in a direction transverse to the line
            transX = random()
            transY = random()
            # now ensure perpendicularity by having dot product be == 0
            dirct = l.Pt2()-l.Pt1()
            transZ = (-transX*dirct[0]-transY*dirct[1])/dirct[2]
            vectTranslate = geoalgo.Vector(transX,transY,transZ)
            # generate point away starting from p1
            pt = geoalgo.Vector(p1+vectTranslate)
            # answer will be vectTranslate sq. length
            answer = vectTranslate.SqLength()
            # closest point will be p1
            tim = time()
            a1 = dAlgo.SqDist(pt,l)
            sqdistT += (time()-tim)
            a2 = dAlgo.SqDist(l,pt)
            tim = time()
            pAnswer1 = dAlgo.ClosestPt(pt,l)
            closestT += (time()-tim)
            pAnswer2 = dAlgo.ClosestPt(l,pt)
            if not ( np.abs(answer-a1) < _epsilon ) : success = 0
            if not ( np.abs(answer-a2) < _epsilon ) : success = 0
            for x in range(3):
                if not ( np.abs(p1[x]-pAnswer1[x]) < _epsilon) : success = 0
                if not ( np.abs(p1[x]-pAnswer2[x]) < _epsilon) : success = 0
            totSuccess += success
        if ( float(totSuccess)/tests < 0.999):
            info(NO + "Success: {0}%".format(100*float(totSuccess)/tests) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess)/tests) + ENDC)
        info("Time for SqDist                       : {0:.3f} us".format(1E6*sqdistT/tests))
        info("Time for ClosestPt                    : {0:.3f} us".format(1E6*closestT/tests))

        # test distance from point to segment
        info('Testing Point & LineSegment Distance')
        totSuccess = 0
        sqdistT_out = 0.
        closestT_out = 0.
        sqdistT_in = 0.
        closestT_in = 0.
        for y in range(tests):
            success = 1
            # generate a random segment
            l = geoalgo.LineSegment(random(),random(),random(),random(),random(),random())
            # get the segment length
            lLen = l.Dir().Length()
            # get the segment direction
            d = l.Dir()
            # generate a segment parallel to it
            # to get a line parallel to the first,
            # translate the Start and End points
            # by some amount in a direction perpendicular
            # to that of the original line
            transX = random()
            transY = random()
            # now ensure perpendicularity by having dot product be == 0
            transZ = (-transX*d[0]-transY*d[1])/d[2]
            vectTranslate = geoalgo.Vector(transX,transY,transZ)
            p1 = l.Start()+vectTranslate
            p2 = l.End()+vectTranslate
            # parallel segment:
            lPar = geoalgo.LineSegment(p1,p2)

            # first, test a point that is "before start point"
            # distance to this point should be distance between lines
            # in quadrature with how much further from start point we go
            dist = random()
            # direction vector outwards from line
            dirOut = lPar.Dir()*(-1*dist/lPar.Dir().Length())
            pTest = p1+dirOut
            answer = dirOut.SqLength()+vectTranslate.SqLength()
            tim = time()
            a1 = dAlgo.SqDist(pTest,l)
            sqdistT_out += (time() - tim)
            a2 = dAlgo.SqDist(l,pTest)
            if not ( np.abs(answer-a1) < _epsilon): success = 0
            if not ( np.abs(answer-a2) < _epsilon): success = 0
            tim = time()
            point1 = dAlgo.ClosestPt(pTest,l)
            closestT_out += (time() - tim)
            point2 = dAlgo.ClosestPt(l,pTest)
            for x in range(3):
                if not ( (l.Start()[x]-point1[x]) < _epsilon ) : success = 0
                if not ( (l.Start()[x]-point2[x]) < _epsilon ) : success = 0

            # now, test a point inside the segment.
            # distance between point & segment should be 
            # pre-determined distance between the two segments
            dist = random()
            dirIn = lPar.Dir()*dist #dist ensures < distance of full segment
            pTest = p1+dirIn
            answer = vectTranslate.SqLength()
            tim = time()
            a1 = dAlgo.SqDist(pTest,l)
            sqdistT_in += (time() - tim)
            a2 = dAlgo.SqDist(l,pTest)
            if not (np.abs(answer-a1) < _epsilon): success = 0
            if not (np.abs(answer-a2) < _epsilon): success = 0
            pAns = l.Start()+dirIn
            tim = time()
            point1 = dAlgo.ClosestPt(pTest,l)
            closestT_in += (time() - tim)
            point2 = dAlgo.ClosestPt(l,pTest)
            for x in range(3):
                if not ( (pAns[x]-point1[x]) < _epsilon ) : success = 0
                if not ( (pAns[x]-point2[x]) < _epsilon ) : success = 0

            # now test a point beyond the segment
            dist = random()
            dirOut = lPar.Dir()*(dist/lPar.Dir().Length())
            pTest = p2+dirOut
            answer = dirOut.SqLength()+vectTranslate.SqLength()
            if not ( np.abs(answer-dAlgo.SqDist(pTest,l)) < _epsilon): success = 0
            if not ( np.abs(answer-dAlgo.SqDist(l,pTest)) < _epsilon): success = 0
            point1 = dAlgo.ClosestPt(pTest,l)
            point2 = dAlgo.ClosestPt(pTest,l)
            for x in range(3):
                if not ( (l.End()[x]-point1[x]) < _epsilon ) : success = 0
                if not ( (l.End()[x]-point2[x]) < _epsilon ) : success = 0
                
            if (success == 1) : totSuccess += 1
        if ( float(totSuccess)/tests < 0.999):
            info(NO + "Success: {0}%".format(100*float(totSuccess)/tests) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess)/tests) + ENDC)
        info("Time for SqDist (Pt Out of Segment)   : {0:.3f} us".format(1E6*sqdistT_out/tests))
        info("Time for ClosestPt (Pt Out of Segment): {0:.3f} us".format(1E6*closestT_out/tests))
        info("Time for SqDist (Pt In Segment)       : {0:.3f} us".format(1E6*sqdistT_in/tests))
        info("Time for ClosestPt (Pt In Segment)    : {0:.3f} us".format(1E6*closestT_in/tests))
            
        # test Point to HalfLine distance
        debug('Testing Point & HalfLine Distance')
        success = 1
        totSuccess = 0
        sqdistT_out = 0.
        closestT_out = 0.
        sqdistT_in = 0.
        closestT_in = 0.
        for x in range(tests):
            # generate a random segment
            l = geoalgo.HalfLine(random(),random(),random(),random(),random(),random())
            # generate a segment parallel to it
            # to get a line parallel to the first,
            # translate the Start and End points
            # by some amount in a direction perpendicular
            # to that of the original line
            transX = random()
            transY = random()
            # now ensure perpendicularity by having dot product be == 0
            transZ = (-transX*l.Dir()[0]-transY*l.Dir()[1])/l.Dir()[2]
            vectTranslate = geoalgo.Vector(transX,transY,transZ)
            # create the new translated & parallel hal-fline
            lPar = geoalgo.HalfLine(l.Start()+vectTranslate, l.Dir())
            # first, test a point that is "before start point"
            # distance to this point should be distance between lines
            # in quadrature with how much further from start point we go
            dist = random()
            # direction vector outwards from line
            dirOut = lPar.Dir()*(-1*dist)
            pTest = lPar.Start()+dirOut
            answer = dirOut.SqLength()+vectTranslate.SqLength()
            tim = time()
            a1 = dAlgo.SqDist(pTest,l)
            tim = time() - tim
            sqdistT_out += tim
            a2 = dAlgo.SqDist(l,pTest)
            if not ( np.abs(answer-a1) < _epsilon): success = 0
            if not ( np.abs(answer-a2) < _epsilon): success = 0
            tim = time()
            point1 = dAlgo.ClosestPt(pTest,l)
            tim = time() - tim
            closestT_out += tim
            point2 = dAlgo.ClosestPt(l,pTest)
            for x in range(3):
                if not ( (l.Start()[x]-point1[x]) < _epsilon ) : success = 0
                if not ( (l.Start()[x]-point2[x]) < _epsilon ) : success = 0

            # now, test a point inside the segment.
            # distance between point & segment should be 
            # pre-determined distance between the two segments
            dist = random()
            dirIn = lPar.Dir()*dist #dist ensures < distance of full segment
            pTest = lPar.Start()+dirIn
            answer = vectTranslate.SqLength()
            pAns = l.Start()+dirIn
            tim = time()
            a1 = dAlgo.SqDist(pTest,l)
            tim = time() - tim
            sqdistT_in += tim
            a2 = dAlgo.SqDist(l,pTest)
            if not ( np.abs(answer-a1) < _epsilon): success = 0
            if not ( np.abs(answer-a2) < _epsilon): success = 0
            tim = time()
            point1 = dAlgo.ClosestPt(pTest,l)
            tim = time() - tim
            closestT_in += tim
            point2 = dAlgo.ClosestPt(l,pTest)
            for x in range(3):
                if not ( (pAns[x]-point1[x]) < _epsilon ) : success = 0
                if not ( (pAns[x]-point2[x]) < _epsilon ) : success = 0

            if (success == 1) : totSuccess += 1
        if ( float(totSuccess)/tests < 0.999):
            info(NO + "Success: {0}%".format(100*float(totSuccess)/tests) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess)/tests) + ENDC)
        info("Time for SqDist (Pt Out of Segment)   : {0:.3f} us".format(1E6*sqdistT_out/tests))
        info("Time for ClosestPt (Pt Out of Segment): {0:.3f} us".format(1E6*closestT_out/tests))
        info("Time for SqDist (Pt In Segment)       : {0:.3f} us".format(1E6*sqdistT_in/tests))
        info("Time for ClosestPt (Pt In Segment)    : {0:.3f} us".format(1E6*closestT_in/tests))

        # test Distance between two Infinite Lines
        debug('Testing Inf Line & Inf Line Distance')
        totSuccess = 0
        sqdistT = 0.
        closestT = 0.
        for y in range(tests):
            success = 1
            l1 = geoalgo.Line(random(),random(),random(),random(),random(),random())
            # take a point a fixed distance away
            # generate a random direction in the plane
            # perpendicular to the plane connecting
            # the point and the line
            # the distance between the two lines
            # should be the fixed amount selected previously
            # use half-way point to do calculations
            p1 = (l1.Pt2()+l1.Pt1())/2
            # get line direction
            d1 = (l1.Pt2()-l1.Pt1())
            # move in a random direction perpendicular to line
            dirx = random()
            diry = random()
            dirz = (-dirx*d1[0]-diry*d1[1])/d1[2]
            vectTranslate = geoalgo.Vector(dirx,diry,dirz)
            # need to re-orient in some random direction on this plane
            # this direction has to be perpendicular to both
            # the line's direction as well as to the direction
            # of the segment uniting the two lines
            # use cross-product
            vectRotate = vectTranslate.Cross(d1)
            l2 = geoalgo.Line(p1+vectTranslate,p1+vectTranslate+vectRotate)
            # now calculate distance. Should be vectTranslate.Length()
            answer = vectTranslate.SqLength()
            tim = time()
            a1 = dAlgo.SqDist(l1,l2)
            tim = time() - tim
            sqdistT += tim
            # expect the closest points on both lines to be p1 & l2.Pt1()
            ptL1 = geoalgo.Vector()
            ptL2 = geoalgo.Vector()
            a2 = dAlgo.SqDist(l1,l2,ptL1,ptL2)
            if not (np.abs(answer-a1) < _epsilon): success = 0
            if not (np.abs(answer-a2) < _epsilon) : success = 0
            for x in range(3):
                if not ( np.abs(ptL1[x]-p1[x]) < _epsilon ) : success = 0
                if not ( np.abs(ptL2[x]-l2.Pt1()[x]) < _epsilon ) : success = 0

            if (success == 1) : totSuccess += 1
        if ( float(totSuccess)/tests < 0.999):
            info(NO + "Success: {0}%".format(100*float(totSuccess)/tests) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess)/tests) + ENDC)
        info("Time for SqDist                       : {0:.3f} us".format(1E6*sqdistT/tests))

        # test Distance between two Half-Infinite Lines
        debug('Testing Half-Inf Line & Half-Inf Line Distance')
        totSuccess = 0
        sqdistT_in = 0
        timesIN = 0
        sqdistT_out = 0
        timesOUT = 0
        for y in range(tests):
            success = 1
            l1 = geoalgo.HalfLine(random(),random(),random(),random(),random(),random())
            # take a point a fixed distance away
            # then generate a new half-line starting
            # at the same point (translated) and with
            # the same (or opposite) direction.
            # But rotated outward a bit
            p1 = l1.Start()
            # get line direction
            d1 = l1.Dir()
            # move in a random direction perpendicular to line
            dirx = random()
            diry = random()
            dirz = (-dirx*d1[0]-diry*d1[1])/d1[2]
            vectTranslate = geoalgo.Vector(dirx,diry,dirz)
            # pick some random direction (aligned with 1st or not)
            aligned = -1
            if (random() < 0.5) : aligned = 1
            l2 = geoalgo.HalfLine(p1+vectTranslate,d1+vectTranslate*aligned)
            # now calculate distance.
            # if aligned == -1 distance should be == 0 (lines intersect)
            # if aligned == 1 then the two start points are the closest points
            if (aligned == -1):
                timesIN += 1
                answer = 0.
                tim = time()
                a1 = dAlgo.SqDist(l1,l2)
                tim = time() - tim
                L1 = geoalgo.Vector()
                L2 = geoalgo.Vector()
                a2 = dAlgo.SqDist(l1,l2,L1,L2)
                sqdistT_in += tim
                if not (np.abs(answer-a1) < _epsilon): success = 0
                if not (np.abs(answer-a2) < _epsilon): success = 0
                for x in range(3):
                    if not (L1[x]-L2[x] < _epsilon) : success = 0
            if (aligned == 1):
                timesOUT += 1
                answer = vectTranslate.SqLength()
                tim = time()
                a1 = dAlgo.SqDist(l1,l2)
                tim = time() - tim
                L1 = geoalgo.Vector()
                L2 = geoalgo.Vector()
                a2 = dAlgo.SqDist(l1,l2,L1,L2)
                sqdistT_out += tim
                if not (np.abs(answer-a1) < _epsilon): success = 0
                if not (np.abs(answer-a2) < _epsilon): success = 0
                for x in range(3):
                    if not (L1[x]-l1.Start()[x] < _epsilon) : success = 0
                    if not (L2[x]-l2.Start()[x] < _epsilon) : success = 0

            if (success == 1) : totSuccess += 1
        if ( float(totSuccess)/tests < 0.999):
            info(NO + "Success: {0}%".format(100*float(totSuccess)/tests) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess)/tests) + ENDC)
        info("Time for SqDist (OUT)                 : {0:.3f} us".format(1E6*sqdistT_out/timesOUT))
        info("Time for SqDist (IN)                  : {0:.3f} us".format(1E6*sqdistT_in/timesIN))


        # test Distance between Half-Line and Line Segment
        debug('Testing Half Line & Line Segment Distance')
        totSuccess1 = 0
        sqdistT1 = 0
        totSuccess2 = 0
        sqdistT2 = 0
        totSuccess3 = 0
        sqdistT3 = 0
        totSuccess4 = 0
        sqdistT4 = 0
        for y in range(tests):
            success1 = 1
            success2 = 1
            success3 = 1
            success4 = 1
            l1 = geoalgo.HalfLine(random(),random(),random(),random(),random(),random())
            # take a point a fixed distance away
            # test multiple scenarios:
            # 1) segment "away" from half-line but same direction
            # -> closest points are half-line start and segment end
            # 2) segment "away" from half-line and rotated
            # -> closest points are half-line start and segment rotation pivot
            # 3) segment "close to" half-line and tilted outwards
            # -> closest points are point on half-line and segment start
            # 4) segment "close to" half-line and rotated
            # -> closest points are point on half-line and segment rotation pivot
            p1 = l1.Start()
            # get line direction
            d1 = l1.Dir()
            # move in a random direction perpendicular to line
            dirx = random()
            diry = random()
            dirz = (-dirx*d1[0]-diry*d1[1])/d1[2]
            vectTranslate = geoalgo.Vector(dirx,diry,dirz)
            # 1) generate line-segment "behind" half-line
            # use unit-vector nature of Dir() to go back slightly more
            dist = random()
            seg = geoalgo.LineSegment(p1+vectTranslate-d1*dist,p1+vectTranslate-d1*dist-d1)
            # distance should be dist
            # closest points should be seg.Start() and l1.Start()
            answer = vectTranslate*vectTranslate+dist*dist
            tim = time()
            a1 = dAlgo.SqDist(l1,seg)
            tim = time() - tim
            sqdistT1 += tim
            L1 = geoalgo.Vector()
            L2 = geoalgo.Vector()
            a2 = dAlgo.SqDist(l1,seg,L1,L2)
            if not (np.abs(answer-a1) < _epsilon): success1 = 0
            if not (np.abs(answer-a2) < _epsilon): success1 = 0
            for x in range(3):
                if not (L1[x]-l1.Start()[x] < _epsilon) : success1 = 0
                if not (L2[x]-seg.Start()[x] < _epsilon) : success1 = 0
            if (success1 == 1) : totSuccess1 += 1
            # 2) generate line segment away but rotated
            # rotate in a direction perpendicular to both line direction and translation direction
            vectRotate = vectTranslate.Cross(d1)
            # choose pivot point
            pivot = p1+vectTranslate-d1*dist
            seg = geoalgo.LineSegment(pivot-vectRotate*random(),pivot+vectRotate*random())
            # distance should be pivot distance to half-line start
            # closest points should be pivot and line start point
            answer = vectTranslate*vectTranslate+dist*dist
            tim = time()
            a1 = dAlgo.SqDist(l1,seg,L1,L2)
            sqdistT2 += (time()-tim)
            a2 = dAlgo.SqDist(seg,l1)
            if not (np.abs(answer-a1) < _epsilon): success2 = 0
            if not (np.abs(answer-a2) < _epsilon): success2 = 0
            for x in range(3):
                if not (L1[x]-l1.Start()[x] < _epsilon) : success2 = 0
                if not (L2[x]-pivot[x] < _epsilon) : success2 = 0
            if (success2 == 1) : totSuccess2 += 1
            # 3) generate line segment next to but tilted outwards
            seg = geoalgo.LineSegment(p1+vectTranslate+d1*dist,p1+vectTranslate+d1*dist+vectTranslate*random())
            # distance should be vectTranslate
            # closest point should be segment start and point on line "d1*dist" away from line-start
            answer = vectTranslate*vectTranslate
            tim = time()
            a1 = dAlgo.SqDist(l1,seg,L1,L2)
            sqdistT3 += (time()-tim)
            a2 = dAlgo.SqDist(seg,l1)
            if not (np.abs(answer-a1) < _epsilon): success3 = 0
            if not (np.abs(answer-a2) < _epsilon): success3 = 0
            ptLine = l1.Start()+d1*dist
            for x in range(3):
                if not (L1[x]-ptLine[x] < _epsilon) : success3 = 0
                if not (L2[x]-seg.Start()[x] < _epsilon) : success3 = 0
            if (success3 == 1) : totSuccess3 += 1
            # 4) generate line segment next to line but rotated
            pivot = p1+vectTranslate+d1*dist
            seg = geoalgo.LineSegment(pivot-vectRotate*random(),pivot+vectRotate*random())
            # distance should be vectTranslate
            # closest point should be pivot on segment and d1*dist away from start for line
            answer = vectTranslate*vectTranslate
            tim = time()
            a1 = dAlgo.SqDist(l1,seg,L1,L2)
            sqdistT4 += (time()-tim)
            a2 = dAlgo.SqDist(seg,l1)
            if not (np.abs(answer-a1) < _epsilon): success4 = 0
            if not (np.abs(answer-a2) < _epsilon): success4 = 0
            ptLine = l1.Start()+d1*dist
            for x in range(3):
                if not (L1[x]-ptLine[x] < _epsilon) : success4 = 0
                if not (L2[x]-pivot[x] < _epsilon) : success4 = 0
            if (success4 == 1) : totSuccess4 += 1

        if ( float(totSuccess1)/tests < 0.999):
            info(NO + "Success: {0}%".format(100*float(totSuccess1)/tests) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess1)/tests) + ENDC)
        info("Time for SqDist (Case 1)              : {0:.3f} us".format(1E6*sqdistT1/tests))
        if ( float(totSuccess2)/tests < 0.999):
            info(NO + "Success: {0}%".format(100*float(totSuccess2)/tests) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess2)/tests) + ENDC)
        info("Time for SqDist (Case 2)              : {0:.3f} us".format(1E6*sqdistT2/tests))
        if ( float(totSuccess3)/tests < 0.999):
            info(NO + "Success: {0}%".format(100*float(totSuccess3)/tests) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess3)/tests) + ENDC)
        info("Time for SqDist (Case 3)              : {0:.3f} us".format(1E6*sqdistT3/tests))
        if ( float(totSuccess4)/tests < 0.999):
            info(NO + "Success: {0}%".format(100*float(totSuccess4)/tests) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess4)/tests) + ENDC)
        info("Time for SqDist (Case 4)              : {0:.3f} us".format(1E6*sqdistT4/tests))

        # test Distance between Half-Line and Line Segment
        debug('Testing Line Segment & Line Segment Distance')
        totSuccess1 = 0.
        sqdistT1 = 0
        totSuccess2 = 0.
        sqdistT2 = 0
        totSuccess3 = 0.
        sqdistT3 = 0
        totSuccess4 = 0.
        sqdistT4 = 0
        for y in range(tests):
            success1 = 1
            success2 = 1
            success3 = 1
            success4 = 1
            seg1 = geoalgo.LineSegment(random(),random(),random(),random(),random(),random())
            # take a point a fixed distance away
            # test multiple scenarios:
            # 1) segment "away" from half-line but same direction
            # -> closest points are half-line start and segment end
            # 2) segment "away" from half-line and rotated
            # -> closest points are half-line start and segment rotation pivot
            # 3) segment "close to" half-line and tilted outwards
            # -> closest points are point on half-line and segment start
            # 4) segment "close to" half-line and rotated
            # -> closest points are point on half-line and segment rotation pivot
            p1 = seg1.Start()
            # get line direction
            p2 = seg1.End()
            d  = p2-p1
            length = d.Length()
            # move in a random direction perpendicular to line
            dirx = random()
            diry = random()
            dirz = (-dirx*d[0]-diry*d[1])/d[2]
            vectTranslate = geoalgo.Vector(dirx,diry,dirz)
            # 1) generate line-segment "behind" half-line
            # use unit-vector nature of Dir() to go back slightly more
            dist = random()
            seg = geoalgo.LineSegment(p1+vectTranslate-d*dist,p1+vectTranslate-d*dist-d)
            #seg = geoalgo.LineSegment(p1+vectTranslate,p1+vectTranslate-d)
            # distance should be dist
            # closest points should be seg.Start() and seg1.Start()
            answer = vectTranslate*vectTranslate+dist*dist*d.SqLength()
            tim = time()
            a1 = dAlgo.SqDist(seg1,seg)
            tim = time() - tim
            sqdistT1 += tim
            L1 = geoalgo.Vector()
            L2 = geoalgo.Vector()
            a2 = dAlgo.SqDist(seg1,seg,L1,L2)
            if not (np.abs(answer-a1) < _epsilon): success1 = 0
            if not (np.abs(answer-a2) < _epsilon): success1 = 0
            for x in range(3):
                if not (L1[x]-seg1.Start()[x] < _epsilon) : success1 = 0
                if not (L2[x]-seg.Start()[x] < _epsilon) : success1 = 0
            if (success1 == 1) : totSuccess1 += 1
            # 2) generate line segment away but rotated
            # rotate in a direction perpendicular to both line direction and translation direction
            vectRotate = vectTranslate.Cross(d)
            # choose pivot point
            pivot = p1+vectTranslate-d*dist
            seg = geoalgo.LineSegment(pivot-vectRotate*random(),pivot+vectRotate*random())
            # distance should be pivot distance to half-line start
            # closest points should be pivot and line start point
            answer = vectTranslate*vectTranslate+dist*dist*d.SqLength()
            tim = time()
            a1 = dAlgo.SqDist(seg1,seg,L1,L2)
            sqdistT2 += (time()-tim)
            a2 = dAlgo.SqDist(seg,seg1)
            if not (np.abs(answer-a1) < _epsilon): success2 = 0
            if not (np.abs(answer-a2) < _epsilon): success2 = 0
            for x in range(3):
                if not (L1[x]-seg1.Start()[x] < _epsilon) : success2 = 0
                if not (L2[x]-pivot[x] < _epsilon) : success2 = 0
            if (success2 == 1) : totSuccess2 += 1
            # 3) generate line segment next to but tilted outwards
            distin = length*random() # ensures that we do not pass the segment
            seg = geoalgo.LineSegment(p1+vectTranslate+d*distin,p1+vectTranslate+d*distin+vectTranslate*random())
            # distance should be vectTranslate
            # closest point should be segment start and point on line "d*distin" away from line-start
            answer = vectTranslate*vectTranslate
            tim = time()
            a1 = dAlgo.SqDist(seg1,seg,L1,L2)
            sqdistT3 += (time()-tim)
            a2 = dAlgo.SqDist(seg,seg1)
            if not (np.abs(answer-a1) < _epsilon): success3 = 0
            if not (np.abs(answer-a2) < _epsilon): success3 = 0
            ptLine = seg1.Start()+d*distin
            for x in range(3):
                if not (L1[x]-ptLine[x] < _epsilon) : success3 = 0
                if not (L2[x]-seg.Start()[x] < _epsilon) : success3 = 0
            if (success3 == 1) : totSuccess3 += 1
            # 4) generate line segment next to line but rotated
            pivot = p1+vectTranslate+d*distin
            seg = geoalgo.LineSegment(pivot-vectRotate*random(),pivot+vectRotate*random())
            # distance should be vectTranslate
            # closest point should be pivot on segment and d*distin away from start for line
            answer = vectTranslate*vectTranslate
            tim = time()
            a1 = dAlgo.SqDist(seg1,seg,L1,L2)
            sqdistT4 += (time()-tim)
            a2 = dAlgo.SqDist(seg,seg1)
            if not (np.abs(answer-a1) < _epsilon): success4 = 0
            if not (np.abs(answer-a2) < _epsilon): success4 = 0
            ptLine = seg1.Start()+d*distin
            for x in range(3):
                if not (L1[x]-ptLine[x] < _epsilon) : success4 = 0
                if not (L2[x]-pivot[x] < _epsilon) : success4 = 0
            if (success4 == 1) : totSuccess4 += 1

        if ( float(totSuccess1)/tests < 0.999):
            info(NO + "Success: {0}%".format(100*float(totSuccess1)/tests) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess1)/tests) + ENDC)
        info("Time for SqDist (Case 1)              : {0:.3f} us".format(1E6*sqdistT1/tests))
        if ( float(totSuccess2)/tests < 0.999):
            info(NO + "Success: {0}%".format(100*float(totSuccess2)/tests) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess2)/tests) + ENDC)
        info("Time for SqDist (Case 2)              : {0:.3f} us".format(1E6*sqdistT2/tests))
        if ( float(totSuccess3)/tests < 0.999):
            info(NO + "Success: {0}%".format(100*float(totSuccess3)/tests) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess3)/tests) + ENDC)
        info("Time for SqDist (Case 3)              : {0:.3f} us".format(1E6*sqdistT3/tests))
        if ( float(totSuccess4)/tests < 0.999):
            info(NO + "Success: {0}%".format(100*float(totSuccess4)/float(tests)) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess4)/float(tests)) + ENDC)
        info("Time for SqDist (Case 4)              : {0:.3f} us".format(1E6*sqdistT4/tests))


        # test Distance between Point and AABox
        debug('Testing Point and AABox Distance/Closest Point')
        success1 = 1
        totSuccess1 = 0.
        sqdistT1 = 0
        closestT1 = 0
        success2 = 1
        totSuccess2 = 0.
        sqdistT2 = 0
        closestT2 = 0
        success3 = 1
        totSuccess3 = 0.
        sqdistT3 = 0
        closestT3 = 0
        success4 = 1
        totSuccess4 = 0.
        sqdistT4 = 0
        closestT4 = 0
        for y in range(tests):
            success1 = 1
            success2 = 1
            success3 = 1
            success4 = 1
            # create a simple cubic box from 0,0,0 to 1,1,1
            b = geoalgo.AABox(0,0,0,1,1,1)
            # various cases to test:
            # case 1) Point inside box
            # case 2) Point outside box
            # 1) POINT INSIDE BOX
            p = geoalgo.Vector(random(),random(),random())
            # if point not recognized as inside -> error!
            if not ( b.Contain(p) ) : success1 = 0
            dTop = 1-p[2]
            dBot = p[2]
            dLeft = p[0]
            dRight = 1-p[0]
            dFront = p[1]
            dBack = 1-p[1]
            # find smallest of these
            dMin = dTop
            if ( dBot < dMin )   : dMin = dBot
            if ( dLeft < dMin )  : dMin = dLeft
            if ( dRight < dMin ) : dMin = dRight
            if ( dFront < dMin ) : dMin = dFront
            if ( dBack < dMin )  : dMin = dBack
            answer = dMin*dMin
            tim = time()
            a1 = dAlgo.SqDist(p,b)
            sqdistT1 += (time()-tim)
            a2 = dAlgo.SqDist(b,p)
            # closest point
            tim = time()
            pt1 = dAlgo.ClosestPt(b,p)
            closestT1 += (time()-tim)
            pt2 = dAlgo.ClosestPt(p,b)
            if not (np.abs(answer-a1) < _epsilon): success1 = 0
            if not (np.abs(answer-a2) < _epsilon): success1 = 0
            for x in range(3):
                if not (pt1[x]-p[x] < _epsilon) : success1 = 0
                if not (pt2[x]-p[x] < _epsilon) : success1 = 0
            if (success1 == 1) : totSuccess1 += 1
            # 2) POINT OUT OF BOX
            # pick a side that is should exit on at random
            pick = random()
            side = 0
            if ( (pick > 0.33) and (pick < 0.67) ) : side = 1
            if ( pick > 0.67 ) : side = 2
            # pick a direction to overflow the box by (-1 = back or +1 = front)
            direction = 1
            if ( random() < 0.5 ) : direction = -1
            if ( side == 0 ) :
                p = geoalgo.Vector(random()+direction,random(),random())
            elif ( side == 1 ) :
                p = geoalgo.Vector(random(),random()+direction,random())
            else :
                p = geoalgo.Vector(random(),random(),random()+direction)
            # if point not recognized as outside -> error!
            if ( b.Contain(p) ) : success1 = 0
            # go through cases and find min distance
            if ( (side == 0) and (direction == 1) ):
                #overflow on +x direction
                dMin = p[0]-1
                pMin = geoalgo.Vector(1,p[1],p[2])
            if ( (side == 0) and (direction == -1) ):
                #overflow on -x direction
                dMin = -p[0]
                pMin = geoalgo.Vector(0,p[1],p[2])
            if ( (side == 1) and (direction == 1) ):
                #overflow on +y direction
                dMin = p[1]-1
                pMin = geoalgo.Vector(p[0],1,p[2])
            if ( (side == 1) and (direction == -1) ):
                #overflow on -y direction
                dMin = -p[1]
                pMin = geoalgo.Vector(p[0],0,p[2])
            if ( (side == 2) and (direction == 1) ):
                #overflow on +z direction
                dMin = p[2]-1
                pMin = geoalgo.Vector(p[0],p[1],1)
            if ( (side == 2) and (direction == -1) ):
                #overflow on -z direction
                dMin = -p[2]
                pMin = geoalgo.Vector(p[0],p[1],0)
            answer = dMin*dMin
            tim = time()
            a1 = dAlgo.SqDist(p,b)
            sqdistT2 += (time()-tim)
            a2 = dAlgo.SqDist(b,p)
            # closest point
            tim = time()
            pt1 = dAlgo.ClosestPt(b,p)
            closestT2 += (time()-tim)
            pt2 = dAlgo.ClosestPt(p,b)
            #info("Point: [{0}, {1}, {2}]".format(p[0],p[1],p[2]))
            #info("Expected: {0}. Got: {1}".format(answer,a1))
            #info("Expected: {0}. Got: {1}".format(answer,a2))
            if not (np.abs(answer-a1) < _epsilon): success2 = 0
            if not (np.abs(answer-a2) < _epsilon): success2 = 0
            for x in range(3):
                #info("\tExpected: {0}. Got: {1}".format(p[x],pt1[x]))
                #info("\tExpected: {0}. Got: {1}".format(p[x],pt2[x]))
                if not (pt1[x]-pMin[x] < _epsilon) : success2 = 0
                if not (pt2[x]-pMin[x] < _epsilon) : success2 = 0
            #info("Success: {0}".format(success2))
            if (success2 == 1) : totSuccess2 += 1

        if ( float(totSuccess1)/tests < 0.999):
            info(NO + "Success: {0}%".format(100*float(totSuccess1)/tests) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess1)/tests) + ENDC)
        info("Time for SqDist (Case 1)              : {0:.3f} us".format(1E6*sqdistT1/tests))
        info("Time for ClosestPt (Case 1)           : {0:.3f} us".format(1E6*closestT1/tests))
        if ( float(totSuccess2)/tests < 0.999):
            info(NO + "Success: {0}%".format(100*float(totSuccess2)/tests) + ENDC)
        else:
            info(OK + "Success: {0}%".format(100*float(totSuccess2)/tests) + ENDC)
        info("Time for SqDist (Case 2)              : {0:.3f} us".format(1E6*sqdistT2/tests))
        info("Time for ClosestPt (Case 2)           : {0:.3f} us".format(1E6*closestT2/tests))
            
    except Exception:
        error('geoalgo::DistanceAlgo unit test failed.')
        print(traceback.format_exception(*sys.exc_info())[2])
        return 1
    
    info('geoalgo::DistanceAlgo unit test complete.')
    return 0

if __name__ == '__main__':
    test_dAlgo()

