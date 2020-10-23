from test_msg import debug, info, error, warning
import traceback,sys
from test_import import test_import
test_import()
from ROOT import geoalgo

def test_vector():

    try:
        debug('Testing default constructor')
        k=geoalgo.Vector(5)
        if not k.size() == 5: raise Exception

        debug('Testing copy ctor')
        j=geoalgo.Vector(k)
        for x in range(k.size()):
            if not k[x] == j[x]: raise Exception

        for x in range(5):
            k[x] = 1
            j[x] = 1
        debug('Testing multiplication')
        k *= 2
        for x in range(5):
            if not k[x] == 2: raise Exception

        debug('Testing addition')
        k += j
        for x in range(5):
            if not k[x] == 3: raise Exception

        debug('Testing division')
        k /= 3.
        for x in range(5):
            if not k[x] == 1.: raise Exception

        debug('Testing dot product')
        if not k * j == 5: raise Exception

        debug('Testing length')
        for x in range(5):
            k[x] = 0
        k[0] = 1
        if not k.Length() == 1: raise Exception

        debug('Testing compatibility check')
        error=False
        try:
            k.compat(j)
        except Exception:
            pass
        if error: raise Exception

    except Exception:
        error('geoalgo::Vector unit test failed.')
        print(traceback.format_exception(*sys.exc_info())[2])
        return 1

    info('geoalgo::Vector unit test complete.')
    return 0

if __name__ == '__main__':
    test_vector()

