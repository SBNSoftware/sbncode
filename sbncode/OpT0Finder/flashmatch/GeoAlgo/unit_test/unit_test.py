import sys
import test_msg
from test_msg import debug, info, error, warning

print('\033[93mExecuting a unit test for geoalgo Package\033[00m')

# message level: [0,1,2,3] = [debug,info,warning,error]
MSG_LEVEL=1

test_msg.test_msg.level = MSG_LEVEL

#
# Test import
#
from test_import import test_import
if test_import(): sys.exit(1)

#
# Test Vector
#
from test_vector import test_vector
if test_vector(): sys.exit(1)

#
# Test DistanceAlgo
#
from test_distance import test_dAlgo
if test_dAlgo(): sys.exit(1)
