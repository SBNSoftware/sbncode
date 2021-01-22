from test_msg import debug, info, error, warning

def test_import():

    try:
        debug('Attempting to import geoalgo')
        import ROOT
        ROOT.gSystem.Load("libflashmatch")
        from ROOT import geoalgo    
    except Exception:
        error('Import geoalgo unit test failed.')
        return 1
    info('geoalgo import succeeded.')
    return 0

if __name__ == '__main__':
    test_Import()

