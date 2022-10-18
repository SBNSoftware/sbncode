
find_package(sbnanaobj REQUIRED) 

add_subdirectory(SinglePhotonAnalysis)
add_subdirectory(SBNEventWeight)
add_subdirectory(NuMuSelection)
add_subdirectory(OpDet)
add_subdirectory(OpT0Finder)

add_subdirectory(CAFMaker)
add_subdirectory(Metadata)

add_subdirectory(Calibration)
add_subdirectory(FlashMatch)
add_subdirectory(LArRecoProducer)
add_subdirectory(TPCReco)


add_subdirectory(BeamSpillInfoRetriever)
add_subdirectory(FluxReader)
add_subdirectory(EventGenerator)
add_subdirectory(PID)
add_subdirectory(GeometryTools)
add_subdirectory(CosmicID)
add_subdirectory(LArG4)
