#ifndef _NeutronEDep_h_
#define _NeutronEDep_h_

std::vector< std::pair<double, int> > NeutronEnergyDepForG4Reweight(G4ReweightTraj * theTraj, double res, double mass, bool isElastic){

  std::vector< std::pair<double, int> > result;
  size_t nSteps = theTraj->GetNSteps();

  for(size_t is = 0; is < nSteps; ++is){
      
    auto theStep = theTraj->GetStep(is);

    double Px = theTraj->GetStep(is)->GetPreStepPx();
    double Py = theTraj->GetStep(is)->GetPreStepPy();
    double Pz = theTraj->GetStep(is)->GetPreStepPz();

    // Energy at start of step
    double stepBeginEnergy = sqrt(Px*Px + Py*Py + Pz*Pz + mass*mass);

    std::string theProc = theStep->GetStepChosenProc();
               
    // Assume interaction occurs in last slice
    int slices = round(theStep->GetStepLength()/res);

    for(int i_slice=0;i_slice<slices-1;i_slice++) result.push_back(std::make_pair(stepBeginEnergy,0));

    if ((!isElastic && theProc.find(std::string("Inelastic")) != std::string::npos) ||
	(isElastic && theProc.find(std::string("hadElastic")) != std::string::npos)) {       
      result.push_back(std::make_pair(stepBeginEnergy,1));
    }
    else result.push_back(std::make_pair(stepBeginEnergy,0)); 
  }

  return result;
}

#endif
