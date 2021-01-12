void MakeInput() {

  std::string input = "cd106.txt"; //Cygnus nucleus file, must already exist

  //New files will be created with these names
  std::string prob_output = "cd106_on_ti48.prb"; //contains excitation probabilities
  std::string tens_output = "cd106_on_ti48.ten"; //contains statistical tensors

  NucleusReader* reader = new NucleusReader(input.c_str());
  Nucleus* nuc = reader->GetNucleus();

  //Change these
  int beam_A = 106;
  int beam_Z = 48;
  int target_A = 48;
  int target_Z = 22;
  double beam_energy = 294.0; //MeV
  
  Reaction* reac = new Reaction(beam_A,beam_Z,target_A,target_Z,beam_energy);
  PointCoulEx* poin = new PointCoulEx(nuc,reac);

  //Uncomment these two lines for target excitation
  //poin->SetProjectileExcitation(false);
  //poin->PrepareConnections();

  std::ofstream file;
  file.open(prob_output.c_str(),std::ios::out);

  const int nBins = 25; //this can be changed if desired
  const int nStates = nuc->GetNstates();
  for(int i=0;i<nBins-1;i++) {
    
    double thetaCM = (TMath::Pi()/(double)nBins)*(i+1)*TMath::RadToDeg();
    
    std::cout << "Point " << i << ": ThetaCM = " << thetaCM << " deg" << std::endl;

    poin->CalculatePointProbabilities(thetaCM);
    
    file << thetaCM*TMath::DegToRad();
    for(int j=0;j<nStates;j++) {
      file << " " << poin->GetProbabilitiesVector()(j);
    }
    file << "\n";

    poin->CalculateTensors();
    
    if(!i) {
      poin->WriteTensorsToFile(tens_output.c_str());
    }
    else {
      poin->WriteTensorsToFile(tens_output.c_str(),std::ios_base::app);
    }
    
  }
  file.close();
  
  return;
}
