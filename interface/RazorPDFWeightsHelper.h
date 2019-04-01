#ifndef PhysicsTools_HepMCCandAlgos_RazorPDFWeightsHelper_h
#define PhysicsTools_HepMCCandAlgos_RazorPDFWeightsHelper_h

#include <Eigen/Dense>

#include <iostream>

#include "FWCore/ParameterSet/interface/FileInPath.h"

class RazorPDFWeightsHelper {
  
public:
  
  RazorPDFWeightsHelper();
  
  void Init(unsigned int nreplicas, unsigned int neigenvectors, const edm::FileInPath &incsv);
  void DoMC2Hessian(double *inweights, double *outweights) const;
  
  unsigned int neigenvectors() const { return transformation_.cols(); }
  
protected:
  
  Eigen::MatrixXd transformation_;  

};
#endif
