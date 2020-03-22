#include "../interface/RazorPDFWeightsHelper.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <fstream>

RazorPDFWeightsHelper::RazorPDFWeightsHelper() {
  
  Eigen::initParallel();

}

void RazorPDFWeightsHelper::Init(unsigned int nreplicas, unsigned int neigenvectors, const edm::FileInPath &incsv) {  
  
  transformation_.resize(nreplicas,neigenvectors);
  
  std::ifstream instream(incsv.fullPath());
  if (!instream.is_open()) {
    throw cms::Exception("RazorPDFWeightsHelper")
      << "Could not open csv file" << incsv.relativePath(); 
  }
  
  
  for (unsigned int ireplica = 0; ireplica<nreplicas; ++ireplica) {
    std::string linestr;
    getline(instream,linestr);
    std::istringstream line(linestr);
    for (unsigned int ieigen = 0; ieigen<neigenvectors; ++ieigen) {
      std::string valstr;
      getline(line,valstr,',');
      std::istringstream val(valstr);
      val >> transformation_(ireplica,ieigen);
    }
    
  }
  
}

void RazorPDFWeightsHelper::DoMC2Hessian(double *inweights, double *outweights) const {
  const unsigned int nreplicas = transformation_.rows();
  const unsigned int neigenvectors = transformation_.cols();
  
  Eigen::VectorXd inweightv(nreplicas);
  for (unsigned int irep=0; irep<nreplicas; ++irep) {
    inweightv[irep] = inweights[irep] - 1.;
  }
  
  Eigen::VectorXd outweightv = transformation_.transpose()*inweightv;
  
  for (unsigned int ieig=0; ieig<neigenvectors; ++ieig) {
    outweights[ieig] = outweightv[ieig] + 1.;
  }  
    
}
