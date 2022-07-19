#include "cms_lpc_llp/llp_ntupler/interface/DBSCAN.h"
#include "TMath.h"
#include <iostream>
#include "TVector3.h"
#include "TGraph.h"
#include "TF1.h"
struct largest_nCsc_cluster_
{
  inline bool operator() (const cscCluster& c1, const cscCluster& c2){return c1.nCscSegments > c2.nCscSegments;}
} largest_nCsc_cluster;

struct hits
{
  float time;
  float error;
  bool strip;
};
const double theWireError_ = 8.6;
const double theStripError_ = 7.0;
const double thePruneCut_ = 9.0;

int DBSCAN::run()
{
    int clusterID = 1;
    vector<Point>::iterator iter;
    for(iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
        if ( iter->clusterID == UNCLASSIFIED )
        {
            if ( expandCluster(*iter, clusterID) != FAILURE )
            {
                clusterID += 1;
            }
        }
    }
    nClusters = clusterID-1;
    return (clusterID-1);
}
void DBSCAN::clear_clusters(){
  // nClusters = 0;
  clusterSize.clear();
  cscLabels.clear();
  clusterEta.clear();
  clusterPhi.clear();
  clusterX.clear();
  clusterY.clear();
  clusterZ.clear();
  clusterTime.clear();
  clusterTimeTotal.clear();
  clusterTimeTotal.clear();
  clusterMajorAxis.clear();
  clusterMinorAxis.clear();
  clusterXSpread.clear();
  clusterYSpread.clear();
  clusterXYSpread.clear();
  clusterRSpread.clear();
  clusterZSpread.clear();
  clusterTimeSpread.clear();
  clusterTimeSpreadWeighted.clear();
  clusterTimeSpreadWeightedAll.clear();

  clusterEtaPhiSpread.clear();
  clusterEtaSpread.clear();
  clusterPhiSpread.clear();
  clusterDeltaRSpread.clear();

  CscCluster.clear();
}
int DBSCAN::result(){

  // for (unsigned int i = 0;i < m_pointSize;i++)
  // {
  //   cscLabels.push_back(m_points[i].clusterID);
  //
  // }
  for(int i = 0; i < nClusters; i++)
  {
    float avg_x(0.0), avg_y(0.0), avg_z(0.0), avg_tWire(0.0), avg_tWirePruned(0.0), avg_t(0.0), avg_tTotal(0.0),tTotalSpreadPruned(0.0);
    float avg_x_sl2(0.0), avg_y_sl2(0.0), avg_z_sl2(0.0);
    float avg_eta(0.0), avg_phi(0.0);
    int size(0), size_z(0), size_xy(0);
    vector<float> wireTimes;
    vector<float> stripTimes;
    std::vector<hits> cscHits;

    vector<Point>::iterator iter;
    for(iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
      if ( iter->clusterID == i+1 )
      {
          if (iter->superlayer == 2) //for DT rechits that only have coordinates in Z
          {
            avg_x_sl2 += iter->x;
            avg_y_sl2 += iter->y;
            avg_z_sl2 += iter->z;
            size_z++;
          }
          else if (iter->superlayer == 1 || iter->superlayer == 3)
          {
            avg_x += iter->x;
            avg_y += iter->y;
            avg_z += iter->z;
            avg_t += iter->t;
            size_xy ++;
          }
          else //csc or for DT "wrong" rechit coordinates
          {
            avg_x += iter->x;
            avg_y += iter->y;
            avg_z += iter->z;
            avg_t += iter->t;
            avg_tWire += iter->twire;
            wireTimes.push_back(iter->twire);
            stripTimes.push_back(iter->t);
            hits thisHit;
            thisHit.time = iter->twire;
            thisHit.error = 1./(theWireError_*theWireError_);
	          thisHit.strip = false;
            cscHits.push_back(thisHit);
            thisHit.time = iter->t;
            thisHit.error = 1./(theStripError_*theStripError_);
	          thisHit.strip = true;
            cscHits.push_back(thisHit);

          }
          size ++;

      }
    }
    if (size_xy > 0 && size_z > 0) //for DT correct position, calculate average Z using sl2 and average XY using sl1/3
    {
      avg_x = avg_x/size_xy;
      avg_y = avg_y/size_xy;
      avg_z = avg_z_sl2/size_z;
    }
    else if (size_xy == 0 && size_z == 0) //csc or DT wrong position
    {
      avg_x = avg_x/size;
      avg_y = avg_y/size;
      avg_z = avg_z/size;
      // cout<<avg_x<<","<<avg_y<<","<<avg_z<<endl;
    }
    else if (size_xy > 0 && size_z == 0)
    {
      avg_x = avg_x/size_xy;
      avg_y = avg_y/size_xy;
      avg_z = avg_z/size_xy;

    }
    else
    {
      avg_x = avg_x_sl2/size_z;
      avg_y = avg_y_sl2/size_z;
      avg_z = avg_z_sl2/size_z;

    }
    avg_tTotal = (avg_t + avg_tWire)/(2 * size);
    avg_t = avg_t/size;
    avg_tWire = avg_tWire/size;

    // prune wire time
    //The wire times have a long tail that has to be pruned.  The strip times (tpeak) are fine
    bool modified = true;
    while (modified) {
      modified = false;
      double maxDiff = -1;
      std::vector<float>::iterator maxHit;
      for (std::vector<float>::iterator itWT = wireTimes.begin(); itWT != wireTimes.end(); ++itWT) {
        float diff = fabs(*itWT - avg_tTotal);
        if (diff > maxDiff) {
          maxDiff = diff;
          maxHit = itWT;
        }
      }
      if (maxDiff > 26) {
        int N = size + wireTimes.size();
        avg_tTotal = (avg_tTotal * N - (*maxHit)) / (N - 1);
        wireTimes.erase(maxHit);
        modified = true;
      }
    }

    //new timing calculation, error weighted
    // https://github.com/cms-sw/cmssw/blob/master/RecoMuon/MuonIdentification/src/CSCTimingExtractor.cc
    modified = false;
    double totalWeightTimeVtx = 0;
    double timeVtx = 0;
    double timeSpread = 0;
    do {
      modified = false;
      totalWeightTimeVtx = 0;
      timeVtx = 0;
      timeSpread = 0;
      for (std::vector<hits>::iterator it = cscHits.begin(); it != cscHits.end(); ++it) {
        timeVtx += it->time * it->error;
        totalWeightTimeVtx += it->error;
      }
      timeVtx /= totalWeightTimeVtx;

      // cut away outliers
      double diff_tvtx;
      double chimax = 0.0;
      int tmmax;
      for (unsigned int i = 0; i < cscHits.size(); i++) {
        diff_tvtx = (cscHits[i].time - timeVtx) * (cscHits[i].time - timeVtx) * cscHits[i].error;

        if (diff_tvtx > chimax) {
          tmmax =  i;
          chimax = diff_tvtx;
        }
      }
      // cut away the outliers
      if (chimax > thePruneCut_) {
        cscHits.erase(cscHits.begin()+tmmax);
        modified = true;
      }
    } while (modified);
    int count = 0;
    for (std::vector<hits>::iterator it = cscHits.begin(); it != cscHits.end(); ++it) {
      if (it->strip)
      {
        timeSpread += (it->time - timeVtx)*(it->time - timeVtx);
        count++;
      }

    }
    timeSpread = sqrt(timeSpread/count);


    // calculate cluster eta and phi
    avg_phi = atan(avg_y/avg_x);
    if  (avg_x < 0.0){
      avg_phi = TMath::Pi() + avg_phi;
    }
    avg_phi = deltaPhi(avg_phi,0.0);
    avg_eta = atan(sqrt(pow(avg_x,2)+pow(avg_y,2))/abs(avg_z));
    avg_eta = -1.0*TMath::Sign(1.0, avg_z)*log(tan(avg_eta/2));

    clusterEta.push_back(avg_eta);
    clusterPhi.push_back(avg_phi);
    clusterX.push_back(avg_x);
    clusterY.push_back(avg_y);
    clusterZ.push_back(avg_z);
    clusterTime.push_back(avg_t);
    clusterTimeTotal.push_back(avg_tTotal);
    clusterSize.push_back(size);

    clusterTimeWeighted.push_back(timeVtx);
    clusterTimeSpreadWeighted.push_back(timeSpread);

  }
  return 0;
}
int DBSCAN::clusterMoments()
{

  for(int i = 0; i < nClusters; i++)
  {
    float m11(0.0), m12(0.0), m22(0.0);
    float XSpread(0.0), YSpread(0.0), ZSpread(0.0), TSpread(0.0),  TSpreadAll(0.0), XYSpread(0.0), RSpread(0.0), DeltaRSpread(0.0);



    vector<Point>::iterator iter;
    for(iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
      if ( iter->clusterID == i+1 )
      {

          m11 += (iter->eta-clusterEta[i])*(iter->eta-clusterEta[i]);
          m12 += (iter->eta-clusterEta[i])* deltaPhi(iter->phi,clusterPhi[i]);
          m22 += deltaPhi(iter->phi,clusterPhi[i])*deltaPhi(iter->phi,clusterPhi[i]);
          DeltaRSpread +=  pow(deltaR(clusterEta[i], clusterPhi[i], iter->eta, iter->phi),2);
          XYSpread += (iter->x - clusterX[i])*(iter->y - clusterY[i]);
          XSpread += (iter->x - clusterX[i]) * (iter->x - clusterX[i]);
          YSpread += (iter->y - clusterY[i]) * (iter->y - clusterY[i]);
          ZSpread += (iter->z - clusterZ[i]) * (iter->z - clusterZ[i]);
          TSpread += (iter->t - clusterTime[i]) * (iter->t - clusterTime[i]);
          TSpreadAll += (iter->t - clusterTimeWeighted[i]) * (iter->t - clusterTimeWeighted[i]);
          float radius = sqrt(pow(iter->x, 2) + pow(iter->y, 2));
          RSpread += pow(radius-sqrt(clusterX[i]*clusterX[i]+clusterY[i]*clusterY[i]),2);


      }
    }

    float a = (m11+m22)/2;
    float b = 0.5*sqrt((m11+m22)*(m11+m22)-4*(m11*m22-m12*m12));
    clusterXSpread.push_back(sqrt(XSpread/(float)clusterSize[i]));
    clusterYSpread.push_back(sqrt(YSpread/(float)clusterSize[i]));
    clusterZSpread.push_back(sqrt(ZSpread/(float)clusterSize[i]));
    clusterRSpread.push_back(sqrt(RSpread/(float)clusterSize[i]));
    clusterDeltaRSpread.push_back(sqrt(DeltaRSpread/(float)clusterSize[i]));
    clusterXYSpread.push_back(sqrt(abs(XYSpread)/(float)clusterSize[i]));
    clusterTimeSpread.push_back(sqrt(TSpread/(float)clusterSize[i]));
    clusterTimeSpreadWeightedAll.push_back(sqrt(TSpreadAll/(float)clusterSize[i]));
    clusterEtaSpread.push_back(sqrt(abs(m11)/clusterSize[i]));
    clusterEtaPhiSpread.push_back(sqrt(abs(m12)/clusterSize[i]));
    clusterPhiSpread.push_back(sqrt(abs(m22)/clusterSize[i]));
    clusterMajorAxis.push_back(sqrt((a+b)/clusterSize[i]));
    clusterMinorAxis.push_back(sqrt((a-b)/clusterSize[i]));

  }

  return 0;
}
//input x,y,z of segments in cluster, avgX, avgY, avgZ
void DBSCAN::sort_clusters()
{
  for(int i = 0; i < nClusters; i++){
    vector<int>cscStations;
    vector<int>cscStations_copy;
    vector<int>cscChambers;
    vector<int>cscChambers_copy;
    vector<int>cscLayersPlus11;
    vector<int>cscLayersPlus12;
    vector<int>cscLayersPlus13;
    vector<int>cscLayersPlus21;
    vector<int>cscLayersPlus22;
    vector<int>cscLayersPlus31;
    vector<int>cscLayersPlus32;
    vector<int>cscLayersPlus41;
    vector<int>cscLayersPlus42;

    vector<int>cscLayersMinus11;
    vector<int>cscLayersMinus12;
    vector<int>cscLayersMinus13;
    vector<int>cscLayersMinus21;
    vector<int>cscLayersMinus22;
    vector<int>cscLayersMinus31;
    vector<int>cscLayersMinus32;
    vector<int>cscLayersMinus41;
    vector<int>cscLayersMinus42;

    vector<int>segment_index;
    int nSegments_Me11 = 0;
    int nSegments_Me12 = 0;

    cscCluster tmpCluster;
    tmpCluster.nCscSegmentChamberPlus11 = 0;
    tmpCluster.nCscSegmentChamberPlus12 = 0;
    tmpCluster.nCscSegmentChamberPlus13 = 0;
    tmpCluster.nCscSegmentChamberPlus21 = 0;
    tmpCluster.nCscSegmentChamberPlus22 = 0;
    tmpCluster.nCscSegmentChamberPlus31 = 0;
    tmpCluster.nCscSegmentChamberPlus32 = 0;
    tmpCluster.nCscSegmentChamberPlus41 = 0;
    tmpCluster.nCscSegmentChamberPlus42 = 0;
    tmpCluster.nCscSegmentChamberMinus11 = 0;
    tmpCluster.nCscSegmentChamberMinus12 = 0;
    tmpCluster.nCscSegmentChamberMinus13 = 0;
    tmpCluster.nCscSegmentChamberMinus21 = 0;
    tmpCluster.nCscSegmentChamberMinus22 = 0;
    tmpCluster.nCscSegmentChamberMinus31 = 0;
    tmpCluster.nCscSegmentChamberMinus32 = 0;
    tmpCluster.nCscSegmentChamberMinus41 = 0;
    tmpCluster.nCscSegmentChamberMinus42 = 0;

    tmpCluster.nDtSegmentStation1 = 0;
    tmpCluster.nDtSegmentStation2 = 0;
    tmpCluster.nDtSegmentStation3 = 0;
    tmpCluster.nDtSegmentStation4 = 0;
    for(unsigned int l=0; l < m_pointSize; l++){
      if (m_points[l].clusterID == i+1){
        cscStations.push_back(m_points[l].station);
        cscChambers.push_back(m_points[l].chamber);
        cscStations_copy.push_back(m_points[l].station);
        cscChambers_copy.push_back(m_points[l].chamber);
        segment_index.push_back(l);

        if (m_points[l].chamber == 11) cscLayersPlus11.push_back(m_points[l].layer);
        if (m_points[l].chamber == 12) cscLayersPlus12.push_back(m_points[l].layer);
        if (m_points[l].chamber == 13) cscLayersPlus13.push_back(m_points[l].layer);
        if (m_points[l].chamber == 21) cscLayersPlus21.push_back(m_points[l].layer);
      	if (m_points[l].chamber == 22) cscLayersPlus22.push_back(m_points[l].layer);
        if (m_points[l].chamber == 31) cscLayersPlus31.push_back(m_points[l].layer);
      	if (m_points[l].chamber == 32) cscLayersPlus32.push_back(m_points[l].layer);
        if (m_points[l].chamber == 41) cscLayersPlus41.push_back(m_points[l].layer);
      	if (m_points[l].chamber == 42) cscLayersPlus42.push_back(m_points[l].layer);

        if (m_points[l].chamber == -11) cscLayersMinus11.push_back(m_points[l].layer);
        if (m_points[l].chamber == -12) cscLayersMinus12.push_back(m_points[l].layer);
        if (m_points[l].chamber == -13) cscLayersMinus13.push_back(m_points[l].layer);
        if (m_points[l].chamber == -21) cscLayersMinus21.push_back(m_points[l].layer);
      	if (m_points[l].chamber == -22) cscLayersMinus22.push_back(m_points[l].layer);
        if (m_points[l].chamber == -31) cscLayersMinus31.push_back(m_points[l].layer);
      	if (m_points[l].chamber == -32) cscLayersMinus32.push_back(m_points[l].layer);
        if (m_points[l].chamber == -41) cscLayersMinus41.push_back(m_points[l].layer);
      	if (m_points[l].chamber == -42) cscLayersMinus42.push_back(m_points[l].layer);

        if (abs(m_points[l].chamber) == 11) nSegments_Me11++;
        if (abs(m_points[l].chamber) == 12) nSegments_Me12++;
      	if (m_points[l].chamber == 11) tmpCluster.nCscSegmentChamberPlus11++;
      	if (m_points[l].chamber == 12) tmpCluster.nCscSegmentChamberPlus12++;
      	if (m_points[l].chamber == 13) tmpCluster.nCscSegmentChamberPlus13++;
      	if (m_points[l].chamber == 21) tmpCluster.nCscSegmentChamberPlus21++;
      	if (m_points[l].chamber == 22) tmpCluster.nCscSegmentChamberPlus22++;
      	if (m_points[l].chamber == 31) tmpCluster.nCscSegmentChamberPlus31++;
      	if (m_points[l].chamber == 32) tmpCluster.nCscSegmentChamberPlus32++;
      	if (m_points[l].chamber == 41) tmpCluster.nCscSegmentChamberPlus41++;
      	if (m_points[l].chamber == 42) tmpCluster.nCscSegmentChamberPlus42++;
      	if (m_points[l].chamber == -11) tmpCluster.nCscSegmentChamberMinus11++;
      	if (m_points[l].chamber == -12) tmpCluster.nCscSegmentChamberMinus12++;
      	if (m_points[l].chamber == -13) tmpCluster.nCscSegmentChamberMinus13++;
      	if (m_points[l].chamber == -21) tmpCluster.nCscSegmentChamberMinus21++;
      	if (m_points[l].chamber == -22) tmpCluster.nCscSegmentChamberMinus22++;
      	if (m_points[l].chamber == -31) tmpCluster.nCscSegmentChamberMinus31++;
      	if (m_points[l].chamber == -32) tmpCluster.nCscSegmentChamberMinus32++;
      	if (m_points[l].chamber == -41) tmpCluster.nCscSegmentChamberMinus41++;
      	if (m_points[l].chamber == -42) tmpCluster.nCscSegmentChamberMinus42++;


        if (m_points[l].station == 1) tmpCluster.nDtSegmentStation1++;
        if (m_points[l].station == 2) tmpCluster.nDtSegmentStation2++;
        if (m_points[l].station == 3) tmpCluster.nDtSegmentStation3++;
        if (m_points[l].station == 4) tmpCluster.nDtSegmentStation4++;

      }
    }
    tmpCluster.x = clusterX[i];
    tmpCluster.y = clusterY[i];
    tmpCluster.z = clusterZ[i];
    tmpCluster.eta = clusterEta[i];
    tmpCluster.phi = clusterPhi[i];
    tmpCluster.t = clusterTime[i];
    tmpCluster.tWeighted = clusterTimeWeighted[i];
    tmpCluster.tTotal = clusterTimeTotal[i];

    tmpCluster.MajorAxis = clusterMajorAxis[i];
    tmpCluster.MinorAxis = clusterMinorAxis[i];
    tmpCluster.XSpread = clusterXSpread[i];
    tmpCluster.XYSpread = clusterXYSpread[i];
    tmpCluster.RSpread = clusterRSpread[i];


    tmpCluster.YSpread = clusterYSpread[i];
    tmpCluster.ZSpread = clusterZSpread[i];
    tmpCluster.TSpread = clusterTimeSpread[i];
    tmpCluster.TSpreadWeighted = clusterTimeSpreadWeighted[i];
    tmpCluster.TSpreadWeightedAll = clusterTimeSpreadWeightedAll[i];

    tmpCluster.EtaPhiSpread = clusterEtaPhiSpread[i];
    tmpCluster.EtaSpread = clusterEtaSpread[i];
    tmpCluster.DeltaRSpread = clusterDeltaRSpread[i];
    tmpCluster.PhiSpread = clusterPhiSpread[i];
    tmpCluster.nCscSegments = clusterSize[i];
    tmpCluster.Me11Ratio = 1.0*nSegments_Me11/clusterSize[i];
    tmpCluster.Me12Ratio = 1.0*nSegments_Me12/clusterSize[i];

    // count number of cscLayers
    std::sort(cscLayersPlus11.begin(), cscLayersPlus11.end());
    std::sort(cscLayersPlus12.begin(), cscLayersPlus12.end());
    std::sort(cscLayersPlus13.begin(), cscLayersPlus13.end());
    std::sort(cscLayersPlus21.begin(), cscLayersPlus21.end());
    std::sort(cscLayersPlus22.begin(), cscLayersPlus22.end());
    std::sort(cscLayersPlus31.begin(), cscLayersPlus31.end());
    std::sort(cscLayersPlus32.begin(), cscLayersPlus32.end());
    std::sort(cscLayersPlus41.begin(), cscLayersPlus41.end());
    std::sort(cscLayersPlus42.begin(), cscLayersPlus42.end());

    std::sort(cscLayersMinus11.begin(), cscLayersMinus11.end());
    std::sort(cscLayersMinus12.begin(), cscLayersMinus12.end());
    std::sort(cscLayersMinus13.begin(), cscLayersMinus13.end());
    std::sort(cscLayersMinus21.begin(), cscLayersMinus21.end());
    std::sort(cscLayersMinus22.begin(), cscLayersMinus22.end());
    std::sort(cscLayersMinus31.begin(), cscLayersMinus31.end());
    std::sort(cscLayersMinus32.begin(), cscLayersMinus32.end());
    std::sort(cscLayersMinus41.begin(), cscLayersMinus41.end());
    std::sort(cscLayersMinus42.begin(), cscLayersMinus42.end());


    tmpCluster.nLayersChamberPlus11 = std::unique(cscLayersPlus11.begin(), cscLayersPlus11.end())-cscLayersPlus11.begin();
    tmpCluster.nLayersChamberPlus12 = std::unique(cscLayersPlus12.begin(), cscLayersPlus12.end())-cscLayersPlus12.begin();
    tmpCluster.nLayersChamberPlus13 = std::unique(cscLayersPlus13.begin(), cscLayersPlus13.end())-cscLayersPlus13.begin();
    tmpCluster.nLayersChamberPlus21 = std::unique(cscLayersPlus21.begin(), cscLayersPlus21.end())-cscLayersPlus21.begin();
    tmpCluster.nLayersChamberPlus22 = std::unique(cscLayersPlus22.begin(), cscLayersPlus22.end())-cscLayersPlus22.begin();
    tmpCluster.nLayersChamberPlus31 = std::unique(cscLayersPlus31.begin(), cscLayersPlus31.end())-cscLayersPlus31.begin();
    tmpCluster.nLayersChamberPlus32 = std::unique(cscLayersPlus32.begin(), cscLayersPlus32.end())-cscLayersPlus32.begin();
    tmpCluster.nLayersChamberPlus41 = std::unique(cscLayersPlus41.begin(), cscLayersPlus41.end())-cscLayersPlus41.begin();
    tmpCluster.nLayersChamberPlus42 = std::unique(cscLayersPlus42.begin(), cscLayersPlus42.end())-cscLayersPlus42.begin();
    tmpCluster.nLayersChamberMinus11 = std::unique(cscLayersMinus11.begin(), cscLayersMinus11.end())-cscLayersMinus11.begin();
    tmpCluster.nLayersChamberMinus12 = std::unique(cscLayersMinus12.begin(), cscLayersMinus12.end())-cscLayersMinus12.begin();
    tmpCluster.nLayersChamberMinus13 = std::unique(cscLayersMinus13.begin(), cscLayersMinus13.end())-cscLayersMinus13.begin();
    tmpCluster.nLayersChamberMinus21 = std::unique(cscLayersMinus21.begin(), cscLayersMinus21.end())-cscLayersMinus21.begin();
    tmpCluster.nLayersChamberMinus22 = std::unique(cscLayersMinus22.begin(), cscLayersMinus22.end())-cscLayersMinus22.begin();
    tmpCluster.nLayersChamberMinus31 = std::unique(cscLayersMinus31.begin(), cscLayersMinus31.end())-cscLayersMinus31.begin();
    tmpCluster.nLayersChamberMinus32 = std::unique(cscLayersMinus32.begin(), cscLayersMinus32.end())-cscLayersMinus32.begin();
    tmpCluster.nLayersChamberMinus41 = std::unique(cscLayersMinus41.begin(), cscLayersMinus41.end())-cscLayersMinus41.begin();
    tmpCluster.nLayersChamberMinus42 = std::unique(cscLayersMinus42.begin(), cscLayersMinus42.end())-cscLayersMinus42.begin();
    tmpCluster.segment_id = segment_index;

    // count the number of chambers and max chamber segments
    std::vector<int>::iterator chamber_it;
    std::sort(cscChambers.begin(), cscChambers.end());
    chamber_it = std::unique(cscChambers.begin(), cscChambers.end());
    cscChambers.resize( std::distance(cscChambers.begin(),chamber_it) );
    int max_chamber = 999; // station with the maximum number of cscsegment in this cluster
    int max_chamber_segment = 0; // station with the maximum number of cscsegment in this cluster
    tmpCluster.nChamber = 0;
    for (unsigned int l = 0; l < cscChambers.size(); l++)
    {
      int counter = 0;
      for(unsigned int j = 0; j < cscChambers_copy.size(); j ++)
      {
        if (cscChambers_copy[j] == cscChambers[l]) counter++;
      }
      if (counter>max_chamber_segment)
      {
        max_chamber_segment = counter;
        max_chamber = cscChambers[l];
      }
      if(counter>5)tmpCluster.nChamber++;
    }
    tmpCluster.maxChamber = max_chamber;
    tmpCluster.maxChamberSegment = max_chamber_segment;
    // count the number of chambers and max chamber segments
    std::vector<int>::iterator station_it;
    std::sort(cscStations.begin(), cscStations.end());
    station_it = std::unique(cscStations.begin(), cscStations.end());
    cscStations.resize( std::distance(cscStations.begin(),station_it) );//list of unique stations
    int max_station = 999; // station with the maximum number of cscsegment in this cluster
    int max_station_segment = 0; // station with the maximum number of cscsegment in this cluster
    tmpCluster.nStation = 0;
    tmpCluster.nStation5 = 0;
    tmpCluster.nStation10 = 0;
    tmpCluster.nStation10perc = 0;
    tmpCluster.avgStation = 0.0;
    tmpCluster.avgStation5 = 0.0;
    tmpCluster.avgStation10 = 0.0;
    tmpCluster.avgStation10perc = 0.0;
    int nSeg10perc = 0;
    int nSeg5 = 0;
    int nSeg10 = 0;
    for (unsigned int l = 0; l < cscStations.size(); l++)
    {
      int counter = 0;
      for(unsigned int j = 0; j < cscStations_copy.size(); j ++)
      {
        if (cscStations_copy[j] == cscStations[l]) counter++;
      }
      if (counter>max_station_segment)
      {
        max_station_segment = counter;
        max_station = cscStations[l];
      }
      tmpCluster.nStation++;
      tmpCluster.avgStation += counter * cscStations[l];

      if(counter>=5.0){
        tmpCluster.avgStation5 += counter * cscStations[l];
        tmpCluster.nStation5++;
        nSeg5 += counter;

      }
      if(counter>=10.0){
        tmpCluster.avgStation10 += counter * cscStations[l];
        tmpCluster.nStation10++;
        nSeg10 += counter;
      }
      if(1.0*counter/clusterSize[i] > 0.1)
      {
        tmpCluster.avgStation10perc += counter * cscStations[l];
        tmpCluster.nStation10perc++;
        nSeg10perc += counter;
      }
    }
    tmpCluster.avgStation10perc = 1.0* tmpCluster.avgStation10perc/nSeg10perc;
    tmpCluster.avgStation5 = 1.0* tmpCluster.avgStation5/nSeg5;
    tmpCluster.avgStation10 = 1.0* tmpCluster.avgStation10/nSeg10;
    tmpCluster.avgStation = 1.0* tmpCluster.avgStation/clusterSize[i];

    tmpCluster.maxStation = max_station;
    tmpCluster.maxStationSegment = max_station_segment;

    CscCluster.push_back(tmpCluster);


  }

  //sort the clusters by size
  sort(CscCluster.begin(), CscCluster.end(), largest_nCsc_cluster);

}


void DBSCAN::merge_clusters()
{
  // clear all the cluster variables
  //change cluster ID of points
  bool modified = true;
  while(modified){
    modified = false;
    float mindR = 15;
    int cluster1 = 999;
    int cluster2 = 999;

    for(unsigned int i = 0; i < clusterEta.size(); i++){
      for(unsigned int j = i+1; j < clusterEta.size(); j++){
        float current_dR = deltaR(clusterEta[i], clusterPhi[i], clusterEta[j], clusterPhi[j]);
        if(current_dR<mindR)
        {
          mindR = current_dR;
          cluster1 = i;
          cluster2 = j;

        }
      }
    }
    if (mindR < MERGE_CLUSTER_DR){
      vector<Point>::iterator iter;
      int count = 0;
      for(iter = m_points.begin(); iter != m_points.end(); ++iter)
      {
        if ( iter->clusterID == cluster2+1 ){
          iter->clusterID = cluster1+1;
          count++;
        }
        if ( iter->clusterID > cluster2+1 )iter->clusterID = iter->clusterID-1;
      }
      clusterEta.erase(clusterEta.begin() + cluster2);
      clusterPhi.erase(clusterPhi.begin() + cluster2);
      nClusters--;
      modified = true;
      // can't use cscCluster, because its sorted, but the other vectors and clusterID are not sorted.
    }
  }
  clear_clusters(); // clear all the vectors, but nClusters is kept to keep track of the number of clusters.
}



int DBSCAN::expandCluster(Point point, int clusterID)
{
    vector<int> clusterSeeds = calculateCluster(point);//neighbors, including itself

    if ( clusterSeeds.size() < m_minPoints )
    {
        point.clusterID = NOISE;
        return FAILURE;
    }
    else//core points
    {
        int index = 0, indexCorePoint = 0;
        vector<int>::iterator iterSeeds;
        //loop through neighbors of core points
        for( iterSeeds = clusterSeeds.begin(); iterSeeds != clusterSeeds.end(); ++iterSeeds)
        {
            m_points.at(*iterSeeds).clusterID = clusterID;//setting all the neighbors to the same clusterID
            //get the index of the core point itself

            // if (m_points.at(*iterSeeds).x == point.x && m_points.at(*iterSeeds).y == point.y && m_points.at(*iterSeeds).z == point.z )
            if (m_points.at(*iterSeeds).eta == point.eta && m_points.at(*iterSeeds).phi == point.phi )
            {
                indexCorePoint = index;
            }
            ++index;
        }

        clusterSeeds.erase(clusterSeeds.begin()+indexCorePoint);
        //now clusterSeeds only contains neighbors, not itself; loop through the neighbors
        for( vector<int>::size_type i = 0, n = clusterSeeds.size(); i < n; ++i )
        {
            vector<int> clusterNeighors = calculateCluster(m_points.at(clusterSeeds[i]));
            //if this neighbor point is a core point
            if ( clusterNeighors.size() >= m_minPoints )
            {
                vector<int>::iterator iterNeighors;
                for ( iterNeighors = clusterNeighors.begin(); iterNeighors != clusterNeighors.end(); ++iterNeighors )
                {
                    if ( m_points.at(*iterNeighors).clusterID == UNCLASSIFIED || m_points.at(*iterNeighors).clusterID == NOISE )
                    {
                        if ( m_points.at(*iterNeighors).clusterID == UNCLASSIFIED )
                        {
                            clusterSeeds.push_back(*iterNeighors);
                            n = clusterSeeds.size();
                        }
                        m_points.at(*iterNeighors).clusterID = clusterID;
                    }
                }
            }
        }

        return SUCCESS;
    }
}

vector<int> DBSCAN::calculateCluster(Point point)
{
    int index = 0;
    vector<Point>::iterator iter;
    vector<int> clusterIndex;
    Point minimum;
    float min_dis = 1000000.0;
    for( iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
        if (min_dis > calculateDistance(point, *iter)) minimum  = *iter;
        if ( calculateDistance(point, *iter) <= m_epsilon )
        {
            clusterIndex.push_back(index);
        }
        index++;
    }
    // if (index > 10) std::cout<<calculateDistance(point, minimum)<<std::endl;
    return clusterIndex;
}

inline double DBSCAN::calculateDistance( Point pointCore, Point pointTarget )
{
    // return sqrt(pow(pointCore.x - pointTarget.x,2)+pow(pointCore.y - pointTarget.y,2)+pow(pointCore.z - pointTarget.z,2));
    // return sqrt(pow(pointCore.eta - pointTarget.eta,2)+pow(deltaPhi(pointCore.phi, pointTarget.phi),2)+pow((pointCore.t - pointTarget.t)/100.0,2));
    return sqrt(pow(pointCore.eta - pointTarget.eta,2)+pow(deltaPhi(pointCore.phi, pointTarget.phi),2));

}
double DBSCAN::deltaPhi(double phi1, double phi2)
{
  double dphi = phi1-phi2;
  while (dphi > TMath::Pi())
  {
    dphi -= TMath::TwoPi();
  }
  while (dphi <= -TMath::Pi())
  {
    dphi += TMath::TwoPi();
  }
  return dphi;
};
double DBSCAN::deltaR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = deltaPhi(phi1,phi2);
  double deta = eta1 - eta2;
  return sqrt( dphi*dphi + deta*deta);
}
