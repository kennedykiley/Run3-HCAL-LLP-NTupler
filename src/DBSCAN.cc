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

int DBSCAN::result(){

  // for (unsigned int i = 0;i < m_pointSize;i++)
  // {
  //   cscLabels.push_back(m_points[i].clusterID);
  //
  // }
  for(int i = 0; i < nClusters; i++)
  {
    float avg_x(0.0), avg_y(0.0), avg_z(0.0), avg_t(0.0);
    float avg_eta(0.0), avg_phi(0.0);
    int size(0);
    vector<Point>::iterator iter;
    for(iter = m_points.begin(); iter != m_points.end(); ++iter)
    {

      if ( iter->clusterID == i+1 )
      {

          avg_x += iter->x;
          avg_y += iter->y;
          avg_z += iter->z;
          avg_t += iter->t;
          size ++;
      }
    }
    avg_x = avg_x/size;
    avg_y = avg_y/size;
    avg_z = avg_z/size;
    avg_t = avg_t/size;
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
    clusterSize.push_back(size);


  }
  return 0;
}
int DBSCAN::clusterMoments()
{


  for(int i = 0; i < nClusters; i++)
  {
    float m11(0.0), m12(0.0), m22(0.0);
    float XSpread(0.0), YSpread(0.0), ZSpread(0.0), TSpread(0.0);
    vector<Point>::iterator iter;
    for(iter = m_points.begin(); iter != m_points.end(); ++iter)
    {


      if ( iter->clusterID == i+1 )
      {
          // float phi = atan(iter->y/iter->x);
          // if  (iter->x < 0.0){
          //   phi = TMath::Pi() + phi;
          // }
          // phi = deltaPhi(phi,0.0);
          // float eta = atan(sqrt(pow(iter->x,2)+pow(iter->y,2))/abs(iter->z));
          // eta = -1.0*TMath::Sign(1.0, iter->z)*log(tan(eta/2));
          m11 += (iter->eta-clusterEta[i])*(iter->eta-clusterEta[i]);
          m12 += (iter->eta-clusterEta[i])* deltaPhi(iter->phi,clusterPhi[i]);
          m22 += deltaPhi(iter->phi,clusterPhi[i])*deltaPhi(iter->phi,clusterPhi[i]);
          XSpread += (iter->x - clusterX[i]) * (iter->x - clusterX[i]);
          YSpread += (iter->y - clusterY[i]) * (iter->y - clusterY[i]);
          ZSpread += (iter->z - clusterZ[i]) * (iter->z - clusterZ[i]);
          TSpread += (iter->t - clusterTime[i]) * (iter->t - clusterTime[i]);

      }
    }
    float a = (m11+m22)/2;
    float b = 0.5*sqrt((m11+m22)*(m11+m22)-4*(m11*m22-m12*m12));
    clusterXSpread.push_back(sqrt(XSpread/(float)clusterSize[i]));
    clusterYSpread.push_back(sqrt(XSpread/(float)clusterSize[i]));
    clusterZSpread.push_back(sqrt(XSpread/(float)clusterSize[i]));
    clusterTimeSpread.push_back(sqrt(TSpread/(float)clusterSize[i]));
    clusterEtaSpread.push_back(sqrt(m11/clusterSize[i]));
    clusterEtaPhiSpread.push_back(sqrt(abs(m12)/clusterSize[i]));
    clusterPhiSpread.push_back(sqrt(m22/clusterSize[i]));
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
    for(unsigned int l=0; l < m_pointSize; l++){
      if (m_points[l].clusterID == i+1){
        cscStations.push_back(m_points[l].station);
        cscChambers.push_back(m_points[l].chamber);
        cscStations_copy.push_back(m_points[l].station);
        cscChambers_copy.push_back(m_points[l].chamber);
        segment_index.push_back(l);
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
      }
    }
    tmpCluster.x = clusterX[i];
    tmpCluster.y = clusterY[i];
    tmpCluster.z = clusterZ[i];
    tmpCluster.eta = clusterEta[i];
    tmpCluster.phi = clusterPhi[i];
    tmpCluster.t = clusterTime[i];
    tmpCluster.MajorAxis = clusterMajorAxis[i];
    tmpCluster.MinorAxis = clusterMinorAxis[i];
    tmpCluster.XSpread = clusterXSpread[i];
    tmpCluster.YSpread = clusterYSpread[i];
    tmpCluster.ZSpread = clusterZSpread[i];
    tmpCluster.TSpread = clusterTimeSpread[i];
    tmpCluster.EtaPhiSpread = clusterEtaPhiSpread[i];
    tmpCluster.EtaSpread = clusterEtaSpread[i];
    tmpCluster.PhiSpread = clusterPhiSpread[i];
    tmpCluster.nCscSegments = clusterSize[i];
    tmpCluster.Me11Ratio = 1.0*nSegments_Me11/clusterSize[i];
    tmpCluster.Me12Ratio = 1.0*nSegments_Me12/clusterSize[i];
    tmpCluster.vertex_r = clusterVertexR[i];
    tmpCluster.vertex_z = clusterVertexZ[i];
    tmpCluster.vertex_dis = clusterVertexDis[i];
    tmpCluster.vertex_chi2 = clusterVertexChi2[i];
    tmpCluster.vertex_n = clusterVertexN[i];
    tmpCluster.vertex_n1 = clusterVertexN1cm[i];
    tmpCluster.vertex_n5 = clusterVertexN5cm[i];
    tmpCluster.vertex_n10 = clusterVertexN10cm[i];
    tmpCluster.vertex_n15 = clusterVertexN15cm[i];
    tmpCluster.vertex_n20 = clusterVertexN20cm[i];
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
    cscStations.resize( std::distance(cscStations.begin(),station_it) );
    int max_station = 999; // station with the maximum number of cscsegment in this cluster
    int max_station_segment = 0; // station with the maximum number of cscsegment in this cluster
    tmpCluster.nStation = 0;
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
      if(counter>5)tmpCluster.nStation++;
    }
    tmpCluster.maxStation = max_station;
    tmpCluster.maxStationSegment = max_station_segment;

    CscCluster.push_back(tmpCluster);


  }

  //sort the clusters by size
  sort(CscCluster.begin(), CscCluster.end(), largest_nCsc_cluster);

}

int DBSCAN::vertexing()
{
  for(int i = 0; i < nClusters; i++){
    TVector3 vecDir;
    TVector3 vecCsc;
    TF1 *fit = new TF1();

    TGraph *gr = new TGraph(clusterSize[i]);

    for(unsigned int j = 0;j < m_pointSize;j++)
    {
      if ( m_points[j].clusterID == i+1 ){
        vecCsc.SetXYZ(m_points[j].x, m_points[j].y, m_points[j].z);
        vecDir.SetXYZ(m_points[j].dirX, m_points[j].dirY, m_points[j].dirZ);
        double slope = (vecCsc.Perp() - (vecCsc+vecDir).Perp()) / (vecCsc.Z() - (vecCsc+vecDir).Z());
        double beta = -1.0*vecCsc.Z() * slope + vecCsc.Perp();
        gr->SetPoint(j,slope,beta);
      }

    }
    // Here's where the fit happens
    // Process is repeated until a good vertex is found or fewer than 3 segments left
    // The segments form a line in RZ plane
    // Ideally, the vertex (z,r) is a point on every line in the cluster
    // So using all the measured segment slopes (a_i) and a vertex guess (z,r) gives predicted intercept (b_i)
    // r = a_i*z + b_i --> b_i = r - a_i*z
    // Fitting (a_i, b_i) points to "b_i = p[0] + p[1]*a_i" by minimizing chi-squared between measured andß
    // predicted intercepts gives best (r,z) estimate
    // Then calculate distance of closest approach of all segment lines and vertex
    // If max distance > 30cm throw out farthest point and refit until all segments within 30cm or < 3 segments left

    bool goodVtx = false;
    double distance = 0.0;
    double maxDistance = 0.0;
    double farSegment = -1;
    int vertexN1cm(0), vertexN5cm(0), vertexN10cm(0), vertexN15cm(0), vertexN20cm(0), vertexN(0);
    float vertexChi2(0.0), vertexDis(0.0), vertexZ(0.0), vertexR(0.0);


    // cout << "starting with " << gr->GetN() << " segments" << endl;
    while (!goodVtx && gr->GetN()>=3){
    	maxDistance = 0.0;
      // gr->Fit("pol 1");
      // fit = gr->GetFunction("pol 1");
    	gr->Fit("1++x","Q");
    	fit = gr->GetFunction("1++x");
    	for (int j=0; j<gr->GetN(); j++){
    	  // Distance from point (x_0, y_0) to line ax + by + c = 0 is |ax_0 + by_0 + c| / sqrt(a^2 + b^2)
    	  // Here: a = slope = "X"; b = -1,;c = intercept = "Y"; (x_0, y_0) = (z, r) = (-p[1], p[0]);
    	  distance = abs(gr->GetX()[j]*fit->GetParameter(1)*-1.0 + -1.0*fit->GetParameter(0) + gr->GetY()[j]);
    	  distance = distance / sqrt(pow(gr->GetX()[j],2)+pow(-1.0,2));
    	  if (distance > maxDistance){
    	    maxDistance = distance;
    	    farSegment = j;
    	  }
    	}
    	if (maxDistance < 30.0){
    	  goodVtx = true;
    	}
    	else {
    	  gr->RemovePoint(farSegment);
    	}
    }
    // count number of segment in distance 1,5,10,20cm
    for (int j=0; j<gr->GetN(); j++){
      // Distance from point (x_0, y_0) to line ax + by + c = 0 is |ax_0 + by_0 + c| / sqrt(a^2 + b^2)
      // Here: a = slope = "X"; b = -1,;c = intercept = "Y"; (x_0, y_0) = (z, r) = (-p[1], p[0]);
      distance = abs(gr->GetX()[j]*fit->GetParameter(1)*-1.0 + -1.0*fit->GetParameter(0) + gr->GetY()[j]);
      distance = distance / sqrt(pow(gr->GetX()[j],2)+pow(-1.0,2));
      if (distance < 1.0) vertexN1cm ++;
      if (distance < 5.0) vertexN5cm ++;
      if (distance < 10.0) vertexN10cm ++;
      if (distance < 15.0) vertexN15cm ++;
      if (distance < 20.0) vertexN20cm ++;

    }
    if (goodVtx && gr->GetN()>=3){
      vertexR = fit->GetParameter(0);
      vertexZ = -1.0*fit->GetParameter(1);
      vertexN = gr->GetN();
      vertexDis = maxDistance;
      vertexChi2 = fit->GetChisquare() ;
    }
    else{
      vertexR = 0.0;
      vertexZ = 0.0;
      vertexN = 0;
      vertexDis =-999;
      vertexChi2 = -999.;

    }
    clusterVertexR.push_back(vertexR);
    clusterVertexZ.push_back(vertexZ);
    clusterVertexN.push_back(vertexN);
    clusterVertexDis.push_back(vertexDis);
    clusterVertexChi2.push_back(vertexChi2);

    clusterVertexN1cm.push_back(vertexN1cm);
    clusterVertexN5cm.push_back(vertexN5cm);
    clusterVertexN10cm.push_back(vertexN10cm);
    clusterVertexN15cm.push_back(vertexN15cm);
    clusterVertexN20cm.push_back(vertexN20cm);
  }
  // cout << "vertex? " << goodVtx << ", # segments = " << gr->GetN() << ", maxDistance = " << maxDistance << endl;
  // cout << "vertex R, Z, N:  " << clusterVertexR <<", " << clusterVertexZ << ",  " << clusterVertexN << endl;
  // cout << "maxDistance = " << maxDistance << endl;

  return 0;

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
