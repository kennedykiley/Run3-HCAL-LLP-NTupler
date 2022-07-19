#ifndef DBSCAN_H
#define DBSCAN_H

#include <vector>
#include <cmath>
using namespace std;

#define UNCLASSIFIED -1
// #define CORE_POINT 1
// #define BORDER_POINT 2
#define NOISE -2
#define SUCCESS 0
#define FAILURE -3
// #define MERGE_DR
const double MERGE_CLUSTER_DR = 0.6;

struct cscCluster
{
  float x, y, z, t, tTotal, tWeighted, eta, phi;//t is t_strip, tWire, tTotal is sum
  int nCscSegments;
  float jetVeto, calojetVeto, muonVeto;
  int maxChamber, maxChamberSegment, nChamber;
  // vector<int> cscChambers;
  int maxStation, maxStationSegment, nStation, nStation5, nStation10, nStation10perc;
  float avgStation, avgStation5, avgStation10, avgStation10perc;
  int nCscSegmentChamberPlus11, nCscSegmentChamberPlus12, nCscSegmentChamberPlus13, nCscSegmentChamberPlus21, nCscSegmentChamberPlus22, nCscSegmentChamberPlus31, nCscSegmentChamberPlus32, nCscSegmentChamberPlus41, nCscSegmentChamberPlus42;
  int nCscSegmentChamberMinus11, nCscSegmentChamberMinus12, nCscSegmentChamberMinus13, nCscSegmentChamberMinus21, nCscSegmentChamberMinus22, nCscSegmentChamberMinus31, nCscSegmentChamberMinus32, nCscSegmentChamberMinus41, nCscSegmentChamberMinus42;

  int nLayersChamberPlus11, nLayersChamberPlus12, nLayersChamberPlus13, nLayersChamberPlus21, nLayersChamberPlus22, nLayersChamberPlus31, nLayersChamberPlus32, nLayersChamberPlus41, nLayersChamberPlus42;
  int nLayersChamberMinus11, nLayersChamberMinus12, nLayersChamberMinus13, nLayersChamberMinus21, nLayersChamberMinus22, nLayersChamberMinus31, nLayersChamberMinus32, nLayersChamberMinus41, nLayersChamberMinus42;
  int nDtSegmentStation1,nDtSegmentStation2,nDtSegmentStation3,nDtSegmentStation4;


  float Me11Ratio, Me12Ratio;
  float MajorAxis, MinorAxis, EtaSpread, PhiSpread, EtaPhiSpread, XYSpread, DeltaRSpread;
  float XSpread, YSpread, ZSpread, TSpread, RSpread, TSpreadWeighted, TSpreadWeightedAll;


  vector<int>segment_id;
};

// struct largest_nCsc_cluster_
// {
//   inline bool operator() (const cscCluster& c1, const cscCluster& c2){return c1.nCscSegments > c2.nCscSegments;}
// } largest_nCsc_cluster;

typedef struct Point_
{
    float x, y, z, t, twire;  // X, Y, Z position
    float eta,phi;
    float dirX, dirY, dirZ;
    int station, chamber, layer, superlayer; //superlayer exists only for DT
    int clusterID;  // clustered ID
}Point;

class DBSCAN {
public:
    DBSCAN(unsigned int minPts, float eps, vector<Point> points){
        m_minPoints = minPts;
        m_epsilon = eps;
        m_points = points;
        m_pointSize = points.size();
    }
    ~DBSCAN(){}
    int nClusters;
    vector<int>clusterSize;
    vector<int>cscLabels;
    vector<float>clusterEta;
    vector<float>clusterPhi;
    vector<float>clusterX;
    vector<float>clusterY;
    vector<float>clusterZ;
    vector<float>clusterTime;
    vector<float>clusterTimeWeighted;
    vector<float>clusterTimeTotal;
    vector<float>clusterMajorAxis;
    vector<float>clusterMinorAxis;
    vector<float>clusterXSpread;
    vector<float>clusterYSpread;
    vector<float>clusterXYSpread;
    vector<float>clusterRSpread;

    vector<float>clusterZSpread;
    vector<float>clusterTimeSpread;
    vector<float>clusterTimeSpreadWeighted;
    vector<float>clusterTimeSpreadWeightedAll;
    vector<float>clusterEtaPhiSpread;
    vector<float>clusterEtaSpread;
    vector<float>clusterPhiSpread;
    vector<float>clusterDeltaRSpread;




    vector<cscCluster> CscCluster;

    int run();
    double deltaPhi(double phi1, double phi2);
    double deltaR(double eta1, double phi1, double eta2, double phi2);

    int result();
    int clusterMoments();
    void clear_clusters();
    void sort_clusters();
    void merge_clusters();//change cluster id of the points if the clusters are close, and clear all the vectors at the end.
    vector<int> calculateCluster(Point point);
    int expandCluster(Point point, int clusterID);
    inline double calculateDistance(Point pointCore, Point pointTarget);

    int getTotalPointSize() {return m_pointSize;}
    int getMinimumClusterSize() {return m_minPoints;}
    int getEpsilonSize() {return m_epsilon;}
private:
    vector<Point> m_points;
    unsigned int m_pointSize;
    unsigned int m_minPoints;
    float m_epsilon;
};

#endif // DBSCAN_H
