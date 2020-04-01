#ifndef DBSCAN_H
#define DBSCAN_H

#include <vector>
#include <cmath>

#define UNCLASSIFIED -1
// #define CORE_POINT 1
// #define BORDER_POINT 2
#define NOISE -2
#define SUCCESS 0
#define FAILURE -3

using namespace std;
struct cscCluster
{
  float x, y, z, t, tTotal, eta, phi;
  int nCscSegments;
  float jetVeto, calojetVeto, muonVeto;
  int maxChamber, maxChamberSegment, nChamber;
  int maxStation, maxStationSegment, nStation;
  int nCscSegmentChamberPlus11, nCscSegmentChamberPlus12, nCscSegmentChamberPlus13, nCscSegmentChamberPlus21, nCscSegmentChamberPlus22, nCscSegmentChamberPlus31, nCscSegmentChamberPlus32, nCscSegmentChamberPlus41, nCscSegmentChamberPlus42;
  int nCscSegmentChamberMinus11, nCscSegmentChamberMinus12, nCscSegmentChamberMinus13, nCscSegmentChamberMinus21, nCscSegmentChamberMinus22, nCscSegmentChamberMinus31, nCscSegmentChamberMinus32, nCscSegmentChamberMinus41, nCscSegmentChamberMinus42;
  int nDtSegmentStation1,nDtSegmentStation2,nDtSegmentStation3,nDtSegmentStation4;

  float Me11Ratio, Me12Ratio;
  float MajorAxis, MinorAxis, EtaSpread, PhiSpread, EtaPhiSpread;
  float XSpread, YSpread, XYSpread, ZSpread, TSpread;
  float vertex_r, vertex_z, vertex_dis, vertex_chi2;
  int vertex_n, vertex_n1, vertex_n5, vertex_n10, vertex_n15, vertex_n20;
  vector<int>segment_id;
};

// struct largest_nCsc_cluster_
// {
//   inline bool operator() (const cscCluster& c1, const cscCluster& c2){return c1.nCscSegments > c2.nCscSegments;}
// } largest_nCsc_cluster;

typedef struct Point_
{
    float x, y, z, t;  // X, Y, Z position
    float eta,phi;
    float dirX, dirY, dirZ;
    int station, chamber;
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
    vector<float>clusterTimeTotal;
    vector<float>clusterMajorAxis;
    vector<float>clusterMinorAxis;
    vector<float>clusterXSpread;
    vector<float>clusterYSpread;
    vector<float>clusterXYSpread;
    vector<float>clusterZSpread;
    vector<float>clusterTimeSpread;
    vector<float>clusterEtaPhiSpread;
    vector<float>clusterEtaSpread;
    vector<float>clusterPhiSpread;

    vector<float>clusterVertexR;
    vector<float>clusterVertexZ;
    vector<float>clusterVertexDis;
    vector<float>clusterVertexChi2;
    vector<int>clusterVertexN;

    vector<int>clusterVertexN1cm;
    vector<int>clusterVertexN5cm;
    vector<int>clusterVertexN10cm;
    vector<int>clusterVertexN15cm;
    vector<int>clusterVertexN20cm;


    vector<cscCluster> CscCluster;

    int run();
    double deltaPhi(double phi1, double phi2);
    int result();
    int clusterMoments();
    void sort_clusters();
    int vertexing();
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
