#ifndef _PATHS_H
#define _PATHS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "Param.h"
#include "Distribution.h"
#include "bid.h"
#include "BidSet.h"
//#include "normal.h"

#define max(x,y) (((x) > (y)) ? (x) : (y))
#define sqr(x) ((x)*(x))

class Paths : public Distribution {
public:
    // constructor
    Paths(Param &p, char *argv[], int argc);

    // destructor
    virtual ~Paths();

    static void usage();
    BidSet *generate(Param &p);
    void randomizeParams(Param &p);

protected:
    // constants
    static const int Z_MAX_BID_SET_SIZE; // = 100
    static const int Z_NUMBIDS; // = 20000
    static const double INF; // = 500000.0

    // parameters
    int maxCities;
    int initialNumCities;
    int initial_connections;
    double edge_density;
    double building_penalty;
    double shipping_cost_factor;
    int max_bid_set_size;

    //   FOR NORMAL CITY PLACEMENT
    //   CityPlacementType city_placement_type;
    //   double cluster_std_dev;
    //   int num_clusters;

    bool keepDisconnected;
    bool countUnusedEdges;
    double goodError;
    bool display_graph;
    int numcities;

    // local structs
    struct DistanceTo {
        int city;
        double distance;
    };

    struct Path {
        int * pathArr, pathlen;
        double cost;

        Path(int maxCities) {
            pathArr = new int[maxCities];
        }

        Path(const Path &orig, int maxCities) {
            pathlen = orig.pathlen;
            cost = orig.cost;
            pathArr = new int[maxCities];

            for (int i = 0; i < pathlen; i++) {
                pathArr[i] = orig.pathArr[i];
            }
        }

        ~Path() {
            delete[] pathArr;
        }
    };

    char *output_buffer;

    struct Point {
        double x, y;
    };

    Point * cities;
    double ** edges;
    // edges is a matrix representing the distance between cities
    // 0 means unknown, neg means known, but no direct connection
    // pos means known with direct connection

    int ** edgesgood;
    int *goodedges1;
    int *goodedges2;
    int **adjacentTo;
    int *numAdjacentTo;
    int *q, *pi;
    double *d;

    Path **best_paths;
    int num_best_paths;
    double maxcost;
    BidSet *bids;
    int num_edges;

    //  Normal *norm;
    //  helper functions

    //  L2 distance between the two cities caching of already computed values
    inline double RealDistance(int city1, int city2) {
        if (!edges[city1][city2]) {
            edges[city1][city2] = edges[city2][city1] = -sqrt(sqr(cities[city2].x - cities[city1].x) +
                sqr(cities[city2].y - cities[city1].y));
        }

        return fabs(edges[city1][city2]);
    }

    // distance between two points taking into account building penalty if the cities aren't connected
    inline double BuildingDistance(int u, int v) {
        return (RealDistance(u, v) * building_penalty);
    }

    void DisplayGraph(char * edgefile, char * cityfile);
    void DisplayBid(char * edgefile, Bid b);
    void Relax(int u, int v, double w, double d[]);
    int ExtractMin(int q[], int &qsize, double d[]);
    void ConstructPath(int s, int t);
    void InitGraph();
    bool Build(int city1, int city2);
    void BuildShortEdges();
    bool BuildViaPaths();
    void BuildSmallWorldEdges();
    int BuildEdgeGoodList();
    void BuildGraph(int num_goods);
    void BidOnPath(Path &p);
    void FindBids(Path &p, int dest);
    bool generateBids(Param &p, BidSet *bids);
    inline bool isConnected();
    int numComponents();
    char *outputSettings(bool to_file);
    static double param_mins[];
    static double param_maxs[];
};

#endif // _PATHS_H

