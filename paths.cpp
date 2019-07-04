#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <queue>

#include "Param.h"
#include "Distribution.h"
#include "bid.h"
#include "BidSet.h"
#include "paths.h"

const int Paths::Z_MAX_BID_SET_SIZE = 100;
const double Paths::INF = 500000.0;
double Paths::param_mins[] = {1, 1.5, 1.2, 2, 1};
double Paths::param_maxs[] = {5, 5, 2, 12, 2};

// constructor
Paths::Paths(Param &p, char *argv[], int argc)
{
    output_buffer = NULL;
    bids = NULL;

    // set initial values
    display_graph = false;
    initial_connections = 2;
    edge_density = 3;
    building_penalty = 1.7;
    shipping_cost_factor = 1.5;
    max_bid_set_size = 5;
    initialNumCities = -1;

    //  FOR NORMAL CITY PLACEMENT
    //  city_placement_type = UNIFORM;
    //  cluster_std_dev = 0.2;
    //  num_clusters = 2;
    //  if (city_placement_type == NORMAL) {
    //      norm = new Normal(p);
    //  }
    //  else
    //      norm = NULL;
    //  }

    int i;

    for (i = 1; i < argc; i++) {
        if (!strcmp("-cities", argv[i])) {
            p.argRecognized[i] = true;

            if (i + 1 >= argc) {
                usage();
            }

            initialNumCities = atoi(argv[i + 1]);
            i++;
            p.argRecognized[i] = true;
        }
        else if (!strcmp("-conn", argv[i])) {
            p.argRecognized[i] = true;

            if (i + 1 >= argc) {
                usage();
            }

            initial_connections = atoi(argv[i + 1]);
            i++;
            p.argRecognized[i] = true;
        }
        else if (!strcmp("-shcost", argv[i])) {
            p.argRecognized[i] = true;

            if (i + 1 >= argc) {
                usage();
            }

            shipping_cost_factor = atof(argv[i + 1]);
            i++;
            p.argRecognized[i] = true;
        }
        else if (!strcmp("-maxbid", argv[i])) {
            p.argRecognized[i] = true;

            if (i + 1 >= argc) {
                usage();
            }

            max_bid_set_size = atoi(argv[i + 1]);

            if (max_bid_set_size > Z_MAX_BID_SET_SIZE) {
                max_bid_set_size = Z_MAX_BID_SET_SIZE;
            }

            i++;
            p.argRecognized[i] = true;
        }

        else if (!strcmp("-bp", argv[i])) {
            p.argRecognized[i] = true;

            if (i + 1 >= argc) {
                usage();
            }

            building_penalty = atof(argv[i + 1]);
            i++;
            p.argRecognized[i] = true;
        }
        else if (!strcmp("-display", argv[i])) {
            p.argRecognized[i] = true;
            display_graph = true;
        }
        else if (!stricmp("-ed", argv[i])) {
            p.argRecognized[i] = true;

            if (i + 1 >= argc) {
                usage();
            }

            edge_density = atof(argv[i + 1]);
            i++;
            p.argRecognized[i] = true;
        }

        //  FOR NORMAL CITY PLACEMENT
        //  else if (!stricmp("-normal", argv[i])) {
        // 	    p.argRecognized[i] = true;
        // 	    city_placement_type = NORMAL;
        //  }
        //  else if (!stricmp("-clusters", argv[i])) {
        // 	    p.argRecognized[i] = true;
        // 	
        // 	    if (i+1 >= argc) {
        // 	        usage();
        // 	    }
        //
        // 	    num_clusters = atoi(argv[i+1]);
        // 	    i++;
        // 	    p.argRecognized[i] = true;
        //  }
        //  else if (!stricmp("-std_dev",argv[i])) {
        // 	    p.argRecognized[i] = true;
        //
        // 	    if (i+1 >= argc) {
        // 	        usage();
        //  	}
        //
        // 	    cluster_std_dev = atof(argv[i+1]);
        // 	    i++;
        // 	    p.argRecognized[i] = true;
        //  }
    }

    // allocate memory
    best_paths = new Path*[Z_MAX_BID_SET_SIZE];
    maxCities = p.num_goods * 5;
    cities = new Point[maxCities];
    goodedges1 = new int[maxCities * 2];
    goodedges2 = new int[maxCities * 2];
    edges = new double*[maxCities];
    edgesgood = new int*[maxCities];
    adjacentTo = new int*[maxCities];
    numAdjacentTo = new int[maxCities];

    for (i = 0; i < maxCities; i++) {
        edges[i] = new double[maxCities];
        edgesgood[i] = new int[maxCities];
        adjacentTo[i] = new int[maxCities];
    }

    q = new int[maxCities];
    pi = new int[maxCities];
    d = new double[maxCities];
}

// destructor
Paths::~Paths()
{
    int i;

    delete[] best_paths;
    delete[] cities;

    for (i = 0; i < maxCities; i++) {
        delete[] edges[i];
        delete[] edgesgood[i];
        delete[] adjacentTo[i];
    }

    delete[] edges;
    delete[] edgesgood;
    delete[] adjacentTo;
    delete[] numAdjacentTo;
    delete[] goodedges1;
    delete[] goodedges2;
    delete[] q;
    delete[] pi;
    delete[] d;

    if (output_buffer) {
        delete[] output_buffer;
    }

    //  if (norm) {
    //      delete norm;
    //  }
}

// output the graph 
void Paths::DisplayGraph(char * edgefile, char * cityfile)
{
    int i, j;

    FILE * edgep, *cityp;
    edgep = fopen(edgefile, "w");

    if (edgep == NULL) {
        printf("error opening display file for writing\n");
        exit(2);
    }

    cityp = fopen(cityfile, "w");

    if (cityp == NULL) {
        printf("error opening display file for writing\n");
        exit(2);
    }

    for (i = 0; i < numcities; i++) {
        for (j = 0; j < numAdjacentTo[i]; j++) {
            fprintf(edgep, "%d\t%d\n", i, adjacentTo[i][j]);
        }

        fprintf(cityp, "%f %f\n", cities[i].x, cities[i].y);
    }

    fclose(edgep);
    fclose(cityp);
}

void Paths::DisplayBid(char * edgefile, Bid b)
{
    int i;

    FILE * fp;
    fp = fopen(edgefile, "w");

    if (fp == NULL) {
        printf("error opening display file for writing\n");
        exit(2);
    }

    for (i = 0; i < (int) b.num_goods; i++) {
        fprintf(fp, "%f %f\n%f %f\n\n", cities[goodedges1[i]].x, cities[goodedges1[i]].y,
                cities[goodedges2[i]].x, cities[goodedges2[i]].y);
    }

    fclose(fp);
}

void Paths::Relax(int u, int v, double w, double d[])
{
    if (d[v] > d[u] + w) {
        d[v] = d[u] + w;
        pi[v] = u;
    }
}

// note: this could be made faster if q were represented with a heap.
int Paths::ExtractMin(int q[], int &qsize, double d[])
{
    int i, besti = 0, b;
    double bestval = d[q[0]];

    for (i = 1; i < qsize; i++) {
        if (d[q[i]] < bestval) {
            besti = i;
            bestval = d[q[i]];
        }
    }

    b = q[besti];
    q[besti] = q[--qsize];
    return b;
}

void Paths::ConstructPath(int s, int t)
{
    int u, i, j, qsize = numcities;

    for (i = 0; i < numcities; i++) {
        q[i] = i;
        d[i] = INF;
        pi[i] = -1;
    }

    d[s] = 0;

    while (qsize != 0) {
        u = ExtractMin(q, qsize, d);
        /* printf("Extracted %3d; target is %3d\n", u, t); */

        if (u == t) {
            return;
        }

        for (j = 0; j < qsize; j++) {
            if (q[j] != u) {
                Relax(u, q[j], BuildingDistance(u, q[j]), d);
            }
        }
    }
}

void Paths::InitGraph()
{
    int i, j;

    //  if (city_placement_type == UNIFORM) {
    for (i = 0; i < numcities; i++) {
        cities[i].x = Param::DRand();
        cities[i].y = Param::DRand();

        for (j = i; j < numcities; j++) {
            edges[i][j] = edges[j][i] = 0;
            edgesgood[i][j] = edgesgood[j][i] = -1;
        }

        numAdjacentTo[i] = 0;
    }

    //  FOR NORMAL CITY PLACEMENT
    //  } else {
    //      assert (city_placement_type == NORMAL);
    //      double *x_means = new double[num_clusters];
    //      double *y_means = new double[num_clusters];
    //
    //      for (i = 0; i<num_clusters; i++) {
    //          x_means[i] = Param::DRand();
    //          y_means[i] = Param::DRand();
    //      }
    //
    //      for (i = 0; i < numcities; i++) {      
    //          bool goodDraw = false;
    //
    //          while (!goodDraw) {
    // 	            int whichCluster = Param::LRand(0,num_clusters-1);
    // 	            cities[i].x = norm->draw(x_means[whichCluster], cluster_std_dev);
    // 	            cities[i].y = norm->draw(y_means[whichCluster], cluster_std_dev);
    //
    // 	            if (cities[i].x >= 0 && cities[i].x <= 1 &&
    // 	                cities[i].y >= 0 && cities[i].y <= 1) {
    // 	                goodDraw = true;
    // 	            }
    //          }
    //
    //          for (j = i+1; j < numcities; j++) {
    // 	            edges[i][j] = edges[j][i] = 0;
    // 	            edgesgood[i][j] = edgesgood[j][i] = -1;
    //          }
    //      }
    //
    //      delete[] x_means;
    //      delete[] y_means;
    //  }

    num_edges = 0;
}

bool Paths::Build(int city1, int city2)
{
    ConstructPath(city1, city2);
    bool madeRoad = false;

    while (city2 != city1) {
        if (edges[city2][pi[city2]] < 0) {
            madeRoad = true;

            edges[pi[city2]][city2] = edges[city2][pi[city2]] *= -1;
            adjacentTo[city2][numAdjacentTo[city2]++] = pi[city2];
            adjacentTo[pi[city2]][numAdjacentTo[pi[city2]]++] = city2;
            goodedges1[num_edges] = city2;
            goodedges2[num_edges] = pi[city2];
            edgesgood[pi[city2]][city2] = edgesgood[city2][pi[city2]] = num_edges++;
        }

        city2 = pi[city2];
    }

    return madeRoad;
}

void Paths::BuildShortEdges()
{
    DistanceTo *shortest = new DistanceTo[initial_connections];
    int i, j, k;

    for (i = 0; i < numcities; i++) {
        for (j = 0; j < initial_connections; j++) {
            shortest[j].distance = INF;
        }

        for (j = 0; j < numcities; j++) {
            if (i == j) {
                continue;
            }

            double dist = RealDistance(i, j);

            for (k = 0; k < initial_connections; k++) {
                if (dist < shortest[k].distance) {
                    break;
                }
            }

            if (k < initial_connections) {
                for (int l = initial_connections - 1; l > k; l--) {
                    shortest[l] = shortest[l - 1];
                }

                shortest[k].distance = dist;
                shortest[k].city = j;
            }
        }

        for (j = 0; j < initial_connections; j++) {
            Build(i, shortest[j].city);
        }
    }

    delete[] shortest;
}

bool Paths::BuildViaPaths()
{
    int city1, city2;
    int timeSinceLastBuild = 0;

    while (num_edges < numcities * edge_density) {
        city1 = Param::LRand() % numcities;

        do {
            city2 = Param::LRand() % numcities;
        }
        while (city1 == city2);

        if (Build(city1, city2)) {
            timeSinceLastBuild = 0;
        }
        else {
            timeSinceLastBuild++;
        }

        if (timeSinceLastBuild > numcities) {
            return false;
        }

    }

    return true;
}

void Paths::BidOnPath(Path &path)
{
    int j;

    if (num_best_paths < max_bid_set_size) {
        j = num_best_paths++;
    }
    else {
        delete best_paths[j = num_best_paths - 1];
    }

    while ((j > 0) && (best_paths[j - 1]->cost > path.cost)) {
        best_paths[j] = best_paths[j - 1];
        j--;
    }

    best_paths[j] = new Path(path, numcities);

    if (num_best_paths == max_bid_set_size) {
        maxcost = best_paths[max_bid_set_size - 1]->cost;
    }
}

void Paths::FindBids(Path &path, int dest)
{
    int lastInPath = path.pathArr[path.pathlen - 1];

    for (int i = 0; i < numAdjacentTo[lastInPath]; i++) {
        int nextCity = adjacentTo[lastInPath][i];
        bool looping = false;

        for (int j = path.pathlen - 1; j >= 0; j--) {
            if (path.pathArr[j] == nextCity) {
                looping = true;
                break;
            }
        }

        if (!looping &&
            (path.cost + edges[lastInPath][nextCity] +
            RealDistance(nextCity, dest) < maxcost)) {
            double origCost = path.cost;
            path.cost += edges[lastInPath][nextCity];
            path.pathArr[path.pathlen++] = nextCity;

            if (path.pathArr[path.pathlen - 1] == dest) {
                BidOnPath(path);
            }
            else {
                FindBids(path, dest);
            }

            path.cost = origCost;
            path.pathlen--;
        }
    }
}

bool Paths::isConnected()
{
    return (numComponents() == 1);
}

int Paths::numComponents()
{
    bool *citiesLocated = new bool[numcities];

    for (int i = 0; i < numcities; i++) {
        citiesLocated[i] = false;
    }

    int numCitiesLocated = 0;
    int compsFound = 0;
    int *tree = new int[numcities];
    int depth;

    while (numCitiesLocated < numcities) {
        // find first unfound city as beginning point of comp.
        tree[0] = -1;

        for (int i = 0; i < numcities; i++) {
            if (!citiesLocated[i]) {
                tree[0] = i;
                break;
            }
        }

        // if all cities are found, we shouldn't have started this loop
        assert(tree[0] >= 0);
        depth = 0;

        // mark this city
        citiesLocated[tree[0]] = true;
        numCitiesLocated++;
        //  printf("In component %d:\t", compsFound+1);

        while (depth >= 0) {
            // find adjacent unmarked city
            int next = -1;

            for (int i = 0; i < numAdjacentTo[tree[depth]]; i++) {
                if (!citiesLocated[adjacentTo[tree[depth]][i]]) {
                    next = adjacentTo[tree[depth]][i];
                    break;
                }
            }

            if (next == -1) {
                //  printf("%d\t", tree[depth]);
                depth--; // no adjacent unmarked cities.  backtrack
            }
            else {
                // move to next city and mark it
                tree[++depth] = next;
                citiesLocated[tree[depth]] = true;
                numCitiesLocated++;
            }
        }

        //  printf("\n");

        compsFound++;
    }

    delete[] citiesLocated;
    delete[] tree;
    return compsFound;
}

bool Paths::generateBids(Param &p, BidSet *bids)
{
    int city1, city2;
    Path path(numcities);
    double d, dist;
    int pass = 0;
    bool *edgeInBid = new bool[num_edges];

    for (int i = 0; i < num_edges; i++) {
        edgeInBid[i] = false;
    }

    int numEdgesInBids = 0;
    int numUnusedEdges = 0;
    int badBids = 0;
    bool done = false, prunedBadEdges = false, tooManyBids = false;

    // if max_bid_set_size > 1, there cannot be dominated bids, so why check?

    bool removeDominatedBids = (p.remove_dominated_bids && max_bid_set_size == 1);

    while (!done) {
        //  printf("edges: %d, edgesInBids: %d, unusedEdges: %d, badBids: %d, totalBids: %d\n", 
        //  num_edges, numEdgesInBids, numUnusedEdges, badBids, bids->numBids());

        city1 = Param::LRand() % numcities;

        do {
            city2 = Param::LRand() % numcities;
        }
        while (city2 == city1);

        path.pathArr[0] = city1;
        d = (Param::DRand() * (shipping_cost_factor - 1)) + 1;
        dist = RealDistance(city1, city2);
        maxcost = d * dist;
        path.pathlen = 1;
        path.cost = 0;
        num_best_paths = 0;

        // this sorts the cities adj lists to explore
        // the ones closer to the dest first

        for (int h = 0; h < numcities; h++) {
            for (int i = numAdjacentTo[h] - 1; i > 0; i--) {
                for (int j = 0; j < i; j++) {
                    if (RealDistance(city2, adjacentTo[h][j]) >
                        RealDistance(city2, adjacentTo[h][j + 1])) {
                        int temp = adjacentTo[h][j + 1];
                        adjacentTo[h][j + 1] = adjacentTo[h][j];
                        adjacentTo[h][j] = temp;
                    }
                }
            }

            //  if (pass == 0 && h == 0) {
            // 	    for (int i=0; i<numAdjacentTo[h]; i++) {
            // 	        printf("%f\n", RealDistance(city2, adjacentTo[h][i]));
            // 	    }
            //  }
        }

        FindBids(path, city2);

        // construct all the substitutes for the bidder
        bool onlyOneGenerated = false;

        for (int i = 0; i < num_best_paths; i++) {
            Bid *bid = new Bid(maxcost - path.cost);

            for (int j = 1; j < best_paths[i]->pathlen; j++) {
                int good = edgesgood[best_paths[i]->pathArr[j - 1]][best_paths[i]->pathArr[j]];
                bid->addGood(good);

                if (!edgeInBid[good]) {
                    // once bad edges have been pruned they shouldn't be bid on anymore
                    assert(!prunedBadEdges);
                    edgeInBid[good] = true;
                    numEdgesInBids++;
                }
            }

            // is it time to stop?
            if (num_best_paths == 1 ||
                (!prunedBadEdges && numEdgesInBids >= p.num_goods)) {
                if (i == 0) {
                    // only one was generated -- don't use Xor
                    bids->add(bid, removeDominatedBids);
                    onlyOneGenerated = true;
                }
                else {
                    bids->addXor(bid);
                }

                break;
            }

            bids->addXor(bid);
        }

        if (!onlyOneGenerated) {
            bids->doneAddingXor(removeDominatedBids);
        }

        for (int i = 0; i < num_best_paths; i++) {
            delete best_paths[i];
        }

        //  int nc = numcities;

        if (!prunedBadEdges && numEdgesInBids >= p.num_goods) {
            if (!bids->isConnected()) {
                //  printf("Removing bad bids.  Starting with %d edges in bids.\n", numEdgesInBids);
                assert(bids->disconnectedBids && bids->numDisconnectedBids > 0);
                badBids += bids->numDisconnectedBids;

                for (int i = 0; i < bids->numDisconnectedBids; i++) {
                    int bidNum = bids->disconnectedBids[i];

                    for (int j = 0; j < bids->getBid(bidNum)->num_goods; j++) {
                        int good = bids->getBid(bidNum)->getGood(j);

                        // don't count dummy goods or recount goods
                        if (good < num_edges && edgeInBid[good]) {
                            numEdgesInBids--;
                            edgeInBid[good] = false;
                        }
                    }

                    bids->removeBid(bidNum);

                    for (int j = i + 1; j < bids->numDisconnectedBids; j++) {
                        if (bids->disconnectedBids[j] == bids->numBids()) {
                            bids->disconnectedBids[j] = bidNum;
                        }
                    }
                }

                //	printf("Removed %d bad bids.  %d edges left\n", bids->numDisconnectedBids, 
                //	numEdgesInBids);

                assert(bids->isConnected());
            }

            // if there are still enough goods, remove unused edges
            if (numEdgesInBids >= p.num_goods) {
                int origNumEdges = num_edges;

                for (int i = 0; i < origNumEdges; i++) {
                    if (!edgeInBid[i]) {
                        edges[goodedges1[i]][goodedges2[i]] = edges[goodedges2[i]][goodedges1[i]] *= -1;

                        for (int j = 0; j < numAdjacentTo[goodedges1[i]]; j++) {
                            if (adjacentTo[goodedges1[i]][j] == goodedges2[i]) {
                                adjacentTo[goodedges1[i]][j] = adjacentTo[goodedges1[i]][--numAdjacentTo[goodedges1[i]]];
                                break;
                            }

                            assert(j != numAdjacentTo[goodedges1[i]] - 1);
                        }

                        for (int j = 0; j < numAdjacentTo[goodedges2[i]]; j++) {
                            if (adjacentTo[goodedges2[i]][j] == goodedges1[i]) {
                                adjacentTo[goodedges2[i]][j] = adjacentTo[goodedges2[i]][--numAdjacentTo[goodedges2[i]]];
                                break;
                            }

                            assert(j != numAdjacentTo[goodedges2[i]] - 1);
                        }

                        num_edges--;
                        numUnusedEdges++;
                    }
                }

                prunedBadEdges = true;
                //	printf("Pruned %d unused edges.\n", numUnusedEdges);
            }
        }

        if (bids->numBids() >= p.num_bids) {
            if (!prunedBadEdges) {
                tooManyBids = true;
            }

            done = true;
        }

        // output progress
        if (p.output_parameter_settings &&
            (++pass % p.output_frequency == 0 || done)) {
            int count = printf("%d %sbids using %d goods created out of %d.",
                               bids->numBids(), p.remove_dominated_bids ? "non-dominated " : "", 
                               numEdgesInBids, p.num_bids);

            if (!done) {
                for (int counter = 0; counter < count; counter++) {
                    printf("\b");
                }
            }
            else {
                printf("\n");
            }
        }
    }

    //  printf("%d unused edges pruned including %d in bad bids.  %d remain.  %d in bids.\n", 
    //  numUnusedEdges, badBids, num_edges, numEdgesInBids);

    delete[] edgeInBid;
    return !tooManyBids;
}

// this is the function that generates the bids
BidSet *Paths::generate(Param &p)
{
    if (initialNumCities == -1) {
        initialNumCities = (int) (1.2 * (double) p.num_goods / edge_density);
    }

    numcities = initialNumCities;
    bool done = false;

    while (!done) {
        if (p.output_parameter_settings) {
            printf("Building city graph with %d cities and %d edges.\n", numcities, 
                    (int) (numcities * edge_density));
        }

        if (numcities > maxCities) {
            printf("Failed to produce valid graph with a reasonable number of cities.\n");
            printf("Try lowering building penalty or raising shipping cost.\n");
            exit(1);
        }

        InitGraph();
        BuildShortEdges();

        if (!BuildViaPaths()) {
            numcities++;
            continue;
        }

        if (!isConnected()) {
            continue;
        }

        // allow for all edges, even though some will probably not end up in a bid
        bids = new BidSet(num_edges, p.num_bids + max_bid_set_size); // maximum number of bids

        if (generateBids(p, bids)) {
            done = true;
        }
        else {
            numcities++;
            delete bids;
        }
    }

    //   if (p.output_parameter_settings) {
    //     printf("Built graph using %d cities.\n", numcities);
    //   }

    if (display_graph) {
        DisplayGraph("edges", "cities");
    }

    bids->compact();
    return bids;
}

void Paths::randomizeParams(Param &p)
{
    if (p.param_dists[(int) Param::PATHS] != NULL) {
        double params[5];
        randomizeParamsWithDistribution(p.param_dists[(int) Param::PATHS], 5, 
                param_mins, param_maxs, params, p);

        initial_connections = (int) params[0];
        edge_density = params[1];
        shipping_cost_factor = params[2];
        max_bid_set_size = (int) params[3];
        building_penalty = params[4];

        if (initial_connections >= initialNumCities) {
            initial_connections = initialNumCities - 1;
        }
    }
    else {
        initial_connections = Param::LRand((int) param_mins[0], (int) param_maxs[0]);
        edge_density = Param::DRand(param_mins[1], param_maxs[1]);
        shipping_cost_factor = Param::DRand(param_mins[2], param_maxs[2]);
        max_bid_set_size = Param::LRand((int) param_mins[3], (int) param_maxs[3]);
        building_penalty = Param::DRand(param_mins[4], param_maxs[4]);
    }
}

void Paths::usage()
{
    int count = printf("Usage for Paths Distribution:\n");
    for (int t = 0; t < count; t++) {
        printf("=");
    }

    printf(" -cities    : sets the number of initial nodes in the graph.\n");
    printf(" -conn      : sets the number of initial connections.\n");
    printf(" -ed        : sets the edge density (number of edges per city).\n");
    printf(" -shcost    : sets the shipping cost factor.\n");
    printf(" -maxbid    : sets the maximum number of bids each bidder can make.\n");
    printf(" -bp        : sets the building penalty.\n");
    printf(" -bpalph    : sets alpha for building penalty.\n");
    //  FOR NORMAL CITY PLACEMENT
    //  printf(" -normal    : uses mixture of normal distributions for city placement (instead of uniform).\n");
    //  printf(" -clusters  : sets number of clusters for normal city placement type.\n");
    //  printf(" -std_dev   : sets standard deviation for normal distributions of city placement.\n");
    printf(" -display   : output the connectivity graph\n");
}

// output distribution parameters, either to the screen or to a file
char *Paths::outputSettings(bool tofile)
{
    // allocate variables
    char comment[3];

    if (tofile) {
        comment[0] = '%';
        comment[1] = ' ';
        comment[2] = 0;
    }

    else {
        comment[0] = 0;
    }

    char buffer1[80], buffer2[80], buffer3[80], buffer4[80], buffer5[80];
    char buffer6[80];

    // generate output
    sprintf(buffer1, "%sPaths Distribution Parameters:\n", comment);
    sprintf(buffer2, "%sInitial connections:       = %d\n", comment, initial_connections);
    sprintf(buffer3, "%sEdge density:              = %f\n", comment, edge_density);
    sprintf(buffer4, "%sShipping cost factor:      = %f\n", comment, shipping_cost_factor);
    sprintf(buffer5, "%sMax bid set size:          = %d\n", comment, max_bid_set_size);
    sprintf(buffer6, "%sBuilding penalty:          = %f\n", comment, building_penalty);

    // prepare the string and return it
    int length = strlen(buffer1) + strlen(buffer2) + strlen(buffer3) + strlen(buffer4) +
            strlen(buffer5) + strlen(buffer6) + 20;

    if (output_buffer) {
        delete[] output_buffer;
    }

    output_buffer = new char[length];
    sprintf(output_buffer, "%s%s%s%s%s%s\n", buffer1, buffer2, buffer3, buffer4, buffer5, buffer6);

    return output_buffer;
}

