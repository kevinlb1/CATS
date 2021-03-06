/* INPUT FILE FORMAT:

   ==================

% comments follow percentage symbols
% blank lines are ignored
% files are not case-sensitive
% [CR]'s are optional, all whitespace is ignored

goods 15
dummy 10
bids 800

% bidnum, price, goodnum, goodnum, ..., goodnum, #
0	20.75	2	6	8	#
1	30		1	2	5	7	13	#
2	5.20	4	#
  .
  .
  .

% implicitly, all goods numbered 15 or above (in this example) are dummies
% note that all counting starts from zero
% bidnum's do not need to be sequential between 0..799.  They just need to be unique.
*/

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#ifndef BidSet_H
#define BidSet_H

#include <stdio.h>
#include <string.h>
#include <set>
#include <map>

using namespace std;

#include "bid.h"
#include "Param.h"
#include "Distribution.h"

class BidSet {
public:
    void refine();
    void compact();
    bool parse(char *filename);
    void parse();

    // functions
    void print(int maxrec = -1);
    void writeCPLEXFile(Param &p, int index = 0);
    void writeOutput(Param &p, Distribution *d, int index = 0);

    bool isConnected();
    int numComponents();
    int *disconnectedBids;
    int numDisconnectedBids;
    void makeBidGraph(void);

    // constructors/destructor
    BidSet(const int num_goods = 0, int init_array_size = 100, bool deleteRemovedBids = true);
    BidSet(const BidSet& orig);
    BidSet(const char *filename, bool deleteRemovedBids = true);
    ~BidSet();

    // VCC apparently doesn't use copy constructor for assignment
    BidSet& operator=(const BidSet& orig);

    // bid access functions
    void removeBid(int index);
    void add(Bid *added_bid, bool doRefine);
    void addXor(Bid *added_bid);
    void doneAddingXor(bool run_refine);

    inline int numBids() {
        return num_bids;
    }

    inline int numGoods() {
        return num_goods;
    }

    inline int numDummyGoods() {
        return num_dummy_items;
    }

    inline Bid *getBid(int index) {
        return bid[index];
    }

    map<int, set<int> > *getGoodMap(void) {
        return &goodMap;
    }

    int **bidGraph; // matrix representation of conflicts

private:
    // file functions
    void skipComments();
    void fileRead(int &int_input);
    void fileRead(double &double_input);
    void fileRead(const char *label, int &int_input);
    signed int fileReadGood(int bidnum);

    bool isAllWhitespace(char *input, int len);
    void initialize();
    bool removeSingletonDummyGoods();

    // processing
    int deleteDominatedBids(int t, bool singleton_cleared);

    // variables
    int num_goods, num_dummy_items, num_bids;

    Bid **bid;
    FILE *fp;
    int num_killed;
    char buffer[1000];
    int first_unchecked_bid, bid_array_size;
    bool delRemovedBids;
    map<int, set<int> > goodMap; // -- map goods to bids
};

#endif

