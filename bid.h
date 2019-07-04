///////////////////////////////////////////////
/* Stores a full or partial (bid) allocation */
///////////////////////////////////////////////

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#ifndef BID_H
#define BID_H

#include <set>

using namespace std;

class Bid {
public:

    Bid(double a, int i = 0);
    ~Bid();
    double amount;
    int bid_num;
    int indexInBidSet;
    int num_goods;
    bool visited;

    //set<Bid*> conflicts; // -- Set of conflicting bids

    int conflicts; //Number of conflicts - o.w. runs out of memory
    int firstDummyGood(unsigned total_goods);
    void addGood(unsigned g);

    inline void removeGood(unsigned g) {
        if (num_goods > 0) {
            array[indexOf(g)] = array[--num_goods]; // it's OK if these are the same
            sorted = false;
        }
    }

    inline unsigned getGood(int index) {
        return array[index];
    }

    inline bool contains(unsigned g) {
        return (indexOf(g) >= 0);
    }

    bool subsetEqualOf(Bid *other);
    bool conflictsWith(Bid *other);
    int largestGood();
    void renumber(int *newNums);

protected:
    unsigned *array;
    int array_size;
    int indexOf(unsigned g);
    void sort();
    bool sorted;
    void sortFrom(int begin, int end, unsigned *tempArray);
    void merge(int leftBegin, int rightBegin, int rightEnd, unsigned *tempArray);
};

#endif

