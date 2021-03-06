#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#ifdef _WIN32
#pragma warning (disable : 4786)
#endif

#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <fstream>
#include <string.h>
#include <strings.h>
#include <time.h>

#include "BidSet.h"
#include "Param.h"

#include <map>
#include <set>
#include <algorithm>

using namespace std;

/* constructor */
BidSet::BidSet(const int n_g, int init_array_size, bool deleteRemovedBids)
{
    delRemovedBids = false; //deleteRemovedBids;
    num_goods = n_g;
    num_dummy_items = 0;
    first_unchecked_bid = 1;
    bid_array_size = init_array_size;

    bid = new Bid*[bid_array_size];
    num_bids = 0;
    num_killed = 0;

    goodMap.clear();
    bidGraph = NULL;

    disconnectedBids = NULL;
    numDisconnectedBids = 0;
}

BidSet::BidSet(const char *filename, bool deleteRemovedBids)
{
    delRemovedBids = deleteRemovedBids;
    num_dummy_items = 0;
    first_unchecked_bid = 1;
    num_killed = 0;

    goodMap.clear();
    bidGraph = NULL;

    disconnectedBids = NULL;
    numDisconnectedBids = 0;

    fp = fopen(filename, "rt");
    if (fp == NULL) {
        printf("Error opening file \"%s\"\n", filename);
        exit(1);
    }

    int temp, t;
    signed int foundGood;

    // read number of goods
    fileRead("goods", num_goods);

    // read number of bids
    fileRead("bids", num_bids);
    bid_array_size = num_bids;
    bid = new Bid*[bid_array_size];

    // read number of dummy goods
    fileRead("dummy", num_dummy_items);

    // make sure it's not a multi-unit file
    long pos = ftell(fp);

    skipComments();
    fscanf(fp, "%s", buffer);

    if (!strcasecmp(buffer, "maximums")) {
        printf("Error reading input file.\nThis is a multi-unit file, and CATS currently only supports single-unit mode.\n");
        exit(1);
    }

    fseek(fp, pos, SEEK_SET);

    // read bids
    for (t = 0; t < num_bids; t++) {
        // create a new bid
        bid[t] = new Bid(0, t);

        // read bid number
        fileRead(temp);
        bid[t]->bid_num = temp;

        // read bid price
        fileRead(bid[t]->amount);

        // read goods until end-of-record encountered
        foundGood = fileReadGood(bid[t]->bid_num);

        while (foundGood != -1) {
#ifdef MULTI_UNIT
#ifdef ASSUME_SINGLE_UNITS
            temp = 1;
#else
            fileRead(&temp);
#endif
            bid[t]->vec->set(foundGood, temp);
#endif

            bid[t]->addGood(foundGood);
            foundGood = fileReadGood(bid[t]->bid_num);

        }

        // update goodMap
        //for (int i=0; i<bid[t]->num_goods; i++)
        //	goodMap[bid[t]->getGood(i)].insert(bid[t]);
        // no longer incremental

    }

    // zero terminate the array
    if (num_bids < bid_array_size) {
        bid[num_bids] = NULL;
    }

    fclose(fp);
}

// VC++ won't create operator= from copy constructor, so we have to do it the ugly way
BidSet::BidSet(const BidSet& orig)
{
    *this = orig;
}

BidSet& BidSet::operator=(const BidSet& orig)
{
    delRemovedBids = orig.delRemovedBids;
    num_goods = orig.num_goods;
    num_dummy_items = orig.num_dummy_items;
    first_unchecked_bid = orig.first_unchecked_bid;
    bid_array_size = orig.bid_array_size;
    num_bids = orig.num_bids;
    bid_array_size = orig.bid_array_size;

    bid = new Bid*[bid_array_size];
    assert(num_bids <= bid_array_size);

    int i;
    for (i = 0; i < num_bids; i++) {
        bid[i] = orig.bid[i];
    }

    if (num_bids < bid_array_size) {
        bid[i] = NULL;
    }

    num_killed = orig.num_killed;
    goodMap = orig.goodMap;
    bidGraph = orig.bidGraph;

    return *this;
}

/* destructor */
BidSet::~BidSet()
{
    if (delRemovedBids) {
        for (int t = 0; t < bid_array_size && bid[t]; t++) {
            delete bid[t];
        }
    }

    delete[] bid;

    if (bidGraph) {
        for (int i = 0; i < num_bids; i++) {
            delete[] bidGraph[i];
        }

        delete[] bidGraph;
    }

    if (disconnectedBids) {
        delete[] disconnectedBids;
    }

}

/* skip over comments */
void BidSet::skipComments() 
{
    int pos;
    int t;
    char c;

    do {
        pos = ftell(fp);
        fgets(buffer, 1000, fp);

        for (t = 0, c = ' '; t < (int) strlen(buffer) && isspace(c); t++) {
            c = buffer[t]; // allow comments to follow whitespace
        }
    }
    while (!feof(fp) && (c == '%' || isAllWhitespace(buffer, strlen(buffer))));

    fseek(fp, pos, SEEK_SET);

}

bool BidSet::isAllWhitespace(char *input, int len)
{
    for (int t = 0; t < len; t++) {
        if (!isspace(input[t])) {
            return false;
        }
    }

    return true;
}

/* read an integer from the file */

void BidSet::fileRead(int &int_input)
{
    skipComments();
    fscanf(fp, "%d", &int_input);
}

/* read a double from the file */
void BidSet::fileRead(double &double_input)
{
    skipComments();
    fscanf(fp, "%lf", &double_input);
}

/* read a label followed by an integer from the file */
void BidSet::fileRead(const char *label, int &int_input)
{
    skipComments();
    fscanf(fp, "%s", buffer);

    if (strcasecmp(buffer, label)) {
        printf("Error reading input file.\nExpected label %s, found %s.\n", label, buffer);
        exit(1);
    }

    fscanf(fp, "%d", &int_input);

}

/* read a good, unless '#' is the first thing found */
signed int BidSet::fileReadGood(int bidnum)
{
    int goodnum;
    skipComments();

    // read input
    fscanf(fp, "%s", buffer);

    // we found a pound sign -- end of record
    if (buffer[0] == '#') {
        return -1;
    }

    // we found a good number
    sscanf(buffer, "%d", &goodnum);

    if (goodnum >= num_goods + num_dummy_items) {
        printf("Good #%d was named in bid %d, but %d is the largest possible good #.\n",
               goodnum, bidnum, num_goods);

        exit(1);
    }

    return goodnum;
}

/* output what has been read, after the read is complete */
void BidSet::print(int maxrec)
{
    printf("Goods: %d, DGoods: %d, Bids: %d\n", num_goods, num_dummy_items, num_bids);

    if (maxrec == -1) {
        maxrec = num_bids;
    }

    for (int t = 0; t < maxrec; t++) {
        printf("Bid#: %d, Amount: %f, Goods: ", bid[t]->bid_num, bid[t]->amount);

        for (int tt = 0; tt < bid[t]->num_goods; tt++) {
            printf("%d, ", bid[t]->getGood(tt));
        }

        printf("\b\b  \n");
    }
}

int BidSet::deleteDominatedBids(int t, bool singleton_cleared)
{
    // check for dominated bids
    for (int tt = (singleton_cleared ? 0 : t + 1); tt < num_bids; tt++) {
        // skip deleted bids
#ifdef MULTI_UNIT
        if (bid[tt]->numGoods == 1) {
            continue;
        }
#endif

        // don't check a bid against itself
        if (tt == t) {
            continue;
        }

        // amount(tt) is more, tt subset_of t
        if (bid[tt]->amount >= bid[t]->amount && bid[tt]->subsetEqualOf(bid[t])) {
#ifdef MULTI_UNIT
            appendDominatedBid(tt, t);
            num_killed++;
#else
            removeBid(t);
            //t--;
            num_killed++;
            //goto next_bid;  // because we just deleted bid t--now we have to compare it to everything.
            return -1; // instead of the commented goto
#endif
        }

        // amount(t) is more, t subset_of tt
        if (bid[t]->amount >= bid[tt]->amount && bid[t]->subsetEqualOf(bid[tt])) {
#ifdef MULTI_UNIT
            appendDominatedBid(t, tt);
#else
            removeBid(tt);
            if (t == num_bids) {
                return deleteDominatedBids(tt, true); 
                // the bid we just deleted was the one we were comparing against
            }

            else tt--; // check the new bid that has been substituted in for tt
#endif
            num_killed++;
        }
    }

    return 0;

}

// go through the list of bids and clear dummy goods that are only referenced in a single bid.
// if such bids are found, delete any bids that this bid now dominates, 
// or delete this bid if it is dominated.

bool BidSet::removeSingletonDummyGoods()
{
    int dummy;
    bool pair, got_one = false;

    // don't waste time if there aren't any dummies anyway
    if (num_dummy_items == 0) {
        return false;
    }

    // look at each bid
    for (int t = 0; t < num_bids - 1; t++) {
        dummy = bid[t]->firstDummyGood(this->num_goods); // assumes that bids only contain one dummy good ?
        pair = false;

        if (dummy != -1) {
            //	for (int tt = t+1;!pair && tt<num_bids;tt++)
            for (int tt = 0; !pair && tt < num_bids; tt++) {
                if (tt == t) {
                    continue; // don't compare a bid to itself
                }

                int dummy2 = bid[tt]->firstDummyGood(this->num_goods);
                if (dummy2 == dummy) {
                    pair = true;
                }
            }

            if (!pair) {
                bid[t]->removeGood(dummy);
                int old_num_killed = num_killed;
                deleteDominatedBids(t, true); // do we have to decrement t here?

                if (num_killed > old_num_killed) {
                    got_one = true;
                }
            }

            if (!pair) {
                t--;
            }
        }
    }

    return got_one;
}

// make good numbers consecutive
void BidSet::compact()
{
    bool *goodExists = new bool[num_goods + num_dummy_items];

    for (int i = 0; i < num_goods + num_dummy_items; i++) {
        goodExists[i] = false;
    }

    for (int i = 0; i < num_bids; i++) {
        for (int j = 0; j < bid[i]->num_goods; j++) {   
            goodExists[bid[i]->getGood(j)] = true;
        }
    }

    int numSkipped = 0;
    int *newNums = new int[num_goods + num_dummy_items];

    for (int i = 0; i < num_goods + num_dummy_items; i++) {
        if (goodExists[i]) {
            newNums[i] = i - numSkipped;
        } 
        else {
            numSkipped++;
        }
    }

    for (int i = 0; i < num_bids; i++) {
        bid[i]->renumber(newNums);
    }

    num_goods -= numSkipped;
    delete[] goodExists;
    delete[] newNums;
}

// remove dominated bids
void BidSet::refine()
{
    int t;

    // first, delete dominated bids once
    // consider each bid
    for (t = first_unchecked_bid; t < num_bids; t++) {
        t += deleteDominatedBids(first_unchecked_bid, true);
    }

    // delete bids over and over as long as dummies are removed
    while (removeSingletonDummyGoods());

    // update the bids that have now been checked
    first_unchecked_bid = num_bids;
}

// remove the bid 'index'
void BidSet::removeBid(int index)
{
    assert(num_bids > 0);

    // update goodMap
    //	for (int i=0; i<bid[index]->num_goods; i++) {
    //		goodMap[bid[index]->getGood(i)].erase(bid[index]);
    //	}

    // no longer incremental
    // remove the bid

    num_bids--;

    if (delRemovedBids) {
        delete bid[index];
    }

    if (num_bids > 0) {
        bid[index] = bid[num_bids]; // if these are equal, so what?
        bid[index]->indexInBidSet = index;
    }

    bid[num_bids] = NULL;
}

// add a bid
void BidSet::add(Bid *added_bid, bool doRefine)
{
    // make sure there's room in the new array
    if (num_bids == bid_array_size) {
        bid_array_size *= 2;
        Bid** newbid = new Bid*[bid_array_size];

        // copy the old bids into the new array
        int t;
        for (t = 0; t < num_bids; t++) {
            newbid[t] = bid[t];
        }

        if (bid) {
            delete[] bid;
        }

        bid = newbid;
    }

    bid[num_bids++] = added_bid;

    // zero terminate the array
    if (num_bids < bid_array_size) {
        bid[num_bids] = NULL;
    }

    //  update goodMap
    //  for (int i=0; i<added_bid->num_goods; i++) {
    //      goodMap[added_bid->getGood(i)].insert(added_bid);
    //  }

    //  no longer incremental
    if (doRefine) {
        refine();
    }
}

// add a set of mutex bids
void BidSet::addXor(Bid *added_bid)
{
    added_bid->addGood(num_goods + num_dummy_items);
    add(added_bid, false);
}

// done adding a set of mutex bids: now refine can be done.
void BidSet::doneAddingXor(bool run_refine)
{
    num_dummy_items++;
    if (run_refine) refine();
}

// write out the output file
void BidSet::writeOutput(Param &p, Distribution *d, int index)
{
    char filename[strlen(p.filename) + 10];
    const char* file_ext = strstr(p.filename, ".txt") ? "" : ".txt";

    if (p.verbatim) {
        sprintf(filename, "%s%s", p.filename, file_ext);
    }
    else {
        sprintf(filename, "%s%04d%s", p.filename, index, file_ext);
    }

    ofstream outfile(filename);

    // comments
    long time_buf = time(NULL);
    char* time_str = ctime(&time_buf);
    outfile << "%% File generated by CATS v." << Param::CATS_VERSION_STRING << time_str << endl;
    outfile << "%% The CATS webpage is http://robotics.stanford.edu/CATS\n\n";
    outfile << "%% PARAMETER SETTINGS:\n";
    outfile << p.outputSettings(true) << endl;
    outfile << d->outputSettings(true) << endl << endl;

    // header
    outfile << "goods " << p.num_goods << "\nbids " << num_bids << "\ndummy " << num_dummy_items << "\n\n";

    // bids
    for (int t = 0; t < numBids(); t++) {
        // bid number
        outfile << t << "\t";

        // prices
        if (p.integer_prices) {
            outfile << (int) getBid(t)->amount << "\t";
        }
        else {
            outfile << getBid(t)->amount << "\t";
        }

        // all the goods
        for (int tt = 0; tt < getBid(t)->num_goods; tt++) {
            outfile << getBid(t)->getGood(tt) << "\t";
        }

        outfile << "#\n";
    }

    outfile.close();
    if (p.output_parameter_settings) {
        printf("\nWrote CATS output file \"%s\".\n", filename);
    }
}

/* write a CPLEX input file */
// assumes that bid numbers are unique
void BidSet::writeCPLEXFile(Param &p, int index)
{
    FILE *fp2;
    int i;
    char filename[strlen(p.cplex_filename) + 10];
    const char* file_ext = strstr(p.cplex_filename, ".lp") ? "" : ".lp";

    if (p.verbatim) {
        sprintf(filename, "%s%s", p.cplex_filename, file_ext);
    }
    else {
        sprintf(filename, "%s%04d%s", p.cplex_filename, index, file_ext);
    }

    fp2 = fopen(filename, "wt");

    if (fp2 == NULL) {
        printf("Error writing CPLEX file \"%s\"\n", filename);
    }

    // first section: maximize revenue
    fprintf(fp2, "max\n\n");

    for (i = 0; i < num_bids; i++) {
        fprintf(fp2, "%lgx%d ", bid[i]->amount, i);
        //  fprintf (fp2, "%g%d ", bid[i]->amount, i);
        if (i < num_bids - 1) {
            fprintf(fp2, "+ ");
        }
        if (i % 10 == 9) {
            fprintf(fp2, "\n");
        }
    }

    // second section: constraints
    fprintf(fp2, "\n\nst\n\n");

    for (i = 0; i < num_goods + num_dummy_items; i++) {
        bool make_a_plus = false;
        int printed = 0;

        for (int j = 0; j < num_bids; j++) {
            // list each bid that names the current good
            if (bid[j]->contains(i)) {
                if (make_a_plus) {
                    fprintf(fp2, "+ ");
                }

                //	#ifdef MULTI_UNIT
                //	    fprintf(fp2,"%dx%d ",bid[j]->vec->readArray(i),bid[j]->bid_num);
                //	#else
                fprintf(fp2, "x%d ", j);
                //	#endif

                make_a_plus = true;
                if (printed++ % 10 == 9) {
                    fprintf(fp2, "\n");
                }
            }
        }

        // finish off the constraint
        //	#ifdef MULTI_UNIT
        //  	fprintf (fp2, " <= %d\n\n",max_array[i]);
        //	#else
        fprintf(fp2, " <= 1\n\n");
        //	#endif
    }

    // third section: integrality
    fprintf(fp2, "integer\n\n");

    for (i = 0; i < num_bids; i++) {
        fprintf(fp2, "x%d ", i);
        if (i % 10 == 9) {
            fprintf(fp2, "\n");
        }
    }

    fprintf(fp2, "\n\nend\n");
    fclose(fp2);

    if (p.output_parameter_settings) {
        printf("Wrote CPLEX output file \"%s\".\n", filename);
    }
    //	delete[] filename;
}

bool BidSet::isConnected()
{
    return (numComponents() == 1);
}

// counts components w/o conflict sets
int BidSet::numComponents()
{
    //  printf("Counting components...\n");
    bool *bidsLocated = new bool[num_bids];

    for (int i = 0; i < num_bids; i++) {
        bidsLocated[i] = false;
    }

    int numBidsLocated = 0;
    int compsFound = 0;
    int *tree = new int[num_bids];
    int depth;
    int *largestComp = new int[num_bids];
    int sizeOfLargestComp = 0;
    int *thisComp = new int[num_bids];

    while (numBidsLocated < num_bids) {
        // find first unfound bid as beginning point of component
        tree[0] = -1;

        for (int i = 0; i < num_bids; i++) {
            if (!bidsLocated[i]) {
                tree[0] = i;
                break;
            }
        }

        // if all bids are found, we shouldn't have started this loop
        assert(tree[0] >= 0);
        depth = 0;

        // mark this bid
        bidsLocated[tree[0]] = true;
        numBidsLocated++;
        thisComp[0] = tree[0];
        int sizeOfThisComp = 1;

        while (depth >= 0) {
            // find adjacent unmarked bid
            int next = -1;

            for (int i = 0; i < num_bids; i++) {
                if (!bidsLocated[i] && bid[tree[depth]]->conflictsWith(bid[i])) {
                    next = i;
                    break;
                }
            }

            if (next == -1) {
                depth--; // no adjacent unmarked bids.  backtrack
            }
            else {
                // move to next bid and mark it
                tree[++depth] = next;
                bidsLocated[next] = true;
                numBidsLocated++;
                thisComp[sizeOfThisComp++] = next;
            }
        }

        compsFound++;
        //  printf("In component %d, %d bids:\t", compsFound, sizeOfThisComp);
        //  for (int i=0; i<sizeOfThisComp; i++) {
        //      printf("%d\t", thisComp[i]);      
        //  }
        //  printf("\n");

        if (sizeOfThisComp > sizeOfLargestComp) {
            int *temp = largestComp;
            largestComp = thisComp;
            thisComp = temp;
            sizeOfLargestComp = sizeOfThisComp;
        }
    }

    if (compsFound > 1) {
        assert(sizeOfLargestComp < num_bids);

        if (disconnectedBids) {
            delete[] disconnectedBids;
        }

        disconnectedBids = new int[num_bids - sizeOfLargestComp];
        numDisconnectedBids = 0;
        bool *inLargestComp = new bool[num_bids];

        for (int i = 0; i < num_bids; i++) {
            inLargestComp[i] = false;
        }

        for (int i = 0; i < sizeOfLargestComp; i++) {
            inLargestComp[largestComp[i]] = true;
        }

        for (int i = 0; i < num_bids; i++) {
            if (!inLargestComp[i]) {
                disconnectedBids[numDisconnectedBids++] = i;
            }
        }

        delete[] inLargestComp;
        assert(numDisconnectedBids == num_bids - sizeOfLargestComp);
    }

    delete[] largestComp;
    delete[] thisComp;
    delete[] bidsLocated;
    delete[] tree;
    return compsFound;
}

// -- Create a bid graph
void BidSet::makeBidGraph(void)
{
    bidGraph = new int *[num_bids];
    goodMap.clear();
    int i, j; // - Stupid MSVC doesn't support ansi

    for (i = 0; i < num_bids; i++) {
        Bid* b = bid[i];
        bidGraph[i] = new int[num_bids];

        for (j = 0; j < num_bids; j++) {
            bidGraph[i][j] = 0;
        }

        for (j = 0; j < b->num_goods; j++) {
            goodMap[b->getGood(j)].insert(i);
        }
    }

    //-- make conflict sets (aka adjacency lists)
    for (i = 0; i < num_bids; i++) {
        Bid* b = bid[i];

        for (j = 0; j < b->num_goods; j++) {
            for (set<int>::iterator it = goodMap[b->getGood(j)].begin(); 
                    it != goodMap[b->getGood(j)].end(); it++) {
                if ((*it) != i && bidGraph[i][*it] == 0) {
                    b->conflicts++;
                    bid[*it]->conflicts++;
                    bidGraph[i][*it] = 1;
                    bidGraph[*it][i] = 1;
                }
            }
        }

        //b->conflicts.erase(b);
    }
}

