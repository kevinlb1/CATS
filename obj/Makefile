VPATH = ..

all : ${OBJ}# DO NOT DELETE
# DO NOT DELETE

testPoly.o: polyModel.h
arbitrary.o: Distribution.h Param.h polyModel.h regions.h normal.h bid.h
arbitrary.o: BidSet.h arbitrary.h
bid.o: bid.h
BidSet.o: BidSet.h bid.h Param.h polyModel.h Distribution.h
Distribution.o: Distribution.h Param.h polyModel.h
featureCalc.o: featureCalc.h BidSet.h bid.h Param.h polyModel.h
featureCalc.o: Distribution.h
Legacy.o: Legacy.h bid.h normal.h Param.h polyModel.h BidSet.h Distribution.h
main.o: BidSet.h bid.h Param.h polyModel.h Distribution.h Legacy.h normal.h
main.o: regions.h arbitrary.h scheduling.h matching.h paths.h featureCalc.h
matching.o: bid.h BidSet.h Param.h polyModel.h Distribution.h matching.h
normal.o: normal.h Param.h polyModel.h
Param.o: Param.h polyModel.h Legacy.h bid.h normal.h BidSet.h Distribution.h
Param.o: featureCalc.h
paths.o: Param.h polyModel.h Distribution.h bid.h BidSet.h paths.h
polyModel.o: polyModel.h Param.h
regions.o: normal.h Param.h polyModel.h bid.h BidSet.h Distribution.h
regions.o: regions.h
scheduling.o: Distribution.h Param.h polyModel.h BidSet.h bid.h scheduling.h
