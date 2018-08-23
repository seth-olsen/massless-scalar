CC = g++
CXX = g++
DEBUG = -g
OPTFLAGS = -O2
LDFLAGS = -Wall $(DEBUG) $(OPTFLAGS)
LDLIBS = -lbbhutil
CFLAGS = -std=c++11 -c $(LDFLAGS) $(LDLIBS)
CXXFLAGS = -std=c++11 -c $(LDFLAGS) $(LDLIBS)

p1 : p1.o

p1m : p1m.o

p1-ctest : p1-ctest.o

p1-new : p1-new.o

p1-new.o : fda-fns.h fda-io.h

.PHONY : clean
clean :
	rm -f *.o *~
