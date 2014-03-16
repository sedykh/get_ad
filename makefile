CC=g++
CFLAGS=-c -m32 -O3 -Wall -static
all: get_ad

get_ad: get_ad.o apstring.o set.o core.o qsar.o dataset.o
	$(CC) get_ad.o apstring.o set.o core.o qsar.o dataset.o -m32 -o get_ad

apstring.o: apstring.cpp
	$(CC) $(CFLAGS) apstring.cpp
set.o:
	$(CC) $(CFLAGS) set.cpp

core.o: core.cpp
	$(CC) $(CFLAGS) core.cpp

qsar.o: qsar.cpp
	$(CC) $(CFLAGS) qsar.cpp

dataset.o: dataset.cpp
	$(CC) $(CFLAGS) dataset.cpp

get_ad.o: get_ad.cpp
	$(CC) $(CFLAGS) get_ad.cpp

clean:
	rm -rf *.o
