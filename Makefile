CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
CFLAGS = -Wall -I$(CFITSIO) $(shell root-config --cflags) -O3
LIBS = -L$(CFITSIO) -lcfitsio -lm $(shell root-config --glibs) -lTreeViewer
GLIBS = 
GLIBS += 
OBJECTS = qNet2root6000_GPS.o 
HEADERS = 

ALL : qNet2root6000_GPS.exe
	echo "Listo!"

qNet2root6000_GPS.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o qNet2root6000_GPS.exe $(LIBS) $(GLIBS) $(CFLAGS)

qNet2root6000_GPS.o : qNet2root6000_GPS.cc $(HEADERS)
	$(CPP) -c qNet2root6000_GPS.cc -o qNet2root6000_GPS.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe
