# @(#)makefile	2-99
# $Revision:$

.SUFFIXES: .o .C .c .cpp

# FLAGS=  -O2
FLAGS= -g

# NOTE: Set path to bamtools here
TOOLS=/usr/bin

HAPFLAGS= -Wl,-rpath,$(TOOLS)/bamtools-master/lib/

SRCDIR = .
BINDIR = .
OBJDIR = obj/

HEADER = read.h  snp.h nb.h

OBJ = read.o snp.o nb.o

OBJ3 = hapmut.o

#===============================================================================

hapmut : $(HEADER) $(OBJ) $(OBJ3)
	g++ -g $(FLAGS) $(HAPFLAGS) -o $(BINDIR)/hapmut $(OBJ) $(OBJ3) -L$(TOOLS)/bamtools-master/lib/ -lm -lbamtools

#===============================================================================

#===============================================================================

snp.o: snp.h

read.o: read.h

nb.o: snp.h read.h nb.h

hapmut.o: nb.h read.h snp.h

#===============================================================================

#===============================================================================

.cpp.o: $(HEADERI)
	g++ -c -g -I$(TOOLS)/bamtools-master/include $(FLAGS) $<
.C.o: $(HEADERI)
	g++ -c -g $(FLAGS) $<

clean:
	rm -f *.o core

tags: $(SRC)
	ctags *.c

lint:	rm -f lint.out
	lint *.c > lint.out

tar.all: $(SRC) $(DOC) $(HEADERS) makefile
	tar cvf tar.all $(SRC) $(DOC) $(HEADERS) makefile

#===============================================================================
