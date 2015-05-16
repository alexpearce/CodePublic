ROOTLIBS=$(shell root-config --libs) -lMLP -lXMLIO -lTreePlayer
ROOTCFLAGS=$(shell root-config --cflags)

LDFLAGS= -L./lib -fPIC  -Wl,--no-as-needed

CC=g++
ODIR=./obj
SRCDIR=./src

CFLAGS = -I./include -g

RCXX= $(CFLAGS) $(ROOTCFLAGS)


_Objects_Run=EfficiencyBin.o EfficiencyBinCollection.o Run.o
Objects_Run=$(patsubst %,$(ODIR)/%,$(_Objects_Run))


all: Run

Run: $(Objects_Run)
	$(CC) $(LDFLAGS) $(ROOTLIBS) $(Objects_Run) -o ./bin/Run
	@echo '-> ./bin/Run created!'

$(ODIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) -c -o $@ $< $(RCXX)


clean:
		@rm -f *~
		@rm -f $(ODIR)/*.o 
		@rm -f lib/*.so 
		@rm -f ./bin/*.x
		@rm -f ./dict/*.cxx
		@rm -f ./*/*.pcm
