OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR) -I$(HEEDDIR)

# Debug flags
# CFLAGS += -g

LDFLAGS = `root-config --glibs` -lGeom -lgfortran -lm
LDFLAGS += -L$(LIBDIR) -lGarfield
# LDFLAGS += -g

PPAC: PPAC.C 
	$(CXX) $(CFLAGS) PPAC.C
	$(CXX) -o PPAC PPAC.o $(LDFLAGS)
	rm PPAC.o

gasfile: gasfile.C
	$(CXX) $(CFLAGS) gasfile.C
	$(CXX) -o gasfile gasfile.o $(LDFLAGS)
	rm gasfile.o

