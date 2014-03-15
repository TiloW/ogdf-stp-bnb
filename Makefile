include Makefile.config
CC = g++
ifeq ($(DEBUG), 1)
  OGDF_LIBRARY_DIR = $(OGDF_LIBRARY_DIR_DEBUG)
else
  OGDF_LIBRARY_DIR = $(OGDF_LIBRARY_DIR_RELEASE)
endif
INCLUDES = -I$(OGDF_INCLUDE_DIR)
CXXFLAGS = -std=c++11 $(INCLUDES) -Wall -Wno-unused-parameter
LDFLAGS = -pthread -lm -L$(OGDF_LIBRARY_DIR)
ifeq ($(DEBUG),1)
  CXXFLAGS += -O0 -g -ggdb -W -DOGDF_DEBUG
else
  ifeq ($(DEBUG), 2)
    # profiler mode but release
    CXXFLAGS += -O3 -g -ggdb -pg
    LDFLAGS += -pg
  else
    CXXFLAGS += -O3
  endif
endif
LDFLAGS = -pthread -lm -L$(OGDF_LIBRARY_DIR)
LDLIBS = -lOGDF -lCOIN

BINARIES = stp-solver
OBJS = 

run: stp-solver
	./stp-solver data/example.stp

run-ext: stp-solver
	./stp-solver data/B/b01.stp

all: $(BINARIES)

stp-solver: stp-solver.o

%o: %cc

Makefile.config:
	$(error Need to know OGDF, CPLEX and CONCERT paths. Please rename Makefile.config.default to Makefile.config and edit appropriately)

.PHONY: clean distclean
clean:
	$(RM) $(BINARIES) $(BINARIES:%=%.o) $(OBJS)

distclean: clean
	$(RM) Makefile.config
