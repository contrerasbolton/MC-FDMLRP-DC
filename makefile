CPLEX    = /opt/ibm/ILOG/CPLEX_Studio128
CPP      = g++
#DEFINES  = -DIL_STD -DNDEBUG -std=c++14 -g -ggdb
DEFINES  = -DIL_STD -std=c++14
CXXFLAGS = -O3 -Wall -fPIC -fno-strict-aliasing $(DEFINES) -Wno-reorder -Wno-ignored-attributes
INCPATH  = -I$(CPLEX)/cplex/include -I$(CPLEX)/concert/include
OBJDIR   = obj
LFLAGS   = -Wl,-O3
OBJ      = $(addprefix $(OBJDIR)/, $(patsubst %.cpp, %.o, $(wildcard *.cpp)))
TARGET   = drone
LIBS1    = -L$(CPLEX)/concert/lib/x86-64_linux/static_pic -lconcert
LIBS2    = -L$(CPLEX)/cplex/lib/x86-64_linux/static_pic -lilocplex -lcplex
LIBS3    = -lm -lpthread -lgomp -fopenmp -ldl
LIBS     = $(LIBS1) $(LIBS2) $(LIBS3)

.PHONY: all clean

all: $(OBJDIR) $(TARGET)

$(OBJDIR):
	mkdir $(OBJDIR)

$(OBJDIR)/%.o: %.cpp
	$(CPP) $(CXXFLAGS) $(INCPATH) -c $< -o $@

$(TARGET): $(OBJ)
	$(CPP) $(CPPFLAGS) -o $@ $^ $(LFLAGS) $(LIBS)

clean:
	@rm -f $(TARGET) $(wildcard *~) $(wildcard *#)
	@rm -rf $(OBJDIR)
