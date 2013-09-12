EXEC = om
OBJS = $(EXEC).o \
	MSYReferencePoints.o \
	OperatingModel.o     \
	HarvestControlRule.o \
	Scenario.o\
	EstimatorClass.o\
	LRGS.o
SRCS = $(OBJS:.o=.cpp)  #!!!
TPLS = $(EXEC).tpl

# identify some extra file name suffixes
.SUFFIXES: .tpl .cpp .o .obj

# tell make not to delete these intermediate targets
.PRECIOUS: %.cpp %.o %.obj

# make some special targets
.PHONY: help rules clean depend

# does this if no target is specified on command line

export ADMB_HOME=/Users/stevenmartell1/admb-trunk/build/dist

# establish the C++ compiler
CC=g++
# and linker
LL=$(CC)

ifeq ($(OPT),TRUE)
  LFLAGS = $(ADMB_HOME)/lib/libadmbo.a $(ADMB_HOME)/contrib/lib/libcontribo.a
else
  LFLAGS = $(ADMB_HOME)/lib/libadmb.a $(ADMB_HOME)/contrib/lib/libcontrib.a
endif

CFLAGS = -c -O3 -D__GNUDOS__ -Dlinux -DSAFE_ALL -DUSE_LAPLACE -I. -I$(ADMB_HOME)/include -I$(ADMB_HOME)/contrib/include


all: $(EXEC)


$(EXEC): $(EXEC).cpp .depend $(OBJS)
	$(LL) $(LFLAGS) -o $@ $(OBJS)
	

depend: .depend

.depend: cmd = $(CC) -MM -MF depend $(var); cat depend >> .depend;
.depend: 
	@echo "Generating dependencies ..."
	@$(foreach var, $(SRCS), $(cmd))
	@rm -f depend

-include .depend

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%: %.cpp
	$(CC) $(CFLAGS) -o $@ $<

$(EXEC).cpp: $(TPLS)
	@echo ------ Making tpl2cpp: $(TPLS)    ------
	$(ADMB_HOME)/bin/tpl2cpp $(TPLS:.tpl=)

# $(EXEC): $(OBJS);
# 	$(LL) $(LFLAGS) -o $@ $^ 

# $(OBJS): $(SRCS);
# 	@echo 
# 	@echo ------ Compiling: $@      ------
# 	@echo
# 	$(CC) $(CFLAGS) -o $@ $<
# 	# $(CC) $(CFLAGS) -o $(SRCS) $(OBJS)

# $(SRCS): $(TPLS);
# 	@echo 
# 	@echo ------ Making tpl2cpp: $(TPLS)    ------
# 	@echo
# 	$(ADMB_HOME)/bin/tpl2cpp $(TPLS:.tpl=)

# *** adcomp -s hello.cpp
# c++ -c -O3 -D__GNUDOS__ -Dlinux -DSAFE_ALL -DUSE_LAPLACE -I. -I"/Users/stevenmartell1/admb-trunk/build/dist"/include -I"/Users/stevenmartell1/admb-trunk/build/dist"/contrib/include hello.cpp

# *** adlink -s hello
# c++ -o hello hello.o "/Users/stevenmartell1/admb-trunk/build/dist"/contrib/lib/libcontrib.a "/Users/stevenmartell1/admb-trunk/build/dist"/lib/libadmb.a

# link the object file into the executable 
# %:%.o;
# 	$(LL) -o $@ $^ $(LFLAGS)

# # compile the C++ into an object file
# %.o:%.cpp;
# 	$(CC) $(CFLAGS)  -o $@ $<

# # translate the ADMB template into C++ code
# %.cpp:%.tpl;
# 	@echo 
# 	@echo ------ Making tpl2cpp: $(TPLS)    ------
# 	@echo
# 	$(ADMB_HOME)/bin/tpl2cpp $*



help:
	@echo No target specified.
	@echo Usage: make target, e. g. make simple
	@echo        specify OPT=TRUE to build optimized version

# generate some information about what your are doing
rules:
	@echo EXEC = $(EXEC)
	@echo OBJS = $(OBJS)
	@echo SRCS = $(SRCS)
	@echo TPLS = $(TPLS)
	@echo OPT = $(OPT)
	@echo CC_OPT = $(CC_OPT)
	@echo PWD = $(PWD)
	@echo OSTYPE = $(OSTYPE)
	@echo OS     = $(OS)
	@echo ADMB_HOME = $(ADMB_HOME)
	@echo CC = $(CC)
	$(CC) --version
	@echo LL = $(LL)
	@echo CFLAGS = $(CFLAGS)
	@echo LFLAGS = $(LFLAGS)

clean:
	rm -f $(EXEC).[bcehlps]* $(EXEC) $(OBJS) 
