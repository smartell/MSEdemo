# | Makefile for building ADMB models including additional objects.


EXEC  = om
TPL   = $(EXEC).tpl
SRCS  = MSYReferencePoints.cpp \
		OperatingModel.cpp     \
		HarvestControlRule.cpp \
		Scenario.cpp           
OBJS  = ${SRCS:.cpp=.o}
HDRS  = ${SRCS:.cpp=.h}
MCA   = -mcmc 500000 -mcsave 200 -nosdmcmc -mcscale 100000
CC    = g++
LL    = $(CC)
NOSIM = 8


all: $(EXEC)


$(EXEC): $(TPL) $(SRCS) $(HDRS)
	admb $(EXEC) $(SRCS)

run: $(EXEC)
	./$(EXEC) -ind $(DAT) $(ARG)

runmse: 
	./$(EXEC) -ind $(DAT) $(ARG)

mse: $(EXEC)
	./$(EXEC) -ind $(DAT) $(ARG) -mse 12345

mcmc: runmc mceval

runmc: $(EXEC) 
	./$(EXEC) -ind $(DAT) $(MCA) $(ARG)

mceval: $(EXEC).psv
	./$(EXEC) -ind $(DAT) $(ARG) -mceval

clean:
	rm -f $(EXEC) *.eva *.log *.rpt *.htp *.cor *.par *.std *.r* *.p0* *.b* *.o
	rm -f *.rep *.bar *.psv *.std $(EXEC).cpp admodel.* variance *.hst *.ecm *.mc2 *.mcm
	rm -f hybrid_seed sims NUL simdone

.PHONY: help

help:
	@echo
	@echo USAGE 
	@echo 'make        : -compiles EXEC'
	@echo 'make run    : -compiles and executes code'
	@echo 'make mcmc   : -compiles ane executes code in mcmc and mceval modes.'
	@echo 'make mse    : -compiles and runs Management Strategy Evaluation.'
	@echo 'make clean  : -remove EXEC, intermediate files, and model output.'
	@echo


# |------------------------------------------------------------------------------------|
# | BUILD DIRECTORIES AND COPY FILES INTO EACH DIRECTORY (target = data)
# |------------------------------------------------------------------------------------|
# |
simdirs := $(shell echo 'cat(formatC(1:$(NOSIM), digits=3, flag="0"))' | R --vanilla --slave)
datadone:= $(foreach dir,$(simdirs),$(dir)/datadone)
runsims := $(foreach dir,$(simdirs),$(dir)/runsims)

$(datadone): 
	mkdir $(@D);
	cp $(EXEC) $(EXEC).dat Makefile $(@D)
	cd $(@D); touch datadone

data: $(EXEC) $(datadone)

$(runsims):
	cd $(@D); make runmse ARG='-mse $(@D)'

sims: data $(runsims)
	touch simdone


cleansims:
	rm -r 0* *.Rdata


# |------------------------------------------------------------------------------------|
# | R SCRIPTS FOR READING REPORT & PARAMETER FILES AND CREATING AN .Rdata FILE.
# |------------------------------------------------------------------------------------|
# |

COLLECTALL='dn<-dir(pattern="^[[:digit:]]");          \
	        sims <- lapply(dn,function(d){            \
		                       require(Riscam);       \
		                       setwd(d);              \
							   A<-read.admb("OM");    \
							   setwd("..");           \
							   return(A)     });      \
	        save(sims,file="allSims.Rdata")'

allSims.Rdata: simdone
	@echo $(COLLECTALL) | R --vanilla --slave

collect: allSims.Rdata


