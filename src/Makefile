# | Makefile for building ADMB models including additional objects.
EXEC = om
TPL  = $(EXEC).tpl
OBJS = MSYReferencePoints.cpp OperatingModel.cpp HarvestControlRule.cpp
MCA  = -mcmc 500000 -mcsave 200 -nosdmcmc -mcscale 100000
CC   = g++
LL   = $(CC)


all: $(EXEC)

$(EXEC): $(TPL) $(OBJS)
	admb $(EXEC) $(OBJS)

run: $(EXEC)
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
	rm -f *.rep *.bar *.psv *.std $(EXEC).cpp admodel.* variance *.hst *.ecm *.mc2
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


