# | Makefile for building ADMB models including additional objects.
EXEC   = om
TPL    = $(EXEC).tpl
OBJS   = MSYReferencePoints.cpp OperatingModel.cpp HarvestControlRule.cpp
CC     = g++
LL     = $(CC)


all: $(EXEC)

$(EXEC): $(TPL) $(OBJS)
	admb $(EXEC) $(OBJS)

clean:
	rm -f $(EXEC) *.eva *.log *.rpt *.htp *.cor *.par *.std *.r* *.p0* *.b* *.o

.PHONY: help

help:
	@echo
	@echo USAGE 
	@echo 'make        : -compiles EXEC'
	@echo 'make clean  : -remove EXEC, intermediate files, and model output.'
	@echo


