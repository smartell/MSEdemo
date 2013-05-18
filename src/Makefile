#Makefile for building and running the LSMR model on a Linux or Mac OSX
EXEC=OM
TPL=$(EXEC).tpl
DAT=$(EXEC).dat
MCA=-mcmc 500000 -mcsave 200 -nosdmcmc -mcscale 100000
ARG=
NOSIM=8

all: $(EXEC) $(EXEC).rep dust

$(EXEC): $(TPL)
	admb -O $(EXEC)

$(EXEC).rep: $(TPL) $(DAT)
	./$(EXEC) -ind $(DAT) $(ARG)

mse: $(EXEC)
	./$(EXEC) -ind $(DAT) $(ARG) -mse 12345

mcmc: runmc mceval dust

runmc: $(EXEC) $(EXEC).rep
	./$(EXEC) -ind $(DAT) $(MCA) $(ARG)

mceval: $(EXEC).psv
	./$(EXEC) -ind $(DAT) $(ARG) -mceval



dust:
	rm -f *.log *.rpt *.htp admodel.* variance *.bar *.mcm *.[bpr][0123456789]*

clean:
	rm -f $(EXEC) *.eva *.log *.rpt *.htp *.cor *.par *.r* *.p0* *.b*
	rm -f *.rep *.bar *.psv *.std $(EXEC).cpp admodel.* variance *.hst *.ecm *.mc2
	rm -f hybrid_seed sims NUL


# |------------------------------------------------------------------------------------|
# | BUILD DIRECTORIES AND COPY FILES INTO EACH DIRECTORY (target = data)
# |------------------------------------------------------------------------------------|
# |
simdirs := $(shell echo 'cat(formatC(1:$(NOSIM), digits=3, flag="0"))' | R --vanilla --slave)
datadone:= $(foreach dir,$(simdirs),$(dir)/datadone)

$(datadone): 
	mkdir $(@D);
	cp $(EXEC) $(EXEC).dat Makefile $(@D)

data: $(EXEC) $(datadone)

cleansims:
	rm -r 0* *.Rdata