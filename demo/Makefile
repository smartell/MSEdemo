EXEC  = om
DATA  = om.dat
DIST  = ../dist/
OBJ   = ../dist/objects/
NOSIM = 8

all:  $(EXEC) $(EXEC).rep

$(EXEC): $(DIST)
	cp $(DIST)$(EXEC) $@
	# cp $(DIST)$(DATA) $(DATA)

$(EXEC).rep:
	./$(EXEC)

<<<<<<< HEAD
runmse:
=======
mse: $(EXEC)
	./$(EXEC) -ind $(DAT) -mse 12345

runmse: 
>>>>>>> 38392d7268fe339bae0e3e020bb49e7ed3718fc8
	./$(EXEC) -ind $(DAT) $(ARG)

clean:
	rm -f $(EXEC) $(EXEC).rep *.[bcelprs]* *.hes variance admodel.*
	rm -f simdone


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

COLLECTALL='source(file.path("..","src","r","collectAll.r"))'

allSims.Rdata: simdone
	@echo $(COLLECTALL) | R --vanilla --slave

collect: allSims.Rdata



