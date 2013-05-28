EXEC  = OM
DATA  = OM.dat
DIST  = ../dist/
NOSIM = 8



all:  $(EXEC) $(EXEC).rep


$(EXEC): $(DIST) 
	cp $(DIST)$(EXEC) $@
	cp $(DIST)$(DATA) $(DATA)

$(EXEC).rep:
	./$(EXEC)

runmse: 
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

COLLECTALL='dn<-dir(pattern="^[[:digit:]]");                    \
			nf   <- paste(basename(getwd()),".Rdata",sep="");   \
	        sims <- lapply(dn,function(d){                      \
		                       require(Riscam);                 \
		                       setwd(d);                        \
							   A<-read.admb("OM");              \
							   setwd("..");                     \
							   return(A)     });                \
	        save(sims,file=nf)'

allSims.Rdata: simdone
	@echo $(COLLECTALL) | R --vanilla --slave

collect: allSims.Rdata



