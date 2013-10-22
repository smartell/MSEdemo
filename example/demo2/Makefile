# Makefile for running scenarios and management procedures.
# How this works:
# 	1) Loop over SCEN and copy SCEN file to Simulation Directory
#	2) Loop over PROC and copy PROC file to Simulation Directory
#   3) Run Realizations
#   4) Read R files
#   5) Concatenate SCEN & PROC to Rdatabase
SCEN = S1.scn S2.scn
PROC = MP1.mp MP2.mp MP3.mp

scenarios = $(wildcard SIMDIR/*.scn)
procedure = $(wildcard SIMDIR/*.mp)

bSCEN := $(notdir $(basename $(scenarios)))
bPROC := $(notdir $(basename $(procedure)))

sdir := $(foreach dir, $(bSCEN),$(dir)/sdir)
pdir := $(foreach dir, $(bSCEN),$(dir)/pdir)
cdir := $(foreach dir, $(bSCEN),$(foreach sdir, $(bPROC), $(dir)/$(sdir)/cdir))
ddir := $(foreach dir, $(notdir $(scenarios)), \
        $(foreach sdir, $(notdir $(procedure)), job-$(dir)-$(sdir)/ddir))

runam := $(foreach dir, $(bSCEN),$(foreach sdir, $(bPROC), $(dir)/$(sdir)/runam))
runmse:= $(foreach dir, $(bSCEN),$(foreach sdir, $(bPROC), $(dir)/$(sdir)/runmse))

.PHONY: run configure

# Run all of the assessment models.
run: $(configure) $(runam)

$(runam):
	@echo $(@D)
	@cd $(@D); make

# Run a Candidate Management Procedure
TARGET = sims
runcmp: $(runmse)
$(runmse):
	cd $(@D) && $(MAKE) $(TARGET) 
	# cd $(@D) && $(MAKE) collect

# Configure script to set up directories
configure: $(sdir) $(pdir) $(cdir) $(ddir)
	

$(sdir): 
	@mkdir -p $(@D)

$(pdir): 
	@cd $(@D) && mkdir -p $(bPROC)

$(cdir): $(scenarios)
	@cp SIMDIR/Makefile $(@D)/Makefile
	@cp SIMDIR/OM.dat   $(@D)/OM.dat
	

s  = $(firstword $(subst -, ,$*))
p  = $(lastword  $(subst -, ,$*))
$(ddir): job-%/ddir:
	@echo $(basename $s)/$(basename $p)
	@cp -f SIMDIR/$s $(basename $s)/$(basename $p)/Scenario.scn
	@cp -f SIMDIR/$p $(basename $s)/$(basename $p)/ManagementProcedure.mp





# |	HELP FOR USERS | #
.PHONY: help
help:
	@echo '___________________________________________________________'                        
	@echo  HELP HELP HELP HELP HELP HELP HELP HELP HELP HELP HELP HELP                        
	@echo '-----------------------------------------------------------'
	@echo 1. run 'make configure' to set up simulation directories.
	@echo 2. run 'make run' to run the assessment models only.
	@echo
	@echo Scenarios:            
	@echo $(bSCEN)
	@echo
	@echo Management Procedures: 
	@echo $(bPROC)
	@echo
	@echo '___________________________________________________________'                        
	

clean:
	rm -fr $(bSCEN) configure