# Makefile for running scenarios and management procedures.
# How this works:
# 	1) Loop over SCEN and copy SCEN file to Simulation Directory
#	2) Loop over PROC and copy PROC file to Simulation Directory
#   3) Run Realizations
#   4) Read R files
#   5) Concatenate SCEN & PROC to Rdatabase
SCEN = S1.scn S2.scn
PROC = MP1.mp MP2.mp MP3.mp

scendone:=$(foreach dir,$(SCEN),$(dir)/scendone)
procdone:=$(foreach dir,$(SCEN),$(foreach sdir,$(PROC),job-$(dir)-$(sdir)))
dir = $(firstword $(subst -, ,$*))
sdir = $(lastword $(subst -, ,$*))


all: $(scendone) $(procdone)



$(scendone):
	-mkdir $(@D);
	cd $(@D); touch scendone

$(procdone): job-%:
	@echo The directory is $(dir) and sdir is $(sdir);
	-cd $(dir) && mkdir $(sdir);
	-cp SIMDIR/om.dat SIMDIR/Makefile $(dir)/$(sdir)/
	-cp SIMDIR/$(dir) $(dir)/$(sdir)/Scenario.scn
	-cp SIMDIR/$(sdir) $(dir)/$(sdir)/ManagementProcedure.mp
	cd $(dir)/$(sdir); touch procdone




clean:
	rm -fr *.scn *.mp