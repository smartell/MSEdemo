# SUBDIR= HCR_NOTHEFT/fortyten   			    \
# 		HCR_NOTHEFT/fixedEsc    			\
# 		HCR_NOTHEFT/fixedEscCap 			\
# 		HCR_NOTHEFT/fixedHarvestRate 		\
# 		HCR_NOTHEFT/fixedCCC 				\
# 		HCR_NOTHEFT/thirtytwenty  			\
# 		HCR_NOTHEFT/thirtytwentyfloor		\
# 		HCR_NOTHEFT/fixedHRdelta  			\
# 		HCR_THEFT/fortyten 					\
# 		HCR_THEFT/fixedEsc 					\
# 		HCR_THEFT/fixedEscCap 				\
# 		HCR_THEFT/fixedHarvestRate 			\
# 		HCR_THEFT/fixedCCC 					\
# 		HCR_THEFT/thirtytwenty  			\
# 		HCR_THEFT/thirtytwentyfloor			\
# 		HCR_THEFT/fixedHRdelta  			\
# 		PDO_NOTHEFT/fortyten 				\
# 		PDO_NOTHEFT/fixedEsc 				\
# 		PDO_NOTHEFT/fixedEscCap 			\
# 		PDO_NOTHEFT/fixedHarvestRate 		\
# 		PDO_NOTHEFT/fixedCCC 				\
# 		PDO_NOTHEFT/thirtytwenty 			\
# 		PDO_NOTHEFT/thirtytwentyfloor		\
# 		PDO_NOTHEFT/fixedHRdelta  			\
# 		PDO_THEFT/fortyten 					\
# 		PDO_THEFT/fixedEsc 					\
# 		PDO_THEFT/fixedEscCap 				\
# 		PDO_THEFT/fixedHarvestRate 			\
# 		PDO_THEFT/fixedCCC 					\
# 		PDO_THEFT/thirtytwenty				\
# 		PDO_THEFT/thirtytwentyfloor			\
# 		PDO_THEFT/fixedHRdelta  

SUBDIR= \
		HCR_NOTHEFT/fixedHRdelta		\
		HCR_THEFT/fixedHRdelta			\
		PDO_NOTHEFT/fixedHRdelta		\
		PDO_THEFT/fixedHRdelta			

TARGET = sims 
ifndef DISK
  DISK = dist
endif

all:
	mkdir -p $(DISK)
	make --directory=./src/
	cp -r ./src/om ./src/om.dat ./src/Makefile $(DISK)

.PHONY: runcmp $(SUBDIR)

runcmp: $(all) $(SUBDIR)
$(SUBDIR):
	cd $@ && $(MAKE) $(TARGET) -j8
	cd $@ && $(MAKE) collect

.PHONY: clean	
clean_files := $(foreach dir,$(SUBDIR),$(dir)/clean)

clean: $(clean_files)
	rm -rvf $(DISK)
	# cd src && $(MAKE) clean

$(clean_files):
	-cd $(@D) && $(MAKE) clean	
	-cd $(@D) && $(MAKE) cleansims


copy_files := $(foreach dir,$(SUBDIR), cp -r OM_Makefile $(dir)/Makefile~)

updatemake:
	$(foreach dir,$(SUBDIR), cp -r OM_Makefile $(dir)/Makefile;)
