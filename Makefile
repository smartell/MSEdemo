SUBDIR = fortyten fixedEsc fixedEscCap fixedHarvestRate fixedCCC
TARGET = sims 
ifndef DISK
  DISK = dist
endif

all: 
	mkdir -p $(DISK)
	make --directory=./src/
	cp -r ./src/om ./src/om.dat ./src/Makefile $(DISK)

.PHONY: runcmp $(SUBDIR)

runcmp: $(SUBDIR)
$(SUBDIR):
	cd $@ && $(MAKE) $(TARGET) -j8
	cd $@ && $(MAKE) collect

.PHONY: clean	
clean_files := $(foreach dir,$(SUBDIR),$(dir)/clean)

clean: $(clean_files)
	rm -rvf $(DISK)

$(clean_files):
	cd $(@D) && $(MAKE) clean	
	cd $(@D) && $(MAKE) cleansims


