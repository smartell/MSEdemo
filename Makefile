SUBDIR = fortyten fixedEsc fixedEscCap fixedHarvestRate fixedCCC
TARGET = sims 
ifndef DISK
  DISK = dist
endif

all: 
	mkdir -p $(DISK)
	make --directory=./src/
	cp -r ./src/OM ./src/OM.dat ./src/Makefile $(DISK)

.PHONY: runcmp $(SUBDIR)

runcmp: $(SUBDIR)
$(SUBDIR):
	cd $@ && $(MAKE) $(TARGET) -j8
	cd $@ && $(MAKE) collect

.PHONY: clean
clean_files := $(foreach dir,$(SUBDIR),$(dir)/clean)

clean: $(clean_files)
	rm -rf $(DISK)

$(clean_files):
	cd $(@D) && $(MAKE) clean	
	cd $(@D) && $(MAKE) cleansims


