CFILES += $(shell find ./*/ -name '*.cpp')
HFILES += $(shell find ./*/ -name '*.h')

all:
	find ./*/ -name Makefile -type f -execdir make all \;

clean:
	find ./*/ -name Makefile -type f -execdir make clean \;

stylefix:
	astyle --style=allman --indent=force-tab=4 --indent-switches --pad-oper --pad-header --unpad-paren $(CFILES) $(HFILES)
