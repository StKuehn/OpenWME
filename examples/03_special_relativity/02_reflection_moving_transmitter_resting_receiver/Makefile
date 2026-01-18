# target
BIN = reflection_moving_transmitter_resting_receiver

OPENWME_LIBRARYDIR = "../../../library"

CFILES += $(wildcard *.cpp)
HFILES += $(wildcard *.h)
CFILES += $(shell find $(OPENWME_LIBRARYDIR) -name '*.cpp')
HFILES += $(shell find $(OPENWME_LIBRARYDIR) -name '*.h')

OBJDIR = obj
OBJPATTERNS = $(notdir $(CFILES:.cpp=))
OBJS = $(addprefix $(OBJDIR)/,$(notdir $(CFILES:.cpp=.o)))

# compiler
CPP = g++

# flags
LIBS = `pkg-config --libs x11 cairo` -lm
INCS = `pkg-config --cflags x11 cairo` -I$(OPENWME_LIBRARYDIR)
FLAGS = $(INCS) -Wall -fverbose-asm

ifdef DEBUG
	FLAGS += -g2
else
	FLAGS += -O3 -fopenmp
	LDFLAGS += -fopenmp
endif

# create all new
all: clean target 

target: $(BIN)

$(BIN): $(OBJPATTERNS)
	$(CPP) $(OBJS) -o $(BIN) $(LDFLAGS) $(LIBS)

$(OBJPATTERNS): $(CFILES)
	@mkdir -p $(OBJDIR)
	$(CPP) $(FLAGS) -c $(filter %$@.cpp,$(CFILES)) -o $(OBJDIR)/$@.o

stylefix:
	astyle --style=allman --indent=force-tab=4 --indent-switches --pad-oper --pad-header --unpad-paren $(CFILES) $(HFILES)

# clean up all
clean:
	rm -rf obj
	rm -rf *.orig
	rm -f $(BIN)

