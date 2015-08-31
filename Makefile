SUFFIXES += .d

# Which libraries to include
USERLIBS= $(shell root-config --libs) -lRooFitCore -lRooFit -lFoam -lMinuit

USERLIBDIRS = $(root-config --libdir) 

# Where to search for the include files.
INCDIRS = -I$(shell root-config --incdir) -I../external/include/ 

CXX = g++
LD = g++
CFLAGS = -fPIC $(shell root-config --cflags)
LDFLAGS = $(shell root-config --ldflags) -L../external/lib -Lobjs/ 

# Here's where you'd want to define any #define with -DVariableName
# (which means you can now use '#ifdef VariableName' preprocessor
# lines in your source code).
CFLAGS += $(INCDIRS) -std=c++0x -O2 -g -Wall -Wpointer-arith -Wunused-variable -Wextra

# Everything it needs to know about libraries
DYNAMIC_LIBS = $(USERLIBDIRS) $(USERLIBS) 
LIBS =  $(DYNAMIC_LIBS)

# Define what to compile as well as the necessary dependencies and
# object names.
SRCS =  $(wildcard *.cc) 
DEPS =  $(patsubst %.cc,objs/%.d,$(SRCS)) 
OBJS =  $(patsubst %.cc,objs/%.o,$(SRCS))

link = echo "Linking $@"; $(LD) $(LDFLAGS) $< $(LIBS) -o $@

all: tau_fit 

# Make the executable
.SECONDEXPANSION:

OBJ = objs/$(@).o #objs/libUtils.a

tau_fit : $$(OBJ)
	@$(link)


# Make the objects
objs/%.o: %.cc
	@echo "Compiling $<"
	@$(CXX) $(CFLAGS) -c -o $@ $<

# Make the dependencies
objs/%.d: %.cc
	@echo "Generating dependencies for $<"
	@$(CXX) $(CFLAGS) -MM -MT '$(patsubst %.cc,objs/%.o,$<)' $< -MF $@

# Clean everythingg
clean:
	@rm -f core* ${DEPS} ${OBJS} objs/*.a

# If we aren't doing a 'gmake clean' then make the dependencies
ifneq ($(MAKECMDGOALS), clean)
-include $(DEPS)
endif

