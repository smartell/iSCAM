## Make file for building iSCAM
## Author: Steven Martell
# SRCS = msy.cpp OpMod.cpp Selex.cpp baranov.cpp
EXEC = iscam
OBJS = $(EXEC).o msy.o OpMod.o Selex.o baranov.o
SRCS = $(OBJS:.o=.cpp)  #!!!
TPLS = $(EXEC).tpl


# identify some extra file name suffixes
.SUFFIXES: .tpl .cpp .o .obj

# tell make not to delete these intermediate targets
.PRECIOUS: %.cpp %.o %.obj

# make some special targets
.PHONY: help info clean depend dust

# Export the path to your ADMB dist directory
export ADMB_HOME=/Users/stevenmartell1/admb-trunk/build/dist


# establish the C++ compiler
CC=g++
# and linker
LL=$(CC)

ifeq ($(OPT),TRUE)
  CC_OPT = -O3 -DOPT_LIB
  LFLAGS = $(ADMB_HOME)/lib/libadmbo.a $(ADMB_HOME)/contrib/lib/libcontribo.a
else
  CC_OPT = -Wall -ggdb -DSAFE_ALL
  LFLAGS = $(ADMB_HOME)/lib/libadmb.a $(ADMB_HOME)/contrib/lib/libcontrib.a
endif

#CFLAGS = -c -O3 -D__GNUDOS__ -Dlinux -DSAFE_ALL -DUSE_LAPLACE -I. -I$(ADMB_HOME)/include -I$(ADMB_HOME)/contrib/include
CFLAGS = -c $(CC_OPT) -D__GNUDOS__ -Dlinux -DUSE_LAPLACE -I. -I$(ADMB_HOME)/include -I$(ADMB_HOME)/contrib/include

all: $(EXEC)

run: $(EXEC)
	./$(EXEC)	


$(EXEC): $(EXEC).cpp .depend $(OBJS)
	$(LL) $(LFLAGS) -o $@ $(OBJS)
	

depend: .depend

.depend: cmd = $(CC) -MM -MF depend $(var); cat depend >> .depend;
.depend: 
	@echo "Generating dependencies ..."
	@$(foreach var, $(SRCS), $(cmd))
	@rm -f depend

-include .depend

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%: %.cpp
	$(CC) $(CFLAGS) -o $@ $<

$(EXEC).cpp: $(TPLS)
	@echo ------ Making tpl2cpp: $(TPLS)    ------
	$(ADMB_HOME)/bin/tpl2cpp $(TPLS:.tpl=)
	gsed -i '/friend class model_parameters;/ a\
    \ \ friend class OperatingModel;' $(EXEC).htp





help:
	@echo No target specified.
	@echo Usage: make target, e. g. make simple
	@echo        specify OPT=TRUE to build optimized version

# generate some information about what your are doing
rules:
	@echo EXEC = $(EXEC)
	@echo OBJS = $(OBJS)
	@echo SRCS = $(SRCS)
	@echo TPLS = $(TPLS)
	@echo OPT = $(OPT)
	@echo CC_OPT = $(CC_OPT)
	@echo PWD = $(PWD)
	@echo OSTYPE = $(OSTYPE)
	@echo OS     = $(OS)
	@echo ADMB_HOME = $(ADMB_HOME)
	@echo CC = $(CC)
	$(CC) --version
	@echo LL = $(LL)
	@echo CFLAGS = $(CFLAGS)
	@echo LFLAGS = $(LFLAGS)

clean:
	rm -f $(EXEC).[bcehlps]* $(EXEC) $(OBJS) 

dust: 
	rm -f *.r* *.p0* *.b* *.cor *.b0* *.eva *.std *.par *.log

