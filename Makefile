LIBSRC=src
LIBRARY=$(LIBDIR)/libff.a

CMDSRC=$(wildcard cmd/*)

BINDIR:=bin
OBJDIR:=obj
INCDIR:=include
LIBDIR:=lib

ARCHIVER:=ar rcs
COMPILER:=g++ -c
LINKER:=g++

LIBS+=-lm
FLAGS+=-Wall -pedantic

# Derived Library Variables
LIBSOURCES=$(wildcard $(addsuffix /*.cpp,$(LIBSRC)))
LIBOBJ=$(subst src,$(OBJDIR),$(LIBSOURCES:%.cpp=%.o))

# Derived Executable Variables
CMDSOURCES=$(addsuffix /main.cpp,$(CMDSRC))
CMDOBJ=$(CMDSOURCES:%.cpp=%.o)
CMDBIN=$(subst cmd,$(BINDIR),$(CMDSRC))

all: FLAGS:=-Wall -pedantic -g -pg
all: $(LIBRARY)

release: FLAGS:=-Wall -pedantic -O3
release: $(CMDBIN)

exec: FLAGS:=-Wall -pedantic -g -pg
exec: $(CMDBIN)

# Archiving, Compiling and Cleaning the library

$(LIBRARY): $(LIBOBJ)
	@echo "Creating the FF Framework Library"
	@mkdir -p $(@D)
	@$(ARCHIVER) $(LIBRARY) $(LIBOBJ)

obj/%.o: src/%.cpp
	@echo "Compiling ($(FLAGS)) $@"
	@mkdir -p $(@D)
	@$(COMPILER) -I$(INCDIR) $(FLAGS) -o $@ $<

clean-lib:
	@echo "Cleaning FF Framework objects and static libraries"
	@rm -f $(LIBRARY) $(LIBOBJ)

# Linking, Compiling and Cleaning the executables

$(BINDIR)/%: LIBS+=-lff
$(BINDIR)/%: cmd/%/main.cpp $(LIBRARY)
	@echo "Compiling and linking ($(FLAGS)) $@"
	@mkdir -p $(BINDIR)
	@$(LINKER) -I$(INCDIR) -L$(LIBDIR) $(FLAGS) -o $@ $< $(LIBS)

clean-exec:
	@echo "Cleaning command object and binary files"
	@rm -f $(CMDBIN)

clean: clean-lib clean-exec

# Making Documentation
doc:
	@doxygen
	@make -C share/latex

clean-doc:
	@echo "Cleaning documentation"
	@rm -rf share

.PHONY: clean clean-lib clean-exec clean-doc
