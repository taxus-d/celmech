####### Files

SRCDIR  = src
# TESTDIR = test
#DEPDIR  = .d
#$(shell mkdir -p $(DEPDIR)/$(SRCDIR) >/dev/null)

SOURCES = \
	  Prec.f90 \
	  IO.f90 \
	  IO_array.f90 \
	  Celmech.f90 \
	  Inival.f90 \
	  FuncIfaces.f90 \
	  Poly.f90 \
	  Utils.f90 \
	  LseSolvers.f90 \
	  Deriv.f90 \
	  NewtSolve.f90 \
	  Integrators.f90 \
	  ODEsolve.f90 \
	  Poincare.f90 \
	  cauchy.f90
SOURCES_RAW = $(addprefix $(SRCDIR)/, $(SOURCES))
OBJECTS  = ${SOURCES_RAW:.f90=.o}
INCPATH   = $(SRCDIR)


# TESTS_RAW   = testmul.f95
# TESTS    = $(addprefix $(TESTDIR)/, $(TESTS_RAW))
# TESTS_EX = $(basename $(TESTS_RAW))

TARGET = cauchy

####### Compiler, tools and options
ifndef CN
CN = gnu
endif

ifeq ($(CN),gnu)
CC       = gfortran
# Flags, splitted for flexibility
CSTD     = -std=f2008 -pedantic
WARN  = -Wall \
	-Wimplicit-interface \
	-Wcompare-reals \
#         -Warray-temporaries

PARALLEL  = #-fopenmp
DEBUG     = -pg
RELEASE   = -O2
MODE      = $(DEBUG)
CFLAGS    = -fmax-errors=10 -fcheck=all -fbacktrace -fall-intrinsics \
	    -ftree-vectorize  -march=native  -ffast-math -cpp \
	    $(CSTD) $(MODE) $(WARN_GF) $(PARALLEL)
LFLAGS    = $(PARALLEL) $(MODE)
DEPFLAGS  = -MT $@ -MMD -MP -MF $*.Td
MODOUT    = -J$(SRCDIR)
LIBS      =
POSTCOMPILE = @mv -f $*.Td $*.d && touch $@
endif

ifeq ($(CN), intel)
CC        = ~/opt/intel/bin/ifort
WARN   =
PARALLEL  = #-qopenmp
SANITY =  -no-wrap-margin
CFLAGS    = -free -p -g -stand f08  -fpp $(PARALLEL) $(SANITY) -xcore-avx2  #-ffast-math
DEPFLAGS  = -gen-dep=$*.Td
LFLAGS    = $(PARALLEL)
MODOUT    = -module $(SRCDIR)
POSTCOMPILE = @mv -f $*.Td $*.d && touch $@

endif
ifndef CC
$(error "Unknown compiler")
endif

COPY      = cp -f
COPY_FILE = $(COPY) -p
COPY_DIR  = $(COPY) -pR
DEL_FILE  = rm -f
####### Implicit rules

# % -- общая часть
# $@ -- target
# $< -- первый элемент в списке зависимостей
# $+ -- список всех зависимостей
# $* -- то, что сопоставилось в %

#.SUFFIXES: .c .cpp .cc .cxx .C

%.d: ;
%.o: %.f90 %.d
	$(CC) -c $(DEPFLAGS) $(CFLAGS) -I$(INCPATH) $< -o $@ $(MODOUT)
	$(POSTCOMPILE)

####### Build rules

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) -o $@ $(LFLAGS) $+


test: all
	./$(TARGET)

clean:
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) $(SRCDIR)/*.mod
	-$(DEL_FILE) $(SRCDIR)/*.d
	-$(DEL_FILE) *.mod *.o

distclean: clean
	-$(DEL_FILE) $(TARGET)

backup: distclean
	easybackup -f $(SOURCES_RAW) Makefile -

dummy:
	@echo $(OBJECTS)
view:
	@gnuplot viz/orbit.plt&

.PHONY: test clean distclean backup dummy
.PRECIOUS: %.d

include $(wildcard $(patsubst %, %.d, $(basename $(OBJECTS))))
