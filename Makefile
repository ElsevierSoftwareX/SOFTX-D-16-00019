LIBS = -lm -lc -lSDL2 -lgomp -lSDL2_gfx
LIBSPATH =
LDLIBS = $(LIBS)
DEPLIBS =
CC = gcc
C++ = g++
CPP = cc -E
LD = ld
CFLAGS = -O3 -fopenmp
CCLINK = $(CC)
CXXLINK = $(CXX)
OBJS = segy-change.o

default:
	@$(MAKE) segy-change \
	"CFLAGS = -O3 -fopenmp -DWITH_SDL"

nosdl:
	@$(MAKE) segy-change \
	"CFLAGS = -O3 -fopenmp" \
	"LIBS =  -lm -lc -lgomp"

debug:
	@$(MAKE) segy-change \
	"CFLAGS = -g -fopenmp -DWITH_SDL"

segy-change: $(OBJS) $(DEPLIBS)
	$(RM) $@
	$(CCLINK) -fopenmp -o $@ $(LDOPTIONS) $(OBJS) $(LOCAL_LIBRARIES) $(LDLIBS) $(LIBS) $(EXTRA_LOAD_FLAGS)

clean:
	$(RM) *.CKP *.ln *.BAK *.bak *.o core errs *~ *.a segy-change

.depend:
	for i in *.c;do $(CPP) -M $$i;done >.depend
