#
#
# Makefile for g2 library
#
#

G2_VERSION = 0.70

#
# g2 installation directories
#
LIBDIR = /usr/local/lib
INCDIR = /usr/local/include


#LIBDIR = $(HOME)/local/lib
#INCDIR = $(HOME)/local/include

#
#------------------------ do not edit ------------------------
#
SHELL = /bin/sh

CC           = gcc
CFLAGS       = -I./src -g -O2  -I/usr/X11R6/include -I/usr/local/include -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DLINUX=1 -DDO_PS=1 -DDO_FIG=1 -DDO_X11=1 -DDO_GD=1 -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIMITS_H=1 
INSTALL      = /usr/bin/install -c
INSTALL_DATA = ${INSTALL} -m 644
FIND         = find
MAKEDEPEND   = makedepend
AR           = ar
ARFLAGS      = -cr
RANLIB       = ranlib
LD           = ld
LDFLAGS      =   -L/usr/X11R6/lib64 -L/usr/local/lib -lm -lX11 -lgd


BASE_DIR = ./src
BASE_SRC = $(BASE_DIR)/g2_device.c         $(BASE_DIR)/g2_ui_control.c \
           $(BASE_DIR)/g2_util.c           $(BASE_DIR)/g2_fif.c \
           $(BASE_DIR)/g2_virtual_device.c $(BASE_DIR)/g2_physical_device.c \
           $(BASE_DIR)/g2_graphic_pd.c     $(BASE_DIR)/g2_control_pd.c \
           $(BASE_DIR)/g2_ui_graphic.c     $(BASE_DIR)/g2_ui_virtual_device.c \
           $(BASE_DIR)/g2_ui_device.c      $(BASE_DIR)/g2_splines.c
BASE_INS = $(BASE_DIR)/g2.h


PS_DIR = ./src/PS
PS_SRC = $(PS_DIR)/g2_PS.c
PS_INS = $(PS_DIR)/g2_PS.h

FIG_DIR = ./src/FIG
FIG_SRC = $(FIG_DIR)/g2_FIG.c
FIG_INS = $(FIG_DIR)/g2_FIG.h

X11_DIR = ./src/X11
X11_SRC = $(X11_DIR)/g2_X11.c
X11_INS = $(X11_DIR)/g2_X11.h

#WIN32_DIR = ./src/WIN32
#WIN32_SRC = $(WIN32_DIR)/g2_win32.c $(WIN32_DIR)/g2_win32_thread.c 
#WIN32_INS = $(WIN32_DIR)/g2_win32.h

GD_DIR = ./src/GD
GD_SRC = $(GD_DIR)/g2_gd.c
GD_INS = $(GD_DIR)/g2_gd.h


SRC = $(BASE_SRC) $(PS_SRC) $(FIG_SRC) $(X11_SRC) $(WIN32_SRC) $(GD_SRC)
OBJ = $(SRC:.c=.o)

INS =  $(BASE_INS) $(PS_INS) $(FIG_INS) $(X11_INS) $(WIN32_INS) $(GD_INS)

#
# define some phony targets
#
.PHONY: all shared clean doc release demo

#
# major rule
#
.c.o:  
	$(CC) $(CFLAGS) -c $< -o $@


all: libg2.a
	test -d ./include || mkdir ./include
	cp $(INS) ./include/

shared: libg2.so.0.$(G2_VERSION)
	@echo


libg2.a: $(OBJ)
	$(AR) $(ARFLAGS) libg2.a $(OBJ)
	test ! -n "$(RANLIB)" || $(RANLIB) $@

libg2.so.0.$(G2_VERSION): $(OBJ)
	ld -shared -soname libg2.so.0 -o $@ $(OBJ) 

install: libg2.a
	test -d $(LIBDIR) || mkdir -p $(LIBDIR)
	test -d $(INCDIR) || mkdir -p $(INCDIR)
	$(INSTALL_DATA) libg2.a $(LIBDIR)
	for IFILE in $(INS); do \
		$(INSTALL_DATA) $$IFILE  $(INCDIR); \
	done

release: clean doc
	cp ./doc/latex/refman.ps  ./doc/g2_refman.ps
	cp ./doc/latex/refman.pdf ./doc/g2_refman.pdf
	rm -r ./doc/latex

doc:
	(cd ./doc ; doxygen Doxyfile)
	(cd ./doc/latex ; make ps ; make pdf)


clean:
	-(cd ./demo ; make clean)
	-rm -f $(OBJ)
	-rm -f libg2.a config.cache config.log Makefile.bak config.status
	-rm -f ./include/*.h
	-$(FIND) . -name "*~" -exec rm -f {} \;
	-(cd ./g2_perl ; make clean)
	-rm -f ./g2_perl/test.ps
	-rm -f ./g2_perl/test.png
	-rm -f ./g2_perl/test.jpg
	-rm -f ./g2_perl/test.fig
	-rm -f ./g2_perl/Makefile.old
	-rm -f ./libg2.so.0.$(G2_VERSION)
	-rm -f libg2.$(G2_VERSION).a
	-rm -f a.out
	-rm -fr doc/html/ doc/latex/ doc/g2_refman.ps doc/g2_refman.pdf

demo:	libg2.a
	(cd ./demo ; make)

depend:
	$(MAKEDEPEND) -- $(CFLAGS) -- $(SRC)
	@echo
	@echo "Run 'make' to build g2 library."
	@echo

distdir:
	cp -r * $(distdir)


# DO NOT DELETE THIS LINE -- make  depend  depends  on it.

