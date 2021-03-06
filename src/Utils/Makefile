# H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# H0 X
# H0 X   libAtoms+QUIP: atomistic simulation library
# H0 X
# H0 X   Portions of this code were written by
# H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
# H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
# H0 X
# H0 X   Copyright 2006-2010.
# H0 X
# H0 X   These portions of the source code are released under the GNU General
# H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# H0 X
# H0 X   If you would like to license the source code under different terms,
# H0 X   please contact Gabor Csanyi, gabor@csanyi.net
# H0 X
# H0 X   Portions of this code were written by Noam Bernstein as part of
# H0 X   his employment for the U.S. Government, and are not subject
# H0 X   to copyright in the USA.
# H0 X
# H0 X
# H0 X   When using this software, please cite the following reference:
# H0 X
# H0 X   http://www.libatoms.org
# H0 X
# H0 X  Additional contributions by
# H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
# H0 X
# H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

all: libquiputils.a


ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif
include Makefile.inc
include Makefile.rules


QUIP_UTILS_F95_FILES =  ts_params transition_state elasticity phonons restraints_constraints_xml crackparams cracktools force_machine_learning_module real_space_covariance structure_analysis_traj_routines
QUIP_UTILS_F95_SOURCES = ${addsuffix .f95, ${QUIP_UTILS_F95_FILES}}
QUIP_UTILS_F95_OBJS = ${addsuffix .o, ${QUIP_UTILS_F95_FILES}}

QUIP_UTILS_CXX_FILES =
QUIP_UTILS_CXX_SOURCES = ${addsuffix .cpp ${QUIP_UTILS_CXX_FILES}}
QUIP_UTILS_CXX_OBJS = ${addsuffix .o, ${QUIP_UTILS_CXX_FILES}}

QUIP_UTILS_OBJS = ${QUIP_UTILS_F95_OBJS} ${QUIP_UTILS_CXX_OBJS}

ifeq (${HAVE_CGAL},1)
  QUIP_UTILS_CXX_FILES += alphashape
endif

.PHONY : clean allclean 


libquiputils.a:  ${QUIP_UTILS_OBJS}
ifneq (${LIBTOOL},)
	${LIBTOOL} -o libquiputils.a ${QUIP_UTILS_OBJS}
else
	${AR} ${AR_ADD} libquiputils.a $?
endif

clean:
	rm -f *.o *.mod *.mod.save libquiputils.a Utils.depend

allclean: clean

ifeq (${USE_MAKEDEP},1)
QUIP_UTILS_F95_FPP_FILES = ${addsuffix .fpp, ${QUIP_UTILS_F95_FILES}}

Utils.depend: ${QUIP_UTILS_F95_FPP_FILES}
	${SCRIPT_PATH}/${MAKEDEP} ${MAKEDEP_ARGS} -- ${addprefix ../../src/Utils/,${QUIP_UTILS_F95_SOURCES}} > Utils.depend

-include Utils.depend
endif

install: 
