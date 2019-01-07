# $Id: Makefile 2501 2007-11-20 02:33:29Z benkirk $


# The location of the mesh library
#meshdir := . 

# include the library options determined by configure.  This will
# set the variables INCLUDE and LIBS that we will need to build and
# link with the library.
include /scratch/libs/LIBMESHSUBDIV/Make.common
 
#
###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#
# source files
#srcfiles 	:= $(wildcard *.C)
srcfiles       := biharmonic.C 

#
# object files
bihobj		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))

###############################################################################



.PHONY: clean clobber distclean

###############################################################################
# Target:
#
biharmonic	   := ./biharmonic-$(METHOD)

bh:: $(biharmonic)

$(bihobj): biharmonic.C
	@echo "Compiling bharmonic.C..."
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $(libmesh_INCLUDE) -c biharmonic.C  -o $(bihobj)


# Production rules:  how to make the target - depends on library configuration
$(biharmonic): $(bihobj)
	@echo "Linking "$@"..."
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $(bihobj) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS) 


# Useful rules.
clean:
	@rm -f $(objects) *~

clobber:
	@$(MAKE) clean
	@rm -f $(target) out.gmv

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.g.o *.pg.o

# Warning, the 3D problem may be extremely slow if you are running in debug mode.
run: $(target)
	@echo "***************************************************************"
	@echo "* Running Example " $(LIBMESHRUN) $(target) $(LIBMESHOPTIONS)
	@echo "***************************************************************"
	@echo " "
	@$(LIBMESHRUN) $(target) -d 1 -n 20 $(LIBMESHOPTIONS)
	@$(LIBMESHRUN) $(target) -d 2 -n 15 $(LIBMESHOPTIONS)
	@$(LIBMESHRUN) $(target) -d 3 -n 6 $(LIBMESHOPTIONS)
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Example " $(LIBMESHRUN) $(target) $(LIBMESHOPTIONS)
	@echo "***************************************************************"


# include the dependency list
include .depend


#
# Dependencies
#
.depend:
	@$(perl) /scratch/libs/LIBMESHSUBDIV/contrib/bin/make_dependencies.pl -I. $(foreach i, $(wildcard $(meshdir)/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend
	@echo "Updated .depend"

###############################################################################
