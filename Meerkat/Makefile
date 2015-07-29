 ########################################################################## 
 # Copyright:                                                             #
 # University of Warwick, UK                                              #
 # Budker Institute of Nuclear Physics, Novosibirsk, Russia               #
 # 2014                                                                   #
 #                                                                        #
 # Distributed under the Boost Software License, Version 1.0.             #
 # (See accompanying file LICENSE_1_0.txt or copy at                      #
 # http://www.boost.org/LICENSE_1_0.txt)                                  #
 #                                                                        #
 # Author:                                                                #
 #   Anton Poluektov                                                      #
 #                                                                        #
 # ---------------------------------------                                #
 # Standalone Makefile for Meerkat package                                #
 # ---------------------------------------                                #
 #                                                                        #
 # Instructions                                                           #
 #     - Review 'external configuration' section below                    #
 #       to match systems compilers setup                                 #
 #                                                                        #
 #     - Make sure the ROOTSYS environment variable is set and points     #
 #       to your ROOT release or the root-config script is in your PATH   #
 #                                                                        #
 #     - run 'make <target>'                                              #
 #                                                                        #
 # Build targets                                                          #
 #   lib   - make libMeerkat.a                                            #
 #   shlib - make libMeerkat.so (default)                                 #
 #   clean - delete all intermediate and final build objects              #
 #                                                                        #
 ########################################################################## 


# --- External configuration ----------------------------------

# first check that ROOTSYS is defined
ifndef ROOTSYS
  ROOTSYS := $(shell root-config --prefix)
  ROOTBINDIR := $(shell root-config --bindir)
  ifeq ($(ROOTSYS), )
    $(error running of root-config failed or reported null value)
  endif 
else
  ROOTBINDIR := $(ROOTSYS)/bin
endif

ROOTCONFIG := $(ROOTBINDIR)/root-config
ARCH       := $(shell $(ROOTCONFIG) --arch)
PLATFORM   := $(shell $(ROOTCONFIG) --platform)

INCLUDES = 
SRCDIR   = src
INCDIR   = inc
LIBDIR   = lib
WORKDIR  = tmp

ifeq ($(findstring linux, $(ARCH)),linux)
# This set here should work for Linux.
CXX      = g++
LD       = g++
CXXFLAGS = -g -O3 -Wall -Wextra -Wshadow -Woverloaded-virtual -Werror -fPIC -std=c++11
MFLAGS   = -MM
SOFLAGS  = -shared
endif

ifeq ($(ARCH),macosx64)
# For Mac OS X you may need to put -m64 in CXXFLAGS and SOFLAGS.
CXX      = g++
LD       = g++
CXXFLAGS = -g -O3 -Wall -Wextra -Wshadow -Woverloaded-virtual -Werror -fPIC -m64 -std=c++11
MFLAGS   = -MM
SOFLAGS  = -m64 -dynamiclib -single_module -undefined dynamic_lookup
endif

# --- Internal configuration ----------------------------------
PACKAGE=Meerkat
DEPDIR=$(WORKDIR)/dependencies
OBJDIR=$(WORKDIR)/objects

INCLUDES += -I$(INCDIR) -I$(shell $(ROOTBINDIR)/root-config --incdir)
CXXFLAGS += $(INCLUDES)
SKIPLIST  = test.cpp
CINTFILE  = $(WORKDIR)/$(PACKAGE)Cint.cc
CINTOBJ   = $(OBJDIR)/$(PACKAGE)Cint.o
LIBFILE   = $(LIBDIR)/lib$(PACKAGE).a
SHLIBFILE = $(LIBDIR)/lib$(PACKAGE).so
ROOTMAPFILE := $(patsubst %.so,%.rootmap,$(SHLIBFILE))

default: shlib

# List of all header files
HHLIST:=$(wildcard $(INCDIR)/*.hh)

# List of all source files to build
CCLIST:=$(filter-out $(SKIPLIST),$(wildcard $(SRCDIR)/*.cpp))

# List of all object files to build
OLIST:=$(patsubst %.cpp,%.o,$(addprefix $(OBJDIR)/,$(notdir $(CCLIST))))

# List of all dependency files to make
DLIST:=$(patsubst %.cpp,%.d,$(addprefix $(DEPDIR)/,$(notdir $(CCLIST))))

# Implicit rule making all dependency Makefiles included at the end of this makefile
$(DEPDIR)/%.d: $(SRCDIR)/%.cpp
	@echo "Making $@"
	@mkdir -p $(DEPDIR)
	@set -e; $(CXX) $(MFLAGS) $(CXXFLAGS) $< \
	          | sed 's#\($(notdir $*)\)\.o[ :]*#$(OBJDIR)/\1.o $@ : #g' > $@; \
	        [ -s $@ ] || rm -f $@

# Implicit rule to compile all classes
$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	@$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to make ROOTCINT output file
$(CINTOBJ): $(HHLIST) $(INCDIR)/$(PACKAGE)_LinkDef.h
	@mkdir -p $(OBJDIR)
	@echo "Running rootcint"
	@$(ROOTBINDIR)/rootcint -f $(CINTFILE) -c -p $(INCLUDES) $(notdir $(HHLIST)) $(INCDIR)/$(PACKAGE)_LinkDef.h
	@echo "Compiling $(CINTFILE)"
	@$(CXX) $(CXXFLAGS) -c $(CINTFILE) -o $(CINTOBJ)

# Rule to combine objects into a library
$(LIBFILE): $(OLIST)
	@echo "Making $(LIBFILE)"
	@mkdir -p $(LIBDIR)
	@rm -f $(LIBFILE)
	@ar rcs $(LIBFILE) $(OLIST)

# Rule to combine objects into a shared library
$(SHLIBFILE): $(OLIST) $(CINTOBJ)
	@echo "Making $(SHLIBFILE)"
	@mkdir -p $(LIBDIR)
	@rm -f $(SHLIBFILE)
	@$(CXX) $(OLIST) $(CINTOBJ) $(SOFLAGS) -o $(SHLIBFILE)

# Rule to create rootmap file
$(ROOTMAPFILE): $(SHLIBFILE)
	@echo "Making $(ROOTMAPFILE)"
	@mkdir -p $(LIBDIR)
	@rm -f $(ROOTMAPFILE)
	@rlibmap -f -o $(ROOTMAPFILE) -l $(SHLIBFILE) -d libCore.so libEG.so libHist.so libMathCore.so libMatrix.so libNet.so libRIO.so libTree.so -c $(INCDIR)/$(PACKAGE)_LinkDef.h

# Useful build targets
lib: $(LIBFILE) 
shlib: $(SHLIBFILE)  $(ROOTMAPFILE)
#$(ROOTMAPFILE)
clean:
	rm -rf $(WORKDIR)
	rm -f $(LIBFILE)
	rm -f $(SHLIBFILE)

.PHONY : shlib lib default clean

-include $(DLIST)
