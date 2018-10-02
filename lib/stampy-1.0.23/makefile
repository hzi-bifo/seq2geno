#
# - the linux x86-64 platform is supported, with Python 2.6 or 2.7.  Both 
#   2- and 4-byte 4-byte Unicode representations (UCS2 and UCS4) are supported.
#   
# - the Mac OS-X 10.6 with x86-64 platform is supported, with Python 2.6, and 2-byte 
#   Unicode representation.  No other combinations are currently supported.
#

SHELL=/bin/bash
python?=python

#############################
# User targets
#############################

all:
	@$(MAKE) ucs=`$(python) -c "import sys;print sys.maxunicode > 65536 and 'ucs4' or 'ucs2'"` platform=`$(python) -c "import sys;print {'linux2':'linux-x86_64','darwin':'macosx-10.6-universal'}.get(sys.platform,'unknown')"` build-stampy

clean:
	-rm maptools.so
	-rm build/*/*.pyc

#############################
# End of user targets
#############################

build-stampy:
	@if [[ `$(python) --version 2>&1 |awk -F "." '{print $$2}'` -eq 6 ]]; then \
	  $(MAKE) bdir=build/$(platform)-2.6-$(ucs) python=python2.6 platform=$(platform) buildall ; \
	elif [[ `$(python) --version 2>&1 |awk -F "." '{print $$2}'` -eq 7 ]]; then \
	  $(MAKE) bdir=build/$(platform)-2.7-$(ucs) python=python2.7 platform=$(platform) buildall ; \
	else \
	  echo Python version 2.6 or 2.7 is required; \
	fi

objs=$(bdir)/pyx/maptools.o $(bdir)/c/maputils.o $(bdir)/c/alignutils.o $(bdir)/readalign.o $(bdir)/algebras.o $(bdir)/frontend.o

$(bdir)/%.o:
	gcc -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -fPIC `$(python)-config --includes` -c ext/readalign/$*.cc -o $@

buildall: $(bdir)/readalign.o $(bdir)/frontend.o $(bdir)/algebras.o
	cp -r build/$(python)/* .
ifeq ($(platform),linux-x86_64)
	g++ `$(python)-config --ldflags` -pthread -shared -Wl $(objs) -o maptools.so
else
	g++ `$(python)-config --ldflags` -pthread -dynamiclib -Wl $(objs) -o maptools.so
endif

