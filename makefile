# Directory specific to the chosen BCH code
BCHDIR   := inc/codes/bch_796_256_60

override BCHDIR := $(patsubst %/,%,$(BCHDIR))
BCH      := $(subst bch_,,$(notdir $(BCHDIR)))
BCHPARAM := $(BCHDIR)/bch_$(BCH).h
BCHPOLY  := $(BCHDIR)/bch_$(BCH)_poly.h

INCDIR   := inc
LUTDIR   := $(INCDIR)/lut
SRCDIR   := src
OBJDIR   := obj
BINDIR   := bin
DOCDIR   := doc

CC       := gcc
CPPFLAGS := -include $(BCHPARAM) -iquote $(INCDIR) -iquote $(LUTDIR)
LTOFLAGS := -O3 -mavx -mavx2 -flto
CFLAGS   := -std=c99 -pedantic -Wall -Wno-parentheses $(LTOFLAGS)
LDFLAGS  := $(LTOFLAGS)

vpath %.c $(SRCDIR)
vpath %.h $(INCDIR) $(LUTDIR)


# Objects list

OBJECTS := gf fft bch test time
OBJECTS := $(addsuffix _lutmul, $(OBJECTS)) $(addsuffix _pclmul, $(OBJECTS))
OBJECTS := $(addsuffix _$(BCH).o, $(OBJECTS))
OBJECTS += time_lutmul_$(BCH)_once.o time_pclmul_$(BCH)_once.o
OBJECTS += generate.o
OBJECTS := $(addprefix $(OBJDIR)/, $(OBJECTS))


# Binary targets

generate: $(addprefix $(OBJDIR)/, generate.o bch_lutmul_$(BCH).o gf_lutmul_$(BCH).o) | $(BINDIR)
	$(CC) $(LDFLAGS) $^ -o $(BINDIR)/$@

.SECONDEXPANSION:

lutmul pclmul: test_$$@_$(BCH) time_$$@_$(BCH) time_$$@_$(BCH)_once

%_pclmul: override LDFLAGS += -mpclmul

test_%: $(addprefix $(OBJDIR)/, test_%.o bch_%.o fft_%.o gf_%.o) | $(BINDIR)
	$(CC) $(LDFLAGS) $^ -o $(BINDIR)/$@

time_%: $(addprefix $(OBJDIR)/, time_%.o bch_%.o fft_%.o gf_%.o) | $(BINDIR)
	$(CC) $(LDFLAGS) $^ -o $(BINDIR)/$@

time_%_once: $(addprefix $(OBJDIR)/, time_%_once.o bch_%.o fft_%.o gf_%.o) | $(BINDIR)
	$(CC) $(LDFLAGS) $^ -o $(BINDIR)/$@


# Object targets

$(OBJECTS): optimizations.h $(BCHPARAM) | $(OBJDIR)

$(OBJDIR)/%_pclmul_$(BCH).o: override CFLAGS += -mpclmul

$(OBJDIR)/generate.o: generate.c bch.h gf.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJDIR)/test_%mul_$(BCH).o: test.c bch.h fft.h gf.h $(BCHPOLY)
	$(CC) -include $(BCHPOLY) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJDIR)/time_%mul_$(BCH).o: time.c bch.h $(BCHPOLY)
	$(CC) -D_GNU_SOURCE -include $(BCHPOLY) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJDIR)/time_%mul_$(BCH)_once.o: time_once.c bch.c bch.h fft.h gf.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJDIR)/bch_%mul_$(BCH).o: bch.c bch.h fft.h gf.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJDIR)/fft_%mul_$(BCH).o: fft.c fft.h fft_lut_1024_64.h gf.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJDIR)/gf_%mul_$(BCH).o: gf_%mul.c gf.h gf_lut_1024.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@


# Other targets

$(BINDIR) $(OBJDIR):
	@mkdir -p $@

doxygen:
	@cd $(DOCDIR) && doxygen doxygen.conf

clean:
	rm -rf $(BINDIR) $(OBJDIR) $(DOCDIR)/html

.PHONY: clean doxygen lutmul pclmul
