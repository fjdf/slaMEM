CC        = gcc
EXEC      = slaMEM
CFLAGS    = -Wall -Wextra -Wunused -mpopcnt
CDEBUG    = -g -ggdb -fno-inline -dH -DGDB
COPTIMIZE = -Wuninitialized -O9 -fomit-frame-pointer
CLIBS     = -lm

CSRCS     = $(wildcard *.c)
CHDRS     = $(wildcard *.h)
TXTS      = $(wildcard *.txt README* LICENSE*)
SCRIPTS   = $(wildcard Makefile* *.sh *.py)

NAME	:= "slaMEM"
VERSION	:= $(shell sed -n 's/.*VERSION \"\(.*\)\".*/\1/p' slamem.c)
CPUARCH	:= $(shell uname -m)

ifeq ($(MAKECMDGOALS),debug)
	CFLAGS += $(CDEBUG)
	EXEC   := $(addsuffix -debug, $(EXEC))
else ifeq ($(MAKECMDGOALS),pack)
	EXEC   := $(addsuffix -v$(VERSION), $(EXEC))
else
	CFLAGS += $(COPTIMIZE)
endif

.PHONY: all clean pack

all: clean bin

debug: all

bin:
	@echo :: Compiling \"$(NAME) v$(VERSION)\" \($(CPUARCH)\) ...
	$(CC) $(CFLAGS) $(CSRCS) -o $(EXEC) $(CLIBS)
	@echo :: Done

clean:
	@echo :: Cleaning up ...
	@rm -f $(EXEC) $(EXEC)-debug $(EXEC)-v$(VERSION).tar.gz

pack:
	@echo :: Packing files ...
	tar -cvzhf $(EXEC).tar.gz $(CSRCS) $(CHDRS) $(TXTS) $(SCRIPTS)
	@echo :: Done
