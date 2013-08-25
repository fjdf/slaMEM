CC        = gcc
EXEC      = slaMEM
CFLAGS    = -Wall -Wextra -Wunused -mpopcnt
CDEBUG    = -g -ggdb
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
	@rm -f $(EXEC) $(EXEC).tar.gz

pack:
	@echo :: Packing files ...
	tar -cvzhf $(EXEC).tar.gz $(CSRCS) $(CHDRS) $(TXTS) $(SCRIPTS)
	@echo :: Done
