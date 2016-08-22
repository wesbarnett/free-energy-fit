all:
PKGNAME = fefit
PREFIX ?= /usr/local

.PHONY: all install clean

all:
	mkdir -p bin
	g++ -O2 -o bin/fefit -I./include src/* -llapack -lblas -Wall

install:
	install -Dm755 bin/* -t $(DESTDIR)$(PREFIX)/bin/
	install -Dm755 include/* -t $(DESTDIR)$(PREFIX)/include/

clean:
	rm -r bin
