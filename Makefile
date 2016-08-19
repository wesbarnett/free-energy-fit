all:
PKGNAME = fefit
PREFIX ?= /usr/local

.PHONY: all install clean

all:
	mkdir -p bin
	g++ -O2 -o bin/fefit -I./include src/* -llapack -lblas

install:
	install -Dm755 bin/* -t $(DESTDIR)$(PREFIX)/bin/

clean:
	rm -r bin
