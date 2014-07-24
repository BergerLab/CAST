prefix=/usr/local

all: cablast-compress cablast-decompress cablast-search cablat
install:
	make
	install -m 0755 cablast-compress cablast-decompress cablast-search cablat $(prefix)/bin

cablast-compress: src/cablast-compress
	cp src/cablast-compress .

src/cablast-compress:
	(cd src && make compress)

cablast-decompress: src/cablast-decompress
	cp src/cablast-decompress .

src/cablast-decompress:
	(cd src && make decompress)

cablast-search: src/cablast-search
	cp src/cablast-search .

src/cablast-search:
	(cd src && make search)

cablat: src/cablat
	cp src/cablat .

src/cablat:
	(cd src && make cablat)


clean:
	(cd src && make clean)
	rm -f cablast-*

push:
	git push origin master
	git push github master

