all: lib include lib/libds.a lib/libopt.a include/ds.h include/opt.h

lib:
	mkdir -p lib

include:
	mkdir -p include

lib/libds.a: libds/libds.a
	cp libds/libds.a lib/libds.a

libds/libds.a:
	(cd libds && make libds.a)

lib/libopt.a: libopt/libopt.a
	cp libopt/libopt.a lib/libopt.a

libopt/libopt.a:
	(cd libopt && make libopt.a)

include/ds.h: libds/ds.h
	cp libds/ds.h include/

libds/ds.h:
	(cd libds && make ds.h)

include/opt.h: libopt/opt.h
	cp libopt/opt.h include/

clean:
	(cd libds && make clean)
	(cd libopt && make clean)
	rm -rf lib/libds.a lib/libopt.a
	rm -rf include/ds.h include/opt.h

push:
	git push origin master
	git push github master

