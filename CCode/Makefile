CC = gcc
CLIB = -lm
CINCLUDES = 
LDFLAGS=

OTHERS=

OBJS = makellist.o utils.o
SRCS = makellist.c utils.c
HDRS = utils.h 


all: makellist exfasta include_mf

include_mf: utils.o include_mf.o
	${CC} ${CFLAGS} $? ${CLIB} ${LDFLAGS} -o $@

makellist: utils.o makellist.o
	${CC} ${CFLAGS} $? ${CLIB} ${LDFLAGS} -o $@

include_mf.o: include_mf.c
	${CC} ${CFLAGS} ${CINCLUDES} -c -O2 $?

makellist.o: makellist.c
	${CC} ${CFLAGS} ${CINCLUDES} -c -O2 $?

exfasta: exfasta.o
	${CC} ${CFLAGS} $? ${CLIB} ${LDFLAGS} -o $@
	
utils.o: utils.c
	${CC} ${CFLAGS} ${CINCLUDES} -c -O2 $?	

clean:
	rm *.o
	rm makellist
	rm target/makellist

install:
	cp makellist target/

# DO NOT DELETE
