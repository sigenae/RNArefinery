CC = gcc
CLIB = -lm
CINCLUDES = 
LDFLAGS=

OTHERS=

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
	rm exfasta
	rm include_mf

install:
	

# DO NOT DELETE
