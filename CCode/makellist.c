#define _ISOC99_SOURCE 
#include <features.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
//
#include "utils.h"

//
/// This reads the FASTA file and outpust a list of sequences with their length.
//
//
//
int main(int argc, char *argv[]) {

  FILE *inf;
  char tmp[1024], *infn, name[1024];
  int len, counter;
  register int j;

  if (argc != 2) {
    printf("%s\n  %s\n\n\n", "Prints multi-FATSA file sequences lengths.",
        "usage: makellist <in fasta>");
    return 10;
  } else {
    infn = argv[1];
  }

  inf = fopen(infn, "r");

  bool read = true;
  bool done = false;
  counter = 0;

  while (!done) {

    //get line from fasta file
    if (read) {
      if (fgets(tmp, 1024, inf) == NULL) {
        done = true;
      }
    }
    //assuming that we have something in buffer start working

    if (!done) {

      chomp(tmp);

      if (tmp[0] == '>') {

        for (j = 1; tmp[j] && !isspace(tmp[j]); j++)
          ;
        tmp[j] = '\0';
        strcpy(name, tmp + 1);

        int k = 0;
        while ((fgets(tmp, 1024, inf) != NULL)) {
          if (tmp[0] == '>') {
            break;
          }
          chomp(tmp);
          int tmpl = strlen(tmp);
          k = k + tmpl;
        }
        len = k;

        printf("%s %d\n", name, len);

      } //if == '>'

    } // if !done

    if (!done) {
      if (tmp[0] == '>') {
        read = false;
      } else {
        read = true;
      }
    }

  } //main while( !done )

  fclose(inf);
}
