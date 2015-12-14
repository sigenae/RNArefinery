#define _ISOC99_SOURCE
#include <features.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

int main(int argc, char *argv[]) {

  FILE *inf;
  char tmpf[256], tmp[256], *infn, *s;
  register int j;

  if (argc != 3) {
    printf(
        "Extracts a single fasta sequence from multifasta file\n usage: exfasta <seq name> <in fasta>\n");
    return 10;
  } else {
    s = argv[1];
    infn = argv[2];
  }

  inf = fopen(infn, "r");

  while (fgets(tmpf, 256, inf) != NULL) {

    if (tmpf[0] == '>') {
      strcpy(tmp, tmpf);
      for (j = 1; tmpf[j] && !isspace(tmpf[j]); j++)
        ;
      tmpf[j] = '\0';
      if (strcmp(&tmpf[1], s) == 0) {
        printf("%s", tmp);
        while (fgets(tmpf, 256, inf) && (tmpf[0] != '>')) {
          printf("%s", tmpf);
        }
        break;
      }
    }

  } //while

  fclose(inf);

}
