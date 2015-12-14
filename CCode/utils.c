#define _ISOC99_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "utils.h"

void chomp(char *str) {
	int l = strlen(str) - 1;
	while ((str[l] == '\n') || (str[l] == ' ')) {
		str[l] = '\0';
		l--;
	}
}

// prints the FASTA sequence into the file
//
void print_fasta(FASTASEQ *f, FILE *out) {

  // line buffer
  char tmp[60];

  // print the name
  fprintf(out, ">%s\n", f->real_name);

  // while sequence 'reminder' greater than 50 symbols print chunks to the file
  register int c = 0;
  while (c + 50 < f->len) {
    strncpy(tmp, (char *) (f->fasta) + c, 50);
    tmp[50] = '\0';
    fprintf(out, "%s\n", tmp);
    c = c + 50;
  }

  // if there is a leftover, print it too
  if (c < f->len) {
    strncpy(tmp, (char *) (f->fasta) + c, f->len - c);
    tmp[f->len - c] = '\0';
    fprintf(out, "%s\n", tmp);
  }
}

int cmp_seq_header(const void* _a, const void* _b) {
  const char* a = (const char*) (((const SEQHEADER*) _a)->name);
  const char* b = (const char*) (((const SEQHEADER*) _b)->name);
  //printf(" --DEBUG-- comparing %s and %s\n", a, b);
  return strcmp(a, b);
}
