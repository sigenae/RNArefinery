#ifndef UTILS_H_
#define UTILS_H_

#include <stdbool.h>

// memory structure used for extracting the sequence, MAX size of seq - 8388607
//
struct fasta_sequence {
  char name[256];
  char real_name[256];
  int start, stop, len;
  int fasta[8388608];
};
typedef struct fasta_sequence FASTASEQ;

// memory structure for keeping sequence names, max seq name length is 256 symbols
//
struct sequence_header {
  char name[257];
  bool extracted; // just keep track of duplicates - don't write them into output
};
typedef struct sequence_header SEQHEADER;

void chomp(char *str);

void print_fasta(FASTASEQ *f, FILE *out);

int cmp_seq_header(const void* _a, const void* _b);

#endif /* UTILS_H_ */
