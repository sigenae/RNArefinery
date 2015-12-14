#define _ISOC99_SOURCE
#include <features.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
# include "utils.h"

char* const skipped_suffix = ".skipped";

// ------------------------------------------------------------------------------------------------
// back to business...
//
int main(int argc, char *argv[]) {

  // input file,include list, output file and the 'leftovers' file - not included sequences
  FILE *in_file, *include_file, *out_file, *skipped_file;

  // reading buffer, file names
  char buffer[256], *in_file_name, *include_list_name, *out_file_name,
      *skipped_list_name;

  // number of sequences to include
  int include_count;

  // pointer to the include list and couple generic pointers
  SEQHEADER *include_list, *current_seq, include_seq;

  // the counter
  int register i;

  // pointer onto the full sequence buffer
  FASTASEQ *fasta_sequence;

  // check for proper args, i know it may be fancier but whatever :)
  if (argc != 4) {
    //print help
    printf("%s\n  %s\n  %s\n\n\n", "Filters the FASTA file according to the list",
        "usage: include_mf <in fasta> <out fasta> <include list>",
        "produces also <include list>.skipped just for your convenience");
    return 10;
  } else {
    //got the arguments, no validation...
    in_file_name = argv[1];
    out_file_name = argv[2];
    include_list_name = argv[3];
    // makeup the skipped sequences filename
    skipped_list_name = calloc(
        strlen(include_list_name) + strlen(skipped_suffix) + 1, 1);
    strncpy(skipped_list_name, include_list_name, strlen(include_list_name));
    strncat(skipped_list_name, skipped_suffix, strlen(skipped_suffix));
  }

  // malloc memory for sequence
  //
  printf("phase 0 - getting buffers & file handlers ... \n");
  fasta_sequence = (FASTASEQ*) malloc(sizeof(FASTASEQ));

  // get the file handlers
  //
  in_file = fopen(in_file_name, "r");
  printf("   . input fasta: %s \n", in_file_name);

  out_file = fopen(out_file_name, "w");
  printf("   . output fasta: %s \n", out_file_name);

  include_file = fopen(include_list_name, "r");
  printf("   . include list: %s \n", include_list_name);

  skipped_file = fopen(skipped_list_name, "w");
  printf("   . skipped sequences list: %s \n", skipped_list_name);

  // back to the business... start counting
  //
  printf("phase 1 - counting sequences in the include list ... \n");
  i = 0;
  while (fgets(buffer, sizeof(buffer), include_file) != NULL) {
    i++;
    if ((i % 10000) == 0) {
      fprintf(stderr, "%d\r", i);
    }
  }
  include_count = i;

  printf("   ... %d sequences is about to be included, reading the list...\n",
      include_count);
  include_list = malloc(sizeof(SEQHEADER) * include_count);
  rewind(include_file);
  i = 0;
  while (fgets(buffer, sizeof(buffer), include_file) != NULL) {
    chomp(buffer);
    strcpy(include_list[i].name, buffer);
    include_list[i].extracted = false;
    i++;
    if ((i % 10000) == 0) {
      fprintf(stderr, "%d\r", i);
    }
  }

  printf("   ... %d sequences loaded, sorting ...\n", include_count);
  qsort((void*) include_list, include_count, sizeof(SEQHEADER), cmp_seq_header);

  // start filtering
  //
  printf("phase 2 - filtering input ... \n");

  bool should_read = true;
  bool done = false;
  i = 0;

  while (!done) {

    //get line from fasta file if read flag is true
    if (should_read) {
      if (fgets(buffer, sizeof(buffer), in_file) == NULL) {
        done = true;
      }
    }

    //assuming that we have something in buffer start working
    //
    if (!done) {

      // clean the buffer ending
      chomp(buffer);

      // check if the buffer contains the sequence name
      if (buffer[0] == '>') {
        // save the full name with all the numbers and other junk, just in case of LUCY
        strcpy(fasta_sequence->real_name, buffer + 1); // (tmp1 + 1) to skip the '>'
        char *sn = strtok(buffer, " "); // split by spaces
        strcpy(fasta_sequence->name, sn + 1);

      } else {
        // okay if '>' not found in the very first line of the block -
        // we hit something wrong, exit from here, i.e we are done
        done = true;
      }

      // read the sequence block
      int k = 0;
      while ((fgets(buffer, sizeof(buffer), in_file) != NULL)) {
        if (buffer[0] == '>') {
          break;
        }
        chomp(buffer);
        int str_length = strlen(buffer);
        strncpy((char *) fasta_sequence->fasta + k, buffer, str_length);
        k = k + str_length;
      }
      fasta_sequence->fasta[k] = '\0';
      fasta_sequence->len = k;

      //printf(" -DEBUG- %s readed from fasta file \n", fs->name);
      strcpy(include_seq.name, fasta_sequence->name);

      //fprintf(stderr," - checking - %s\n", fs.name );

      // at this point we got the sequence read, let's see if we ned to include it
      //
      current_seq = bsearch((const void*) &include_seq,
          (const void*) include_list, include_count, sizeof(SEQHEADER),
          cmp_seq_header);

      // if found - print it into the output
      //
      if (current_seq != NULL) {
        if (current_seq->extracted != true) {
          current_seq->extracted = true;
          print_fasta(fasta_sequence, out_file);
          i++;
          if ((i % 10000) == 0) {
            fprintf(stderr, "%d\r", i);
          }
          //					printf(" -DEBUG- included: %s\n", fs->name);
        } else {
          //					printf(" -DEBUG- oops! dupe found: %s\n", fs->name);
        }
      } else {
        fprintf(skipped_file, "%s\n", fasta_sequence->real_name);
        //    printf(" --DEBUG-- excluded %s\n", fs.name);
      }

    }			// if( !done )

    //fprintf(stderr, "sequences processed so far:%d\r", counter);

    if (!done) {
      //printf("%s\n", tmpf);
      if (buffer[0] == '>') {
        should_read = false;
      } else {
        should_read = true;
        done = true;
      }
    }

    //done = true;

  }			//main while( !done )

  printf("finished - included %d sequences\n", i);

  fclose(in_file);
  fclose(out_file);
  fclose(include_file);
  fclose(skipped_file);

  return 0;

}
