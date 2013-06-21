/* 
   Copyright (c) 2013 by Margus Lukk <margus.lukk@cruk.cam.ac.uk>
   
   This software is distributed under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
   The licence text is available in the file COPYING.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "redblackSN.h"

#ifndef VERSION
#define VERSION "0.1beta"
#endif

void visit(char *key,char *value, void *data) {
  fprintf(stdout,"%s",value);
}

int samSort(char *fname) {

  rbsnTree tree;
  FILE *infd;
  int bsize = 1024;
  char buffer[bsize];
  char *chr;
  char *chrpos;
  int pos;
  int postimes;
  int i;
  // int data=0;
  char *keyp;
  char key[bsize];

  struct timeval t0, t1;
  gettimeofday(&t0, 0);

  tree = rbsn_new();
  
  infd = fopen(fname,"r");

  while (fgets(buffer,bsize,infd) != NULL) {
    // process only if not header
    if(buffer[0] != '@') {      
      int blen = strlen(buffer);
      pos = 0;
      postimes = 0;
      key[0]='\0';
      for(i=0; i<blen; i++) {
	if(buffer[i] == '\t') {
	  if(postimes == 2) {
	    int len = i-pos;
	    chr = (char *) malloc (len+1);
	    strncpy(chr,buffer+pos,len);
	    chr[len] = '\0';
	    strcpy(key,chr);
	    free(chr);
	  }
	  if(postimes == 3) {
	    int len = i-pos;
	    chrpos = (char *) malloc (len+1);
	    strncpy(chrpos,buffer+pos,len);
	    chrpos[len] = '\0';
	    strcat(key,chrpos);
	    free(chrpos);
	    break;
	  }
	  pos = i+1;
	  postimes++;
	}
      }
      keyp = strdup(key);
      char *data;
      data = (char *) malloc(strlen(buffer)+1);
      strcpy(data,buffer);
      rbsn_insert(tree,keyp,data);
    }
    // print header
    else {
      fprintf(stdout,"%s",buffer);
    }
    if (feof(infd)) {
      break;
    }
  }
  fclose(infd);

  int tsize = rbsn_count(tree);
 
  gettimeofday(&t1, 0);
  long elapsed = (t1.tv_sec-t0.tv_sec);
  fprintf(stderr,"[fast_sam_sort] Tree of size %d constructed in %lu sec.\n",tsize,elapsed);

  // traverse tree
  rbsn_traverse(tree,visit,NULL,0);

  gettimeofday(&t1, 0);
  elapsed = (t1.tv_sec-t0.tv_sec);
  fprintf(stderr,"[fast_sam_sort] Tree constructed and trasversal finished in %lu sec.\n",elapsed);

  return 0;
}

int main(int argc,char **argv) {

  if (argc != 2) {
    fprintf(stderr, "\nProgram: samSort (sorts Sam file in RB tree)\n");
    fprintf(stderr, "Version: %s\n", VERSION);
    fprintf(stderr, "Usage:   samSort <filename.sam>\n");
    fprintf(stderr, "Description: Sorted sam file is written to stdout.\n");
    fprintf(stderr, "             NB! The program is in the proof of concept stage, has not been optimised, consumes 1.5x size of the memory of the input file; and may not produce meaningful output.\n\n");
    return 1;
  }
  
  return samSort(argv[1]);

}
