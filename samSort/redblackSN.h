/* 
  Copyright (c) 2013 by Gordon D Brown <gordon.brown@cruk.cam.ac.uk>

  Modified by Margus Lukk <margus.lukk@cruk.cam.ac.uk>

  This software is distributed under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007                           
  The licence text is available in file COPYING.
*/

#ifndef REDBLACKSN_H
#define REDBLACKSN_H

typedef struct rbsnNode {
  char *key; /* node key value */
  char isRed; /* 1 for red, 0 for black */
  struct rbsnNode *left;
  struct rbsnNode *right;
  struct rbsnNode *parent;
  // int data;
  char *data;
} *rbsnNode;

typedef struct rbsnTree {
  int count; /* nodes currently in tree */
  rbsnNode root;
} *rbsnTree;

rbsnTree rbsn_new();
void rbsn_insert(rbsnTree rbt,char *key,char *data);
int rbsn_count(rbsnTree rbt);
// int rbsn_search(rbsnTree rbt,char *key);
// void rbsn_update(rbsnTree rbt,char *key,int value);
// void rbsn_traverse(rbsnTree rbt,void(*visit)(char *key,int value,void *data),void *data,int reverse);
char *rbsn_search(rbsnTree rbt,char *key);
void rbsn_update(rbsnTree rbt,char *key,char *value);
void rbsn_traverse(rbsnTree rbt,void(*visit)(char *key,char *value,void *data),void *data,int reverse);

int rbsn_isNull(rbsnNode rbn);

#endif
