/* 
   Copyright (c) 2013 by Gordon D Brown <gordon.brown@cruk.cam.ac.uk>

   Modified by Margus Lukk <margus.lukk@cruk.cam.ac.uk>

   This software is distributed under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
   The licence text is available in file COPYING.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "redblackSN.h"

#define ISRED(x) ((x)->isRed==1)
#define SETRED(x) ((x)->isRed = 1)
#define SETBLACK(x) ((x)->isRed = 0)

rbsnTree rbsn_new() {
  rbsnTree x = (rbsnTree) malloc(sizeof(struct rbsnTree));
  x->root = NULL;
  x->count = 0;
  return x;
}

// static rbsnNode node_new(char *k,int data) {
static rbsnNode node_new(char *k,char *data) {
  rbsnNode x = (rbsnNode) malloc(sizeof(struct rbsnNode));
  x->key = k;
  x->isRed = 1;
  x->left = NULL;
  x->right = NULL;
  x->parent = NULL;
  x->data = data;
  return x;
}

static void insert(rbsnTree rbt,rbsnNode rbn) {
  rbsnNode y = NULL;
  rbsnNode x = rbt->root;
  while (x != NULL) {
    y = x;
    if (strcmp(rbn->key,x->key) < 0) {
      x = x->left;
    } else {
      x = x->right;
    }
  }
  rbn->parent = y;
  if (y == NULL) {
    rbt->root = rbn;
  } else if (strcmp(rbn->key,y->key) < 0) {
    y->left = rbn;
  } else {
    y->right = rbn;
  }
}

static void leftRotate(rbsnTree rbt,rbsnNode rbn) {
  rbsnNode y = rbn->right;
  rbn->right = y->left;
  if (y->left != NULL) {
    y->left->parent = rbn;
  }
  y->parent = rbn->parent;
  if (rbn->parent == NULL) {
    rbt->root = y;
  } else if (rbn == rbn->parent->left) {
    rbn->parent->left = y;
  } else {
    rbn->parent->right = y;
  }
  y->left = rbn;
  rbn->parent = y;
}

static void rightRotate(rbsnTree rbt,rbsnNode rbn) {
  rbsnNode y = rbn->left;
  rbn->left = y->right;
  if (y->right != NULL) {
    y->right->parent = rbn;
  }
  y->parent = rbn->parent;
  if (rbn->parent == NULL) {
    rbt->root = y;
  } else if (rbn == rbn->parent->right) {
    rbn->parent->right = y;
  } else {
    rbn->parent->left = y;
  }
  y->right = rbn;
  rbn->parent = y;
}

// void rbsn_insert(rbsnTree rbt,char *key,int data) {
void rbsn_insert(rbsnTree rbt,char *key,char *data) {
  rbsnNode x = node_new(key,data);
  insert(rbt,x);
  while (x!=rbt->root && ISRED(x->parent)) {
    if (x->parent->parent->left == x->parent) {
      rbsnNode y = x->parent->parent->right;
      if (y != NULL && ISRED(y)) {
        SETBLACK(x->parent);
        SETBLACK(y);
        SETRED(x->parent->parent);
        x = x->parent->parent;
      } else {
        if (x == x->parent->right) {
          x = x->parent;
          leftRotate(rbt,x);
        }
        SETBLACK(x->parent);
        SETRED(x->parent->parent);
        rightRotate(rbt,x->parent->parent);
      }
    } else {
      rbsnNode y = x->parent->parent->left;
      if (y != NULL && ISRED(y)) {
        SETBLACK(x->parent);
        SETBLACK(y);
        SETRED(x->parent->parent);
        x = x->parent->parent;
      } else {
        if (x == x->parent->left) {
          x = x->parent;
          rightRotate(rbt,x);
        }
        SETBLACK(x->parent);
        SETRED(x->parent->parent);
        leftRotate(rbt,x->parent->parent);
      }
    }
  }
  SETBLACK(rbt->root);
  rbt->count++;
}

static rbsnNode search(rbsnNode rbn,char *key) {
  int cmp;
  if (rbn == NULL) {
    return rbn;
  }
  cmp = strcmp(rbn->key,key);
  if (cmp == 0) {
    return rbn;
  }
  if (cmp > 0) {
    return search(rbn->left,key);
  } else {
    return search(rbn->right,key);
  }
}

//int rbsn_search(rbsnTree rbt,char *key) {
char *rbsn_search(rbsnTree rbt,char *key) {
  rbsnNode rbn;
  // int result = 0;
  rbn = search(rbt->root,key);
  if (rbn != NULL) {
    // result = rbn->data;
    return rbn->data;
  }
  // return result;
  return NULL;
}

// void rbsn_update(rbsnTree rbt,char *key,int value) {
void rbsn_update(rbsnTree rbt,char *key,char *value) {
  rbsnNode x = search(rbt->root,key);
  if (x != NULL) {
    x->data = value;
  } else {
    // fprintf(stderr,"error: attempt to update non-existent node (%s,%d)\n",key,value);
    fprintf(stderr,"error: attempt to update non-existent node (%s,%s)\n",key,value);
  }
}

// static void traverse(rbsnNode rbn,void (*visit)(char *key,int value,void *data),void *data) {
static void traverse(rbsnNode rbn,void (*visit)(char *key,char *value,void *data),void *data) {
  if (rbn == NULL) {
    return;
  } else {
    traverse(rbn->left,visit,data);
    visit(rbn->key,rbn->data,data);
    traverse(rbn->right,visit,data);
  }
}

//static void traverseRL(rbsnNode rbn,void (*visit)(char *key,int value,void *data),void *data) {
static void traverseRL(rbsnNode rbn,void (*visit)(char *key,char *value,void *data),void *data) {
  if (rbn == NULL) {
    return;
  } else {
    traverseRL(rbn->right,visit,data);
    visit(rbn->key,rbn->data,data);
    traverseRL(rbn->left,visit,data);
  }
}

//void rbsn_traverse(rbsnTree r,void(*visit)(char *key,int value,void *data),void *data,int reverse) {
void rbsn_traverse(rbsnTree r,void(*visit)(char *key,char *value,void *data),void *data,int reverse) {
  if (!reverse) {
    traverse(r->root,visit,data);
  } else {
    traverseRL(r->root,visit,data);
  }
}

int rbsn_count(rbsnTree rbt) {
  return rbt->count;
}

int rbsn_isNull(rbsnNode rbn) {
  return rbn == NULL;
}

