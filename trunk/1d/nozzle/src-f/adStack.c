static char adSid[]="$Id: adStack.c,v 1.6 10/.0/.0 .0:.5:.0 llh Exp $";

#include <stdlib.h>

#define ONE_BLOCK_SIZE 16384
#define CHUNK_SIZE 4096

/* The main stack is a double-chain of DoubleChainedBlock objects.
 * Each DoubleChainedBlock holds an array[ONE_BLOCK_SIZE] of char. */
typedef struct _doubleChainedBlock{
  struct _doubleChainedBlock *prev ;
  char                       *contents ;
  struct _doubleChainedBlock *next ;
} DoubleChainedBlock ;

/* Globals that define the current position in the stack: */
static DoubleChainedBlock *curStack = NULL ;
static char               *curStackTop    = NULL ;
/* Globals that define the current LOOKing position in the stack: */
static DoubleChainedBlock *lookStack = NULL ;
static char               *lookStackTop    = NULL ;

/* PUSHes "nbChars" consecutive chars from a location starting at address "x".
 * Resets the LOOKing position if it was active.
 * Checks that there is enough space left to hold "nbChars" chars.
 * Otherwise, allocates the necessary space. */
void pushNarray(char *x, int nbChars) {
  int nbmax = (curStack)?ONE_BLOCK_SIZE-(curStackTop-(curStack->contents)):0 ;
  lookStack = NULL ;
  if (nbChars <= nbmax) {
    memcpy(curStackTop,x,nbChars);
    curStackTop+=nbChars ;
  } else {
    char *inx = x+(nbChars-nbmax) ;
    if (nbmax>0) memcpy(curStackTop,inx,nbmax) ;
    while (inx>x) {
      /* Create new block: */
      if ((curStack == NULL) || (curStack->next == NULL)) {
	DoubleChainedBlock *newStack ;
	char *contents = (char*)malloc(ONE_BLOCK_SIZE*sizeof(char)) ;
	newStack = (DoubleChainedBlock*)malloc(sizeof(DoubleChainedBlock)) ;
	if ((contents == NULL) || (newStack == NULL)) {
	  DoubleChainedBlock *stack = curStack ;
	  int nbBlocks = (stack?-1:0) ;
	  while(stack) {
	      stack = stack->prev ;
	      nbBlocks++ ;
	  }
	  printf("Out of memory (allocated %i blocks of %i bytes)\n",
		 nbBlocks, ONE_BLOCK_SIZE) ;
          exit(0);
	}
	if (curStack != NULL) curStack->next = newStack ;
	newStack->prev = curStack ;
	newStack->next = NULL ;
	newStack->contents = contents ;
	curStack = newStack ;
      } else
	curStack = curStack->next ;
      /* new block created! */
      inx -= ONE_BLOCK_SIZE ;
      if(inx>x)
	memcpy(curStack->contents,inx,ONE_BLOCK_SIZE) ;
      else {
	int nbhead = (inx-x)+ONE_BLOCK_SIZE ;
	curStackTop = curStack->contents ;
	memcpy(curStackTop,x,nbhead) ;
	curStackTop += nbhead ;
      }
    }
  }
}

/* POPs "nbChars" consecutive chars to a location starting at address "x".
 * Resets the LOOKing position if it was active.
 * Checks that there is enough data to fill "nbChars" chars.
 * Otherwise, pops as many blocks as necessary. */
void popNarray(char *x, int nbChars) {
  int nbmax = curStackTop-(curStack->contents) ;
  lookStack = NULL ;
  if (nbChars <= nbmax) {
    curStackTop-=nbChars ;
    memcpy(x,curStackTop,nbChars);
  } else {
    char *tlx = x+nbChars ;
    if (nbmax>0) memcpy(x,curStack->contents,nbmax) ;
    x+=nbmax ;
    while (x<tlx) {
      curStack = curStack->prev ;
      if (x+ONE_BLOCK_SIZE<tlx) {
	memcpy(x,curStack->contents,ONE_BLOCK_SIZE) ;
	x += ONE_BLOCK_SIZE ;
      } else {
	int nbtail = tlx-x ;
	curStackTop=(curStack->contents)+ONE_BLOCK_SIZE-nbtail ;
	memcpy(x,curStackTop,nbtail) ;
	x = tlx ;
      }
    }
  }
}

/* LOOKs "nbChars" consecutive chars to a location starting at address "x".
 * Activates the LOOKing position if it was reset.
 * LOOKing is just like POPping, except that the main pointer
 * remains in place, so that the value is not POPped.
 * Further PUSHs or POPs will start from the same place as if
 * no LOOK had been made. */
void lookNarray(char *x, int nbChars) {
  int nbmax ;
  if (lookStack == NULL) {
    lookStack = curStack ;
    lookStackTop = curStackTop ;
  }
  nbmax = lookStackTop-(lookStack->contents) ;
  if (nbChars <= nbmax) {
    lookStackTop-=nbChars ;
    memcpy(x,lookStackTop,nbChars);
  } else {
    char *tlx = x+nbChars ;
    if (nbmax>0) memcpy(x,lookStack->contents,nbmax) ;
    x+=nbmax ;
    while (x<tlx) {
      lookStack = lookStack->prev ;
      if (x+ONE_BLOCK_SIZE<tlx) {
	memcpy(x,lookStack->contents,ONE_BLOCK_SIZE) ;
	x += ONE_BLOCK_SIZE ;
      } else {
	int nbtail = tlx-x ;
	lookStackTop=(lookStack->contents)+ONE_BLOCK_SIZE-nbtail ;
	memcpy(x,lookStackTop,nbtail) ;
	x = tlx ;
      }
    }
  }
}

/****** Exported PUSH/POP/LOOK functions for ARRAYS: ******/

void pushcharacterarray_(char *x, int *n) {
  pushNarray(x,*n) ;
}
void popcharacterarray_(char *x, int *n) {
  popNarray(x,*n) ;
}
void lookcharacterarray_(char *x, int *n) {
  lookNarray(x,*n) ;
}

void pushbooleanarray_(char *x, int *n) {
  pushNarray(x,(*n*4)) ;
}
void popbooleanarray_(char *x, int *n) {
  popNarray(x,(*n*4)) ;
}
void lookbooleanarray_(char *x, int *n) {
  lookNarray(x,(*n*4)) ;
}

void pushinteger4array_(char *x, int *n) {
  pushNarray(x,(*n*4)) ;
}
void popinteger4array_(char *x, int *n) {
  popNarray(x,(*n*4)) ;
}
void lookinteger4array_(char *x, int *n) {
  lookNarray(x,(*n*4)) ;
}

void pushinteger8array_(char *x, int *n) {
  pushNarray(x,(*n*8)) ;
}
void popinteger8array_(char *x, int *n) {
  popNarray(x,(*n*8)) ;
}
void lookinteger8array_(char *x, int *n) {
  lookNarray(x,(*n*8)) ;
}

void pushinteger16array_(char *x, int *n) {
  pushNarray(x,(*n*16)) ;
}
void popinteger16array_(char *x, int *n) {
  popNarray(x,(*n*16)) ;
}
void lookinteger16array_(char *x, int *n) {
  lookNarray(x,(*n*16)) ;
}

void pushreal4array_(char *x, int *n) {
  pushNarray(x,(*n*4)) ;
}
void popreal4array_(char *x, int *n) {
  popNarray(x,(*n*4)) ;
}
void lookreal4array_(char *x, int *n) {
  lookNarray(x,(*n*4)) ;
}

void pushreal8array_(char *x, int *n) {
  pushNarray(x,(*n*8)) ;
}
void popreal8array_(char *x, int *n) {
  popNarray(x,(*n*8)) ;
}
void lookreal8array_(char *x, int *n) {
  lookNarray(x,(*n*8)) ;
}

void pushreal16array_(char *x, int *n) {
  pushNarray(x,(*n*16)) ;
}
void popreal16array_(char *x, int *n) {
  popNarray(x,(*n*16)) ;
}
void lookreal16array_(char *x, int *n) {
  lookNarray(x,(*n*16)) ;
}

void pushreal32array_(char *x, int *n) {
  pushNarray(x,(*n*32)) ;
}
void popreal32array_(char *x, int *n) {
  popNarray(x,(*n*32)) ;
}
void lookreal32array_(char *x, int *n) {
  lookNarray(x,(*n*32)) ;
}

void pushcomplex4array_(char *x, int *n) {
  pushNarray(x,(*n*4)) ;
}
void popcomplex4array_(char *x, int *n) {
  popNarray(x,(*n*4)) ;
}
void lookcomplex4array_(char *x, int *n) {
  lookNarray(x,(*n*4)) ;
}

void pushcomplex8array_(char *x, int *n) {
  pushNarray(x,(*n*8)) ;
}
void popcomplex8array_(char *x, int *n) {
  popNarray(x,(*n*8)) ;
}
void lookcomplex8array_(char *x, int *n) {
  lookNarray(x,(*n*8)) ;
}

void pushcomplex16array_(char *x, int *n) {
  pushNarray(x,(*n*16)) ;
}
void popcomplex16array_(char *x, int *n) {
  popNarray(x,(*n*16)) ;
}
void lookcomplex16array_(char *x, int *n) {
  lookNarray(x,(*n*16)) ;
}

void pushcomplex32array_(char *x, int *n) {
  pushNarray(x,(*n*32)) ;
}
void popcomplex32array_(char *x, int *n) {
  popNarray(x,(*n*32)) ;
}
void lookcomplex32array_(char *x, int *n) {
  lookNarray(x,(*n*32)) ;
}

/************* Debug displays of the state of the stack: ***********/

void printtopplace_() {
    DoubleChainedBlock *stack = curStack ;
    int nbBlocks = (stack?-1:0) ;
    int remainder = 0;
    while(stack) {
	stack = stack->prev ;
	nbBlocks++ ;
    }
    if (curStack && curStackTop) remainder = curStackTop-(curStack->contents) ;
    /*printf("Stack  top: %i*%i+%i\n",nbBlocks,ONE_BLOCK_SIZE,remainder) ;*/
    printf("Stack size: %f Kbytes\n",nbBlocks*(ONE_BLOCK_SIZE/1024)+((float)remainder)/1024) ;
}

void printlookingplace_() {
    if (lookStack == NULL)
	printtopplace_() ;
    else {
	DoubleChainedBlock *stack = lookStack ;
	int nbBlocks = (stack?-1:0) ;
	while(stack) {
	    stack = stack->prev ;
	    nbBlocks++ ;
	}
	printf("Stack look: %i*%i+%i\n",nbBlocks,ONE_BLOCK_SIZE,lookStackTop-(lookStack->contents)) ;
    }
}
