/*
 * SMART: string matching algorithms research tool.
 * Copyright (C) 2012  Simone Faro and Thierry Lecroq
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 * 
 * contact the authors at: faro@dmi.unict.it, thierry.lecroq@univ-rouen.fr
 * download the tool at: http://www.dmi.unict.it/~faro/smart/
 *
 * This is an implementation of the Optimal Hash algorithm
 */

#include "include/define.h"
#include "include/main.h"



#define ASIZE 256
#define DSIGMA 65536
#define MAX(a,b) ((a) > (b) ? (a) : (b))


int HS(unsigned char *x, int i, int q) {
  int j;

  unsigned int res = x[i];
  for (j = 1; j < q; ++j)
    res += (x[i+j] << j);
  return res%DSIGMA;
}

// Structure to store information of a suffix
struct suffix {
  int index;
  unsigned char *suff;
};

// A comparison function used by sort() to compare two suffixes
int cmp(struct suffix *a, struct suffix *b) {
  return strcmp(a->suff, b->suff) < 0 ? 0 : 1;
}

// This is the main function that takes a string 'txt' of size n as an
// argument, builds and return the suffix array for the given string
int *buildSuffixArray(unsigned char *txt, int n) {
  // A structure to store suffixes and their indexes
  struct suffix suffixes[n];
  int i;

  // Store suffixes and their indexes in an array of structures.
  // The structure is needed to sort the suffixes alphabetically
  // and maintain their old indexes while sorting
  for (i = 0; i < n; i++) {
    suffixes[i].index = i;
    suffixes[i].suff = (txt+i);
  }

  // Sort the suffixes using the comparison function
  // defined above.
  qsort(suffixes, n, sizeof(struct suffix), cmp);

  // Store indexes of all sorted suffixes in the suffix array
  int *suffixArr = (int *)malloc(n*sizeof(int));
  for (i = 0; i < n; i++)
    suffixArr[i] = suffixes[i].index;

  // Return the suffix array
  return suffixArr;
}


int computeMaxLCP(unsigned char *y, int n, int *SA) {
  //int *ISA;
  int ISA[n];
  int j, r;

  for (r = 0; r < n; r++) {
    ISA[SA[r]] = r;
  }
  int ell = 0;
  int res = 0;
  for (j = 0; j < n; j++) {
    ell = MAX(0, ell-1);
    if (ISA[j] > 0) {
      while (MAX(j+ell, SA[ISA[j]-1]+ell) < n && y[j+ell] == y[SA[ISA[j]-1]+ell]) {
        ell++;
      }
    }
    else {
      ell = 0;
    }
    if (ell > res) res = ell;
  }
  return res;
}


int search(unsigned char *x, int m, unsigned char *y, int n) {
  int i, count, j;
  int z[DSIGMA], hs;

  BEGIN_PREPROCESSING
  int *SA = buildSuffixArray(x, m);
  int q = computeMaxLCP(x, m, SA);
  free(SA);
  for (i = 0; i <= m; ++i) {
    y[n+i] = x[i];
  }
  ++q;
  switch (q) {
    case 1 :
      return ohash1(x, m, y, n);
    case 2 :
      return ohash2(x, m, y, n);
    case 3 :
      return ohash3(x, m, y, n);
    case 4 :
      return ohash4(x, m, y, n);
    case 5 :
      return ohash5(x, m, y, n);
    case 6 :
      return ohash6(x, m, y, n);
    case 7 :
      return ohash7(x, m, y, n);
    case 8 :
      return ohash8(x, m, y, n);
    case 9 :
      return ohash9(x, m, y, n);
    case 10 :
      return ohash10(x, m, y, n);
    default :
      return hash8(x, m, y, n);
  }
}

#define RANK8 8

int hash8(unsigned char *x, int m, unsigned char *y, int n) {
   int i, j, sh, shift[WSIZE], sh1, mMinus1, mMinus7, count;
   unsigned int h;
   if(m<8) return -1;

   count = 0;
   mMinus1 = m-1;
   mMinus7 = m-7;
   for (i = 0; i < WSIZE; ++i)
      shift[i] = mMinus7;

   h = x[0];
   h = ((h<<1) + x[1]);
   h = ((h<<1) + x[2]);
   h = ((h<<1) + x[3]);
   h = ((h<<1) + x[4]);
   h = ((h<<1) + x[5]);
   h = ((h<<1) + x[6]);
   h = ((h<<1) + x[7]);
   shift[h%WSIZE] = m-RANK8;
   for (i=RANK8; i < mMinus1; ++i) {
      h = x[i-7];
      h = ((h<<1) + x[i-6]);
      h = ((h<<1) + x[i-5]);
      h = ((h<<1) + x[i-4]);
      h = ((h<<1) + x[i-3]);
      h = ((h<<1) + x[i-2]);
      h = ((h<<1) + x[i-1]);
      h = ((h<<1) + x[i]);
      shift[h%WSIZE] = mMinus1-i;
   }
   h = x[i-7];
   h = ((h<<1) + x[i-6]);
   h = ((h<<1) + x[i-5]);
   h = ((h<<1) + x[i-4]);
   h = ((h<<1) + x[i-3]);
   h = ((h<<1) + x[i-2]);
   h = ((h<<1) + x[i-1]);
   h = ((h<<1) + x[i]);
   sh1 = shift[h%WSIZE];
   shift[h%WSIZE] = 0;
   if(sh1==0) sh1=1;
   END_PREPROCESSING   

   /* Searching */
   BEGIN_SEARCHING
   i = mMinus1;
   memcpy(y+n, x, m);
   while (1) {
      sh = 1;
      while (sh != 0) {
         h = y[i-7];
         h = ((h<<1) + y[i-6]);
         h = ((h<<1) + y[i-5]);
         h = ((h<<1) + y[i-4]);
         h = ((h<<1) + y[i-3]);
         h = ((h<<1) + y[i-2]);
         h = ((h<<1) + y[i-1]);
         h = ((h<<1) + y[i]);
         sh = shift[h%WSIZE];
         i+=sh;
      }
      if (i < n) {
         j=0;
         while(j<m && x[j]==y[i-mMinus1+j]) j++;
         if (j>=m) {
            OUTPUT(i-mMinus1);
         }
         i+=sh1;
      }
      else {
   	END_SEARCHING
      	return count;
      }
   }
}


int ohash1(unsigned char *x, int m, unsigned char *y, int n) {
   int count, j, i, sh, sh1, mMinus1, mMinus2, shift[ASIZE];
   unsigned char h;
   if(m<2) return -1;
   count = 0;
   mMinus1 = m-1;

   for (i = 0; i < ASIZE; ++i)
      shift[i] = m;

   h = x[0];
   shift[h] = mMinus1;
   for (i=1; i < mMinus1; ++i) {
      shift[x[i]] = mMinus1-i;
   }
   h = x[i];
   sh1 = shift[h];
   shift[h] = 0;
   if(sh1==0) sh1=1;
   END_PREPROCESSING

   BEGIN_SEARCHING
   i = mMinus1;
   while (1) {
      sh = 1;
      while (sh != 0) {
         sh = shift[y[i]];
         i+=sh;
      }
      if (i < n) {
         j=0;
         while(j<mMinus1 && x[j]==y[i-mMinus1+j]) j++;
         if (j>=mMinus1) {
            ++count;
         }
         i+=sh1;
      }
      else {
        END_SEARCHING
        return count;
      }
   }
}



#define RANK2 2

int ohash2(unsigned char *x, int m, unsigned char *y, int n) {
   int count, j, i, sh, sh1, mMinus1, mMinus2, shift[DSIGMA];
   unsigned int h;
   if(m<2) return -1; 
   count = 0;
   mMinus1 = m-1;
   mMinus2 = m-2;

   sh = mMinus2;
   if (sh == 0) sh=1;
   for (i = 0; i < DSIGMA; ++i)
      shift[i] = sh;

   h = x[0];
   h = ((h<<8) + x[1]);
   shift[h] = m-RANK2;
   for (i=RANK2; i < mMinus1; ++i) {
      h = x[i-1];
      h = ((h<<8) + x[i]);
      shift[h] = mMinus1-i;
   }   
   h = x[mMinus2];
   h = ((h<<8) + x[mMinus1]);
   sh1 = shift[h];
   shift[h] = 0;
   if(sh1==0) sh1=1;
   END_PREPROCESSING

   BEGIN_SEARCHING
   i = mMinus1;
   while (1) {
      sh = 1;
      while (sh != 0) {
         h = y[i-1];
         h = ((h<<8) + y[i]);
         sh = shift[h];
         i+=sh;
      }
      if (i < n) {
         j=0;
         while(j<mMinus2 && x[j]==y[i-mMinus1+j]) j++;
         if (j>=mMinus2) {
            ++count;
         }
         i+=sh1;
      }
      else {
        END_SEARCHING
        return count;
      }
   }
}


#define RANK3 3

int ohash3(unsigned char *x, int m, unsigned char *y, int n) {
   int count, j, i, sh, sh1, mMinus1, mMinus3, shift[DSIGMA];
   unsigned int h;
   if(m<3) return -1; 
   count = 0;
   mMinus1 = m-1;
   mMinus3 = m-3;

   sh = mMinus3;
   if (sh == 0) sh=1;
   for (i = 0; i < DSIGMA; ++i)
      shift[i] = sh;

   h = x[0];
   h = ((h<<1) + x[1]);
   h = ((h<<1) + x[2]);
   shift[h] = m-RANK3;
   for (i=RANK3; i < mMinus1; ++i) {
      h = x[i-2];
      h = ((h<<1) + x[i-1]);
      h = ((h<<1) + x[i]);
      shift[h] = mMinus1-i;
   }   
   h = x[mMinus3];
   h = ((h<<1) + x[mMinus1-1]);
   h = ((h<<1) + x[mMinus1]);
   sh1 = shift[h];
   shift[h] = 0;
   if(sh1==0) sh1=1;
   END_PREPROCESSING

   BEGIN_SEARCHING
   i = mMinus1;
   while (1) {
      sh = 1;
      while (sh != 0) {
         h = y[i-2];
         h = ((h<<1) + y[i-1]);
         h = ((h<<1) + y[i]);
         sh = shift[h];
         i+=sh;
      }
      if (i < n) {
         j=0;
         while(j<m && x[j]==y[i-mMinus1+j]) j++;
         if (j>=m) {
            ++count;
         }
         i+=sh1;
      }
      else {
        END_SEARCHING
        return count;
      }
   }
}


#define RANK4 4

int ohash4(unsigned char *x, int m, unsigned char *y, int n) {
   int count, j, i, sh, sh1, mMinus1, mMinus4, shift[DSIGMA];
   unsigned int h;
   if(m<4) return -1; 
   count = 0;
   mMinus1 = m-1;
   mMinus4 = m-4;

   sh = mMinus4;
   if (sh == 0) sh=1;
   for (i = 0; i < DSIGMA; ++i)
      shift[i] = sh;

   h = x[0];
   h = ((h<<1) + x[1]);
   h = ((h<<1) + x[2]);
   h = ((h<<1) + x[3]);
   shift[h] = m-RANK4;
   for (i=RANK4; i < mMinus1; ++i) {
      h = x[i-3];
      h = ((h<<1) + x[i-2]);
      h = ((h<<1) + x[i-1]);
      h = ((h<<1) + x[i]);
      shift[h] = mMinus1-i;
   }   
   h = x[mMinus4];
   h = ((h<<1) + x[mMinus1-2]);
   h = ((h<<1) + x[mMinus1-1]);
   h = ((h<<1) + x[mMinus1]);
   sh1 = shift[h];
   shift[h] = 0;
   if(sh1==0) sh1=1;
   END_PREPROCESSING

   BEGIN_SEARCHING
   i = mMinus1;
   while (1) {
      sh = 1;
      while (sh != 0) {
         h = y[i-3];
         h = ((h<<1) + y[i-2]);
         h = ((h<<1) + y[i-1]);
         h = ((h<<1) + y[i]);
         sh = shift[h];
         i+=sh;
      }
      if (i < n) {
         j=0;
         while(j<m && x[j]==y[i-mMinus1+j]) j++;
         if (j>=m) {
            ++count;
         }
         i+=sh1;
      }
      else {
        END_SEARCHING
        return count;
      }
   }
}


#define RANK5 5

int ohash5(unsigned char *x, int m, unsigned char *y, int n) {
   int count, j, i, sh, sh1, mMinus1, mMinus5, shift[DSIGMA];
   unsigned int h;
   if(m<5) return -1; 
   count = 0;
   mMinus1 = m-1;
   mMinus5 = m-5;

   sh = mMinus5;
   if (sh == 0) sh=1;
   for (i = 0; i < DSIGMA; ++i)
      shift[i] = sh;

   h = x[0];
   h = ((h<<1) + x[1]);
   h = ((h<<1) + x[2]);
   h = ((h<<1) + x[3]);
   h = ((h<<1) + x[4]);
   shift[h] = m-RANK5;
   for (i=RANK5; i < mMinus1; ++i) {
      h = x[i-4];
      h = ((h<<1) + x[i-3]);
      h = ((h<<1) + x[i-2]);
      h = ((h<<1) + x[i-1]);
      h = ((h<<1) + x[i]);
      shift[h] = mMinus1-i;
   }   
   h = x[mMinus5];
   h = ((h<<1) + x[mMinus1-3]);
   h = ((h<<1) + x[mMinus1-2]);
   h = ((h<<1) + x[mMinus1-1]);
   h = ((h<<1) + x[mMinus1]);
   sh1 = shift[h];
   shift[h] = 0;
   if(sh1==0) sh1=1;
   END_PREPROCESSING

   BEGIN_SEARCHING
   i = mMinus1;
   while (1) {
      sh = 1;
      while (sh != 0) {
         h = y[i-4];
         h = ((h<<1) + y[i-3]);
         h = ((h<<1) + y[i-2]);
         h = ((h<<1) + y[i-1]);
         h = ((h<<1) + y[i]);
         sh = shift[h];
         i+=sh;
      }
      if (i < n) {
         j=0;
         while(j<m && x[j]==y[i-mMinus1+j]) j++;
         if (j>=m) {
            ++count;
         }
         i+=sh1;
      }
      else {
        END_SEARCHING
        return count;
      }
   }
}


#define RANK6 6

int ohash6(unsigned char *x, int m, unsigned char *y, int n) {
   int count, j, i, sh, sh1, mMinus1, mMinus6, shift[DSIGMA];
   unsigned int h;
   if(m<6) return -1; 
   count = 0;
   mMinus1 = m-1;
   mMinus6 = m-6;

   sh = mMinus6;
   if (sh == 0) sh=1;
   for (i = 0; i < DSIGMA; ++i)
      shift[i] = sh;

   h = x[0];
   h = ((h<<1) + x[1]);
   h = ((h<<1) + x[2]);
   h = ((h<<1) + x[3]);
   h = ((h<<1) + x[4]);
   h = ((h<<1) + x[5]);
   shift[h] = m-RANK6;
   for (i=RANK6; i < mMinus1; ++i) {
      h = x[i-5];
      h = ((h<<1) + x[i-4]);
      h = ((h<<1) + x[i-3]);
      h = ((h<<1) + x[i-2]);
      h = ((h<<1) + x[i-1]);
      h = ((h<<1) + x[i]);
      shift[h] = mMinus1-i;
   }   
   h = x[mMinus6];
   h = ((h<<1) + x[mMinus1-4]);
   h = ((h<<1) + x[mMinus1-3]);
   h = ((h<<1) + x[mMinus1-2]);
   h = ((h<<1) + x[mMinus1-1]);
   h = ((h<<1) + x[mMinus1]);
   sh1 = shift[h];
   shift[h] = 0;
   if(sh1==0) sh1=1;
   END_PREPROCESSING

   BEGIN_SEARCHING
   i = mMinus1;
   while (1) {
      sh = 1;
      while (sh != 0) {
         h = y[i-5];
         h = ((h<<1) + y[i-4]);
         h = ((h<<1) + y[i-3]);
         h = ((h<<1) + y[i-2]);
         h = ((h<<1) + y[i-1]);
         h = ((h<<1) + y[i]);
         sh = shift[h];
         i+=sh;
      }
      if (i < n) {
         j=0;
         while(j<m && x[j]==y[i-mMinus1+j]) j++;
         if (j>=m) {
            ++count;
         }
         i+=sh1;
      }
      else {
        END_SEARCHING
        return count;
      }
   }
}


#define RANK7 7

int ohash7(unsigned char *x, int m, unsigned char *y, int n) {
   int count, j, i, sh, sh1, mMinus1, mMinus7, shift[DSIGMA];
   unsigned int h;
   if(m<7) return -1; 
   count = 0;
   mMinus1 = m-1;
   mMinus7 = m-7;

  sh = mMinus7;
  if (sh == 0) sh=1;
   for (i = 0; i < DSIGMA; ++i)
      shift[i] = sh;

   h = x[0];
   h = ((h<<1) + x[1]);
   h = ((h<<1) + x[2]);
   h = ((h<<1) + x[3]);
   h = ((h<<1) + x[4]);
   h = ((h<<1) + x[5]);
   h = ((h<<1) + x[6]);
   shift[h] = m-RANK7;
   for (i=RANK7; i < mMinus1; ++i) {
      h = x[i-6];
      h = ((h<<1) + x[i-5]);
      h = ((h<<1) + x[i-4]);
      h = ((h<<1) + x[i-3]);
      h = ((h<<1) + x[i-2]);
      h = ((h<<1) + x[i-1]);
      h = ((h<<1) + x[i]);
      shift[h] = mMinus1-i;
   }   
   h = x[mMinus7];
   h = ((h<<1) + x[mMinus1-5]);
   h = ((h<<1) + x[mMinus1-4]);
   h = ((h<<1) + x[mMinus1-3]);
   h = ((h<<1) + x[mMinus1-2]);
   h = ((h<<1) + x[mMinus1-1]);
   h = ((h<<1) + x[mMinus1]);
   sh1 = shift[h];
   shift[h] = 0;
   if(sh1==0) sh1=1;
   END_PREPROCESSING

   BEGIN_SEARCHING
   i = mMinus1;
   while (1) {
      sh = 1;
      while (sh != 0) {
         h = y[i-6];
         h = ((h<<1) + y[i-5]);
         h = ((h<<1) + y[i-4]);
         h = ((h<<1) + y[i-3]);
         h = ((h<<1) + y[i-2]);
         h = ((h<<1) + y[i-1]);
         h = ((h<<1) + y[i]);
         sh = shift[h];
         i+=sh;
      }
      if (i < n) {
         j=0;
         while(j<m && x[j]==y[i-mMinus1+j]) j++;
         if (j>=m) {
            ++count;
         }
         i+=sh1;
      }
      else {
        END_SEARCHING
        return count;
      }
   }
}


#define RANK8 8

int ohash8(unsigned char *x, int m, unsigned char *y, int n) {
   int count, j, i, sh, sh1, mMinus1, mMinus8, shift[DSIGMA];
   unsigned int h;
   if(m<8) return -1; 
   count = 0;
   mMinus1 = m-1;
   mMinus8 = m-8;

   sh = mMinus8;
   if (sh == 0) sh=1;
   for (i = 0; i < DSIGMA; ++i)
      shift[i] = sh;

   h = x[0];
   h = ((h<<1) + x[1]);
   h = ((h<<1) + x[2]);
   h = ((h<<1) + x[3]);
   h = ((h<<1) + x[4]);
   h = ((h<<1) + x[5]);
   h = ((h<<1) + x[6]);
   h = ((h<<1) + x[7]);
   shift[h] = m-RANK8;
   for (i=RANK8; i < mMinus1; ++i) {
      h = x[i-7];
      h = ((h<<1) + x[i-6]);
      h = ((h<<1) + x[i-5]);
      h = ((h<<1) + x[i-4]);
      h = ((h<<1) + x[i-3]);
      h = ((h<<1) + x[i-2]);
      h = ((h<<1) + x[i-1]);
      h = ((h<<1) + x[i]);
      shift[h] = mMinus1-i;
   }   
   h = x[mMinus8];
   h = ((h<<1) + x[mMinus1-6]);
   h = ((h<<1) + x[mMinus1-5]);
   h = ((h<<1) + x[mMinus1-4]);
   h = ((h<<1) + x[mMinus1-3]);
   h = ((h<<1) + x[mMinus1-2]);
   h = ((h<<1) + x[mMinus1-1]);
   h = ((h<<1) + x[mMinus1]);
   sh1 = shift[h];
   shift[h] = 0;
   if(sh1==0) sh1=1;
   END_PREPROCESSING

   BEGIN_SEARCHING
   i = mMinus1;
   while (1) {
      sh = 1;
      while (sh != 0) {
         h = y[i-7];
         h = ((h<<1) + y[i-6]);
         h = ((h<<1) + y[i-5]);
         h = ((h<<1) + y[i-4]);
         h = ((h<<1) + y[i-3]);
         h = ((h<<1) + y[i-2]);
         h = ((h<<1) + y[i-1]);
         h = ((h<<1) + y[i]);
         sh = shift[h];
         i+=sh;
      }
      if (i < n) {
         j=0;
         while(j<m && x[j]==y[i-mMinus1+j]) j++;
         if (j>=m) {
            //OUTPUT(i-mMinus1);
            ++count;
         }
         i+=sh1;
      }
      else {
        END_SEARCHING
        return count;
      }
   }
}


#define RANK9 9

int ohash9(unsigned char *x, int m, unsigned char *y, int n) {
   int count, j, i, sh, sh1, mMinus1, mMinus9, shift[DSIGMA];
   unsigned int h;
   if(m<9) return -1; 
   count = 0;
   mMinus1 = m-1;
   mMinus9 = m-9;

  sh = mMinus9;
  if (sh == 0) sh=1;
   for (i = 0; i < DSIGMA; ++i)
      shift[i] = sh;

   h = x[0];
   h = ((h<<1) + x[1]);
   h = ((h<<1) + x[2]);
   h = ((h<<1) + x[3]);
   h = ((h<<1) + x[4]);
   h = ((h<<1) + x[5]);
   h = ((h<<1) + x[6]);
   h = ((h<<1) + x[7]);
   h = ((h<<1) + x[8]);
   shift[h] = m-RANK9;
   for (i=RANK9; i < mMinus1; ++i) {
      h = x[i-8];
      h = ((h<<1) + x[i-7]);
      h = ((h<<1) + x[i-6]);
      h = ((h<<1) + x[i-5]);
      h = ((h<<1) + x[i-4]);
      h = ((h<<1) + x[i-3]);
      h = ((h<<1) + x[i-2]);
      h = ((h<<1) + x[i-1]);
      h = ((h<<1) + x[i]);
      shift[h] = mMinus1-i;
   }   
   h = x[mMinus9];
   h = ((h<<1) + x[mMinus1-7]);
   h = ((h<<1) + x[mMinus1-6]);
   h = ((h<<1) + x[mMinus1-5]);
   h = ((h<<1) + x[mMinus1-4]);
   h = ((h<<1) + x[mMinus1-3]);
   h = ((h<<1) + x[mMinus1-2]);
   h = ((h<<1) + x[mMinus1-1]);
   h = ((h<<1) + x[mMinus1]);
   sh1 = shift[h];
   shift[h] = 0;
   if(sh1==0) sh1=1;
   END_PREPROCESSING

   BEGIN_SEARCHING
   i = mMinus1;
   while (1) {
      sh = 1;
      while (sh != 0) {
         h = y[i-8];
         h = ((h<<1) + y[i-7]);
         h = ((h<<1) + y[i-6]);
         h = ((h<<1) + y[i-5]);
         h = ((h<<1) + y[i-4]);
         h = ((h<<1) + y[i-3]);
         h = ((h<<1) + y[i-2]);
         h = ((h<<1) + y[i-1]);
         h = ((h<<1) + y[i]);
         sh = shift[h];
         i+=sh;
      }
      if (i < n) {
         j=0;
         while(j<m && x[j]==y[i-mMinus1+j]) j++;
         if (j>=m) {
            ++count;
         }
         i+=sh1;
      }
      else {
        END_SEARCHING
        return count;
      }
   }
}


#define RANK10 10

int ohash10(unsigned char *x, int m, unsigned char *y, int n) {
   int count, j, i, sh, sh1, mMinus1, mMinus10, shift[DSIGMA];
   unsigned int h;
   if(m<10) return -1; 
   count = 0;
   mMinus1 = m-1;
   mMinus10 = m-10;

   sh = mMinus10;
   if (sh == 0) sh=1;
   for (i = 0; i < DSIGMA; ++i)
      shift[i] = sh;

   h = x[0];
   h = ((h<<1) + x[1]);
   h = ((h<<1) + x[2]);
   h = ((h<<1) + x[3]);
   h = ((h<<1) + x[4]);
   h = ((h<<1) + x[5]);
   h = ((h<<1) + x[6]);
   h = ((h<<1) + x[7]);
   h = ((h<<1) + x[8]);
   h = ((h<<1) + x[9]);
   shift[h%DSIGMA] = m-RANK10;
   for (i=RANK10; i < mMinus1; ++i) {
      h = x[i-9];
      h = ((h<<1) + x[i-8]);
      h = ((h<<1) + x[i-7]);
      h = ((h<<1) + x[i-6]);
      h = ((h<<1) + x[i-5]);
      h = ((h<<1) + x[i-4]);
      h = ((h<<1) + x[i-3]);
      h = ((h<<1) + x[i-2]);
      h = ((h<<1) + x[i-1]);
      h = ((h<<1) + x[i]);
      shift[h%DSIGMA] = mMinus1-i;
   }   
   h = x[mMinus10];
   h = ((h<<1) + x[mMinus1-8]);
   h = ((h<<1) + x[mMinus1-7]);
   h = ((h<<1) + x[mMinus1-6]);
   h = ((h<<1) + x[mMinus1-5]);
   h = ((h<<1) + x[mMinus1-4]);
   h = ((h<<1) + x[mMinus1-3]);
   h = ((h<<1) + x[mMinus1-2]);
   h = ((h<<1) + x[mMinus1-1]);
   h = ((h<<1) + x[mMinus1]);
   sh1 = shift[h%DSIGMA];
   shift[h%DSIGMA] = 0;
   if(sh1==0) sh1=1;
   END_PREPROCESSING

   BEGIN_SEARCHING
   i = mMinus1;
   while (1) {
      sh = 1;
      while (sh != 0) {
         h = y[i-9];
         h = ((h<<1) + y[i-8]);
         h = ((h<<1) + y[i-7]);
         h = ((h<<1) + y[i-6]);
         h = ((h<<1) + y[i-5]);
         h = ((h<<1) + y[i-4]);
         h = ((h<<1) + y[i-3]);
         h = ((h<<1) + y[i-2]);
         h = ((h<<1) + y[i-1]);
         h = ((h<<1) + y[i]);
         sh = shift[h%DSIGMA];
         i+=sh;
      }
      if (i < n) {
         j=0;
         while(j<m && x[j]==y[i-mMinus1+j]) j++;
         if (j>=m) {
            ++count;
         }
         i+=sh1;
      }
      else {
        END_SEARCHING
        return count;
      }
   }
}
