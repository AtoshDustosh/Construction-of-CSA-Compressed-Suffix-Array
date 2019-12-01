#ifndef SABUILDFUNC_H_INCLUDED
#define SABUILDFUNC_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "HelperFunction.h"

/*
 * Test functions.
 */

void _quickSortTest();
void _suffixArrayQuickSortTest();
void _compareSuffixTest();
void _inverseSAWholeTest();
void _psiArrayBuildWholeTest();
void _binarySearchBoundTest();
void _CSABinaryBoundSearchTest();
void _CSABinarySearchOrderValueTest();    // this test is unnecessary ... perhaps
void _convertPsiToBWTTest();


/*
 * Important functions.
 */
void convertPsiToBWT(char* T, int* Psi, char* BWT);
void CSABinarySearchOrderValue(int* SA, int* Psi, int lc, int rc, int prevOrderValue,
                               int* max_b);
void CSABinaryBoundSearch(char* T, int* SA, char c, int* left, int* right);
void directBinarySearchBound(char* chArray, char c, int* left, int* right);
void psiArrayBuildWhole(int SA[], int SA_inverse[], int Psi[], int length);
void inverseSAWhole(int SA[], int SA_inverse[], int length);
void quickSort(char *str[], int left, int right);
void suffixArrayQuickSort(int SA[], char T[], int left, int right);
int compareSuffix(int i, int j, char T[]);


/**
 * Quick sort an array that stores strings.
 *
 * @param str[] an array of strings
 * @param left start position (index-start)
 * @param right end position (index-end)
 */
void quickSort(char *str[], int left, int right) {
    if(left >= right) {
        return;
    }
    int i = left;
    int j = right;
    char *key = str[left];
//  printf("key: %s ****", key);
//    int index = 0;
//    for(index = 0; index < right-left+1; index++){
//        printf(" %s", str[index]);
//    }
//    printf("\n");


    while(i < j) {
        // find the str that is samller than key in right part
        while(i < j && strcmp(key, str[j]) <= 0) {
            j--;
        }
        str[i] = str[j];    // update that str to left part
        // find the str that is larger than key in left part
        while(i < j && strcmp(str[i], key) <= 0) {
            i++;
        }
        str[j] = str[i];
    }
    // finally update the key str
    str[i] = key;

    /*
     * iteration part
     */
    quickSort(str, left, i);
    quickSort(str, i + 1, right);
}

/**
 * Use quick sort to sort suffix arrays of T[] and store the lex-order in
 * SA[].
 *
 * @param SA[] suffix array
 * @param T[] char array responding to SA[] - a string
 * @param left left bounder index
 * @param right right bounder index
 */
void suffixArrayQuickSort(int SA[], char T[], int left, int right) {
    if(left >= right) {
        return;
    }
    int i = left;
    int j = right;
    int key = SA[left]; // the key-suffix starts at SA[left]

    while(i < j) {
        // find the suffix that is smaller than key-suffix in right part
        while(i < j && compareSuffix(key, SA[j], T) <= 0) {
            j--;
        }
        SA[i] = SA[j];  // update that SA to left part
        // find the suffix that is larger than key-suffix in left part
        while(i < j && compareSuffix(SA[i], key, T) <= 0) {
            i++;
        }
        SA[j] = SA[i];  // update that SA to right part
    }
    // finally update the key SA
    SA[i] = key;

    /*
     * iteration part
     */
    suffixArrayQuickSort(SA, T, left, i);
    suffixArrayQuickSort(SA, T, i + 1, right);
}

/**
 * Compare the 2 suffices in a String.
 *
 * @param i start position of suffix 1
 * @param j start position of suffix 2
 * @param T the whole string (char array)
 * @return 1 when suffix[i] > suffix[j]; -1 when suffix[i] < suffix[j];
 *      0 when suffix[i] equals suffix[j]
 */
int compareSuffix(int i, int j, char T[]) {
//    printf("comparing suffix(%d, %c, %d) and suffix(%d, %c, %d)\n", i, T[i], T[i], j, T[j], T[j]);
    while(1) {
        if((T[i] == '$' && T[j] == '$') || (T[i] == '\0' && T[j] == '\0')) {
            // if suffix[i] and suffix[j] both end, suffix[i] = suffix[j]
            return 0;
        } else if((T[i] == '$' && T[j] != '$') || (T[i] == '\0' && T[j] != '\0')) {
            // if suffix[i] ends first, suffix[i] < suffix[j]
            return -1;
        } else if((T[i] != '$' && T[j] == '$') || (T[i] != '\0' && T[j] == '\0')) {
            // if suffix[j] ends first, suffix[i] > suffix[j]
            return 1;
        }
        if(T[i] == T[j]) {
            // if T[i] == T[j], continue to compare the next character
            i++;
            j++;
        } else if(T[i] < T[j]) {
            return -1;
        } else if(T[i] > T[j]) {
            return 1;
        }
    }
}

/**
 * Build psi function array of the whole array.
 *
 * @param SA[] suffix array
 * @param SA_inverse[] inverse of suffix array
 * @param Psi[] psi function array
 * @param length length of these arrays (all the same)
 */
void psiArrayBuildWhole(int SA[], int SA_inverse[], int Psi[], int length) {
    int i = 0;
    Psi[0] = SA_inverse[0];
    for(i = 1; i < length; i++) {
        Psi[i] = SA_inverse[SA[i] + 1];
    }
}

/**
 * Create an inverse of suffix array SA[] of the whole array.
 *
 * @param SA[] suffix array
 * @param SA_inverse[] inverse of suffix array (to be created)
 * @param length length of suffix array
 */
void inverseSAWhole(int SA[], int SA_inverse[], int length) {
    int i = 0;
    for(i = 0; i < length; i++) {
        SA_inverse[SA[i]] = i;
    }
}

/**
 * Binary search chArray to find the left and right bounds of a character c.
 * Note that chArray is a sorted array by lex-order.
 *
 * @param chArray a char array sorted by lex-order
 * @param c character whose bounds are to be searched for
 * @param / @return left left bound
 * @param / @return right right bound
 */
void directBinarySearchBound(char* chArray, char c, int* left, int* right) {
    int leftBorder = *left;
    int rightBorder = *right;
    int lc = *left;
    int rc = *right;

    // find the left bound
    leftBorder = *left;
    rightBorder = *right;
    while(leftBorder < rightBorder) {
        int mid = leftBorder + (rightBorder - leftBorder) / 2;
        char midCh = chArray[mid];
        if(mid == leftBorder) {
            if(chArray[leftBorder] == c) {
                break;
            } else {
                leftBorder++;
            }
        } else if(midCh < c) {
            leftBorder = mid;
        } else if(midCh >= c) {
            rightBorder = mid;
        }
    }
    lc = leftBorder;

    // find the right bound
    leftBorder = *left;
    rightBorder = *right;
    while(leftBorder < rightBorder) {
        int mid = leftBorder + (rightBorder - leftBorder) / 2;
        char midCh = chArray[mid];
        if(mid == leftBorder) {
            if(chArray[rightBorder] == c) {
                break;
            } else {
                rightBorder--;
            }
        } else if(midCh <= c) {
            leftBorder = mid;
        } else if(midCh > c) {
            rightBorder = mid;
        }
    }
    rc = rightBorder;

    // assign the results to pointers
    *left = lc;
    *right = rc;
}

/**
 * Binary search T_SA_inverse[] to find the left and right bounds of a character c.
 * Note that T[SA[i]]] is like a map from i to the i-th (by lex-order) suffix of the
 * designated part in T (from \arg left to \arg right).
 *
 * @param T DNA sequence plus a '$'
 * @param SA SA of T in the designated interval
 * @param c character whose bounds are to be searched for
 * @param / @return left left bound
 * @param / @return right right bound
 */
void CSABinaryBoundSearch(char* T, int* SA, char c, int* left, int* right) {
    int leftBorder = *left;
    int rightBorder = *right;
    int lc = *left;
    int rc = *right;

    // find the left bound
    leftBorder = *left;
    rightBorder = *right;
    while(leftBorder < rightBorder) {
        int mid = leftBorder + (rightBorder - leftBorder) / 2;
        char midCh = T[SA[mid]];
        if(mid == leftBorder) {
            if(T[SA[leftBorder]] == c) {
                break;
            } else {
                leftBorder++;
            }
        } else if(midCh < c) {
            leftBorder = mid;
        } else if(midCh >= c) {
            rightBorder = mid;
        }
    }
    if(T[SA[rightBorder]] < c) {
        leftBorder++;
    }
    lc = leftBorder;

    // find the right bound
    leftBorder = *left;
    rightBorder = *right;
    while(leftBorder < rightBorder) {
        int mid = leftBorder + (rightBorder - leftBorder) / 2;
        char midCh = T[SA[mid]];
        if(mid == leftBorder) {
            if(T[SA[rightBorder]] <= c) {
                break;
            } else {
                rightBorder--;
            }
        } else if(midCh <= c) {
            leftBorder = mid;
        } else if(midCh > c) {
            rightBorder = mid;
        }
    }
    rc = rightBorder;

    // assign the results to pointers
    *left = lc;
    *right = rc;
}

/**
 * Find the maximum b that is within [lc, rc] and satisfies condition that for any
 * order(cX, T'), Psi_T'[SA[b]] <= order(X, T').
 *
 * @param SA SA of T in the designated interval
 * @param Psi Psi of T in the designated interval
 * @param lc left bound of the interval
 * @param rc right bound of the interval
 * @param prevOrderValue value of previous order func, which actually is order(X, T')
 * @param / @return max_b a pointer to maximum b; returns -1 if there is no such b
 */
void CSABinarySearchOrderValue(int* SA, int* Psi, int lc, int rc, int prevOrderValue,
                               int* max_b) {
//    printf("\n");
//     find the right bound
    while(lc <= rc) {
        int mid = lc + (rc - lc) / 2;
        int midValue = Psi[mid];
//        printf("midValue(%d), lcValue(%d), rcValue(%d)\n", midValue, Psi[lc], Psi[rc]);
        if(mid == lc) {
            // process the last cases
            if(midValue <= prevOrderValue) {
                if(Psi[rc] > prevOrderValue) {
                    *max_b = mid;
                } else {
                    *max_b = rc;
                }
            } else {
                *max_b = -1;
            }
            return;
        } else if(midValue <= prevOrderValue) {
            lc = mid;
        } else if(midValue > prevOrderValue) {
            rc = mid;
        }
    }
}

/**
 * Convert Psi of T to get BWT.
 *
 * @param T DNA sequence plus a '$'
 * @param Psi Psi func of T
 * @param BWT Burrows-Wheeler Transform of T
 */
void convertPsiToBWT(char* T, int* Psi, char* BWT) {
    int i = 0;
    int arrayLength = strlen(T);

    // convert Psi of T to get SA of T
    int* SA = (int*)malloc(sizeof(int) * arrayLength);
    int x = 0;

    SA[Psi[x]] = 0;
    x = Psi[0];
    for(i = 0; i < arrayLength; i++) {
//        printf("%d -> x: %d, Psi[x]: %d, SA[Psi[x]]: %d\n",
//               i, x, Psi[x], SA[Psi[x]]);
        if(x == 0) {
            continue;   // this is necessary - cannot be deleted
        }
        SA[Psi[x]] = SA[x] + 1;
        x = Psi[x];
    }

    // convert SA of T to get BWT of T
    for(i = 0; i < arrayLength; i++) {
        if(SA[i] == 0) {
            BWT[i] = T[arrayLength - 1];
            continue;
        }
        BWT[i] = T[SA[i] - 1];
    }

    free(SA);
}




//////////////////////////////////////// the following funcs will not be used ///////////////////////////////


/**
 * Test converting Psi func to BWT.
 */
void _convertPsiToBWTTest() {
    printf("\n ******* _convertPsiToBWTTest *********\n");
    int i = 0;
    char* T = "acacgt$";
    int Psi[7] = {1, 3, 4, 2, 5, 6, 0};
    char* BWT = (char*)malloc(sizeof(char) * strlen(T));

    convertPsiToBWT(T, Psi, BWT);

    printf("i\tT[]\t");
    printf("Psi[]\t");
    printf("BWT[]\t");
    printf("\n");
    for(i = 0; i < strlen(T); i++){
        printf("%d\t%c\t", i, T[i]);
        printf("%d\t", Psi[i]);
        printf("%c\t", BWT[i]);
        printf("\n");
    }

    free(BWT);
}

/**
 * Test binary search order value.
 */
void _CSABinarySearchOrderValueTest() {
    printf("\n ******* _CSABinarySearchOrderValueTest *********\n");
    int i = 0;
    char* T_apostrophe = "cagac$";
    char* T_i = "gca";
    int SA_apostrophe[6] = {5, 3, 1, 4, 0, 2};
    int Psi_apostrophe[6] = {4, 3, 5, 0, 2, 1};
    int iLength = strlen(T_i);
    int apostropheLength = strlen(T_apostrophe);
    int* order = (int*)malloc(sizeof(int) * iLength);

    int prevOrderValue = 4;
    for(i = iLength - 1; i >= 0; i--) {
        char c = T_i[i];  // the character in the formula
        int lc = 0;
        int rc = apostropheLength - 1;
        int orderValue = 0;
        CSABinaryBoundSearch(T_apostrophe, SA_apostrophe, c, &lc, &rc);
        // T[SA[lc]] ~ T[SA[rc]] represents the field of c
        // implement of condition ¡Á[b] -> ¡Á[SA[b]], lc <= b <= rc
        printf("(%c) -> lc: %d, rc: %d\t", c, lc, rc);

        if(lc > rc) {
            orderValue = lc - 1;   // this is modified ... on my own will
        } else {
            // find the max b that satisfies condition that for any order(cX, T'), Psi_T'[b] <= order(X, T')
            int max_b = 0;
            CSABinarySearchOrderValue(SA_apostrophe, Psi_apostrophe, lc, rc, prevOrderValue, &max_b);
            if(max_b == -1) {
                orderValue = lc - 1;
            } else {
                orderValue = max_b;
            }
        }

        printf("prevOrderValue: %d\torderValue(%d): %d\t", prevOrderValue, i, orderValue);

        order[i] = orderValue;
        prevOrderValue = orderValue;

        printf("\n");
    }

    free(order);
}

/**
 * Test binary search bound.
 */
void _CSABinaryBoundSearchTest() {
    printf("\n ******* _CSABinaryBoundSearchTest *********\n");
    char* T = "ttaaccttaaata$";
    int SA[] = {13, 12, 8, 2, 9, 3, 10, 4, 5, 11, 7, 1, 6, 0};
    char c = 'g';
    int left = 0;
    int right = strlen(T) - 1;
    int i = 0;
    for(i = 0; i < strlen(T); i++) {
        printf("%c", T[SA[i]]);
    }
    printf("\n");

    printf("length: %d\n", (int)strlen(T));

    c = 'g';
    left = 0;
    right = strlen(T) - 1;
    CSABinaryBoundSearch(T, SA, c, &left, &right);
    printf("%c -> %d, %d\n", c, left, right);

    c = 'c';
    left = 0;
    right = strlen(T) - 1;
    CSABinaryBoundSearch(T, SA, c, &left, &right);
    printf("%c -> %d, %d\n", c, left, right);
}

/**
 * Test binary search bound.
 */
void _binarySearchBoundTest() {
    printf("\n ******* _binarySearchBoundTest *********\n");
    char* chArray = "$aaagggggggttttt";
    char c = 'c';
    int left = 0;
    int right = strlen(chArray) - 1;
    int i = 0;

    for(i = 0; i < strlen(chArray); i++) {
        printf("%c", chArray[i]);
    }
    printf("\n");

    c = 't';
    left = 0;
    right = strlen(chArray) - 1;
    directBinarySearchBound(chArray, c, &left, &right);
    printf("length: %d. -> %d, %d\n", (int)strlen(chArray), left, right);

    c = 'c';
    left = 0;
    right = strlen(chArray) - 1;
    directBinarySearchBound(chArray, c, &left, &right);
    printf("length: %d. -> %d, %d\n", (int)strlen(chArray), left, right);
}

/**
 * Test quick sort.
 */
void _quickSortTest() {
    printf("\n ******* _quickSortTest *********\n");
    char *str[] = {"abb", "aab", "aaa", "bbdsas",
                   "bsass", "brabrabra", "blahblah", "atosh"
                  };
    int strLength = 8;

    quickSort(str, 0, strLength - 1);

    int i = 0;
    printf("sorted strings' array: ");
    for(i = 0; i < strLength; i++) {
        printf(" %s", str[i]);
    }
    printf("\n\n");
}

/**
 * Test quick sort of suffix arrays.
 */
void _suffixArrayQuickSortTest() {
    printf("\n ******* _suffixArrayQuickSortTest *********\n");
    int i = 0;
    char T[] = "acaaccg$";
    int SA[] = {0, 1, 2, 3, 4, 5, 6, 7};
    int sortLength = strlen(T);
    printf("total array: %s\n", T);

    suffixArrayQuickSort(SA, T, 0, sortLength - 1);

    printf("sorted suffix array\n");
    for( i = 0; i < sortLength; i++) {
        printf("%d\t%c\n", SA[i], T[SA[i]]);
    }
    printf("\n");
}

/**
 * Test comparing suffixes.
 */
void _compareSuffixTest() {
    printf("\n ******* _compareSuffixTest *********\n");
    char* T = "aabcdabcdefg";
    int result = compareSuffix(0, 5, T);

    printf("result: %d\n", result);
    printf("result should be %d\n", strcmp("aabcdabcdefg", "abcdefg"));

    T = "cttagctt";
    result = compareSuffix(0, 5, T);
    printf("result: %d\n", result);
    printf("result should be %d\n", strcmp("cttagctt", "ctt"));

}

/**
 * Test inverse the SA[].
 */
void _inverseSAWholeTest() {
    printf("\n ******* _inverseSAWholeTest *********\n");
    char T[] = "acaaccg$";
    int SA[] = {0, 1, 2, 3, 4, 5, 6, 7};
    int i = 0;
    int arrayLength = strlen(T);
    printf("total array: %s\n", T);

    // create suffix array of T
    suffixArrayQuickSort(SA, T, 0, arrayLength - 1);
    printf("sorted suffix array: ");
    for(i = 0; i < arrayLength; i++) {
        printf("%d ", SA[i]);
    }
    printf("\n");

    // get inversed suffix array (x->y => y->x)
    int SA_inverse[arrayLength];
    inverseSAWhole(SA, SA_inverse, arrayLength);
    printf("inversed suffix array: ");
    for(i = 0; i < arrayLength; i++) {
        printf("%d ", SA_inverse[i]);
    }
    printf("\n");

    printf("\n");
}

/**
 * Test building Psi[] array. Psi[i] = SA_inverse[SA[i]].
 */
void _psiArrayBuildWholeTest() {
    printf("\n ******* _psiArrayBuildWholeTest *********\n");
    char T[] = "acaaccg$";
    int SA[] = {0, 1, 2, 3, 4, 5, 6, 7};
    int i = 0;
    int arrayLength = strlen(T);
    printf("total array: %s\n", T);

    // create suffix array of T
    suffixArrayQuickSort(SA, T, 0, arrayLength - 1);
    printf("sorted suffix array: ");
    for(i = 0; i < arrayLength; i++) {
        printf("%d ", SA[i]);
    }
    printf("\n");

    // get inversed suffix array (x->y => y->x)
    int SA_inverse[arrayLength];
    inverseSAWhole(SA, SA_inverse, arrayLength);
    printf("inversed suffix array: ");
    for(i = 0; i < arrayLength; i++) {
        printf("%d ", SA_inverse[i]);
    }
    printf("\n");

    // build Psi array based on SA[] and SA_inverse[]
    int Psi[arrayLength];
    psiArrayBuildWhole(SA, SA_inverse, Psi, arrayLength);
    printf("Psi array: ");
    for(i = 0; i < arrayLength; i++) {
        printf("%d ", Psi[i]);
    }
    printf("\n");

    printf("\n");

}


#endif // SABUILDFUNC_H_INCLUDED
