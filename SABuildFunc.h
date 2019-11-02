#ifndef SABUILDFUNC_H_INCLUDED
#define SABUILDFUNC_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include "HelperFunction.h"

/*
 * Test functions.
 */
void _quickSortTest();
void _suffixArrayQuickSortTest();
void _compareSuffixTest();
void _inverseSATest();

/*
 * Important functions.
 */
void inverseSA(int SA[], int SA_inverse[], int length);
void quickSort(char *str[], int left, int right);
void suffixArrayQuickSort(int SA[], char T[], int left, int right);
int compareSuffix(int i, int j, char T[]);


/**
 * Test quick sort.
 */
void _quickSortTest() {

    printf("section\\\\ Quick Sort test\n");
    char *str[] = {"abb", "aab", "aaa", "bbdsas",
                   "bsass", "brabrabra", "blahblah", "atosh"
                  };
    int strLength = 8;

    quickSort(str, 0, strLength-1);

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
    printf("\n ******* Suffix Array Quick Sort Test *********\n");
    char T[] = "acaaccg$";
    int SA[] = {0,1,2,3,4,5,6,7};
    printf("total array: %s\n", T);

    suffixArrayQuickSort(SA, T, 0, 7);

    printf("sorted suffix array: ");
    int i = 0;;
    for( i = 0; i < 8; i++) {
        printf("%d ", SA[i]);
    }
    printf("\n");
}

/**
 * Test comparing suffixes.
 */
void _compareSuffixTest() {
    char T[] = "aabcdabcdefg";
    int result = compareSuffix(0, 5, T);
    printf("result: %d\n", result);
    printf("result should be -1\n");
}

/**
 * Test inverse the SA[].
 */
void _inverseSATest() {
    printf("\n ****** inverse SA test ******* \n");
    char T[] = "acaaccg$";
    int SA[] = {0,1,2,3,4,5,6,7};
    int i = 0;
    int arrayLength = sizeof(SA)/sizeof(int);
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
    inverseSA(SA, SA_inverse, arrayLength);
    printf("inversed suffix array: ");
    for(i = 0; i < arrayLength; i++) {
        printf("%d ", SA_inverse[i]);
    }
    printf("\n");

    printf("\n");
}

void _psiArrayBuildTest(){
    printf("\n ****** psi Array Build Test ******* \n");
    char T[] = "acaaccg$";
    int SA[] = {0,1,2,3,4,5,6,7};
    int i = 0;
    int arrayLength = sizeof(SA)/sizeof(int);
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
    inverseSA(SA, SA_inverse, arrayLength);
    printf("inversed suffix array: ");
    for(i = 0; i < arrayLength; i++) {
        printf("%d ", SA_inverse[i]);
    }
    printf("\n");

    // build Psi array based on SA[] and SA_inverse[]
    int Psi[arrayLength];
    psiArrayBuild(SA, SA_inverse, Psi, arrayLength);
    printf("Psi array: ");
    for(i = 0; i < arrayLength; i++){
        printf("%d ", Psi[i]);
    }
    printf("\n");

    printf("\n");

}
/**
 * Build psi function array.
 *
 * @param SA[] suffix array
 * @param SA_inverse[] inverse of suffix array
 * @param Psi[] psi function array
 * @param length length of these arrays (all the same)
 */
void psiArrayBuild(int SA[], int SA_inverse[], int Psi[], int length){
    int i = 0;
    Psi[0] = SA_inverse[0];
    for(i = 1; i < length; i++){
        Psi[i] = SA_inverse[SA[i] + 1];
    }
}

/**
 * Create an inverse of suffix array SA[].
 *
 * @param SA[] suffix array
 * @param SA_inverse[] inverse of suffix array (to be created)
 * @param length length of suffix array
 */
void inverseSA(int SA[], int SA_inverse[], int length) {
    int i = 0;
    for(i = 0; i < length; i++) {
        SA_inverse[SA[i]] = i;
    }
}

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
//	printf("key: %s ****", key);
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
        str[i] = str[j];	// update that str to left part
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
    quickSort(str, i+1, right);
}

/**
 * Use quick sort to sort suffix arrays.
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
        // find the suffix that is samller than key-suffix in right part
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
    suffixArrayQuickSort(SA, T, i+1, right);
}

/**
 * Compare the 2 suffices in a String.
 *
 * @param i start position of suffix 1
 * @param j start position of suffix 2
 * @param T the whole string (char array)
 * @return 1 when suffix[i] > suffix[j]; -1 when suffix[i] < suffix[j];
 * 		0 when suffix[i] equals suffix[j]
 */
int compareSuffix(int i, int j, char T[]) {
//    printf("comparing suffix(%d, %c, %d) and suffix(%d, %c, %d)\n", i, T[i], T[i], j, T[j], T[j]);
    while(T[i] != '\0' && T[j] != '\0') {
        if(T[i] == '$' && T[j] == '$') {
            // if suffix[i] and suffix[j] both end, suffix[i] = suffix[j]
            return 0;
        } else if(T[i] == '$' && T[j] != '$') {
            // if suffix[i] ends first, suffix[i] < suffix[j]
            return -1;
        } else if(T[i] != '$' && T[j] == '$') {
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

#endif // SABUILDFUNC_H_INCLUDED
