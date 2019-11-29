#ifndef _HEADFILETEST_H_
#define _HEADFILETEST_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "HelperFunction.h"

int _valueTest = 0;
int* _arrayTest = NULL;

void _CLanguageReview() ;
void _mathematicalFuncsTest();
void _globalVariableInHeaderFile(int arrayTest[]);
void _timeOperationTest();

/**
 * Test time operations of C language.
 */
void _timeOperationTest() {
    int i = 0, j = 0;
    long value = 0;
    value = clock();
    printf("clock: %ld\n", value);

    for(i = 0; i < 10000; i++) {
        for(j = 0; j < 100000; j++) {
            // this loop takes about 4 seconds
        }
    }
    value = clock();
    printf("clock: %ld\n", value);
}

/**
 * Review usage of C language.
 */
void _CLanguageReview() {
    printf("section\\\\ C language review test\n");
    // review string[] array construction
    printf("String operation\n");
    char *str[] = {"abb", "aab", "aaa", "bbdsas",
                   "bsass", "brabrabra", "blahblah", "atosh"
                  };
    int strLength = 8;

    // review string[] array access
    int i = 0;
    printf("String length: %d\n", strLength);
    printf("str[4]: %s\n", str[4]);
    for(i = 0; i < strLength; i++) {
        printf("%s\n", str[i]);
    }
    printf("%s\n", str[1]);
    // review string[].getChar(i) -> integer operation
    char *strCh = str[1];
    int count = 0;
    while(*strCh != '\0') {
        printf("%d: %c%d\n", count, *strCh, *strCh);
        strCh = strCh + 1;
        count ++;
    }
    printf("\n\n");
}

/**
 * Test usages of mathematical functions in C programming.
 */
void _mathematicalFuncsTest() {
    if(_valueTest < 1) {
        printf("sin(pi/3) = %f\n", sin(3.1415926 / 3));
        printf("ln 2.7183 = %f\n", log(2.7183));
        printf("log(2) 7 = %f\n", log10(7) / log10(2));
    }
    printf("_valueTest: %d\n", _valueTest);
    _valueTest++;
}

/**
 * Test usage of global variables in a header file. And observe
 * the addresses of these variables.
 *
 * note: global variables are restricted to their own files.
 */
void _globalVariableInHeaderFile(int arrayTest[]) {
    _arrayTest = arrayTest;
    printf("_arrayTest: %p, arrayTest: %p\n", _arrayTest, arrayTest);
    printf("&_arrayTest[1]: %p, &arrayTest[1]: %p\n", &_arrayTest[1], &arrayTest[1]);
}

#endif
