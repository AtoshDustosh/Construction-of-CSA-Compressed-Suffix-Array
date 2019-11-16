#ifndef _HEADFILETEST_H_
#define _HEADFILETEST_H_

#include <stdio.h>
#include <stdlib.h>
#include "HelperFunction.h"

/*
 * Review usage of C language.
 */
void _CLanguageReview(){
    printf("section\\\\ C language review test\n");
    // review string[] array construction
	printf("String operation\n");
    char *str[] = {"abb", "aab", "aaa", "bbdsas",
        "bsass", "brabrabra", "blahblah", "atosh"};
    int strLength = 8;

    // review string[] array access
	int i = 0;
	printf("String length: %d\n", strLength);
	printf("str[4]: %s\n", str[4]);
	for(i = 0; i < strLength; i++){
		printf("%s\n", str[i]);
	}
	printf("%s\n", str[1]);
	// review string[].getChar(i) -> integer operation
	char *strCh = str[1];
	int count = 0;
	while(*strCh != '\0'){
        printf("%d: %c%d\n", count, *strCh, *strCh);
        strCh = strCh+1;
        count ++;
	}
    printf("\n\n");
}

#endif
