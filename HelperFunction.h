#ifndef HELPERFUNCTION_H_INCLUDED
#define HELPERFUNCTION_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * Test funcs.
 */
void _myStrLengthTest();
void _lowerCaseTest();

/*
 * Important funcs.
 */
int myStrLength(char str[]);
char lowerCase(char ch);


void _myStrLengthTest() {
    printf("section\\\\ myStringLength test\n");
    char str[] = {'1', '2', '3', '\0'};
    char *ch = &str[0];
    int strlength = 0;
    int i = 0;
    for(i = 0; i < 3; i++) {
        printf("%c ", str[i]);
    }

    printf("\n******\n");

    while(*ch != '\0') {
        printf("%c ", *ch);
        strlength++;
        ch = ch + 1;
    }
}

void _lowerCaseTest() {
    printf("%c\t%c\t%c\n", lowerCase('A'), lowerCase('a'), lowerCase('D'));

}

/**
 * Get the length of a string. (including the character '\0')
 *
 * @param str input string
 */
int myStrLength(char str[]) {
    char *ch = &str[0];
    int strlength = 0;

    while(*ch != '\0') {
        strlength++;
        ch = ch + 1;
    }

    return strlength;
}

/**
 * Get the lower case of a letter.
 *
 * @param ch a letter
 * @return lower case of ch; ch if it's not a letter or it is already lower case
 */
char lowerCase(char ch) {
    if(ch >= 'A' && ch <= 'Z') {
        return ch - ('A' - 'a');
    } else {
        return ch;
    }
}

#endif // HELPERFUNCTION_H_INCLUDED
