#ifndef HELPERFUNCTION_H_INCLUDED
#define HELPERFUNCTION_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>

void _myStrLengthTest(){
    printf("section\\\\ myStringLength test\n");
    char str[] = {'1','2','3','\0'};
    char *ch = &str[0];
    int strlength = 0;
    for(int i = 0; i < 3; i++){
        printf("%c ", str[i]);
    }

    printf("\n******\n");

    while(*ch!='\0'){
        printf("%c ", *ch);
        strlength++;
        ch = ch + 1;
    }
}

/*
 Get the length of a string. (including the character '\0')
 */
int myStrLength(char str[]){
    char *ch = &str[0];
    int strlength = 0;

    while(*ch!='\0'){
        strlength++;
        ch = ch + 1;
    }

    return strlength;
}

#endif // HELPERFUNCTION_H_INCLUDED
