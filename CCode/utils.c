#define _ISOC99_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void chomp(char *str) {
	int l = strlen(str) - 1;
	while ((str[l] == '\n') || (str[l] == ' ')) {
		str[l] = '\0';
		l--;
	}
}
