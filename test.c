#include <stdio.h>

#ifdef BUILTIN
   #define K 1
#else
     #define K 2
#endif


int main(){
    printf("K = %d\n", K);
    return 0;
}