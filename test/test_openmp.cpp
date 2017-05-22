//
// Created by 季俊哲 on 21/05/2017.
//

#include <iostream>
#include <ctime>
#include <zconf.h>

using namespace std;

int main()
{
    long a[1000000];
    clock_t t_start = clock();
    for (long  i = 0; i < 100000; i++) {
        a[i] = i;
    }
    clock_t t_end1 = clock();
    for (long  i = 0; i < 100000; i++) {
        if (a[i] == i) {
            a[i] = 0;
        }

    }
    clock_t t_end2 = clock();

    cout << (double)((t_end1 - t_start)/(CLOCKS_PER_SEC)) * 10000000 << endl;
    cout << (double)((t_end2 - t_end1)/(CLOCKS_PER_SEC)) * 1000000 << endl;


    return 0;
}