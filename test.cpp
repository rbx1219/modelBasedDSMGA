
#include <iostream>
#include "myrand.h"
#include "statistics.h"
#include "onedarray.h"

using namespace std;

#define N 10000000

MyRand myRand;

int main()
{

    Statistics st1, st2;

    for (int i=0; i<N; i++) {

        OneDArray<double> a(5);
        for (int j=0; j<5; j++)
            a[j] = myRand.normal();

        double mean1=0.0, mean2=0.0;

        for (int j=0; j<5; j++)
            mean1 += a[j];

        for (int j=1; j<=3; j++)
            mean2 += a.getOS(j);

        st1.record(mean1/5);
        st2.record(mean2/3);


    }
    
    cout << st1.getMean() << ", " << st1.getVariance() << endl;
    cout << st2.getMean() << ", " << st2.getVariance() << endl;

    return 0;
}
