
#ifndef _ONE_D_ARRAY_H_
#define _ONE_D_ARRAY_H_


#include <cstdlib>
#include <cstring>

template <class T>
class OneDArray {

public:
    OneDArray() {
        size = 0;
        data = NULL;
    }

    OneDArray(int nSize) {
        size = nSize;
        data = new T[size];
    }
    
    ~OneDArray() {
        if (data!=NULL) 
            delete []data;
    }

    T& operator[] (int index) {
        return data[index];
    }

    T getOS(int index) {
        T* temp = new T[size];
        memcpy(temp, data, size*sizeof(T));

        T result = getOS(index, temp, size);
        delete []temp;
        return result;
    }

    T getOS(int index, T* temp, int size) const {

        int x = rand() % size;
        int cut = partition(x, temp, size);
        if (index < cut) 
            return getOS(index, temp, cut);
        else if (index > cut) 
            return getOS(index-cut-1, temp+cut+1, size-cut-1);
        else
            return temp[cut];

    }
    
    void swap(T& a, T& b) const {
        T temp = a;
        a = b;
        b = temp;
    }

    int partition(int x, T* temp, int size) const {
        int left = 0;
        int right = size-1;

        swap(temp[right], temp[x]);
        int cut = left;        

        for (int i=left; i<right; i++) {
            if (temp[i] < temp[size-1]) {
                swap(temp[i], temp[cut]);
                cut++;
            }
        }

        swap(temp[right], temp[cut]);
        return cut;
    }

private:
    T* data;
    int size;

};


#endif
