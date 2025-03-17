#ifndef HAVE_GENLIST_H
#define HAVE_GENLIST_H

#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>

#define SELF (*this)
#define FOR(i) for(int i=0; i<Size; i++)

template <typename T>
class genlist {
private:
    int Size;
    T* data; /* C-array */

public:
    const int& size = Size;

    // Non-const version (for modification)
    T& operator[](int i) {
        assert(i >= 0 && i < Size);
        return data[i];
    }

    // Const version (for reading from const objects)
    const T& operator[](int i) const {
        assert(i >= 0 && i < Size);
        return data[i];
    }

    genlist() { /* default constructor */
        data = nullptr;
        Size = 0;
    }

    genlist(const genlist& other) { /* copy constructor */
        Size = other.size;
        data = new T[Size]();
        FOR(i) SELF[i] = other[i];  // Now works with const objects
    }

    genlist& operator=(const genlist& other) { /* copy assignment */
        if (this != &other) {
            delete[] data;
            Size = other.size;
            data = new T[Size]();
            FOR(i) SELF[i] = other[i];  // Now works with const objects
        }
        return *this;
    }

    genlist(genlist&&) = delete; /* move constructor */
    genlist& operator=(genlist&&) = delete; /* move assignment */

    ~genlist() { /* destructor */
        delete[] data;
    }

    void add(T item) {
        T* newdata = new T[Size + 1]();
        FOR(i) newdata[i] = data[i];
        newdata[Size] = item;
        Size++;
        delete[] data;
        data = newdata;
    }
};

#endif
