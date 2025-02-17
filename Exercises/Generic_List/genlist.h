#ifndef HAVE_GENLIST_H
#define HAVE_GENLIST_H
#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#define SELF (*this)
#define FOR(i) for(int i=0;i<Size;i++)

template < typename T > class genlist {

private: int Size; T *data; /* C-array */
public : const int& size=Size;

T& operator[](int i){
	assert(i >= 0 && i < Size);
	return data[i];
	}

genlist() { /* default ctor */
	data = nullptr;
	Size = 0;
	}

genlist(genlist& other){ /* copy ctor */
	Size=other.size;
	data = new T[Size]();
	FOR(i) SELF[i]=other[i];
	}

genlist(genlist&&) = delete; /* move ctor */

~genlist() { /* dtor */
	delete[] data;
	Size = 0;
	}

genlist& operator=(genlist& other){ /* copy assignment */
	SELF.~genlist();
	Size=other.size;
	data = new T[Size]();
	FOR(i) SELF[i]=other[i];
	return SELF;
	}

genlist& operator=(genlist&&) = delete; /* move assignment */

void add(T item) {
	T *newdata = new T[Size+1]();
	FOR(i)newdata[i]=data[i];
	newdata[Size] = item;
	Size++;
	delete[] data;
	data=newdata;
	}

};
#endif
