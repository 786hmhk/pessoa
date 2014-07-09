/**
 * @file
 * @author Athanasios Tasoglou <tasoglou@gmail.com>
 * @version 0.92
 *
 * @section LICENSE
 *
 * Copyright (c) 2013, TU Delft: Delft University of Technology
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of TU Delft: Delft University of Technology nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL TU DELFT: DELFT UNIVERSITY OF TECHNOLOGY BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * @section DESCRIPTION
 *
 * This is an abstract base class for various min-heap implementations.
 *
 * 		No details yet.
 */

#ifndef HEAP_H
#define HEAP_H

/**
 * This is an abstract base class from which specific heap classes can be
 * derived.  Different heaps derived from this abstract base class can be used
 * interchangeably by algorithms that were written using the universal
 * interface it provides.
 *
 * This heap stores integer items, and associates with each item a long integer
 * key.  Any derived heap heap must provide the following methods:
 *
 * - extractMin()   removes the item with the minimum key from the heap, and
 *                  returns it.
 * - insert()       inserts an item 'item' with key 'key' into the heap.
 * - decreaseKey()  decreases the key of item 'item' to the new value newKey.
 * - size()         returns the number of items currently in the heap.
 * - isEmpty()      returns true if the heap is empty, false otherwise.
 * - nComps()       returns the number of key comparison operations. (Experimental)
 */
class Heap {
  public:
    virtual ~Heap() { };
    /**
     * Returns the node/item of minimum value from a min heap
     * after removing it from the heap
     *
     * @return the node/key with the minimum value
     */
    virtual int extractMin() = 0;
    /**
     * Adding a new key to the heap
     *
     * @param item the integer the charecterizes the key
     * @param key the value of the key
     */
    virtual void insert(int item, long key) = 0;
    /**
     * Updating a key within min-heap
     *
     * @param item is the integer the charecterizes the key
     * @param newKey is the new value of the key
     */
    virtual void decreaseKey(int item, long newKey) = 0;
    /**
     * Return the number of items in the heap.
     *
     * @return number of items in the heap
     */
    virtual int size() const = 0;
    /**
     * Checks if the heap is empty
     *
     * @return returns true if the heap is empty, false otherwise
     */
    virtual bool isEmpty() const = 0;
    /**
     * Return the number of comparisons made to achieve the final
     * heap structure. (Experimental)
     *
     * @return the number of comparisons
     */
    virtual long nComps() const = 0;
};

class HeapDesc {
  public:
    virtual ~HeapDesc() { };
    virtual Heap *newInstance(int n) const = 0;
};

template <class T>
class HeapD: public HeapDesc {
  public:
    Heap *newInstance(int n) const { return new T(n); };
};

#endif
