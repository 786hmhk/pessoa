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
 * This the Fibonacci implementation of the min-heap.
 *
 * 		No details yet.
 */

#ifndef FIBONACCIHEAP_HH_
#define FIBONACCIHEAP_HH_

/* Defines the base class for heaps. */
#include "heap.hh"



/**
 * Fibonacci heap node class.
 *
 * A nodes has the following pointers:
 * - parent       a pointer to the nodes parent node (if any).
 * - child        a pointer to a child node (typically the highest rank child).
 * - left, right  sibling pointers which provide a circular doubly linked list
 *                containing all the parents nodes children.
 *
 * The remaining fields are:
 * - rank         the nodes rank, that is, the number of children it has.
 * - key          the nodes key.
 * - item         the number of the item that the node is associated with.
 */
class FHeapNode {
  public:
    FHeapNode *parent;
    FHeapNode *left, *right;
    FHeapNode *child;
    int rank;
    int marked;
    long key;
    int item;
};

/**
 * Fibonacci heap class.
 *
 * - trees     An array of pointers to trees at root level in the heap.  Entry i
 *             in the array points to the root node of a tree that has nodes of
 *             dimension i on the main trunk.
 * - nodes     An array of pointers to nodes in the heap.  Nodes are indexed
 *             according to their vertex number.  This array can then be used to
 *             look up the node for corresponding to a vertex number, and is
 *             useful when freeing space taken up by the heap.
 * - maxNodes  The maximum number of nodes allowed in the heap.
 * - maxTrees  The maximum number of trees allowed in the heap (calculated from
 *             maxNodes).
 * - itemCount The current number of nodes in the heap.
 * - treeSum   The binary value represented by trees in the heap.
 *             By maintaining this it is easy to keep track of the maximum rank
 *             tree in the heap.
 * - compCount Can be used for experimental purposes when counting the number
 *             of key comparisons.
 */
class FibonacciHeap: public Heap {
  public:

	/**
	 * Creates an FibonacciHeap object capable of holding up to <strong class="paramname">n</strong> items.
	 * @param n the maximum number of nodes/item
	 */
    FibonacciHeap(int n);
    ~FibonacciHeap();

    /**
     * Returns the node/item of minimum value from a min heap
     * after removing it from the heap
     *
     * @return the node/key with the minimum value
     */
    int extractMin();

    /**
     * Adding a new key to the heap
     *
     * @param item the integer the charecterizes the key
     * @param key the value of the key
     */
    void insert(int item, long key);

    /**
     * Updating a key within min-heap
     *
     * @param item is the integer the charecterizes the key
     * @param newKey is the new value of the key
     */
    void decreaseKey(int item, long newKey);

    /**
     * Return the number of items in the heap.
     *
     * @return number of items in the heap
     */
    int size() const { return itemCount; }

    /**
     * Checks if the heap is empty
     *
     * @return returns true if the heap is empty, false otherwise
     */
    bool isEmpty() const { return itemCount == 0; }


    /**
     * Return the number of comparisons made to achieve the final
     * heap structure. (Experimental)
     *
     * @return the number of comparisons
     */
    long nComps() const { return compCount; }

  private:
    FHeapNode **trees;
    FHeapNode **nodes;
    int maxNodes, maxTrees, itemCount, treeSum;
    long compCount;

    void meld(FHeapNode *treeList);
};

#endif /* FIBONACCIHEAP_HH_ */
