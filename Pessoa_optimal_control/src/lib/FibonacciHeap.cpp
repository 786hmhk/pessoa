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

#include "FibonacciHeap.hh"


#include <cstdlib>
#include <cmath>
#include <assert.h>
#include <cstdio>
#include "FibonacciHeap.hh"


/**
 * Creates an FibonacciHeap object capable of holding up to <strong class="paramname">n</strong> items.
 * @param n the maximum number of nodes/item
 */
FibonacciHeap::FibonacciHeap(int n) {
	int i;
	maxTrees = 1 + (int) (1.44 * log(n) / log(2.0));
	maxNodes = n;

	trees = new FHeapNode*[maxTrees];
	for (i = 0; i < maxTrees; i++)
		trees[i] = 0;

	nodes = new FHeapNode*[n];
	for (i = 0; i < n; i++)
		nodes[i] = 0;

	itemCount = 0;

	// The treeSum of the heap helps to keep track of the maximum rank while
	// nodes are inserted or deleted.
	treeSum = 0;

	// For experimental purposes, we keep a count of the number of key
	// comparisons.
	compCount = 0;
}

/**
 *
 */
FibonacciHeap::~FibonacciHeap() {
	int i;

	for (i = 0; i < maxNodes; i++)
		delete nodes[i];
	delete[] nodes;
	delete[] trees;
}

/**
 * Adding a new key to the heap
 *
 * @param item the integer the charecterizes the key
 * @param key the value of the key
 */
void FibonacciHeap::insert(int item, long key) {
	FHeapNode *newNode;

	/* create an initialise the new node */
	newNode = new FHeapNode;
	newNode->child = NULL;
	newNode->left = newNode->right = newNode;
	newNode->rank = 0;
	newNode->item = item;
	newNode->key = key;

	assert((maxNodes>item) && "Item value cannot be larger than the maximum number of nodes in the heap... Revise!");

	/* maintain a pointer to $item$'s new node in the heap */
	nodes[item] = newNode;

	/* meld the new node into the heap */
	meld(newNode);

	/* update the heaps node count */
	itemCount++;
}

/**
 * Returns the node/item of minimum value from a min heap
 * after removing it from the heap
 *
 * @return the node/key with the minimum value
 */
int FibonacciHeap::extractMin() {
	FHeapNode *minNode, *child, *next;
	long k, k2;
	int r, v, item;

	/* First we determine the maximum rank in the heap. */
	v = treeSum;
	r = -1;
	while (v) {
		v = v >> 1;
		r++;
	};

	/* Now determine which root node is the minimum. */
	minNode = trees[r];
	k = minNode->key;
	while (r > 0) {
		r--;
		next = trees[r];
		if (next) {
			if ((k2 = next->key) < k) {
				k = k2;
				minNode = next;
			}
			compCount++;
		}
	}

	/* We remove the minimum node from the heap but keep a pointer to it. */
	r = minNode->rank;
	trees[r] = NULL;
	treeSum -= (1 << r);

	child = minNode->child;
	if (child)
		meld(child);

	/* Record the vertex no of the old minimum node before deleting it. */
	item = minNode->item;
	nodes[item] = NULL;
	delete minNode;
	itemCount--;

	return item;
}

/**
 * Updating a key within min-heap
 *
 * @param item is the integer the charecterizes the key
 * @param newKey is the new value of the key
 */
void FibonacciHeap::decreaseKey(int item, long newKey) {
	FHeapNode *cutNode, *parent, *newRoots, *r, *l;
	int prevRank;

	/* Obtain a pointer to the decreased node and its parent then decrease the
	 * nodes key.
	 */
	cutNode = nodes[item];

	// I could implement a ifNotFound procedure, but I think leaving it like this
	// and making sure the client will not decrease a key that is not present, will
	// make the whole process of decreaseKey alot faster, since we do not have to
	// search for existing key...
	assert((cutNode!=NULL) && "Item must be in the heap... Revise!");

	parent = cutNode->parent;
	cutNode->key = newKey;

	/* No reinsertion occurs if the node changed was a root. */
	if (!parent) {
		return;
	}

	/* Update the left and right pointers of cutNode and its two neighbouring
	 * nodes.
	 */
	l = cutNode->left;
	r = cutNode->right;
	l->right = r;
	r->left = l;
	cutNode->left = cutNode->right = cutNode;

	/* Initially the list of new roots contains only one node. */
	newRoots = cutNode;

	/* While there is a parent node that is marked a cascading cut occurs. */
	while (parent && parent->marked) {

		/* Decrease the rank of cutNode's parent and update its child pointer.
		 */
		parent->rank--;
		if (parent->rank) {
			if (parent->child == cutNode)
				parent->child = r;
		} else {
			parent->child = NULL;
		}

		/* Update the cutNode and parent pointers to the parent. */
		cutNode = parent;
		parent = cutNode->parent;

		/* Update the left and right pointers of cutNodes two neighbouring
		 * nodes.
		 */
		l = cutNode->left;
		r = cutNode->right;
		l->right = r;
		r->left = l;

		/* Add cutNode to the list of nodes to be reinserted as new roots. */
		l = newRoots->left;
		newRoots->left = l->right = cutNode;
		cutNode->left = l;
		cutNode->right = newRoots;
		newRoots = cutNode;
	}

	/* If the root node is being relocated then update the trees[] array.
	 * Otherwise mark the parent of the last node cut.
	 */
	if (!parent) {
		prevRank = cutNode->rank + 1;
		trees[prevRank] = NULL;
		treeSum -= (1 << prevRank);
	} else {
		/* Decrease the rank of cutNode's parent an update its child pointer.
		 */
		parent->rank--;
		if (parent->rank) {
			if (parent->child == cutNode)
				parent->child = r;
		} else {
			parent->child = NULL;
		}

		parent->marked = 1;
	}

	/* Meld the new roots into the heap. */
	meld(newRoots);
}




/**
 * melds the linked list of trees pointed to by $treeList$ into the heap.
 */
void FibonacciHeap::meld(FHeapNode *treeList) {
	FHeapNode *first, *next, *nodePtr, *newRoot, *temp, *temp2, *lc, *rc;
	int r;

	/* We meld each tree in the circularly linked list back into the root level
	 * of the heap.  Each node in the linked list is the root node of a tree.
	 * The circularly linked list uses the sibling pointers of nodes.  This
	 *  makes melding of the child nodes from a deleteMin operation simple.
	 */
	nodePtr = first = treeList;

	do {

		/* Keep a pointer to the next node and remove sibling and parent links
		 * from the current node.  nodePtr points to the current node.
		 */
		next = nodePtr->right;
		nodePtr->right = nodePtr->left = nodePtr;
		nodePtr->parent = NULL;

		/* We merge the current node, nodePtr, by inserting it into the
		 * root level of the heap.
		 */
		newRoot = nodePtr;
		r = nodePtr->rank;

		/* This loop inserts the new root into the heap, possibly restructuring
		 * the heap to ensure that only one tree for each degree exists.
		 */
		do {

			/* Check if there is already a tree of degree r in the heap.
			 * If there is then we need to link it with newRoot so it will be
			 * reinserted into a new place in the heap.
			 */
			if ((temp = trees[r])) {

				/* temp will be linked to newRoot and relocated so we no
				 * longer will have a tree of degree r.
				 */
				trees[r] = NULL;
				treeSum -= (1 << r);

				/* Swap temp and newRoot if necessary so that newRoot always
				 * points to the root node which has the smaller key of the
				 * two.
				 */
				if (temp->key < newRoot->key) {
					temp2 = newRoot;
					newRoot = temp;
					temp = temp2;
				}
				compCount++;

				/* Link temp with newRoot, making sure that sibling pointers
				 * get updated if rank is greater than 0.  Also, increase r for
				 * the next pass through the loop since the rank of new has
				 * increased.
				 */
				if (r++ > 0) {
					rc = newRoot->child;
					lc = rc->left;
					temp->left = lc;
					temp->right = rc;
					lc->right = rc->left = temp;
				}
				newRoot->child = temp;
				newRoot->rank = r;
				temp->parent = newRoot;
				temp->marked = 0;
			}
			/* Otherwise if there is not a tree of degree r in the heap we
			 * allow newRoot, which possibly carries moved trees in the heap,
			 * to be a tree of degree r in the heap.
			 */
			else {

				trees[r] = newRoot;
				treeSum += (1 << r);
				;

				/* NOTE:  Because newRoot is now a root we ensure it is
				 *        marked.
				 */
				newRoot->marked = 1;
			}

			/* Note that temp will be NULL if and only if there was not a tree
			 * of degree r.
			 */
		} while (temp);

		nodePtr = next;

	} while (nodePtr != first);
}

