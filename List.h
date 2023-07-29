#pragma once
// C++ program to implement Custom Linked List and 
// Iterator pattern. 
#pragma once
#include <iostream> 
using namespace std;

// Custom class to handle Linked List operations 
// Operations like push_back, push_front, pop_back, 
// pop_front, erase, size can be added here 



template <class T>
class List {
public:
	struct Node {
	public:
		T data;
		Node *pNext, *pPrev;
		double score;
	};

	class Iterator {
	public:
		Iterator() {}
		Iterator(Node* node) {
			currentNode = node;
		}
		Iterator& operator++ () {
			currentNode = currentNode->pNext;
			return *this;
		}
		Iterator& operator++ (int) {
			Iterator it = *this;
			++*this;
			return it;
		}
		Iterator& operator-- () {
			currentNode = currentNode->pPrev;
			return *this;
		}
		Iterator& operator-- (int) {
			Iterator it = *this;
			--*this;
			return it;
		}
		bool operator == (const Iterator &rhs) {
			return currentNode == rhs.currentNode;
		}

		bool operator != (const Iterator &rhs) {
			return currentNode != rhs.currentNode;
		}

		T operator*() {
			return currentNode->data;
		}

		double get_score() {
			return currentNode->score;
		}
		Node *currentNode;
	};

	List() noexcept {
		// caution: static members can't be 
		// initialized by initializer list 
		head = new Node();
		tail = new Node();
		head->pNext = tail;
		tail->pPrev = head;
	}

	bool empty() {
		return head->pNext == tail;
	}

	Iterator begin() {
		return Iterator(head->pNext);
	}

	Iterator end() {
		return Iterator(tail);
	}

	void clear() {
		if (head->pNext == tail) {
			return;
		}
		Node* now = head->pNext;
		Node* pNext = now->pNext;
		while (pNext != tail) {
			delete now;
			now = pNext;
			pNext = pNext->pNext;
		}
		head->pNext = tail;
		tail->pPrev = head;
	}

	void print() {/*
				  Node *cur = head->pNext;
				  while (cur != tail) {
				  cout << cur->data << endl;
				  cur = cur->pNext;
				  }*/
		for (Iterator it = begin(); it != end(); it++) {
			cout << *it << ' ';
		}
		cout << endl;
	}

	Iterator insert_after(Iterator it, T data) {
		Node* pCur = createNode(data);
		Node* pNext = it.currentNode->pNext;
		it.currentNode->pNext = pCur;
		pCur->pPrev = it.currentNode;
		pCur->pNext = pNext;
		pNext->pPrev = pCur;
		pCur->score = (pCur->pPrev->score + pNext->score) / 2;
		return Iterator(pCur);
	}

	Iterator insert_before(Iterator it, T data,bool flag,double score) {
		Node* pPrev = it.currentNode->pPrev;
		Node* pCur = createNode(data);
		pPrev->pNext = pCur;
		pCur->pPrev = pPrev;
		pCur->pNext = it.currentNode;
		it.currentNode->pPrev = pCur;
		Node *pNext = pCur->pNext;

		if (flag == false)
		{
			pCur->score = pCur->pNext->score;
			pCur->pNext->score++;
		}
		else pCur->score = pNext->score - score;

		//cout <<"score"<< pCur->score << ' ' << pNext->score << endl;

		return Iterator(pCur);
	}

	void push_back(T data) {
		Node* pCur = createNode(data);
		tail->pPrev->pNext = pCur;
		pCur->pPrev = tail->pPrev;
		pCur->pNext = tail;
		tail->pPrev = pCur;
		pCur->score = tail->score;
		tail->score++;
	}

	void remove(Iterator &it) {
		it--;
		Node *pPrev = it.currentNode;
		Node *pCur = pPrev->pNext;
		Node *pNext = pCur->pNext;
		pPrev->pNext = pNext;
		pNext->pPrev = pPrev;
		delete pCur;
	}

	// Create a new Node 
	Node* createNode(T data) {
		Node* pNewNode = new Node;
		pNewNode->data = data;
		pNewNode->pNext = nullptr;
		pNewNode->pPrev = nullptr;
		return pNewNode;
	}


	Node *head, *tail;
};





