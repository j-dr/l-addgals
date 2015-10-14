#include "Logger.h"
#include "Logger.h"

#define MAXLENGTH 200

#ifndef LIST_H
#define LIST_H

template <class T>
class ListNode
{

private:
	char* name;
	T* element;
	ListNode* next;
	bool flag;
	
public:
	ListNode(T* element, char* name, ListNode<T>* next);
	char* getName();
	T* getElement();
	bool getFlag();
	void setFlag();
	void setNext(ListNode<T>* next);
	void setElement(T* element);
	ListNode<T>* getNext();
	void deleteName();
	void deleteElement();
};

template <class T>
class List
{
private:
	ListNode<T> * head;
	ListNode<T> * tail;
	ListNode<T> * iterator;

public:
	List();
	bool insert(char* name, T* element);
	T* get(char* name);
	void deleteUnused();
	void emptyList(bool toDeleteName, bool toDeleteElement);
	void printList();
	void startIterating();
	char* getNextName();
	bool fillTable(List<T>* mainParameterTable);
	static int caselessCompare(char* string1, char* string2);
	static char lowerChar(char ch);
};




template <class T> ListNode<T>::ListNode(T* element, char* name, ListNode<T>* next)
{
	this->name = name;
	this->element = element;
	this->next = next;
	flag = false;
}

template <class T> char* ListNode<T>::getName()
{
	return name;
}

template <class T> T* ListNode<T>::getElement()
{
	return element;
}

template <class T> bool ListNode<T>::getFlag()
{
	return flag;
}

template <class T> void ListNode<T>::setFlag()
{
	this->flag = true;
}

template <class T> void ListNode<T>::setNext(ListNode<T>* next)
{
	this->next = next;
}

template <class T> ListNode<T>* ListNode<T>::getNext()
{
	return next;
}

template <class T> void ListNode<T>::setElement(T* element)
{
	this->element = element;
}

template <class T> void ListNode<T>::deleteName()
{
	delete name;
}

template <class T> void ListNode<T>::deleteElement()
{
	delete element;
}



template <class T> List<T>::List()
{
	head = NULL;
	tail = NULL;
	iterator = NULL;
}

template <class T> bool List<T>::insert(char* name, T* element)
{
	if(name == NULL)
	{
		Logger::print(3, "Null Sam Name");
		return 0;
	}
	else if(element == NULL)
	{
		Logger::print(3, "Null t object with name %s", name);
		return 0;
	}
	else
	{
		//char* nameToInsert = new char[MAXLENGTH];
		//memcpy(nameToInsert, name, MAXLENGTH);
		ListNode<T>* newListNode = new ListNode<T>(element, name, NULL);
		if(newListNode == NULL)
		{
			Logger::print(3, "Memory could not be created for newListNode");
			return 0;
		}
		if(head == NULL)
		{
			head = tail = newListNode;
			return 1;
		}
		else
		{
			tail->setNext(newListNode);
			tail = newListNode;
			return 1;
		}
	}
}

template <class T> T* List<T>::get(char* name)
{
	if(name==NULL)
	{
		Logger::print(3, "List::get Input name is null");
		return NULL;
	}
	ListNode<T>* ptr = head;
	while( ptr != NULL)
	{
		if( caselessCompare(ptr->getName(), name) == 0)
		{
			ptr->setFlag();
			return ptr->getElement();
		}
		else
		{
			//printf("\n %s %s", toLower(ptr->getName()), toLower(name));
			ptr = ptr->getNext();
		}
	}
	Logger::print(2, "List::get name=%s not found", name);
	return NULL;
}

template <class T> void List<T>::deleteUnused()
{
	ListNode<T>* pptr = NULL;
	ListNode<T>* ptr = head;
	while( ptr != NULL)
	{
		if( ptr->getFlag() == false)
		{
			if(pptr == NULL)
			{
				head = ptr->getNext();
				//delete ptr;
				ptr = head;
			}
			else
			{
				pptr->setNext(ptr->getNext());
				//delete ptr;
				ptr = pptr->getNext();
			}
		}
		else
		{
			pptr = ptr;
			ptr = ptr->getNext();
		}
	}
}


template <class T> void List<T>::emptyList(bool toDeleteName, bool toDeleteElement)
{
	ListNode<T>* pptr = NULL;
	ListNode<T>* ptr = head;
	while( ptr != NULL)
	{
		pptr = ptr->getNext();
		if(toDeleteName)
			ptr->deleteName();
		if(toDeleteElement)
			ptr->deleteElement();
		delete ptr;
		ptr = pptr;
	}
}

template <class T> void List<T>::printList()
{
	ListNode<T>* ptr = head;
	while( ptr != NULL)
	{
		Logger::print(1,ptr->getName());
		ptr = ptr->getNext();
	}
}

template <class T> void List<T>::startIterating()
{
	iterator = head;
}
	
template <class T> char* List<T>::getNextName()
{
	if(iterator != NULL)
	{
		char* name = iterator->getName();
		iterator = iterator->getNext();
		return name;
	}
	else
	{
		return NULL;
	}
}

template <class T> bool List<T>::fillTable(List<T>* mainParameterTable)
{
	bool areAllFilled = true;
	ListNode<T>* ptr  = head;
	while( ptr != NULL)
	{
		T* element = mainParameterTable->get(ptr->getName());
		if(element!=NULL)
		{
			ptr->setElement(element);
			ptr->setFlag();
		}
		else
		{
			//printf("%s",ptr->getName());
			ptr->setElement(element);
			areAllFilled = false;
		}
		ptr = ptr->getNext();
	}
	return areAllFilled;
}

template <class T> int List<T>::caselessCompare(char* string1, char* string2)
{
	int length1 = strlen(string1);
	int length2 = strlen(string2);
	if(length1 != length2)
		return length1 - length2;
	for(int i=0;i<length1;i++)
	{
		if(lowerChar(string1[i]) != lowerChar(string2[i]))
		{
			return (string1[i] - string2[i]);
		}
	}
	return 0;
}

template <class T> char List<T>::lowerChar(char ch)
{
	if ( ( ch >= 'A' && ch <= 'Z' ) )
	{
		return ch + 'a' - 'A';
	}
	else
		return ch;
}


#endif

