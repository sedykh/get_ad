//list.h: interface for the class, 
//
//Implementation of double- and single-linked circular lists
//
//Each empty class instance takes up 12 bytes. 
//Each list element - 8bytes (4 bytes for single-linked) + size of ListType
//
//A.Sedykh, 2000-2002
//////////////////////////////////////////////////////////////////////

#if !defined(LIST_CLASS)
#define LIST_CLASS

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000



template <class ListType>

struct listItem
{
	ListType  index;
	listItem *fw_ptr;
	listItem *bk_ptr;
};

//double linked list
template <class ListType>

class list  
{
	public:
		list();
		list(const list<ListType> &);
		virtual ~list();		

		ListType& Curr();					//returns current item
		ListType& Next();					//shifts CurrentItem forward  & returns that item
		ListType& Prev();					//shifts CurrentITem backward & returns that item
		void	  Insert(ListType &);			//Inserts an item in the position of CurrentItem
		void	  Delete();					//Delets  an item in the position of CurrentItem
		unsigned short Size();			//returns the size of the list
		list<ListType> & operator =(list<ListType> &);
		listItem<ListType> * GetCurItem();
		void Dump();						//cleans up the list, sets it to nothing
		

	private:
		//listItem<ListType> * LIST;
		listItem<ListType> * CurrentItem;
		
		unsigned short  list_size;
};



template <class ListType>

struct slistItem
{
	ListType  index;
	slistItem *ptr;	
};

//single-linked list
template <class ListType>

class slist
{
	public:
		slist();
		slist(const slist<ListType> &);
		virtual ~slist();		

		ListType& Curr();					//returns current item
		ListType& Next();					//shifts CurrentItem forward  & returns that item

		void	  Insert(ListType &);			//Inserts an item in the position of CurrentItem
		void	  Delete();					//Delets  an item in the position of CurrentItem
		
		unsigned short Size();			//returns the size of the list
		slist<ListType> & operator =(slist<ListType> &);
		slistItem<ListType> * GetCurItem();
		void Dump();						//cleans up the list, sets it to nothing
		

	private:		
		slistItem<ListType> * CurrentItem;		
		unsigned short  list_size;
};



#include "list.cpp"
#endif // !defined(LIST_CLASS)
