// list.cpp: implementation of the circular list class.
// this file should not be included in a project, but should be physically present.
// Only header file is used.
//////////////////////////////////////////////////////////////////////

#include "list.h"


//double linked list

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
template <class ListType>
list<ListType>::list()
{		
	CurrentItem = NULL;
	list_size = ZERO;
}

template <class ListType>
list<ListType>::list(const list<ListType> &L)
{
	CurrentItem = NULL;
	list_size = ZERO;
	
	while (list_size < L.Size())
		Insert(L.Prev());
}

template <class ListType>
list<ListType>::~list()
{
	if (list_size == ZERO)	return;

	listItem<ListType> *T;
		
	//this deleting should be faster than Dump()	
	while (list_size > 1)
	{
		list_size--;
		T = CurrentItem;
		CurrentItem		= CurrentItem->fw_ptr;
		delete T;	//DROP_MEM_BLOCK(T);
	}
	
	delete CurrentItem;	//DROP_MEM_BLOCK(CurrentItem);	//last item	
}

template <class ListType>
ListType& list<ListType>::Curr()
//definition:	returns current item
//should not be called if list is empty!
{	
	return (CurrentItem ->index);
}

template <class ListType>
ListType& list<ListType>::Next()
//shifts CurrentItem forward  & returns that item
//should not be called if list is empty!
{
	CurrentItem = CurrentItem ->fw_ptr;
	return (CurrentItem ->index);
}

template <class ListType>
ListType& list<ListType>::Prev()
//shifts CurrentITem backward & returns that item
//should not be called if list is empty!
{
	CurrentItem = CurrentItem ->bk_ptr;
	return (CurrentItem ->index);
}

template <class ListType>
void	list<ListType>::Insert(ListType &Item)
//Inserts an item in the position just before of CurrentItem
//Inserted element becomes CurrentItem
{
	listItem<ListType> *NewItem;
	NewItem = new listItem<ListType>; //GRAB_MEM_BLOCK(listItem<ListType>);

	NewItem->index  = Item;	
	
	if (CurrentItem == NULL)
	{
		CurrentItem	= NewItem;
		NewItem->bk_ptr	= NewItem;	
		NewItem->fw_ptr	= NewItem;
	}
	else
	{		
		NewItem->bk_ptr			= CurrentItem->bk_ptr;
		CurrentItem->bk_ptr		= NewItem;
		NewItem->bk_ptr->fw_ptr = NewItem;
		NewItem->fw_ptr			= CurrentItem;
		CurrentItem				= NewItem;
	};

	list_size++;
}

template <class ListType>
void	list<ListType>::Delete()
//Deletes  an item in the position of CurrentItem
//next element in a list becomes CurrentItem
{
	listItem<ListType> *T;

	if (CurrentItem != NULL)
	{
		list_size--;

		if (CurrentItem->fw_ptr == CurrentItem)
		{//if only one element in a list, then just delete it
			delete CurrentItem;	//DROP_MEM_BLOCK(CurrentItem);
			CurrentItem = NULL;
		}
		else
		{
			
			CurrentItem->fw_ptr->bk_ptr		= CurrentItem ->bk_ptr;
			CurrentItem->bk_ptr->fw_ptr		= CurrentItem ->fw_ptr;
			T								= CurrentItem ->fw_ptr;
			delete CurrentItem;	//DROP_MEM_BLOCK(CurrentItem);
			CurrentItem						= T;
			
		};
	};
}

template <class ListType>
unsigned short list<ListType>::Size()
{
	return (list_size);
}

template <class ListType>
list<ListType> & list<ListType>::operator =(list<ListType> &L)
{	
	Dump();

	while (list_size < L.Size())
		Insert(L.Prev());

	return (*this);
}

template <class ListType>
listItem<ListType> * list<ListType>::GetCurItem()
{//procedure that returns pointer to the current handl, thus providing free browsing of the list
	return (CurrentItem);
}

template <class ListType>
void list<ListType>::Dump()
{
	while (list_size != ZERO)
		Delete();
};




////single-linked list

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
template <class ListType>
slist<ListType>::slist()
{		
	CurrentItem = NULL;
	list_size = ZERO;
}

template <class ListType>
slist<ListType>::slist(const slist<ListType> &L)
{//note, reverses order in the circular list!
	
	list_size = ZERO;
	CurrentItem = NULL;

	while (list_size < L.Size())
		Insert(L.Next());
}

template <class ListType>
slist<ListType>::~slist()
{
	if (list_size == ZERO)	return;

	slistItem<ListType> *T;
		
	//this deleting should be faster than Dump()	
	while (list_size > 1)
	{
		list_size--;
		T = CurrentItem;
		CurrentItem		= CurrentItem->ptr;
		delete T;	//DROP_MEM_BLOCK(T);
	}
	
	delete CurrentItem;	//DROP_MEM_BLOCK(CurrentItem);	//last item	
}

template <class ListType>
ListType& slist<ListType>::Curr()
//definition:	returns current item
//should not be called if list is empty!
{	
	return (CurrentItem ->index);
}

template <class ListType>
ListType& slist<ListType>::Next()
//shifts CurrentItem forward  & returns that item
//should not be called if list is empty!
{
	CurrentItem = CurrentItem ->ptr;
	return (CurrentItem ->index);
}

template <class ListType>
void	slist<ListType>::Insert(ListType &Item)
//Inserts an item in the position just before of CurrentItem
//Inserted element becomes CurrentItem
{
	slistItem<ListType> *NewItem;
	NewItem = new slistItem<ListType>;//GRAB_MEM_BLOCK(slistItem<ListType>);
	
	if (CurrentItem == NULL)
	{
		NewItem->index  = Item;
		CurrentItem		= NewItem;
		NewItem->ptr	= NewItem;			
	}
	else
	{		
		NewItem->ptr			= CurrentItem->ptr;
		NewItem->index			= CurrentItem->index;
		CurrentItem->index		= Item;
		CurrentItem->ptr		= NewItem;
	};

	list_size++;
}

template <class ListType>
void	slist<ListType>::Delete()
//Deletes  an item in the position of CurrentItem
//next element in a list becomes CurrentItem
{
	slistItem<ListType> *T;

	if (CurrentItem != NULL)
	{
		list_size--;

		if (CurrentItem->ptr == CurrentItem)
		{//if only one element in a list, then just delete it
			delete CurrentItem;	//DROP_MEM_BLOCK(CurrentItem);
			CurrentItem = NULL;
		}
		else
		{
			T								= CurrentItem ->ptr;
			CurrentItem->index				= T->index;
			CurrentItem->ptr				= T->ptr->ptr;			
			delete T;	//DROP_MEM_BLOCK(T);			
		};
	};
}

template <class ListType>
unsigned short slist<ListType>::Size()
{
	return (list_size);
}

template <class ListType>
slist<ListType> & slist<ListType>::operator =(slist<ListType> &L)
{	//reverses order in the circular list!
	Dump();

	while (list_size < L.Size())
		Insert(L.Next());

	return (*this);
}

template <class ListType>
slistItem<ListType> * slist<ListType>::GetCurItem()
{//procedure that returns pointer to the current handl, thus providing free browsing of the list
	return (CurrentItem);
}

template <class ListType>
void slist<ListType>::Dump()
{
	while (list_size != ZERO)
		Delete();
}
