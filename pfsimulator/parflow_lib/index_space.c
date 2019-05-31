
#include "parflow.h"
#include "index_space.h"

#include <limits.h>

void IndexCopy(Index index, Index src)
{
  for(int dim = 0; dim < DIM; dim++)
  {
    index[dim] = src[dim];
  }  
}

int BoxSize(Box *box)
{
  return (box->up[0] - box->lo[0] + 1) * (box->up[1] - box->lo[1] + 1) * (box->up[2] - box->lo[2] + 1);
}

void BoxNumberCells(Box* box, Index* number_cells)
{
  for(int dim = 0; dim < DIM; dim++)
  {
    (*number_cells)[dim] = box->up[dim] - box->lo[dim] + 1;
  }
}

void BoxClear(Box *box)
{
  for(int dim = 0; dim < DIM; dim++)
  {
    box -> lo[dim] = INT_MAX;
    box -> up[dim] = INT_MAX;
  }
}

void BoxSet(Box *box, Index lo, Index up)
{
  for(int dim = 0; dim < DIM; dim++)
  {
    box -> lo[dim] = lo[dim];
    box -> up[dim] = up[dim];
  }
}

void BoxCopy(Box *dst, Box *src)
{
  for(int dim = 0; dim < DIM; dim++)
  {
    dst -> lo[dim] = src -> lo[dim];
    dst -> up[dim] = src -> up[dim];
  }
}

void BoxPrint(Box* box)
{
  printf("[%04d,%04d,%04d][%04d,%04d,%04d]", box->lo[0], box->lo[1], box->lo[2], box->up[0], box->up[1], box->up[2]);
}

BoxArray* NewBoxArray(BoxList *box_list)
{
   BoxArray* box_array = ctalloc(BoxArray,1);

   if(box_list)
   {
      box_array -> size = BoxListSize(box_list);
      box_array -> boxes = ctalloc(Box, box_array -> size);

      int i = 0;
      BoxListElement* element = box_list -> head;
      while(element)
      {
	 BoxCopy(&(box_array -> boxes[i++]), &(element -> box));
	 element = element -> next;
      }

   }
   
   return box_array;
}

void FreeBoxArray(BoxArray* box_array)
{
  if(box_array)
  {
     if(box_array -> boxes)
     {
	tfree(box_array -> boxes);
     }

     tfree(box_array);
  }
}

BoxList* NewBoxList(void)
{
  return ctalloc(BoxList, 1);
}

void FreeBoxList(BoxList *box_list)
{
  if(box_list)
  {
     BoxListElement* element = box_list -> head;
     while(element)
     {
	BoxListElement* next = element -> next;
	tfree(element);
	element = next;
     }
  
     tfree(box_list);
  }
}

int BoxListSize(BoxList *box_list)
{
  return box_list -> size;
}

int BoxListIsEmpty(BoxList *box_list)
{
  return box_list -> size > 0 ? 0 : 1;
}

Box* BoxListFront(BoxList *box_list)
{
  return &(box_list -> head -> box);
}

void BoxListAppend(BoxList* box_list, Box* box)
{
  if(box_list -> size == 0)
  {
    box_list -> head = ctalloc(BoxListElement, 1);
    BoxCopy(&(box_list -> head -> box), box);
    
    box_list -> tail = box_list -> head;
    box_list -> head -> next = NULL;
    box_list -> head -> prev = NULL;
    
    box_list -> size = 1;
  }
  else
  {
    BoxListElement* new_element = ctalloc(BoxListElement, 1);
    BoxCopy(&(new_element -> box), box);
    
    new_element -> next = NULL;
    new_element -> prev = box_list -> tail;
    box_list -> tail -> next = new_element;
    box_list -> tail = new_element;
    
    box_list -> size++;
  }
}

void BoxListConcatenate(BoxList *box_list, BoxList *concatenate_list)
{
  BoxListElement* element = concatenate_list -> head;
  while(element)
  {
    BoxListAppend(box_list, &(element -> box));
    element = element -> next;
  }
}

void BoxListClearItems(BoxList* box_list)
{
  BoxListElement* element = box_list -> head;
  while(element)
  {
    BoxListElement* next = element->next;
    tfree(element);
    element = next;
  }

  box_list -> head = NULL;
  box_list -> tail = NULL;
  box_list -> size = 0;
}

void BoxListPrint(BoxList* box_list)
{
  BoxListElement* element = box_list -> head;
  while(element)
  {
    printf("\t");
    BoxPrint(&(element -> box));
    printf("\n");
    element = element->next;
  }
}
