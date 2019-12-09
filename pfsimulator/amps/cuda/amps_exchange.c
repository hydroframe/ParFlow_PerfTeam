/*BHEADER*********************************************************************
 *
 *  Copyright (c) 1995-2009, Lawrence Livermore National Security,
 *  LLC. Produced at the Lawrence Livermore National Laboratory. Written
 *  by the Parflow Team (see the CONTRIBUTORS file)
 *  <parflow@lists.llnl.gov> CODE-OCEC-08-103. All rights reserved.
 *
 *  This file is part of Parflow. For details, see
 *  http://www.llnl.gov/casc/parflow
 *
 *  Please read the COPYRIGHT file or Our Notice and the LICENSE file
 *  for the GNU Lesser General Public License.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License (as published
 *  by the Free Software Foundation) version 2.1 dated February 1999.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms
 *  and conditions of the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 *  USA
 **********************************************************************EHEADER*/

#include "amps.h"
#include "amps_cuda.h"

#ifdef AMPS_MPI_NOT_USE_PERSISTENT

void _amps_wait_exchange(amps_Handle handle)
{
  int notdone;
  int i;

  MPI_Status *status;

  if (handle->package->num_recv + handle->package->num_send)
  {
    status = (MPI_Status*)calloc((handle->package->num_recv +
                                  handle->package->num_send), sizeof(MPI_Status));

    MPI_Waitall(handle->package->num_recv + handle->package->num_send,
                handle->package->requests,
                status);

    free(status);

    for (i = 0; i < handle->package->num_recv; i++)
    {
      if (handle->package->recv_invoices[i]->mpi_type != MPI_DATATYPE_NULL)
      {
        MPI_Type_free(&(handle->package->recv_invoices[i]->mpi_type));
      }

      MPI_Request_free(&handle->package->requests[i]);
    }

    for (i = 0; i < handle->package->num_send; i++)
    {
      if (handle->package->send_invoices[i]->mpi_type != MPI_DATATYPE_NULL)
      {
        MPI_Type_free(&handle->package->send_invoices[i]->mpi_type);
      }

      MPI_Request_free(&handle->package->requests[handle->package->num_recv + i]);
    }
  }
}

amps_Handle amps_IExchangePackage(amps_Package package)
{
  int i;

  /*--------------------------------------------------------------------
   * post receives for data to get
   *--------------------------------------------------------------------*/
  package->recv_remaining = 0;

  for (i = 0; i < package->num_recv; i++)
  {
    amps_create_mpi_type(MPI_COMM_WORLD, package->recv_invoices[i]);

    MPI_Type_commit(&(package->recv_invoices[i]->mpi_type));

    MPI_Irecv(MPI_BOTTOM, 1, package->recv_invoices[i]->mpi_type,
              package->src[i], 0, MPI_COMM_WORLD,
              &(package->requests[i]));
  }

  /*--------------------------------------------------------------------
   * send out the data we have
   *--------------------------------------------------------------------*/
  for (i = 0; i < package->num_send; i++)
  {
    amps_create_mpi_type(MPI_COMM_WORLD, package->send_invoices[i]);

    MPI_Type_commit(&(package->send_invoices[i]->mpi_type));

    MPI_Isend(MPI_BOTTOM, 1, package->send_invoices[i]->mpi_type,
              package->dest[i], 0, MPI_COMM_WORLD,
              &(package->requests[package->num_recv + i]));
  }

  return(amps_NewHandle(NULL, 0, NULL, package));
}

#else

int _amps_send_sizes(amps_Package package, int **sizes)
{
  int size_acc;

  *sizes = (int*)calloc(package->num_send, sizeof(int));

  size_acc = 0;
  for (int i = 0; i < package->num_send; i++)
  {
    (*sizes)[i] = amps_sizeof_invoice(amps_CommWorld, package->send_invoices[i]);
    size_acc += (*sizes)[i];
  }

  /* check to see if enough space is already allocated            */
  if (amps_combuf_send_size < size_acc)
  {
    if (amps_combuf_send_size != 0) CUDA_ERR(cudaFree(amps_combuf_send));

    CUDA_ERR(cudaMalloc((void**)&amps_combuf_send, size_acc));
    amps_combuf_send_size = size_acc;
  }

  size_acc = 0;
  for (int i = 0; i < package->num_send; i++)
  {
    package->send_invoices[i]->combuf = &amps_combuf_send[size_acc];
    amps_pack(amps_CommWorld, package->send_invoices[i], package->send_invoices[i]->combuf);
    
    MPI_Isend(&((*sizes)[i]), 1, MPI_INT, package->dest[i],
              0, amps_CommWorld,
              &(package->send_requests[i]));
    size_acc += (*sizes)[i];
  }

  return(0);
}
int _amps_recv_sizes(amps_Package package)
{
  int *sizes;
  int size_acc;

  MPI_Status status;
  sizes = (int*)calloc(package->num_recv, sizeof(int));

  size_acc = 0;
  for (int i = 0; i < package->num_recv; i++)
  {
    MPI_Recv(&sizes[i], 1, MPI_INT, package->src[i], 0,
          amps_CommWorld, &status);
    size_acc += sizes[i];
  }

  /* check to see if enough space is already allocated            */
  if (amps_combuf_recv_size < size_acc)
  {
    if (amps_combuf_recv_size != 0) CUDA_ERR(cudaFree(amps_combuf_recv));

    CUDA_ERR(cudaMalloc((void**)&(amps_combuf_recv), size_acc));
    CUDA_ERR(cudaMemset(amps_combuf_recv, 0, size_acc));
    amps_combuf_recv_size = size_acc;
  }

  size_acc = 0;
  for (int i = 0; i < package->num_recv; i++)
  {
    package->recv_invoices[i]->combuf = &amps_combuf_recv[size_acc];
    MPI_Recv_init(package->recv_invoices[i]->combuf, sizes[i],
                  MPI_BYTE, package->src[i], 1, amps_CommWorld,
                  &(package->recv_requests[i]));
    size_acc += sizes[i];
  }

  free(sizes);

  return(0);
}

void _amps_wait_exchange(amps_Handle handle)
{
  int i;
  int num;

  num = handle->package->num_send + handle->package->num_recv;

  MPI_Waitall(num, handle->package->recv_requests,
              handle->package->status);

  if (num)
  {
    if (handle->package->num_recv)
    {
      for (i = 0; i < handle->package->num_recv; i++)
      {
        amps_unpack(amps_CommWorld, handle->package->recv_invoices[i],
                    (char *)handle->package->recv_invoices[i]->combuf);
      }
    }
  }
}

amps_Handle amps_IExchangePackage(amps_Package package)
{

  int i;
  int *send_sizes;

  MPI_Status *status_array;

  if (package->num_send)
    _amps_send_sizes(package, &send_sizes);

  if (package->num_recv)
    _amps_recv_sizes(package);

  /*--------------------------------------------------------------------
   * end the sending of sizes
   *--------------------------------------------------------------------*/

  if (package->num_send)
  {
    status_array = (MPI_Status*)calloc(package->num_send,
                                       sizeof(MPI_Status));
    MPI_Waitall(package->num_send, package->send_requests, status_array);
    free(status_array);
  }

  if (package->num_send)
  {
    for (i = 0; i < package->num_send; i++)
    {
      MPI_Send_init(package->send_invoices[i]->combuf,
                    send_sizes[i], MPI_BYTE, package->dest[i], 1,
                    amps_CommWorld, &(package->send_requests[i]));
    }

    free(send_sizes);
  }

  /*--------------------------------------------------------------------
   * post receives for data to get
   *--------------------------------------------------------------------*/

  if (package->num_recv)
    MPI_Startall(package->num_recv, package->recv_requests);

  /*--------------------------------------------------------------------
   * send out the data we have
   *--------------------------------------------------------------------*/
  if (package->num_send)
    MPI_Startall(package->num_send, package->send_requests);

  return(amps_NewHandle(NULL, 0, NULL, package));
}

#endif
