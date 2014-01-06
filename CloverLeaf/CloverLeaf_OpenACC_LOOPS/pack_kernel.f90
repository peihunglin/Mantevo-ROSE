!Crown Copyright 2012 AWE.
!
! This file is part of CloverLeaf.
!
! CloverLeaf is free software: you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! CloverLeaf is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details.
!
! You should have received a copy of the GNU General Public License along with 
! CloverLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Fortran mpi buffer packing kernel
!>  @author Wayne Gaudin
!>  @details Packs/unpacks mpi send and receive buffers

MODULE pack_kernel_module

CONTAINS

SUBROUTINE pack_left_right_buffers(x_min,x_max,y_min,y_max,              &
                                   chunk_left,chunk_right,external_face, &
                                   x_inc,y_inc,depth,size,               &
                                   field,left_snd_buffer,right_snd_buffer)

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max
  INTEGER      :: chunk_left,chunk_right,external_face
  INTEGER      :: x_inc,y_inc,depth,size

  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
  REAL(KIND=8) :: left_snd_buffer(:),right_snd_buffer(:)

  INTEGER      :: j,k,index

  IF(chunk_left.NE.external_face) THEN
!$ACC DATA &
!$ACC PRESENT(left_snd_buffer,field)
!$ACC PARALLEL LOOP PRIVATE(index) ASYNC(1)
    DO k=y_min-depth,y_max+y_inc+depth
      DO j=1,depth
        index=j+(k+depth-1)*depth
        left_snd_buffer(index)=field(x_min+x_inc-1+j,k)
      ENDDO
    ENDDO
!$ACC END PARALLEL LOOP
!$ACC UPDATE HOST (left_snd_buffer(1:size)) ASYNC(1)
!$ACC END DATA
  ENDIF
  IF(chunk_right.NE.external_face) THEN
!$ACC DATA &
!$ACC PRESENT(right_snd_buffer,field)
!$ACC PARALLEL LOOP PRIVATE(index) ASYNC(2)
    DO k=y_min-depth,y_max+y_inc+depth
      DO j=1,depth
        index=j+(k+depth-1)*depth
        right_snd_buffer(index)=field(x_max+1-j,k)
      ENDDO
    ENDDO
!$ACC END PARALLEL LOOP
!$ACC UPDATE HOST (right_snd_buffer(1:size)) ASYNC(2)
!$ACC END DATA
  ENDIF
!$ACC WAIT

END SUBROUTINE pack_left_right_buffers

SUBROUTINE unpack_left_right_buffers(x_min,x_max,y_min,y_max,              &
                                     chunk_left,chunk_right,external_face, &
                                     x_inc,y_inc,depth,size,               &
                                     field,left_rcv_buffer,right_rcv_buffer)

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max
  INTEGER      :: chunk_left,chunk_right,external_face
  INTEGER      :: x_inc,y_inc,depth,size

  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
  REAL(KIND=8) :: left_rcv_buffer(:),right_rcv_buffer(:)

  INTEGER      :: j,k,index

  IF(chunk_left.NE.external_face) THEN
!$ACC DATA &
!$ACC PRESENT(left_rcv_buffer,field)
!$ACC UPDATE DEVICE (left_rcv_buffer(1:size)) ASYNC(3)
!$ACC PARALLEL LOOP PRIVATE(index) ASYNC(3)
    DO k=y_min-depth,y_max+y_inc+depth
      DO j=1,depth
        index=j+(k+depth-1)*depth
        field(x_min-j,k)=left_rcv_buffer(index)
      ENDDO
    ENDDO
!$ACC END PARALLEL LOOP
!$ACC END DATA
  ENDIF
  IF(chunk_right.NE.external_face) THEN
!$ACC DATA &
!$ACC PRESENT(right_rcv_buffer,field)
!$ACC UPDATE DEVICE (right_rcv_buffer(1:size)) ASYNC(4)
!$ACC PARALLEL LOOP PRIVATE(index) ASYNC(4)
    DO k=y_min-depth,y_max+y_inc+depth
      DO j=1,depth
        index=j+(k+depth-1)*depth
        field(x_max+x_inc+j,k)=right_rcv_buffer(index)
      ENDDO
    ENDDO
!$ACC END PARALLEL LOOP
!$ACC END DATA
  ENDIF
!$ACC WAIT

END SUBROUTINE unpack_left_right_buffers

SUBROUTINE pack_top_bottom_buffers(x_min,x_max,y_min,y_max,              &
                                   chunk_bottom,chunk_top,external_face, &
                                   x_inc,y_inc,depth,size,               &
                                   field,bottom_snd_buffer,top_snd_buffer)

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max
  INTEGER      :: chunk_bottom,chunk_top,external_face
  INTEGER      :: x_inc,y_inc,depth,size

  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
  REAL(KIND=8) :: bottom_snd_buffer(:),top_snd_buffer(:)

  INTEGER      :: j,k,index

  IF(chunk_bottom.NE.external_face) THEN
!$ACC DATA &
!$ACC PRESENT(bottom_snd_buffer,field)
!$ACC PARALLEL LOOP PRIVATE(index) ASYNC(5)
    DO j=x_min-depth,x_max+x_inc+depth
      DO k=1,depth
        index=j+depth+(k-1)*(x_max+x_inc+(2*depth))
        bottom_snd_buffer(index)=field(j,y_min+y_inc-1+k)
      ENDDO
    ENDDO
!$ACC END PARALLEL LOOP
!$ACC UPDATE HOST (bottom_snd_buffer(1:size)) ASYNC(5)
!$ACC END DATA
  ENDIF
  IF(chunk_top.NE.external_face) THEN
!$ACC DATA &
!$ACC PRESENT(top_snd_buffer,field)
!$ACC PARALLEL LOOP PRIVATE(index) ASYNC(6)
    DO j=x_min-depth,x_max+x_inc+depth
      DO k=1,depth
        index=j+depth+(k-1)*(x_max+x_inc+(2*depth))
        top_snd_buffer(index)=field(j,y_max+1-k)
      ENDDO
    ENDDO
!$ACC END PARALLEL LOOP
!$ACC UPDATE HOST (top_snd_buffer(1:size)) ASYNC(6)
!$ACC END DATA
  ENDIF
!$ACC WAIT

END SUBROUTINE pack_top_bottom_buffers

SUBROUTINE unpack_top_bottom_buffers(x_min,x_max,y_min,y_max,             &
                                    chunk_bottom,chunk_top,external_face, &
                                    x_inc,y_inc,depth,size,               &
                                    field,bottom_rcv_buffer,top_rcv_buffer)

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max
  INTEGER      :: chunk_bottom,chunk_top,external_face
  INTEGER      :: x_inc,y_inc,depth,size

  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
  REAL(KIND=8) :: bottom_rcv_buffer(:),top_rcv_buffer(:)

  INTEGER      :: j,k,index

  IF(chunk_bottom.NE.external_face) THEN
!$ACC DATA &
!$ACC PRESENT(bottom_rcv_buffer,field)
!$ACC UPDATE DEVICE (bottom_rcv_buffer(1:size)) ASYNC(7)
!$ACC PARALLEL LOOP PRIVATE(index) ASYNC(7)
    DO j=x_min-depth,x_max+x_inc+depth
      DO k=1,depth
        index=j+depth+(k-1)*(x_max+x_inc+(2*depth))
        field(j,y_min-k)=bottom_rcv_buffer(index)
      ENDDO
    ENDDO
!$ACC END PARALLEL LOOP
!$ACC END DATA
  ENDIF
  IF(chunk_top.NE.external_face) THEN
!$ACC DATA &
!$ACC PRESENT(top_rcv_buffer,field)
!$ACC UPDATE DEVICE (top_rcv_buffer(1:size)) ASYNC(8)
!$ACC PARALLEL LOOP PRIVATE(index) ASYNC(8)
    DO j=x_min-depth,x_max+x_inc+depth
      DO k=1,depth
        index=j+depth+(k-1)*(x_max+x_inc+(2*depth))
        field(j,y_max+y_inc+k)=top_rcv_buffer(index)
      ENDDO
    ENDDO
!$ACC END PARALLEL LOOP
!$ACC END DATA
  ENDIF
!$ACC WAIT

END SUBROUTINE unpack_top_bottom_buffers

END MODULE pack_kernel_module
