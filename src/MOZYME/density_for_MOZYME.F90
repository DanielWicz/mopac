! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

subroutine density_for_MOZYME (p, mode, nclose_loc, partp)
   !***********************************************************************
   !
   !   DENSIT COMPUTES THE DENSITY MATRIX GIVEN THE EIGENVECTOR MATRIX, AND
   !          INFORMATION ABOUT THE M.O. OCCUPANCY.
   !
   !  INPUT:  COCC  = COMPRESSED OCCUPIED EIGENVECTOR MATRIX
   !
   !   ON EXIT: P   = DENSITY MATRIX
   !
   !***********************************************************************
    use molkst_C, only: numat, mpack, keywrd
    use MOZYME_C, only : lijbo, nijbo, ncf, ncocc, &
      nncf, iorbs, cocc, icocc, icocc_dim
    use chanel_C, only: iw
#ifdef _OPENMP
    use omp_lib, only: omp_get_max_threads, omp_get_num_threads, omp_get_thread_num
#endif
    implicit none
    integer, intent (in) :: mode, nclose_loc
    double precision, dimension (mpack), intent (in) :: partp
    double precision, dimension (mpack), intent (inout) :: p
    logical :: first = .true.
    logical, save :: prnt
    integer :: max_threads
    logical :: use_parallel_density
    double precision :: spinfa, sum
    integer, external :: ijbo

    integer, allocatable, save :: atom_lmo_ptr(:), atom_lmo_idx(:), atom_lmo_off(:)
    integer, allocatable :: atom_degree(:)
    integer :: pos, i_lmo, jj_idx, j_atom, off_j, ptr, l_idx
    integer :: i, j, j1, k, k1, k2, ka, kk, l

    if (first) then
      first = .false.
      prnt = Index (keywrd, " DIAG ") /= 0
    end if

    use_parallel_density = .false.
#ifdef _OPENMP
    max_threads = omp_get_max_threads()
    use_parallel_density = (max_threads > 1)
#endif

    if (mode == 0) then
#ifdef _OPENMP
      if (use_parallel_density) then
!$omp parallel do schedule(static)
        do l = 1, mpack
          p(l) = 0.d0
        end do
!$omp end parallel do
      else
#endif
        p(:) = 0.d0
#ifdef _OPENMP
      end if
#endif
    else if (mode ==-1) then
#ifdef _OPENMP
      if (use_parallel_density) then
!$omp parallel do schedule(static)
        do l = 1, mpack
          p(l) = -0.5d0 * partp(l)
        end do
!$omp end parallel do
      else
#endif
        p(:) = -0.5d0 * partp(:)
#ifdef _OPENMP
      end if
#endif
    else
#ifdef _OPENMP
      if (use_parallel_density) then
!$omp parallel do schedule(static)
        do l = 1, mpack
          p(l) = 0.5d0 * partp(l)
        end do
!$omp end parallel do
      else
#endif
        p(:) = 0.5d0 * partp(:)
#ifdef _OPENMP
      end if
#endif
    end if

    ! INSPECTOR PHASE
    if (.not. allocated(atom_lmo_ptr) .or. size(atom_lmo_ptr) < numat + 2) then
        if (allocated(atom_lmo_ptr)) deallocate(atom_lmo_ptr, atom_lmo_idx, atom_lmo_off)
        allocate(atom_lmo_ptr(numat + 2))
        allocate(atom_lmo_idx(icocc_dim), atom_lmo_off(icocc_dim))
    end if

    allocate(atom_degree(numat))
    atom_degree = 0
    do i_lmo = 1, nclose_loc
        do jj_idx = nncf(i_lmo) + 1, nncf(i_lmo) + ncf(i_lmo)
            j_atom = icocc(jj_idx)
            atom_degree(j_atom) = atom_degree(j_atom) + 1
        end do
    end do

    atom_lmo_ptr(1) = 1
    do j_atom = 1, numat
        atom_lmo_ptr(j_atom + 1) = atom_lmo_ptr(j_atom) + atom_degree(j_atom)
    end do
    
    atom_degree(1:numat) = atom_lmo_ptr(1:numat)

    do i_lmo = 1, nclose_loc
        off_j = ncocc(i_lmo)
        do jj_idx = nncf(i_lmo) + 1, nncf(i_lmo) + ncf(i_lmo)
            j_atom = icocc(jj_idx)
            pos = atom_degree(j_atom)
            atom_lmo_idx(pos) = i_lmo
            atom_lmo_off(pos) = off_j
            atom_degree(j_atom) = pos + 1
            off_j = off_j + iorbs(j_atom)
        end do
    end do
    deallocate(atom_degree)

    ! EXECUTOR PHASE
#ifdef _OPENMP
    if (use_parallel_density) then
!$omp parallel do schedule(dynamic, 4) private(j, ptr, i, off_j, ka, kk, k, l, l_idx, j1, k1, k2, sum)
      do j = 1, numat
        do ptr = atom_lmo_ptr(j), atom_lmo_ptr(j+1) - 1
            i = atom_lmo_idx(ptr)
            off_j = atom_lmo_off(ptr)

            ka = ncocc(i)
            do kk = nncf(i) + 1, nncf(i) + ncf(i)
                k = icocc(kk)
                if (j == k) then
                    if (lijbo) then
                        l = nijbo(j, k)
                    else
                        l = ijbo(j, k)
                    end if
                    l_idx = l
                    do j1 = 1, iorbs(j)
                        sum = cocc(off_j + j1)
                        do k1 = 1, j1
                            k2 = ka + k1
                            l_idx = l_idx + 1
                            p(l_idx) = p(l_idx) + cocc(k2) * sum
                        end do
                    end do
                else if (j > k) then
                    if (lijbo) then
                        l = nijbo(j, k)
                    else
                        l = ijbo(j, k)
                    end if
                    if (l >= 0) then
                        l_idx = l
                        do j1 = 1, iorbs(j)
                            sum = cocc(off_j + j1)
                            do k1 = 1, iorbs(k)
                                k2 = ka + k1
                                l_idx = l_idx + 1
                                p(l_idx) = p(l_idx) + cocc(k2) * sum
                            end do
                        end do
                    end if
                end if
                ka = ka + iorbs(k)
            end do
        end do
      end do
!$omp end parallel do
    else
#endif
      do j = 1, numat
        do ptr = atom_lmo_ptr(j), atom_lmo_ptr(j+1) - 1
            i = atom_lmo_idx(ptr)
            off_j = atom_lmo_off(ptr)

            ka = ncocc(i)
            do kk = nncf(i) + 1, nncf(i) + ncf(i)
                k = icocc(kk)
                if (j == k) then
                    if (lijbo) then
                        l = nijbo(j, k)
                    else
                        l = ijbo(j, k)
                    end if
                    l_idx = l
                    do j1 = 1, iorbs(j)
                        sum = cocc(off_j + j1)
                        do k1 = 1, j1
                            k2 = ka + k1
                            l_idx = l_idx + 1
                            p(l_idx) = p(l_idx) + cocc(k2) * sum
                        end do
                    end do
                else if (j > k) then
                    if (lijbo) then
                        l = nijbo(j, k)
                    else
                        l = ijbo(j, k)
                    end if
                    if (l >= 0) then
                        l_idx = l
                        do j1 = 1, iorbs(j)
                            sum = cocc(off_j + j1)
                            do k1 = 1, iorbs(k)
                                k2 = ka + k1
                                l_idx = l_idx + 1
                                p(l_idx) = p(l_idx) + cocc(k2) * sum
                            end do
                        end do
                    end if
                end if
                ka = ka + iorbs(k)
            end do
        end do
      end do
#ifdef _OPENMP
    end if
#endif

    if (mode == 0 .or. mode == 1) then
      spinfa = 2.d0
    else if (mode ==-1) then
      spinfa = -2.d0
    else
      spinfa = 1.d0
    end if
   
    if (Abs(spinfa - 1.d0) > 1.d-10) then
#ifdef _OPENMP
      if (use_parallel_density) then
!$omp parallel do schedule(static)
        do i = 1, mpack
          p(i) = spinfa * p(i)
        end do
!$omp end parallel do
      else
#endif
        p(:) = spinfa * p(:)
#ifdef _OPENMP
      end if
#endif
    end if

    if (prnt) then
      sum = 0.d0
      if (lijbo) then
        do i = 1, numat
          j = nijbo(i, i) + 1
          if (iorbs(i) > 0) then
            sum = sum + p(j)
          end if
          if (iorbs(i) == 4) then
            sum = sum + p(j+2) + p(j+5) + p(j+9)
          end if
        end do
      else
        do i = 1, numat
          j = ijbo (i, i) + 1
          if (iorbs(i) > 0) then
            sum = sum + p(j)
          end if
          if (iorbs(i) == 4) then
            sum = sum + p(j+2) + p(j+5) + p(j+9)
          end if
        end do
      end if
      write (iw,*) " COMPUTED NUMBER OF ELECTRONS:", sum
    end if
end subroutine density_for_MOZYME
