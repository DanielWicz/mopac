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
      nncf, iorbs, cocc, icocc
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
    integer :: i, j, j1, ja, jj, k, k1, k2, ka, kk, l, loop, nj
    integer :: max_threads
    logical :: use_parallel_density
    double precision :: spinfa, sum
#ifdef _OPENMP
    double precision, allocatable, save :: p_thread(:, :)
    integer, save :: p_thread_mpack = 0
    integer, save :: p_thread_nthreads = 0
    integer :: tid, nthreads_used
#endif
    integer, external :: ijbo

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
      !
      !  INITIALIZE P:  FULL DENSITY CALCULATION
      !
      if (use_parallel_density) then
#ifdef _OPENMP
!$omp parallel do schedule(static)
        do l = 1, mpack
          p(l) = 0.d0
        end do
!$omp end parallel do
#endif
      else
        p(:) = 0.d0
      end if
    else if (mode ==-1) then
      !
      !   FULL DENSITY MATRIX IS IN PARTP, REMOVE DENSITY DUE TO
      !   OLD LMO's USED IN SCF
      !
      if (use_parallel_density) then
#ifdef _OPENMP
!$omp parallel do schedule(static)
        do l = 1, mpack
          p(l) = -0.5d0 * partp(l)
        end do
!$omp end parallel do
#endif
      else
        p(:) = -0.5d0 * partp(:)
      end if
    else
      !
      !   PARTIAL DENSITY MATRIX IS IN PARTP, BUILD THE REST OF P
      !
      if (use_parallel_density) then
#ifdef _OPENMP
!$omp parallel do schedule(static)
        do l = 1, mpack
          p(l) = 0.5d0 * partp(l)
        end do
!$omp end parallel do
#endif
      else
        p(:) = 0.5d0 * partp(:)
      end if
    end if

#ifdef _OPENMP
    if (use_parallel_density) then
      if ((.not. allocated(p_thread)) .or. p_thread_mpack /= mpack .or. p_thread_nthreads < max_threads) then
        if (allocated(p_thread)) deallocate (p_thread)
        allocate (p_thread(mpack, max_threads))
        p_thread_mpack = mpack
        p_thread_nthreads = max_threads
      end if

      nthreads_used = 1
!$omp parallel default(shared) private(i, loop, ja, jj, j, nj, ka, kk, k, l, j1, k1, k2, sum, tid)
      tid = omp_get_thread_num() + 1
!$omp single
      nthreads_used = omp_get_num_threads()
!$omp end single
      p_thread(:, tid) = 0.d0
!$omp do schedule(static, 1)
      do i = 1, nclose_loc
        loop = ncocc(i)
        ja = 0
        if (lijbo) then
          do jj = nncf(i) + 1, nncf(i) + ncf(i)
            j = icocc(jj)
            nj = iorbs(j)
            ka = loop
            do kk = nncf(i) + 1, nncf(i) + ncf(i)
              k = icocc(kk)
              if (j == k) then
                l = nijbo (j, k)
                do j1 = 1, nj
                  sum = cocc(ja+j1+loop)
                  do k1 = 1, j1
                    k2 = ka + k1
                    l = l + 1
                    p_thread(l, tid) = p_thread(l, tid) + cocc(k2) * sum
                  end do
                end do
              else if (j > k .and. nijbo (j, k) >= 0) then
                l = nijbo (j, k)
                do j1 = 1, nj
                  sum = cocc(ja+j1+loop)
                  do k1 = 1, iorbs(k)
                    k2 = ka + k1
                    l = l + 1
                    p_thread(l, tid) = p_thread(l, tid) + cocc(k2) * sum
                  end do
                end do
              end if
              ka = ka + iorbs(k)
            end do
            ja = ja + nj
          end do
        else
          do jj = nncf(i) + 1, nncf(i) + ncf(i)
            j = icocc(jj)
            nj = iorbs(j)
            ka = loop
            do kk = nncf(i) + 1, nncf(i) + ncf(i)
              k = icocc(kk)
              l = ijbo (j, k)
              if (j == k) then
                do j1 = 1, nj
                  sum = cocc(ja+j1+loop)
                  do k1 = 1, j1
                    k2 = ka + k1
                    l = l + 1
                    p_thread(l, tid) = p_thread(l, tid) + cocc(k2) * sum
                  end do
                end do
              else if (j > k .and. l >= 0) then
                do j1 = 1, nj
                  sum = cocc(ja+j1+loop)
                  do k1 = 1, iorbs(k)
                    k2 = ka + k1
                    l = l + 1
                    p_thread(l, tid) = p_thread(l, tid) + cocc(k2) * sum
                  end do
                end do
              end if
              ka = ka + iorbs(k)
            end do
            ja = ja + nj
          end do
        end if
      end do
!$omp end do
!$omp end parallel

!$omp parallel do schedule(static) private(tid, sum)
      do l = 1, mpack
        sum = 0.d0
        do tid = 1, nthreads_used
          sum = sum + p_thread(l, tid)
        end do
        p(l) = p(l) + sum
      end do
!$omp end parallel do
    else
#endif
      do i = 1, nclose_loc
        loop = ncocc(i)
        ja = 0
        if (lijbo) then
          do jj = nncf(i) + 1, nncf(i) + ncf(i)
            j = icocc(jj)
            nj = iorbs(j)
            ka = loop
            do kk = nncf(i) + 1, nncf(i) + ncf(i)
              k = icocc(kk)
              if (j == k) then
                l = nijbo (j, k)
                do j1 = 1, nj
                  sum = cocc(ja+j1+loop)
                  do k1 = 1, j1
                    k2 = ka + k1
                    l = l + 1
                    p(l) = p(l) + cocc(k2) * sum
                  end do
                end do
              else if (j > k .and. nijbo (j, k) >= 0) then
                l = nijbo (j, k)
                do j1 = 1, nj
                  sum = cocc(ja+j1+loop)
                  do k1 = 1, iorbs(k)
                    k2 = ka + k1
                    l = l + 1
                    p(l) = p(l) + cocc(k2) * sum
                  end do
                end do
              end if
              ka = ka + iorbs(k)
            end do
            ja = ja + nj
          end do
        else
          do jj = nncf(i) + 1, nncf(i) + ncf(i)
            j = icocc(jj)
            nj = iorbs(j)
            ka = loop
            do kk = nncf(i) + 1, nncf(i) + ncf(i)
              k = icocc(kk)
              l = ijbo (j, k)
              if (j == k) then
                do j1 = 1, nj
                  sum = cocc(ja+j1+loop)
                  do k1 = 1, j1
                    k2 = ka + k1
                    l = l + 1
                    p(l) = p(l) + cocc(k2) * sum
                  end do
                end do
              else if (j > k .and. l >= 0) then
                do j1 = 1, nj
                  sum = cocc(ja+j1+loop)
                  do k1 = 1, iorbs(k)
                    k2 = ka + k1
                    l = l + 1
                    p(l) = p(l) + cocc(k2) * sum
                  end do
                end do
              end if
              ka = ka + iorbs(k)
            end do
            ja = ja + nj
          end do
        end if
      end do
#ifdef _OPENMP
    end if
#endif
    if (mode == 0 .or. mode == 1) then
      !
      !    FULL DENSITY CALCULATION.  MULTIPLY BY 2 FOR SPIN
      !
      spinfa = 2.d0
    else if (mode ==-1) then
      !
      !   MAKING PARTIAL DENSITY MATRIX. REVERSE SIGN ONCE MORE
      !
      spinfa = -2.d0
    else
      spinfa = 1.d0
    end if
   !
    if (Abs(spinfa - 1.d0) > 1.d-10) then
#ifdef _OPENMP
      if (use_parallel_density) then
!$omp parallel do schedule(static)
        do i = 1, mpack
          p(i) = spinfa * p(i)
        end do
!$omp end parallel do
      else
        p(:) = spinfa * p(:)
      end if
#else
      p(:) = spinfa * p(:)
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
