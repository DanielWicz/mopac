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

subroutine diagg1 (fao, nocc, nvir, eigv, ws, latoms, ifmo, fmo, fmo_dim, nij, idiagg,  avir, aocc, aov)
   !**********************************************************************
   !
   !  FAO:      FOCK MATRIX OVER ATOMIC ORBITALS
   !  FMO:      FOCK MATRIX OVER MOLECULAR ORBITALS
   !  VECTOR:   EIGENVECTORS OF M.O.S
   !  FILLED:   OCCUPIED M.O.S
   !  eigs:     EIGENVALUES OF OCCUPIED M.O.S
   !  EMPTY:    VIRTUAL M.O.S
   !  EIGV:     EIGENVALUES OF VIRTUAL M.O.S
   !  NOCC:     NUMBER OF OCCUPIED M.O.S
   !  EIG:      ALL THE EIGENVALUES
   !  N:        NUMBER OF M.O.S = NUMBER OF A.O.S*
   !
   !**********************************************************************
   !
   !  ALL OCCUPIED AND UNOCCUPIED LMO ENERGY LEVELS ARE CALCULATED, AS WELL
   !  AS ALL MATRIX ELEMENTS WHICH CAN BE NON ZERO CONNECTING THE OCCUPIED
   !  AND VIRTUAL SETS OF LMO'S.
   !
   !**********************************************************************
   !
	    use molkst_C, only: numat, norbs, mpack, numcal, keywrd
	    use MOZYME_C, only : nvirtual, icocc_dim, &
	       & nfmo, lijbo, nijbo, &
	       tiny, sumt, ijc, ovmax, ncf, nce, nncf, nnce, ncocc, ncvir, &
	     & iorbs, icocc, icvir, cocc, cvir
	    use common_arrays_C, only : eigs, nfirst, nlast, p
#ifdef _OPENMP
	    use omp_lib, only: omp_get_max_threads, omp_get_num_threads, omp_get_thread_num
#endif
	    implicit none
	    integer, intent (in) :: idiagg, nocc, nvir, fmo_dim
	    integer, intent (inout) :: nij
	    logical, dimension (numat), intent (out) :: latoms
    integer, dimension (2, fmo_dim), intent (inout) :: ifmo
    double precision, dimension (fmo_dim), intent (out) :: fmo
    double precision, dimension (numat), intent (out) :: aov
    double precision, dimension (norbs), intent (out) :: avir, ws
	    double precision, dimension (mpack), intent (in) :: fao
	    double precision, dimension (nvirtual), intent (inout) :: eigv
		    double precision, dimension (icocc_dim), intent (out) :: aocc
		    integer, save :: icalcn = 0
		    integer :: i, i1, i2, i4, ii, j, j1, j2, j4, prev_i, &
			         & ij0, jj, jl, jx, k, k1, kk, kl, l, loopi, loopj, kj, i1j1, &
			         & i1j2, i2j1, i2j2, k1j1
		    logical :: lij
		    logical, save :: times
		    double precision :: cutlim, flim
		    double precision, save :: fref, oldlim, safety
		    double precision :: cutoff, sum, sum1
#ifdef _OPENMP
	    type :: diagg1_thread_buffer
	      integer :: n = 0
	      integer, allocatable :: i(:), j(:)
	      double precision, allocatable :: f(:)
	    end type diagg1_thread_buffer
	    type(diagg1_thread_buffer), allocatable :: buffers(:)
	    integer, allocatable :: owner_tid(:), owner_pos(:)
	    double precision, allocatable :: sumt_thr(:), tiny_thr(:)

	    logical :: use_parallel_diagg1
	    integer :: nthreads_used, tid, nthreads, est, off, ncopy, pos, ijc_total, start_buf
	    integer :: max_threads, nij_max
	    double precision :: sumt_par, tiny_par, sumt_loc, tiny_loc
	    integer, allocatable :: start_idx(:)
	    double precision, allocatable :: ws_loc(:), avir_loc(:), aov_loc(:), fstore(:)
	    logical, allocatable :: latoms_loc(:)
	    integer, allocatable :: jstore(:)
#endif
	    integer, external :: ijbo
	    if (numcal /= icalcn) then
	      icalcn = numcal
	      fref = 10.0d0
	      safety = 1.0d0
	      oldlim = 0.0d0
	      times = (Index (keywrd, " TIMES") /= 0)
	      if (Index (keywrd, " OLDENS") /= 0) then
	        fref = 0.d0
	      end if
	    end if
    !
    !
    aocc(:) = 0.d0
    !
	    !   IF THE CONTRIBUTION OF AN ATOM IN AN OCCUPIED LMO IS VERY SMALL
	    !   THEN DO NOT USE THAT ATOM IN CALCULATING THE OCCUPIED-VIRTUAL
	    !   INTERACTION.  PUT THE CONTRIBUTIONS INTO AN ARRAY 'AOCC'.
	    !
!$omp parallel do default(shared) private(loopj,kl,kk,k1,sum,k) schedule(static)
	    do j = 1, nocc
	        loopj = ncocc(j)
	        kl = 0
	        do kk = nncf(j) + 1, nncf(j) + ncf(j)
	          k1 = icocc(kk)
          sum = 0.d0
          do k = nfirst(k1), nlast(k1)
            kl = kl + 1
            sum = sum + cocc(kl+loopj) ** 2
          end do
          !
          !   AOCC(KK) HOLDS THE SQUARE OF THE CONTRIBUTION OF THE KK'TH ATOM
          !   IN THE OCCUPIED SET.  NOTE:  THIS IS NOT ATOM KK.
          !
          aocc(kk) = sum
          !
	        end do
	      !
	    end do
!$omp end parallel do
	    !
	    !
	    !    CUTLIM    PRECISION OF PL
	    !
    !    1.D-6      0.004
    !    1.D-7      0.00004
    !
    cutlim = 1.d-8
    cutoff = Max (cutlim, tiny*10.d0*cutlim)
    flim = Min (3.d0, fref*0.5d0)
    fref = 0.d0
    if (idiagg <= 5) then
      cutoff = cutlim
    end if
    sumt = 0.d0
    ijc = 0
    tiny = 0.d0
	    !
	    !
#ifdef _OPENMP
	    use_parallel_diagg1 = .false.
	    max_threads = omp_get_max_threads()
	    use_parallel_diagg1 = (max_threads > 1 .and. nvir > 8)
	    if (use_parallel_diagg1) then
	      nij_max = nij
	      ijc = 0
	      sumt = 0.d0
	      tiny = 0.d0
		      if (idiagg <= 5 .or. Mod (idiagg, 2) == 0) then
		        nfmo(1:nvir) = 0
		        allocate (buffers(max_threads), owner_tid(nvir), owner_pos(nvir), &
		             & sumt_thr(max_threads), tiny_thr(max_threads), start_idx(nvir))
		        owner_tid(:) = 0
		        owner_pos(:) = 0
		        sumt_thr(:) = 0.d0
		        tiny_thr(:) = 0.d0
		        nthreads_used = 1
!$omp parallel default(shared) private(ws_loc, latoms_loc, avir_loc, aov_loc, jstore, fstore, sumt_loc, tiny_loc, &
!$omp& sumt_par, tiny_par, l, tid, nthreads, est, i, prev_i, start_buf)
			        tid = omp_get_thread_num() + 1
			        nthreads = omp_get_num_threads()
!$omp single
			        nthreads_used = nthreads
!$omp end single
			        est = Max (1024, (nij_max + nthreads - 1) / Max (1, nthreads))
			        call buffer_reserve (buffers(tid), est)

			        allocate (ws_loc(norbs), latoms_loc(numat), avir_loc(numat), aov_loc(numat), &
			             & jstore(nocc), fstore(nocc))
			        sumt_loc = 0.d0
			        tiny_loc = 0.d0
			        prev_i = 0

!$omp do schedule(dynamic,1)
				        ! Irregular-work balancing: each virtual i has variable cost W_i,
				        ! so dynamic scheduling approximates minimizing max_t(sum_{i in t} W_i).
				        do i = 1, nvir
			          call build_virtual (i, prev_i, ws_loc, latoms_loc, avir_loc, aov_loc, cutoff)
			          prev_i = i
			          call rebuild_couplings (i, ws_loc, latoms_loc, aov_loc, cutoff, flim, oldlim, jstore, fstore, &
			               & sumt_par, tiny_par, l)
			          nfmo(i) = l
			          sumt_loc = sumt_loc + sumt_par
		          tiny_loc = Max (tiny_loc, tiny_par)
		          if (l > 0) then
		            start_buf = buffers(tid)%n + 1
		            call buffer_append (buffers(tid), i, jstore, fstore, l)
		            owner_tid(i) = tid
		            owner_pos(i) = start_buf
		          end if
			        end do
!$omp end do

		        sumt_thr(tid) = sumt_loc
		        tiny_thr(tid) = tiny_loc

		        deallocate (ws_loc, latoms_loc, avir_loc, aov_loc, jstore, fstore)
!$omp end parallel

		        sumt = 0.d0
		        tiny = 0.d0
		        do tid = 1, nthreads_used
		          sumt = sumt + sumt_thr(tid)
		          tiny = Max (tiny, tiny_thr(tid))
		        end do

			        ! Eq. (1) Prefix packing offsets:
			        !   start_idx(1) = 1
			        !   start_idx(i) = 1 + sum_{k=1}^{i-1} nfmo(k), i>1
			        start_idx(1) = 1
			        do i = 2, nvir
			          start_idx(i) = start_idx(i-1) + nfmo(i-1)
			        end do
			        ! Eq. (2) Total couplings after rebuild:
			        !   ijc_total = sum_{k=1}^{nvir} nfmo(k)
			        ijc_total = start_idx(nvir) + nfmo(nvir) - 1

			        if (ijc_total <= nij_max) then
			          ijc = ijc_total
!$omp parallel do default(shared) private(i, l, tid, pos, off) schedule(static)
		          do i = 1, nvir
		            l = nfmo(i)
		            if (l <= 0) cycle
		            tid = owner_tid(i)
		            pos = owner_pos(i)
		            off = start_idx(i)
			            ! Deterministic mapping (virtual-major order):
			            ! off = start_idx(i), block length = nfmo(i)
			            ifmo(1, off:off+l-1) = i
		            ifmo(2, off:off+l-1) = buffers(tid)%j(pos:pos+l-1)
		            fmo(off:off+l-1) = buffers(tid)%f(pos:pos+l-1)
		          end do
!$omp end parallel do
		        else
		          ijc = 0
		          do i = 1, nvir
		            l = nfmo(i)
		            if (l <= 0) cycle
		            if (ijc == nij_max) then
		              nfmo(i:nvir) = 0
		              exit
		            end if
		            ncopy = Min (l, nij_max - ijc)
		            if (ncopy > 0) then
		              tid = owner_tid(i)
		              pos = owner_pos(i)
		              ifmo(1, ijc+1:ijc+ncopy) = i
		              ifmo(2, ijc+1:ijc+ncopy) = buffers(tid)%j(pos:pos+ncopy-1)
		              fmo(ijc+1:ijc+ncopy) = buffers(tid)%f(pos:pos+ncopy-1)
		              ijc = ijc + ncopy
		            end if
		            if (ncopy < l) then
		              nfmo(i) = ncopy
		              if (i < nvir) nfmo(i+1:nvir) = 0
		              exit
		            end if
		          end do
		        end if

		        deallocate (buffers, owner_tid, owner_pos, sumt_thr, tiny_thr, start_idx)
		      else
		        allocate (start_idx(nvir))
		        start_idx(1) = 1
		        do i = 2, nvir
		          start_idx(i) = start_idx(i-1) + nfmo(i-1)
		        end do
		        ijc = start_idx(nvir) + nfmo(nvir) - 1
		        ijc = Min (ijc, nij_max)
		        sumt_par = 0.d0
		        tiny_par = 0.d0
!$omp parallel default(shared) private(ws_loc, latoms_loc, avir_loc, aov_loc, l, sum, loopj, kl, jj, j, k, kk, k1, prev_i) &
!$omp& reduction(+:sumt_par) reduction(max:tiny_par)
		        allocate (ws_loc(norbs), latoms_loc(numat), avir_loc(numat), aov_loc(numat))
		        prev_i = 0
!$omp do schedule(dynamic,1)
		        do i = 1, nvir
		          call build_virtual (i, prev_i, ws_loc, latoms_loc, avir_loc, aov_loc, cutoff)
		          prev_i = i
		          l = nfmo(i)
		          if (l <= 0) cycle
			          do jj = 1, l
			            if (start_idx(i)+jj-1 > nij_max) exit
		            j = ifmo(2, start_idx(i)+jj-1)
		            loopj = ncocc(j)
			            sum = 0.d0
			            kl = 0
			            do kk = nncf(j) + 1, nncf(j) + ncf(j)
			              k1 = icocc(kk)
		              if (.not. latoms_loc(k1)) then
		                kl = kl + iorbs(k1)
		              else if (aov_loc(k1)*aocc(kk) < cutoff) then
		                kl = kl + iorbs(k1)
		              else
		                do k = nfirst(k1), nlast(k1)
		                  kl = kl + 1
	                  sum = sum + ws_loc(k) * cocc(kl+loopj)
	                end do
	              end if
	            end do
		            sumt_par = sumt_par + Abs (sum)
		            tiny_par = Max (tiny_par, Abs (sum))
		            fmo(start_idx(i)+jj-1) = sum
		          end do
		        end do
!$omp end do
		        deallocate (ws_loc, latoms_loc, avir_loc, aov_loc)
!$omp end parallel
	        sumt = sumt_par
	        tiny = tiny_par
	        deallocate (start_idx)
	      end if
	    else
#endif
	      do i = 1, nvir
	        !
	        loopi = ncvir(i)
	        latoms(:) = .false.
	        l = 0
	        do j = nnce(i) + 1, nnce(i) + nce(i)
	          j1 = icvir(j)
	          sum = 0.d0
	          do k = l + 1, l + iorbs(j1)
	            sum = sum + cvir(k+loopi) ** 2
	          end do
	          l = l + iorbs(j1)
	          !
	          !   AVIR(J1) HOLDS THE SQUARE OF THE CONTRIBUTION OF THE ATOM J1
	          !   IN THE VIRTUAL LMO 'I'.
	          !
	          avir(j1) = sum
	          latoms(icvir(j)) = .true.
	        end do
	        !
	        if (lijbo) then
	          do jj = nnce(i) + 1, nnce(i) + nce(i)
	            j1 = icvir(jj)
	            do jx = 1, iorbs(j1)
	              ws(nfirst(j1)+jx-1) = 0.0d00
	            end do
	            !
	            kl = loopi
	            do kk = nnce(i) + 1, nnce(i) + nce(i)
	              k1 = icvir(kk)
	              kj = nijbo(k1, j1)
	              if (kj >= 0) then
	                if (avir(k1)*p(kj+1) > cutoff) then
	                  !
	                  !  EXTRACT THE ATOM-ATOM INTERSECTION OF FAO
	                  !
	                  if (iorbs(k1) .eq. 1 .and. iorbs(j1) .eq. 1) then
	                    ws(nfirst(j1)) = ws(nfirst(j1)) + fao(kj+1) * cvir(kl+1)
	                  else
	                    if (k1 > j1) then
	                      ii = kj
	                      do i4 = 1, iorbs(k1)
	                        do jx = 1, iorbs(j1)
	                          ii = ii + 1
	                          ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) &
	                               & * cvir(kl+i4)
	                        end do
	                      end do
	                    else if (k1 < j1) then
	                      ii = kj
	                      do jx = 1, iorbs(j1)
	                        do i4 = 1, iorbs(k1)
	                          ii = ii + 1
	                          ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) &
	                               & * cvir(kl+i4)
	                        end do
	                      end do
	                    else
	                      do jx = 1, iorbs(j1)
	                        do i4 = 1, iorbs(j1)
	                          if (i4 > jx) then
	                            ii = kj + (i4*(i4-1)) / 2 + jx
	                          else
	                            ii = kj + (jx*(jx-1)) / 2 + i4
	                          end if
	                          ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) &
	                               & * cvir(kl+i4)
	                        end do
	                      end do
	                    end if
	                  end if
	                end if
	              end if
	              kl = kl + iorbs(k1)
	            end do
	          end do
	        else
	          do jj = nnce(i) + 1, nnce(i) + nce(i)
	            j1 = icvir(jj)
	            do jx = 1, iorbs(j1)
	              ws(nfirst(j1)+jx-1) = 0.0d00
	            end do
	            !
	            kl = loopi
	            do kk = nnce(i) + 1, nnce(i) + nce(i)
	              k1 = icvir(kk)
	              kj = ijbo (k1, j1)
	              if (kj >= 0) then
	                if (avir(k1)*p(kj+1) > cutoff) then
	                  !
	                  !  EXTRACT THE ATOM-ATOM INTERSECTION OF FAO
	                  !
	                  if (iorbs(k1) == 1 .and. iorbs(j1) == 1) then
	                    ws(nfirst(j1)) = ws(nfirst(j1)) + fao(kj+1) * cvir(kl+1)
	                  else
	                    if (k1 > j1) then
	                      ii = kj
	                      do i4 = 1, iorbs(k1)
	                        do jx = 1, iorbs(j1)
	                          ii = ii + 1
	                          ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) &
	                               & * cvir(kl+i4)
	                        end do
	                      end do
	                    else if (k1 < j1) then
	                      ii = kj
	                      do jx = 1, iorbs(j1)
	                        do i4 = 1, iorbs(k1)
	                          ii = ii + 1
	                          ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) &
	                               & * cvir(kl+i4)
	                        end do
	                      end do
	                    else
	                      do jx = 1, iorbs(j1)
	                        do i4 = 1, iorbs(j1)
	                          if (i4 > jx) then
	                            ii = kj + (i4*(i4-1)) / 2 + jx
	                          else
	                            ii = kj + (jx*(jx-1)) / 2 + i4
	                          end if
	                          ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) &
	                               & * cvir(kl+i4)
	                        end do
	                      end do
	                    end if
	                  end if
	                end if
	              end if
	              kl = kl + iorbs(k1)
	            end do
	          end do
	        end if
	        !
	        do j = 1, numat
	          if (latoms(j)) then
	            sum = 0.d0
	            do k = nfirst(j), nlast(j)
	              sum = sum + ws(k) ** 2
	            end do
	            !
	            !   AOV(J) HOLDS THE SQUARE OF THE ENERGY CONTRIBUTION OF THE
	            !   ATOM J IN THE VIRTUAL LMO 'I'.
	            !
	            aov(j) = sum
	          else
	            aov(j) = 0.d0
	          end if
	        end do
	        !
	        !   EVALUATE THE VIRTUAL ENERGY LEVELS
	        !
	        sum = 0.d0
	        kl = loopi
	        do kk = nnce(i) + 1, nnce(i) + nce(i)
	          k1 = icvir(kk)
	          if (aov(k1)*avir(k1) > cutoff) then
	            do k = nfirst(k1), nlast(k1)
	              kl = kl + 1
	              sum = sum + ws(k) * cvir(kl)
	            end do
	          else
	            kl = kl + iorbs(k1)
	          end if
	        end do
	        !
	        eigv(i) = sum
	        !
	        !  EVALUATE THE OCCUPIED-VIRTUAL LMO INTERACTION ENERGIES
	        !
	        if (idiagg <= 5 .or. Mod (idiagg, 2) == 0) then
	          l = 0
	          i1 = icvir(nnce(i)+1)
	          if (nce(i) > 1) then
	            i2 = icvir(nnce(i)+2)
	          else
	            i2 = i1
	          end if
	          if (lijbo) then
	            do j = 1, nocc
	              if (ijc /= nij) then
	                !
	                !  FAST TEST TO SEE IF THE INTEGRAL IS WORTH EVALUATING
	                !
	                j1 = icocc(nncf(j)+1)
	                if (ncf(j) > 1) then
	                  j2 = icocc(nncf(j)+2)
	                else
	                  j2 = j1
	                end if
	                !
	                ! nijbo is symmetric; use (j,i) for better memory locality.
	                i1j1 = nijbo(j1, i1)
	                i1j2 = nijbo(j2, i1)
	                i2j1 = nijbo(j1, i2)
	                i2j2 = nijbo(j2, i2)
	                !
	                if (i1j1 >= 0 .or. i1j2 >= 0 .or. i2j1 >= 0 .or. i2j2 >= 0) then
	                  sum = 0.d0
	                  if (i1j1 >= 0) then
	                    sum = Abs (fao(i1j1+1))
	                  end if
	                  if (i1j2 >= 0) then
	                    sum = sum + Abs (fao(i1j2+1))
	                  end if
	                  if (i2j1 >= 0) then
	                    sum = sum + Abs (fao(i2j1+1))
	                  end if
	                  if (i2j2 >= 0) then
	                    sum = sum + Abs (fao(i2j2+1))
	                  end if
	                  if (sum >= flim) then
	                    loopj = ncocc(j)
	                    lij = .false.
	                    sum = 0.d0
	                    kl = 0
	                    do kk = nncf(j) + 1, nncf(j) + ncf(j)
	                      k1 = icocc(kk)
	                      if (aocc(kk)*aov(k1) < cutoff .or. .not. latoms(k1)) then
	                        kl = kl + iorbs(k1)
	                      else
	                        lij = .true.
	                        do k = nfirst(k1), nlast(k1)
	                          kl = kl + 1
	                          sum = sum + ws(k) * cocc(kl+loopj)
	                        end do
	                      end if
	                    end do
	                    sumt = sumt + Abs (sum)
	                    tiny = Max (tiny, Abs (sum))
	                    if (lij) then
	                      if (Abs (sum) > oldlim) then
	                        l = l + 1
	                        ijc = ijc + 1
	                        ifmo(1, ijc) = i
	                        ifmo(2, ijc) = j
	                        fmo(ijc) = sum
	                      end if
	                    end if
	                  end if
	                end if
	              end if
	            end do
	          else
	            do j = 1, nocc
	              if (ijc /= nij) then
	                !
	                !  FAST TEST TO SEE IF THE INTEGRAL IS WORTH EVALUATING
	                !
	                j1 = icocc(nncf(j)+1)
	                if (ncf(j) > 1) then
	                  j2 = icocc(nncf(j)+2)
	                else
	                  j2 = j1
	                end if
	                !
	                i1j1 = ijbo (i1, j1)
	                i1j2 = ijbo (i1, j2)
	                i2j1 = ijbo (i2, j1)
	                i2j2 = ijbo (i2, j2)
	                !
	                if (i1j1 >= 0 .or. i1j2 >= 0 .or. i2j1 >= 0 .or. i2j2 >= 0) &
	                     & then
	                  sum = 0.d0
	                  if (i1j1 >= 0) then
	                    sum = Abs (fao(i1j1+1))
	                  end if
	                  if (i1j2 >= 0) then
	                    sum = sum + Abs (fao(i1j2+1))
	                  end if
	                  if (i2j1 >= 0) then
	                    sum = sum + Abs (fao(i2j1+1))
	                  end if
	                  if (i2j2 >= 0) then
	                    sum = sum + Abs (fao(i2j2+1))
	                  end if
	                  if (sum >= flim) then
	                    loopj = ncocc(j)
	                    lij = .false.
	                    sum = 0.d0
	                    kl = 0
	                    do kk = nncf(j) + 1, nncf(j) + ncf(j)
	                      k1 = icocc(kk)
	                      if (aocc(kk)*aov(k1) < cutoff .or. .not. latoms(k1)) then
	                        kl = kl + iorbs(k1)
	                      else
	                        lij = .true.
	                        do k = nfirst(k1), nlast(k1)
	                          kl = kl + 1
	                          sum = sum + ws(k) * cocc(kl+loopj)
	                        end do
	                      end if
	                    end do
	                    sumt = sumt + Abs (sum)
	                    tiny = Max (tiny, Abs (sum))
	                    if (lij) then
	                      if (Abs (sum) > oldlim) then
	                        l = l + 1
	                        ijc = ijc + 1
	                        ifmo(1, ijc) = i
	                        ifmo(2, ijc) = j
	                        fmo(ijc) = sum
	                      end if
	                    end if
	                  end if
	                end if
	              end if
	            end do
	          end if
	          nfmo(i) = l
	        else
	          ij0 = ijc
	          do jj = 1, nfmo(i)
	            if (ijc == nij) exit
	            j = ifmo(2, ij0+jj)
	            loopj = ncocc(j)
	            sum = 0.d0
	            kl = 0
	            do kk = nncf(j) + 1, nncf(j) + ncf(j)
	              k1 = icocc(kk)
	              if (aov(k1)*aocc(kk) < cutoff .or. .not. latoms(k1)) then
	                kl = kl + iorbs(k1)
	              else
	                do k = nfirst(k1), nlast(k1)
	                  kl = kl + 1
	                  sum = sum + ws(k) * cocc(kl+loopj)
	                end do
	              end if
	            end do
	            sumt = sumt + Abs (sum)
	            tiny = Max (tiny, Abs (sum))
	            ijc = ijc + 1
	            fmo(ijc) = sum
	          end do
	        end if
	      end do
#ifdef _OPENMP
	    end if
#endif
	    !
	    if (ijc == nij) then
	      !
	      !   THERE WAS NOT ENOUGH STORAGE TO HOLD ALL THE INTEGRALS.
      !   THEREFORE, ON THE NEXT ITERATION, CALCULATE FEWER INTEGRALS.
      !
      safety = safety * 2.d0
    else
      !
      !  THERE IS ENOUGH STORAGE FOR ALL THE INTEGRALS.  IF NECESSARY,
      !  CALCULATE MORE INTEGRALS.
      !
      safety = Max (safety*0.5d0, 1.d0)
    end if

    nij = ijc
    if (idiagg > 2 .and. Mod (idiagg, 4) /= 0) then
      fref = tiny ** 4
      if (times) then
        call timer (" AFTER DIAGG1 IN ITER")
      end if
        !
    else
	      !
	      !  EVALUATE THE OCCUPIED ENERGY LEVELS
	      !
	        !
#ifdef _OPENMP
!$omp parallel default(shared) private(avir_loc)
	      allocate (avir_loc(numat))
!$omp do schedule(static)
	      do i = 1, nocc
	        call build_occupied_energy (i, avir_loc, cutoff)
	      end do
!$omp end do
	      deallocate (avir_loc)
!$omp end parallel
#else
	      do i = 1, nocc
	        call build_occupied_energy (i, avir(1:numat), cutoff)
	      end do
#endif
	      oldlim = tiny * safety * 1.d-3
	      fref = tiny ** 4
	      if (times) then
	        call timer (" AFTER DIAGG1 IN ITER")
	      end if
	    end if

		    ovmax = tiny
		  contains
#ifdef _OPENMP
		  subroutine buffer_reserve (buf, cap)
		    type(diagg1_thread_buffer), intent (inout) :: buf
		    integer, intent (in) :: cap
		    integer :: new_cap
		    integer, allocatable :: i_new(:), j_new(:)
		    double precision, allocatable :: f_new(:)
		    if (cap <= 0) return
		    if (.not. allocated (buf%i)) then
		      new_cap = cap
		      allocate (buf%i(new_cap), buf%j(new_cap), buf%f(new_cap))
		      buf%n = 0
		    else if (size (buf%i) < cap) then
		      new_cap = cap
		      allocate (i_new(new_cap), j_new(new_cap), f_new(new_cap))
		      if (buf%n > 0) then
		        i_new(1:buf%n) = buf%i(1:buf%n)
		        j_new(1:buf%n) = buf%j(1:buf%n)
		        f_new(1:buf%n) = buf%f(1:buf%n)
		      end if
		      call move_alloc (i_new, buf%i)
		      call move_alloc (j_new, buf%j)
		      call move_alloc (f_new, buf%f)
		    end if
		  end subroutine buffer_reserve

		  subroutine buffer_append (buf, i, jlist, flist, nadd)
		    type(diagg1_thread_buffer), intent (inout) :: buf
		    integer, intent (in) :: i, nadd
		    integer, dimension (:), intent (in) :: jlist
		    double precision, dimension (:), intent (in) :: flist
		    integer :: needed, cap
		    if (nadd <= 0) return
		    needed = buf%n + nadd
		    if (.not. allocated (buf%i)) then
		      cap = Max (1024, needed)
		      call buffer_reserve (buf, cap)
		    else if (size (buf%i) < needed) then
		      cap = Max (needed, 2 * size (buf%i))
		      call buffer_reserve (buf, cap)
		    end if
		    buf%i(buf%n+1:needed) = i
		    buf%j(buf%n+1:needed) = jlist(1:nadd)
		    buf%f(buf%n+1:needed) = flist(1:nadd)
		    buf%n = needed
		  end subroutine buffer_append
#endif
			  subroutine build_virtual (i, prev_i, ws, latoms, avir, aov, cutoff)
			    integer, intent (in) :: i
			    integer, intent (in) :: prev_i
			    double precision, dimension (norbs), intent (inout) :: ws
			    logical, dimension (numat), intent (inout) :: latoms
			    double precision, dimension (numat), intent (inout) :: avir, aov
		    double precision, intent (in) :: cutoff
		    integer :: loopi, l, j, j1, jx, k, k1, kk, kj, kl, ii, i4
		    double precision :: sum

		    loopi = ncvir(i)
		    if (prev_i > 0) then
		      do j = nnce(prev_i) + 1, nnce(prev_i) + nce(prev_i)
		        j1 = icvir(j)
		        latoms(j1) = .false.
		      end do
		    else
		      latoms(:) = .false.
		    end if
		    l = 0
		    do j = nnce(i) + 1, nnce(i) + nce(i)
		      j1 = icvir(j)
		      sum = 0.d0
	      do k = l + 1, l + iorbs(j1)
	        sum = sum + cvir(k+loopi) ** 2
	      end do
	      l = l + iorbs(j1)
	      avir(j1) = sum
	      latoms(icvir(j)) = .true.
	    end do

	    if (lijbo) then
	      do j = nnce(i) + 1, nnce(i) + nce(i)
	        j1 = icvir(j)
	        do jx = 1, iorbs(j1)
	          ws(nfirst(j1)+jx-1) = 0.0d00
	        end do
	        kl = loopi
	        do kk = nnce(i) + 1, nnce(i) + nce(i)
	          k1 = icvir(kk)
	          kj = nijbo(k1, j1)
	          if (kj >= 0) then
	            if (avir(k1)*p(kj+1) > cutoff) then
	              if (iorbs(k1) == 1 .and. iorbs(j1) == 1) then
	                ws(nfirst(j1)) = ws(nfirst(j1)) + fao(kj+1) * cvir(kl+1)
	              else
	                if (k1 > j1) then
	                  ii = kj
	                  do i4 = 1, iorbs(k1)
	                    do jx = 1, iorbs(j1)
	                      ii = ii + 1
	                      ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) * cvir(kl+i4)
	                    end do
	                  end do
	                else if (k1 < j1) then
	                  ii = kj
	                  do jx = 1, iorbs(j1)
	                    do i4 = 1, iorbs(k1)
	                      ii = ii + 1
	                      ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) * cvir(kl+i4)
	                    end do
	                  end do
	                else
	                  do jx = 1, iorbs(j1)
	                    do i4 = 1, iorbs(j1)
	                      if (i4 > jx) then
	                        ii = kj + (i4*(i4-1)) / 2 + jx
	                      else
	                        ii = kj + (jx*(jx-1)) / 2 + i4
	                      end if
	                      ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) * cvir(kl+i4)
	                    end do
	                  end do
	                end if
	              end if
	            end if
	          end if
	          kl = kl + iorbs(k1)
	        end do
	      end do
	    else
	      do j = nnce(i) + 1, nnce(i) + nce(i)
	        j1 = icvir(j)
	        do jx = 1, iorbs(j1)
	          ws(nfirst(j1)+jx-1) = 0.0d00
	        end do
	        kl = loopi
	        do kk = nnce(i) + 1, nnce(i) + nce(i)
	          k1 = icvir(kk)
	          kj = ijbo (k1, j1)
	          if (kj >= 0) then
	            if (avir(k1)*p(kj+1) > cutoff) then
	              if (iorbs(k1) == 1 .and. iorbs(j1) == 1) then
	                ws(nfirst(j1)) = ws(nfirst(j1)) + fao(kj+1) * cvir(kl+1)
	              else
	                if (k1 > j1) then
	                  ii = kj
	                  do i4 = 1, iorbs(k1)
	                    do jx = 1, iorbs(j1)
	                      ii = ii + 1
	                      ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) * cvir(kl+i4)
	                    end do
	                  end do
	                else if (k1 < j1) then
	                  ii = kj
	                  do jx = 1, iorbs(j1)
	                    do i4 = 1, iorbs(k1)
	                      ii = ii + 1
	                      ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) * cvir(kl+i4)
	                    end do
	                  end do
	                else
	                  do jx = 1, iorbs(j1)
	                    do i4 = 1, iorbs(j1)
	                      if (i4 > jx) then
	                        ii = kj + (i4*(i4-1)) / 2 + jx
	                      else
	                        ii = kj + (jx*(jx-1)) / 2 + i4
	                      end if
	                      ws(nfirst(j1)+jx-1) = ws(nfirst(j1)+jx-1) + fao(ii) * cvir(kl+i4)
	                    end do
	                  end do
	                end if
	              end if
	            end if
	          end if
	          kl = kl + iorbs(k1)
	        end do
	      end do
		    end if
		    do j = nnce(i) + 1, nnce(i) + nce(i)
		      j1 = icvir(j)
		      sum = 0.d0
		      do k = nfirst(j1), nlast(j1)
		        sum = sum + ws(k) ** 2
		      end do
		      aov(j1) = sum
		    end do

		    sum = 0.d0
		    kl = loopi
		    do kk = nnce(i) + 1, nnce(i) + nce(i)
	      k1 = icvir(kk)
	      if (aov(k1)*avir(k1) > cutoff) then
	        do k = nfirst(k1), nlast(k1)
	          kl = kl + 1
	          sum = sum + ws(k) * cvir(kl)
	        end do
	      else
	        kl = kl + iorbs(k1)
	      end if
	    end do
		    eigv(i) = sum
		  end subroutine build_virtual

		  subroutine build_occupied_energy (i, avir_occ, cutoff)
		    integer, intent (in) :: i
		    double precision, dimension (numat), intent (inout) :: avir_occ
		    double precision, intent (in) :: cutoff

		    integer :: loopi, l, j, j1, j4, jx, k, k1, k1j1, kl, jl, i4, ii
		    double precision :: sum, sum1

		    loopi = ncocc(i)
		    l = 0
		    do j = nncf(i) + 1, nncf(i) + ncf(i)
		      j1 = icocc(j)
		      sum = 0.d0
		      do k = l + 1, l + iorbs(j1)
		        sum = sum + cocc(k+loopi) ** 2
		      end do
		      l = l + iorbs(j1)
		      avir_occ(j1) = sum
		    end do

		    sum = 0.d0
		    jl = loopi
		    if (lijbo) then
		      do j = nncf(i) + 1, nncf(i) + ncf(i)
		        j1 = icocc(j)
		        kl = loopi
		        do k = nncf(i) + 1, nncf(i) + ncf(i)
		          k1 = icocc(k)
		          k1j1 = nijbo(k1, j1)
		          if (k1j1 >= 0) then
		            if (avir_occ(k1)*p(k1j1+1)*aocc(k) >= cutoff) then
		              if (k1 > j1) then
		                do jx = 1, iorbs(j1)
		                  sum1 = 0.d0
		                  do i4 = 1, iorbs(k1)
		                    ii = k1j1 + (i4-1) * iorbs(j1) + jx
		                    sum1 = sum1 + fao(ii) * cocc(kl+i4)
		                  end do
		                  sum = sum + cocc(jl+jx) * sum1
		                end do
		              else if (k1 < j1) then
		                do jx = 1, iorbs(j1)
		                  sum1 = 0.d0
		                  do i4 = 1, iorbs(k1)
		                    ii = k1j1 + (jx-1) * iorbs(k1) + i4
		                    sum1 = sum1 + fao(ii) * cocc(kl+i4)
		                  end do
		                  sum = sum + cocc(jl+jx) * sum1
		                end do
		              else
		                do jx = 1, iorbs(j1)
		                  sum1 = 0.d0
		                  do j4 = 1, jx
		                    ii = k1j1 + (jx*(jx-1)) / 2 + j4
		                    sum1 = sum1 + fao(ii) * cocc(kl+j4)
		                  end do
		                  ii = k1j1 + (jx*(jx+1)) / 2
		                  do i4 = jx + 1, iorbs(k1)
		                    ii = k1j1 + (i4*(i4-1)) / 2 + jx
		                    sum1 = sum1 + fao(ii) * cocc(kl+i4)
		                  end do
		                  sum = sum + cocc(jl+jx) * sum1
		                end do
		              end if
		            end if
		          end if
		          kl = kl + iorbs(k1)
		        end do
		        jl = jl + iorbs(j1)
		      end do
		    else
		      do j = nncf(i) + 1, nncf(i) + ncf(i)
		        j1 = icocc(j)
		        kl = loopi
		        do k = nncf(i) + 1, nncf(i) + ncf(i)
		          k1 = icocc(k)
		          k1j1 = ijbo (k1, j1)
		          if (k1j1 >= 0) then
		            if (avir_occ(k1)*p(k1j1+1)*aocc(k) >= cutoff) then
		              if (k1 > j1) then
		                do jx = 1, iorbs(j1)
		                  sum1 = 0.d0
		                  do i4 = 1, iorbs(k1)
		                    ii = k1j1 + (i4-1) * iorbs(j1) + jx
		                    sum1 = sum1 + fao(ii) * cocc(kl+i4)
		                  end do
		                  sum = sum + cocc(jl+jx) * sum1
		                end do
		              else if (k1 < j1) then
		                do jx = 1, iorbs(j1)
		                  sum1 = 0.d0
		                  do i4 = 1, iorbs(k1)
		                    ii = k1j1 + (jx-1) * iorbs(k1) + i4
		                    sum1 = sum1 + fao(ii) * cocc(kl+i4)
		                  end do
		                  sum = sum + cocc(jl+jx) * sum1
		                end do
		              else
		                do jx = 1, iorbs(j1)
		                  sum1 = 0.d0
		                  do j4 = 1, jx
		                    ii = k1j1 + (jx*(jx-1)) / 2 + j4
		                    sum1 = sum1 + fao(ii) * cocc(kl+j4)
		                  end do
		                  ii = k1j1 + (jx*(jx+1)) / 2
		                  do i4 = jx + 1, iorbs(k1)
		                    ii = k1j1 + (i4*(i4-1)) / 2 + jx
		                    sum1 = sum1 + fao(ii) * cocc(kl+i4)
		                  end do
		                  sum = sum + cocc(jl+jx) * sum1
		                end do
		              end if
		            end if
		          end if
		          kl = kl + iorbs(k1)
		        end do
		        jl = jl + iorbs(j1)
		      end do
		    end if

		    eigs(i) = sum
		  end subroutine build_occupied_energy
	
		  subroutine rebuild_couplings (i, ws, latoms, aov, cutoff, flim, oldlim, jstore, fstore, sumt_i, tiny_i, nf)
		    integer, intent (in) :: i
		    double precision, dimension (norbs), intent (in) :: ws
	    logical, dimension (numat), intent (in) :: latoms
	    double precision, dimension (numat), intent (in) :: aov
	    double precision, intent (in) :: cutoff, flim, oldlim
	    integer, dimension (nocc), intent (out) :: jstore
	    double precision, dimension (nocc), intent (out) :: fstore
	    double precision, intent (out) :: sumt_i, tiny_i
	    integer, intent (out) :: nf

	    integer :: i1, i2, i1j1, i1j2, i2j1, i2j2
	    integer :: j, j1, j2, k, k1, kk, kl, loopj
	    double precision :: sum
	    logical :: lij

	    nf = 0
	    sumt_i = 0.d0
	    tiny_i = 0.d0

	    i1 = icvir(nnce(i)+1)
	    if (nce(i) > 1) then
	      i2 = icvir(nnce(i)+2)
	    else
	      i2 = i1
	    end if

	    if (lijbo) then
	      do j = 1, nocc
	        j1 = icocc(nncf(j)+1)
	        if (ncf(j) > 1) then
	          j2 = icocc(nncf(j)+2)
	        else
	          j2 = j1
	        end if

	        ! nijbo is symmetric; use (j,i) for better memory locality.
	        i1j1 = nijbo(j1, i1)
	        i1j2 = nijbo(j2, i1)
	        i2j1 = nijbo(j1, i2)
	        i2j2 = nijbo(j2, i2)

	        if (i1j1 >= 0 .or. i1j2 >= 0 .or. i2j1 >= 0 .or. i2j2 >= 0) then
	          sum = 0.d0
	          if (i1j1 >= 0) sum = Abs (fao(i1j1+1))
	          if (i1j2 >= 0) sum = sum + Abs (fao(i1j2+1))
	          if (i2j1 >= 0) sum = sum + Abs (fao(i2j1+1))
	          if (i2j2 >= 0) sum = sum + Abs (fao(i2j2+1))

	          if (sum >= flim) then
	            loopj = ncocc(j)
	            lij = .false.
	            sum = 0.d0
	            kl = 0
	            do kk = nncf(j) + 1, nncf(j) + ncf(j)
	              k1 = icocc(kk)
	              if (.not. latoms(k1)) then
	                kl = kl + iorbs(k1)
	              else if (aocc(kk)*aov(k1) < cutoff) then
	                kl = kl + iorbs(k1)
	              else
	                lij = .true.
	                do k = nfirst(k1), nlast(k1)
	                  kl = kl + 1
	                  sum = sum + ws(k) * cocc(kl+loopj)
	                end do
	              end if
	            end do

	            sumt_i = sumt_i + Abs (sum)
	            tiny_i = Max (tiny_i, Abs (sum))
	            if (lij) then
	              if (Abs (sum) > oldlim) then
	                nf = nf + 1
	                jstore(nf) = j
	                fstore(nf) = sum
	              end if
	            end if
	          end if
	        end if
	      end do
	    else
	      do j = 1, nocc
	        j1 = icocc(nncf(j)+1)
	        if (ncf(j) > 1) then
	          j2 = icocc(nncf(j)+2)
	        else
	          j2 = j1
	        end if

	        i1j1 = ijbo (i1, j1)
	        i1j2 = ijbo (i1, j2)
	        i2j1 = ijbo (i2, j1)
	        i2j2 = ijbo (i2, j2)

	        if (i1j1 >= 0 .or. i1j2 >= 0 .or. i2j1 >= 0 .or. i2j2 >= 0) then
	          sum = 0.d0
	          if (i1j1 >= 0) sum = Abs (fao(i1j1+1))
	          if (i1j2 >= 0) sum = sum + Abs (fao(i1j2+1))
	          if (i2j1 >= 0) sum = sum + Abs (fao(i2j1+1))
	          if (i2j2 >= 0) sum = sum + Abs (fao(i2j2+1))

	          if (sum >= flim) then
	            loopj = ncocc(j)
	            lij = .false.
	            sum = 0.d0
	            kl = 0
	            do kk = nncf(j) + 1, nncf(j) + ncf(j)
	              k1 = icocc(kk)
	              if (.not. latoms(k1)) then
	                kl = kl + iorbs(k1)
	              else if (aocc(kk)*aov(k1) < cutoff) then
	                kl = kl + iorbs(k1)
	              else
	                lij = .true.
	                do k = nfirst(k1), nlast(k1)
	                  kl = kl + 1
	                  sum = sum + ws(k) * cocc(kl+loopj)
	                end do
	              end if
	            end do

	            sumt_i = sumt_i + Abs (sum)
	            tiny_i = Max (tiny_i, Abs (sum))
	            if (lij) then
	              if (Abs (sum) > oldlim) then
	                nf = nf + 1
	                jstore(nf) = j
	                fstore(nf) = sum
	              end if
	            end if
	          end if
	        end if
	      end do
		    end if
		  end subroutine rebuild_couplings
		end subroutine diagg1
