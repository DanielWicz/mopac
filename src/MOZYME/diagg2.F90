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

  subroutine diagg2 (nocc, nvir,  eigv, iused, latoms, &
     & nij, idiagg, storei, storej)
   !***********************************************************************
   !
   !   DIAGG2  PERFORMS A SIMPLE JACOBIAN ANNIHILATION OF THE ENERGY TERMS
   !   CONNECTING THE OCCUPIED LMOS AND THE VIRTUAL LMOS.  THE ENERGY TERMS
   !   ARE IN THE ARRAY FMO, AND THE INDICES OF THE LMOS ARE IN IFMO.
   !
   !***********************************************************************
    use molkst_C, only: numat, norbs, numcal, keywrd, id
    use MOZYME_C, only : nvirtual, icocc_dim, shift, &
       & icvir_dim, cocc_dim, cvir_dim, ipad2, thresh, tiny, sumb, &
       ncf, nce, nncf, nnce, ncocc, ncvir, iorbs, cocc, cvir, icocc, icvir, &
       ifmo, fmo
!
    use common_arrays_C, only : eigs, nat
    use parameters_C, only: main_group
    use chanel_C, only: iw
#ifdef _OPENMP
	    use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_wtime
#endif
    implicit none
    integer, intent (in) :: idiagg, nij, nocc, nvir
    logical, dimension (numat), intent (out) :: latoms
    integer, dimension (numat), intent (out) :: iused
    double precision, dimension (norbs), intent (out) :: storei, storej
    double precision, dimension (nvirtual), intent (in) :: eigv
    logical :: bug = .false.
    logical :: retry
    logical :: use_parallel_diagg2
    logical, save :: debug, times
    integer :: i, ii, jur, l
    integer, save :: icalcn = 0
    integer :: ij, ilr, incv, iur, j, jlr, jncf, k, le, lf, lij, loopi, loopj, &
   & mie, mle, mlee, mlf, mlff, ncei, ncfj, nrej, njlen, nilen
    double precision :: a, alpha, b, beta, biglim, c, d, e, sum
    double precision, save :: const, eps, eta, bigeps
    double precision, external :: reada
    integer, dimension (2) :: nrejct
    data nrejct / 2 * 0 /
    !
    ! Cluster width for DIAGG2 multi-neighbor updates in the OpenMP level scheduler.
    ! Kept outside of _OPENMP guards so the file also compiles when OpenMP is disabled.
    !
    integer, parameter :: diagg2_block_cap = 4
#ifdef _OPENMP
		    integer, parameter :: mode_level = 1, mode_task = 2
		    double precision, parameter :: tune_alpha = 0.25d0, tune_switch_margin = 1.02d0
		    integer, parameter :: tune_explore_period = 24
		    integer :: alloc_stat, max_threads, tid, diagg2_threads, width_hint
		    integer :: lvl, max_level, max_level_est, nactive, pos
		    integer :: filter_threads, filter_tid, filter_start, filter_end
		    integer, allocatable :: filter_count(:), filter_offset(:)
		    integer :: batch_start, batch_end, batch_target, batch_edges
		    integer :: batch_min, batch_max
		    integer :: i_task, j_task
		    integer :: mode_selected, mode_used, tune_slot
		    integer :: tune_calls_slot
	    logical :: task_capable, force_explore
	    double precision :: t_start, elapsed, edge_cost
		    integer, allocatable :: edge_level(:), active_edges(:), last_occ(:), last_vir(:)
			    integer, allocatable :: rem_edges(:), next_edges(:)
			    integer, allocatable :: level_count(:), level_offset(:), level_pos(:), level_edges(:)
			    integer :: max_deg, max_deg_occ, max_deg_vir, key
			    integer :: nocc_active, nvir_active
			    integer, allocatable :: deg_occ(:), deg_vir(:)
			    integer, allocatable :: bucket_count(:), bucket_pos(:)
			    integer :: block_row_max
			    integer :: block_cap, nblocks, blk
		    integer, allocatable :: occ_block_id(:)
		    integer, allocatable :: occ_count(:)
		    integer, allocatable :: level_block_offset(:)
		    integer, allocatable :: block_occ(:), block_nedges(:)
		    integer, allocatable :: block_edges(:,:)
		    integer, allocatable, save :: iused_ws(:,:)
		    logical, allocatable, save :: latoms_ws(:,:)
		    double precision, allocatable, save :: storei_ws(:,:), storej_ws(:,:)
		    double precision, allocatable, save :: block_c_ws(:,:,:), block_c2_ws(:,:,:)
		    integer, allocatable, save :: block_atoms_ws(:,:)
		    double precision, allocatable :: sumb_thread(:)
		    integer, allocatable :: nrej_thread(:)
		    logical :: need_ws_resize
		    logical :: active
		    integer :: occ_span_max, vir_span_max, nrem, nnext, iedge
		    logical :: use_cluster_level
		    integer :: depth_match_lb, depth_cluster_lb
	    integer, allocatable, save :: tune_calls(:), tune_task_samples(:), tune_level_samples(:)
    integer, allocatable, save :: tune_batch(:), tune_batch_step(:), tune_batch_dir(:)
    double precision, allocatable, save :: tune_task_cost(:), tune_level_cost(:), tune_prev_task_cost(:)
    double precision :: sumb_local
    integer :: nrej_local
#endif
    if (numcal /= icalcn) then
      icalcn = numcal
      times = (Index (keywrd, " TIMES") /= 0)
      debug = (Index (keywrd, " DIAGG2") /= 0)
      !
      !   IF THE SYSTEM IS A SOLID, THEN DAMP ROTATION OF VECTORS,
      !   IN AN ATTEMPT TO PREVENT AUTOREGENERATIVE CHARGE OSCILLATION.
      !
      i = Index (keywrd, " DAMP")
      if (i /= 0) then
        const = reada (keywrd, i+5)
      else if (id == 3) then
        const = 0.5d0
      else
        a = 1.d0               !  If a transition metal, set DAMP to 0.5d0
        do i = 1, numat        !
          if (.not. main_group(nat(i)))  a = 0.5d0
        end do
        const = a
      end if
      !
      !   EPS IS THE SMALLEST NUMBER WHICH, WHEN ADDED TO 1.D0, IS NOT
      !   EQUAL TO 1.D0
      call epseta (eps, eta)
      !
      !   INCREASE EPS TO ALLOW FOR A LOT OF ROUND-OFF
      !
      bigeps = 50.d0 * Sqrt (eps)
    end if
   !
   !  RETRY IS .TRUE. IF THE NUMBER OF REJECTED ANNIHILATIONS IS IN THE
   !         RANGE 1 TO 20  AND THE SAME ON TWO ITERATIONS.  THIS WILL
   !         OCCUR NEAR THE END OF A SCF CALCULATION, WHEN ONLY A FEW
   !         LMOS ARE BADLY BEHAVED.
   !
    retry = (nrejct(1) == nrejct(2) .and. nrejct(1) /= 0 .and. nrejct(1) < 20)
    if (Mod(idiagg, 5) == 0 .or. idiagg <= 5) then
      tiny = -1.d0
      biglim = -1.d0
    else
      tiny = 0.01d0 * tiny
      biglim = bigeps
    end if
   !***********************************************************************
   !
   !   DO A CRUDE 2 BY 2 ROTATION TO "ELIMINATE" SIGNIFICANT ELEMENTS
   !
   !***********************************************************************
    iused(:) = -1
    latoms(:) = .false.
    if (debug) then
      write (iw,*)
      write (iw,*) "            SIZE OF OCCUPIED ARRAYS IN DIAGG2"
      write (iw,*)
      write (iw,*) "    LMO    NNCF     NCF   SPACE  ", " NCOCC    SIZE   SPACE"
      do i = 1, nocc - 1
        l = ncocc(i)
        do j = nncf(i) + 1, nncf(i) + ncf(i)
          l = l + iorbs(icocc(j))
        end do
        write (iw, "(7I8)") i, nncf (i), ncf (i), nncf (i+1) - &
       & nncf(i) - ncf(i), ncocc(i), l - ncocc(i), ncocc(i+1) - l
      end do
      i = nocc
      l = ncocc(i)
      do j = nncf(i) + 1, nncf(i) + ncf(i)
        l = l + iorbs(icocc(j))
      end do
      write (iw, "(7I8)") i, nncf (i), ncf (i), icocc_dim - nncf (i) - ncf &
     & (i), ncocc(i), l - ncocc(i), cocc_dim - l
      write (iw,*)
      write (iw,*) "            SIZE OF VIRTUAL ARRAYS IN DIAGG2"
      write (iw,*)
      write (iw,*) "    LMO    NNCE     NCE   SPACE  ", " NCVIR    SIZE   SPACE"
      do i = 1, nvir - 1
        l = ncvir(i)
        do j = nnce(i) + 1, nnce(i) + nce(i)
          l = l + iorbs(icvir(j))
        end do
        write (iw, "(7I8)") i, nnce (i), nce (i), nnce (i+1) - nnce (i) - nce &
       & (i), ncvir(i), l - ncvir(i), ncvir(i+1) - l
      end do
      i = nvir
      l = ncvir(i)
      do j = nnce(i) + 1, nnce(i) + nce(i)
        l = l + iorbs(icvir(j))
      end do
      write (iw, "(7I8)") i, nnce (i), nce (i), icvir_dim - nnce (i) - nce (i), &
     & ncvir(i), l - ncvir(i), cvir_dim - l
      if (bug) then
        write (iw,*)
        write (iw, "(A,I3,A)") " THIS FAULT CAN PROBABLY BE " // &
                             & "CORRECTED BY USE OF KEYWORD 'NLMO=", &
                             & ipad2 + 50, "'"
        write (iw,*)
        call mopend ("VALUE OF NLMO IS TOO SMALL")
      end if
    end if
    sumb = 0.d0
    nrej = 0
    lij = 0
    use_parallel_diagg2 = .false.
#ifdef _OPENMP
	    max_threads = omp_get_max_threads()
	    use_parallel_diagg2 = (max_threads > 1 .and. nij > 64)
		    if (use_parallel_diagg2) then
		      allocate (edge_level(nij), active_edges(nij), rem_edges(nij), next_edges(nij), &
		           & last_occ(nocc), last_vir(nvir), stat=alloc_stat)
		      if (alloc_stat /= 0) then
		        call mopend("Insufficient memory to run DIAGG2 (parallel scheduling)")
		        use_parallel_diagg2 = .false.
		      else
		        edge_level(:) = 0
		        nactive = 0
		        !
		        ! Filter active pairs using the same tiny/biglim criterion as the rotation kernel.
		        !
		        !
		        ! Parallelize the ij scan (nij may be large) while preserving deterministic
		        ! ij ordering in active_edges(1:nactive).
		        !
		        filter_threads = max_threads
		        allocate (filter_count(filter_threads), filter_offset(filter_threads+1), stat=alloc_stat)
		        if (alloc_stat /= 0) then
		          call mopend("Insufficient memory to run DIAGG2 (active filter)")
		        end if
		        filter_count(:) = 0
		        !
		        ! Eq. (D5): per-thread count of active edges over a static ij partition.
		        !
		        !$omp parallel default(shared) &
		        !$omp& private(filter_tid, filter_start, filter_end, ij, i, j, active, c, d)
		          filter_tid = omp_get_thread_num() + 1
		          filter_start = ( (filter_tid-1) * nij ) / filter_threads + 1
		          filter_end   = ( filter_tid * nij ) / filter_threads
		          do ij = filter_start, filter_end
		            i = ifmo(1, ij)
		            j = ifmo(2, ij)
		            active = .true.
		            if (tiny >= 0.d0) then
		              if (Abs (fmo(ij)) < tiny) active = .false.
		            end if
		            if (active .and. biglim >= 0.d0) then
		              c = fmo(ij) * const
		              d = eigs(j) - eigv(i) - shift
		              if (d == 0.d0) then
		                active = (Abs(c) > 0.d0)
		              else
		                active = (Abs (c/d) >= biglim)
		              end if
		            end if
		            if (active) then
		              edge_level(ij) = 1
		              filter_count(filter_tid) = filter_count(filter_tid) + 1
		            end if
		          end do
		        !$omp end parallel
		        !
		        ! Exclusive prefix sum of counts to get deterministic per-thread output ranges.
		        !
		        filter_offset(1) = 1
		        do filter_tid = 1, filter_threads
		          filter_offset(filter_tid+1) = filter_offset(filter_tid) + filter_count(filter_tid)
		        end do
		        nactive = filter_offset(filter_threads+1) - 1
		        !
		        ! Scatter active ij into active_edges in ascending ij order.
		        !
		        !$omp parallel default(shared) &
		        !$omp& private(filter_tid, filter_start, filter_end, ij, pos)
		          filter_tid = omp_get_thread_num() + 1
		          filter_start = ( (filter_tid-1) * nij ) / filter_threads + 1
		          filter_end   = ( filter_tid * nij ) / filter_threads
		          pos = filter_offset(filter_tid)
		          do ij = filter_start, filter_end
		            if (edge_level(ij) == 1) then
		              active_edges(pos) = ij
		              pos = pos + 1
		            end if
		          end do
		        !$omp end parallel
		        !
		        ! Reset markers so edge_level is available for later level assignment.
		        !
		        do pos = 1, nactive
		          edge_level(active_edges(pos)) = 0
		        end do
		        deallocate (filter_count, filter_offset)

		        if (nactive > 0) then
		          !
		          ! Degree-aware edge ordering (deterministic, stable) to reduce schedule depth.
		          !
		          allocate (deg_occ(nocc), deg_vir(nvir), stat=alloc_stat)
		          if (alloc_stat /= 0) then
		            call mopend("Insufficient memory to run DIAGG2 (edge degrees)")
		          end if
		          deg_occ(:) = 0
		          deg_vir(:) = 0
			          do iedge = 1, nactive
			            ij = active_edges(iedge)
			            i = ifmo(1, ij)
			            j = ifmo(2, ij)
			            deg_occ(j) = deg_occ(j) + 1
			            deg_vir(i) = deg_vir(i) + 1
			          end do
			          max_deg_occ = 0
			          nocc_active = 0
			          do j = 1, nocc
			            if (deg_occ(j) > 0) nocc_active = nocc_active + 1
			            max_deg_occ = Max (max_deg_occ, deg_occ(j))
			          end do
			          max_deg_vir = 0
			          nvir_active = 0
			          do i = 1, nvir
			            if (deg_vir(i) > 0) nvir_active = nvir_active + 1
			            max_deg_vir = Max (max_deg_vir, deg_vir(i))
			          end do
			          max_deg = Max (max_deg_occ, max_deg_vir)
		          !
		          ! Eq. (D6): key(e=(i,j)) = max(deg_occ(j), deg_vir(i)).
		          ! Stable counting sort by descending key clusters high-degree nodes early,
		          ! improving the greedy matching fill and reducing max_level/barriers.
		          !
		          if (max_deg > 1) then
		            allocate (bucket_count(max_deg), bucket_pos(max_deg), stat=alloc_stat)
		            if (alloc_stat /= 0) then
		              call mopend("Insufficient memory to run DIAGG2 (edge ordering)")
		            end if
		            bucket_count(:) = 0
		            do iedge = 1, nactive
		              ij = active_edges(iedge)
		              i = ifmo(1, ij)
		              j = ifmo(2, ij)
		              key = Max (deg_occ(j), deg_vir(i))
		              bucket_count(key) = bucket_count(key) + 1
		            end do
		            bucket_pos(max_deg) = 1
		            do key = max_deg - 1, 1, -1
		              bucket_pos(key) = bucket_pos(key+1) + bucket_count(key+1)
		            end do
		            do iedge = 1, nactive
		              ij = active_edges(iedge)
		              i = ifmo(1, ij)
		              j = ifmo(2, ij)
		              key = Max (deg_occ(j), deg_vir(i))
		              pos = bucket_pos(key)
		              next_edges(pos) = ij
		              bucket_pos(key) = pos + 1
		            end do
		            active_edges(1:nactive) = next_edges(1:nactive)
		            deallocate (bucket_count, bucket_pos)
		          end if
		          deallocate (deg_occ, deg_vir)
		          !
		          ! Eq. (D0): one-pass level-depth estimate
		          !   l(e_k) = 1 + max(last_occ(j_k), last_vir(i_k))
		          ! gives an upper bound for schedule depth at O(nactive) cost.
	          ! Used for thread-width and mode selection before deciding whether
	          ! the full repeated-maximal-matching schedule is required.
	          !
	          last_occ(:) = 0
	          last_vir(:) = 0
	          max_level_est = 0
	          do iedge = 1, nactive
	            ij = active_edges(iedge)
	            i = ifmo(1, ij)
	            j = ifmo(2, ij)
	            lvl = Max (last_occ(j), last_vir(i)) + 1
	            last_occ(j) = lvl
	            last_vir(i) = lvl
	            max_level_est = Max (max_level_est, lvl)
	          end do

	            ! Previous-commit improvement kept here: adaptive team sizing from level width.
	            ! Eq. (0): width_hint = ceil(nactive / max_level_est)
	            !         diagg2_threads = min(max_threads, 2*width_hint, nactive)
	            width_hint = Max (1, (nactive + max_level_est - 1) / max_level_est)
	            diagg2_threads = Min (max_threads, Max (1, 2 * width_hint))
	            diagg2_threads = Min (diagg2_threads, nactive)

	            if (diagg2_threads > 1) then
		              occ_span_max = 1
		              do j = 1, nocc
		                jlr = ncocc(j) + 1
		                if (j /= nocc) then
		                  jur = ncocc(j+1)
		                else
		                  jur = cocc_dim
		                end if
		                jur = Min (jlr+norbs-1, jur)
		                occ_span_max = Max (occ_span_max, jur - jlr + 1)
		              end do
		              vir_span_max = 1
		              do i = 1, nvir
		                ilr = ncvir(i) + 1
		                if (i /= nvir) then
		                  iur = ncvir(i+1)
		                else
		                  iur = cvir_dim
		                end if
		                iur = Min (ilr+norbs-1, iur)
		                vir_span_max = Max (vir_span_max, iur - ilr + 1)
		              end do
		              !
		              ! Block workspace rows bound: union(AO) over one occupied + up to diagg2_block_cap virtuals.
		              ! Eq. (B0): m_max <= min(norbs, span_occ_max + k*span_vir_max), k = diagg2_block_cap
		              !
		              block_row_max = Min (norbs, occ_span_max + diagg2_block_cap * vir_span_max)
		              need_ws_resize = .true.
		              ! Previous-commit improvement kept here: persistent per-thread scratch reuse.
		              ! Allocation is skipped when dimensions already satisfy:
		              !   ws_rows >= required_rows and ws_cols >= diagg2_threads.
		              ! Eq. (D2): backup rows are bounded by max reserved LMO span, not global norbs.
		              if (allocated(iused_ws) .and. allocated(latoms_ws) .and. allocated(storei_ws) .and. allocated(storej_ws)) then
		                if (size(iused_ws,1) >= numat .and. size(iused_ws,2) >= diagg2_threads .and. &
		                   & size(latoms_ws,1) >= numat .and. size(latoms_ws,2) >= diagg2_threads .and. &
		                   & size(storei_ws,1) >= vir_span_max .and. size(storei_ws,2) >= diagg2_threads .and. &
		                   & size(storej_ws,1) >= occ_span_max .and. size(storej_ws,2) >= diagg2_threads .and. &
		                   & allocated(block_c_ws) .and. allocated(block_c2_ws) .and. allocated(block_atoms_ws) .and. &
		                   & size(block_c_ws,1) >= block_row_max .and. size(block_c_ws,2) >= diagg2_block_cap+1 .and. &
		                   & size(block_c_ws,3) >= diagg2_threads .and. &
		                   & size(block_c2_ws,1) >= block_row_max .and. size(block_c2_ws,2) >= diagg2_block_cap+1 .and. &
		                   & size(block_c2_ws,3) >= diagg2_threads .and. &
		                   & size(block_atoms_ws,1) >= numat .and. size(block_atoms_ws,2) >= diagg2_threads) then
		                  need_ws_resize = .false.
		                end if
		              end if

	              if (need_ws_resize) then
	                if (allocated(iused_ws)) deallocate (iused_ws)
	                if (allocated(latoms_ws)) deallocate (latoms_ws)
	                if (allocated(storei_ws)) deallocate (storei_ws)
	                if (allocated(storej_ws)) deallocate (storej_ws)
	                if (allocated(block_c_ws)) deallocate (block_c_ws)
	                if (allocated(block_c2_ws)) deallocate (block_c2_ws)
	                if (allocated(block_atoms_ws)) deallocate (block_atoms_ws)

		                ! Scratch footprint scales as O((numat + span_occ + span_vir)*threads), reused.
	                allocate (iused_ws(numat, diagg2_threads), latoms_ws(numat, diagg2_threads), &
	                     & storei_ws(vir_span_max, diagg2_threads), &
	                     & storej_ws(occ_span_max, diagg2_threads), &
	                     & block_c_ws(block_row_max, diagg2_block_cap+1, diagg2_threads), &
	                     & block_c2_ws(block_row_max, diagg2_block_cap+1, diagg2_threads), &
	                     & block_atoms_ws(numat, diagg2_threads), stat=alloc_stat)
	                if (alloc_stat /= 0) then
	                  call mopend("Insufficient memory to run DIAGG2 (parallel scratch workspace)")
	                end if
	                !
	                ! Eq. (D3): initialize thread-local marker state once per (re)allocation.
	                ! Per-pair sparse resets keep markers clean, so we avoid O(numat*threads)
	                ! full clears at each DIAGG2 call.
	                !
	                iused_ws(:, :) = -1
	                latoms_ws(:, :) = .false.
	              end if

	              if (.not. allocated(tune_calls) .or. size(tune_calls) < max_threads) then
	                if (allocated(tune_calls)) deallocate (tune_calls)
	                if (allocated(tune_task_samples)) deallocate (tune_task_samples)
	                if (allocated(tune_level_samples)) deallocate (tune_level_samples)
	                if (allocated(tune_batch)) deallocate (tune_batch)
	                if (allocated(tune_batch_step)) deallocate (tune_batch_step)
	                if (allocated(tune_batch_dir)) deallocate (tune_batch_dir)
	                if (allocated(tune_task_cost)) deallocate (tune_task_cost)
	                if (allocated(tune_level_cost)) deallocate (tune_level_cost)
	                if (allocated(tune_prev_task_cost)) deallocate (tune_prev_task_cost)

	                allocate (tune_calls(max_threads), tune_task_samples(max_threads), tune_level_samples(max_threads), &
	                     & tune_batch(max_threads), tune_batch_step(max_threads), tune_batch_dir(max_threads), &
	                     & tune_task_cost(max_threads), tune_level_cost(max_threads), tune_prev_task_cost(max_threads), &
	                     & stat=alloc_stat)
	                if (alloc_stat /= 0) then
	                  call mopend("Insufficient memory to run DIAGG2 (autotuner workspace)")
	                end if
	                tune_calls(:) = 0
	                tune_task_samples(:) = 0
	                tune_level_samples(:) = 0
	                tune_batch(:) = 0
	                tune_batch_step(:) = 0
	                tune_batch_dir(:) = 1
	                tune_task_cost(:) = 0.d0
	                tune_level_cost(:) = 0.d0
	                tune_prev_task_cost(:) = 0.d0
	              end if

	              tune_slot = diagg2_threads
	              if (tune_batch(tune_slot) <= 0) then
	                tune_batch(tune_slot) = Max (64, 8*diagg2_threads)
	              end if
	              if (tune_batch_step(tune_slot) <= 0) then
	                tune_batch_step(tune_slot) = Max (16, tune_batch(tune_slot)/4)
	              end if
	              if (tune_batch_dir(tune_slot) == 0) then
	                tune_batch_dir(tune_slot) = 1
	              end if

		              tune_calls(tune_slot) = tune_calls(tune_slot) + 1
		              tune_calls_slot = tune_calls(tune_slot)
		              !
		              ! Avoid task mode when it is unlikely to amortize task/depend overhead.
		              ! Heuristic: require enough work per thread (k tasks/thread).
		              !
		              task_capable = (max_level_est > 1 .and. nactive >= 8*diagg2_threads)
		              mode_selected = mode_level

	              if (task_capable) then
	                if (tune_level_samples(tune_slot) == 0) then
	                  mode_selected = mode_level
	                else if (tune_task_samples(tune_slot) == 0) then
	                  mode_selected = mode_task
	                else
	                  force_explore = (tune_task_samples(tune_slot) >= 3 .and. tune_level_samples(tune_slot) >= 3 .and. &
	                       & Mod(tune_calls_slot, tune_explore_period) == 0)
	                  if (force_explore) then
	                    if (tune_task_cost(tune_slot) <= tune_level_cost(tune_slot)) then
	                      mode_selected = mode_level
	                    else
	                      mode_selected = mode_task
	                    end if
	                  else
	                    ! Eq. (4) Exploit faster mode with hysteresis:
	                    ! choose task mode when c_task <= m * c_level, m = tune_switch_margin.
	                    if (tune_task_cost(tune_slot) <= tune_switch_margin*tune_level_cost(tune_slot)) then
	                      mode_selected = mode_task
	                    else
	                      mode_selected = mode_level
	                    end if
	                  end if
	                end if
	              end if

	              mode_used = mode_selected
	              if (mode_selected == mode_task) then
	                allocate (sumb_thread(diagg2_threads), nrej_thread(diagg2_threads), stat=alloc_stat)
	                if (alloc_stat /= 0) then
	                  ! Robust fallback to level mode if temporary task accumulators cannot be allocated.
	                  mode_used = mode_level
	                end if
	              else
	                mode_used = mode_level
	              end if

		              if (mode_used == mode_task) then
		                sumb_thread(:) = 0.d0
		                nrej_thread(:) = 0
		                batch_min = Max (32, 4*diagg2_threads)
		                batch_max = Max (batch_min, Min (nactive, 8192))
		                batch_target = Min (batch_max, Max (batch_min, tune_batch(tune_slot)))

		                t_start = omp_get_wtime()
!$omp parallel num_threads(diagg2_threads) default(shared) &
!$omp& private(ij, pos, tid, i_task, j_task, batch_start, batch_end)
	                tid = omp_get_thread_num() + 1
!$omp single
		                batch_start = 1
		                do while (batch_start <= nactive)
		                  batch_end = Min (nactive, batch_start + batch_target - 1)
		                  do pos = batch_start, batch_end
		                    ij = active_edges(pos)
		                    i_task = ifmo(1, ij)
		                    j_task = ifmo(2, ij)
		                    ! Two tasks are independent iff they do not share indices:
		                    !   (i1 /= i2) AND (j1 /= j2)
		                    ! Dependence tokens (last_vir(i), last_occ(j)) enforce this relation.
!$omp task default(shared) firstprivate(ij, i_task, j_task) private(tid) &
!$omp& depend(inout:last_occ(j_task), last_vir(i_task))
		                    tid = omp_get_thread_num() + 1
		                    call process_pair (ij, retry, tiny, biglim, sumb_thread(tid), nrej_thread(tid), &
		                         & iused_ws(:,tid), latoms_ws(:,tid), storei_ws(:,tid), storej_ws(:,tid))
!$omp end task
		                  end do
		                  batch_start = batch_end + 1
	                end do
		                !
		                ! Eq. (D4): one synchronization per task sweep.
		                ! Dependencies already serialize conflicting (i,j) updates,
		                ! so waiting once avoids O(nbatch) barrier overhead.
		                !
!$omp taskwait
!$omp end single
!$omp end parallel
	                elapsed = Max (1.d-12, omp_get_wtime() - t_start)
	                ! Eq. (1) Normalized cost per active pair: c = elapsed / nactive
	                edge_cost = elapsed / Dble (Max (1, nactive))
	                if (tune_task_samples(tune_slot) == 0) then
	                  tune_task_cost(tune_slot) = edge_cost
	                else
	                  ! Eq. (2) EWMA update: c_new = (1-alpha)*c_old + alpha*c_obs
	                  tune_task_cost(tune_slot) = (1.d0-tune_alpha)*tune_task_cost(tune_slot) + tune_alpha*edge_cost
	                end if
	                tune_task_samples(tune_slot) = tune_task_samples(tune_slot) + 1

	                if (tune_prev_task_cost(tune_slot) > 0.d0) then
	                  ! If no improvement (c_k >= 0.99*c_{k-1}), reverse search direction and shrink step.
	                  if (edge_cost >= 0.99d0*tune_prev_task_cost(tune_slot)) then
	                    tune_batch_dir(tune_slot) = -tune_batch_dir(tune_slot)
	                    tune_batch_step(tune_slot) = Max (16, tune_batch_step(tune_slot)/2)
	                  end if
	                end if
	                tune_prev_task_cost(tune_slot) = edge_cost
	                ! Eq. (3) Batch walk: b_{k+1} = clamp(b_k + dir*step, b_min, b_max)
	                tune_batch(tune_slot) = batch_target + tune_batch_dir(tune_slot)*tune_batch_step(tune_slot)
	                tune_batch(tune_slot) = Min (batch_max, Max (batch_min, tune_batch(tune_slot)))

	                sumb = 0.d0
	                nrej = 0
	                do tid = 1, diagg2_threads
	                  sumb = sumb + sumb_thread(tid)
	                  nrej = nrej + nrej_thread(tid)
	                end do
		                deallocate (sumb_thread, nrej_thread)
			              else
			                !
		                ! Level mode has two implementations:
		                !   (A) occupied-centered multi-virtual units (up to k virtuals per occupied per level)
		                !   (B) repeated maximal matchings (unique i and j per level)
		                !
		                ! The multi-virtual path can reduce barrier depth when occupied degrees dominate,
		                ! but it also increases per-unit work. Gate it with a cheap deterministic heuristic.
		                !
		                block_cap = diagg2_block_cap
		                depth_match_lb = max_deg
		                depth_cluster_lb = Max (max_deg_vir, (max_deg_occ + block_cap - 1) / block_cap)
		                !
		                ! Eq. (G1): enable multi-virtual path only when:
		                !   - expected depth reduction >= 2x (depth_cluster_lb <= depth_match_lb/2)
		                !   - enough occupied blocks exist to keep the team busy
		                !   - enough work exists to amortize schedule + packing overhead
		                !
		                use_cluster_level = (depth_match_lb >= 4 .and. depth_cluster_lb <= depth_match_lb/2 .and. &
		                     & nocc_active >= Max (2, diagg2_threads/2) .and. nactive >= 8*diagg2_threads)

		                if (use_cluster_level) then
		                  !
		                  ! Occupied-centered multi-virtual schedule:
		                  !   - virtual indices are unique within each level
		                  !   - an occupied index may appear up to k times per level (k = diagg2_block_cap)
		                  ! This reduces depth for high-degree occupied LMOs by processing
		                  ! (occupied + up to k virtuals) inside one unit.
		                  !
		                  ! Eq. (L1) Level constraint:
		                  !   for all levels l:
		                  !     each i in V appears at most once in level l
		                  !     each j in O appears at most k times in level l
		                  !
		                  ! Worst case is one block-level per active edge, so allocate offsets as O(nactive).
		                  allocate (occ_block_id(nocc), occ_count(nocc), level_block_offset(nactive+2), &
		                       & block_occ(nactive), block_nedges(nactive), block_edges(block_cap, nactive), stat=alloc_stat)
		                  if (alloc_stat /= 0) then
		                    call mopend("Insufficient memory to run DIAGG2 (multi-virtual level schedule)")
		                  end if
		                  block_edges(:, :) = 0

		                  rem_edges(1:nactive) = active_edges(1:nactive)
		                  nrem = nactive
		                  max_level = 0
		                  nblocks = 0
		                  level_block_offset(1) = 1
		                  last_occ(:) = 0
		                  last_vir(:) = 0
		                  occ_count(:) = 0
		                  occ_block_id(:) = 0

		                  do while (nrem > 0)
		                    max_level = max_level + 1
		                    nnext = 0

		                    do iedge = 1, nrem
		                      ij = rem_edges(iedge)
		                      i = ifmo(1, ij)
		                      j = ifmo(2, ij)
		                      if (last_vir(i) == max_level) then
		                        nnext = nnext + 1
		                        next_edges(nnext) = ij
		                      else
		                        if (last_occ(j) /= max_level) then
		                          last_occ(j) = max_level
		                          occ_count(j) = 0
		                          occ_block_id(j) = 0
		                        end if
		                        if (occ_count(j) >= block_cap) then
		                          nnext = nnext + 1
		                          next_edges(nnext) = ij
		                        else
		                          last_vir(i) = max_level
		                          occ_count(j) = occ_count(j) + 1
		                          blk = occ_block_id(j)
		                          if (blk == 0) then
		                            nblocks = nblocks + 1
		                            blk = nblocks
		                            occ_block_id(j) = blk
		                            block_occ(blk) = j
		                            block_nedges(blk) = 0
		                          end if
		                          block_nedges(blk) = block_nedges(blk) + 1
		                          block_edges(block_nedges(blk), blk) = ij
		                        end if
		                      end if
		                    end do

		                    level_block_offset(max_level+1) = nblocks + 1
		                    if (nnext > 0) rem_edges(1:nnext) = next_edges(1:nnext)
		                    nrem = nnext
		                  end do

		                  sumb_local = 0.d0
		                  nrej_local = 0
		                  t_start = omp_get_wtime()
!$omp parallel num_threads(diagg2_threads) default(shared) private(lvl, tid, blk) &
!$omp& reduction(+:sumb_local, nrej_local)
		                  tid = omp_get_thread_num() + 1

		                  do lvl = 1, max_level
!$omp do schedule(dynamic,1)
		                    do blk = level_block_offset(lvl), level_block_offset(lvl+1) - 1
		                      call process_occ_block (block_occ(blk), block_nedges(blk), block_edges(:, blk), &
		                           & retry, tiny, biglim, sumb_local, nrej_local, &
		                           & iused_ws(:,tid), latoms_ws(:,tid), storei_ws(:,tid), storej_ws(:,tid), &
		                           & block_c_ws(:,:,tid), block_c2_ws(:,:,tid), block_atoms_ws(:,tid))
		                    end do
!$omp end do
		                  end do
!$omp end parallel

		                  deallocate (block_edges, block_nedges, block_occ, level_block_offset, occ_count, occ_block_id)
		                  elapsed = Max (1.d-12, omp_get_wtime() - t_start)
		                else
		                  !
		                  ! Eq. (D1): repeated maximal matchings M_l satisfy:
		                  !   e=(i,j) in M_l => i and j are unique within level l.
		                  ! This preserves conflict-free per-pair rotations in each level.
		                  !
		                  edge_level(:) = 0
		                  last_occ(:) = 0
		                  last_vir(:) = 0
		                  rem_edges(1:nactive) = active_edges(1:nactive)
		                  nrem = nactive
		                  max_level = 0
		                  do while (nrem > 0)
		                    max_level = max_level + 1
		                    nnext = 0
		                    do iedge = 1, nrem
		                      ij = rem_edges(iedge)
		                      i = ifmo(1, ij)
		                      j = ifmo(2, ij)
		                      if (last_occ(j) == max_level .or. last_vir(i) == max_level) then
		                        nnext = nnext + 1
		                        next_edges(nnext) = ij
		                      else
		                        edge_level(ij) = max_level
		                        last_occ(j) = max_level
		                        last_vir(i) = max_level
		                      end if
		                    end do
		                    if (nnext > 0) rem_edges(1:nnext) = next_edges(1:nnext)
		                    nrem = nnext
		                  end do

		                  if (allocated(level_count)) deallocate (level_count)
		                  if (allocated(level_offset)) deallocate (level_offset)
		                  if (allocated(level_pos)) deallocate (level_pos)
		                  if (allocated(level_edges)) deallocate (level_edges)
		                  allocate (level_count(max_level), level_offset(max_level+1), level_pos(max_level), &
		                       & level_edges(nactive), stat=alloc_stat)
		                  if (alloc_stat /= 0) then
		                    call mopend("Insufficient memory to run DIAGG2 (level schedule)")
		                  end if

		                  level_count(:) = 0
		                  do iedge = 1, nactive
		                    ij = active_edges(iedge)
		                    lvl = edge_level(ij)
		                    level_count(lvl) = level_count(lvl) + 1
		                  end do
		                  level_offset(1) = 1
		                  do lvl = 1, max_level
		                    level_offset(lvl+1) = level_offset(lvl) + level_count(lvl)
		                  end do
		                  level_pos(:) = level_offset(1:max_level)
		                  do iedge = 1, nactive
		                    ij = active_edges(iedge)
		                    lvl = edge_level(ij)
		                    pos = level_pos(lvl)
		                    level_edges(pos) = ij
		                    level_pos(lvl) = pos + 1
		                  end do

		                  sumb_local = 0.d0
		                  nrej_local = 0
		                  t_start = omp_get_wtime()
!$omp parallel num_threads(diagg2_threads) default(shared) private(ij, lvl, pos, tid) &
!$omp& reduction(+:sumb_local, nrej_local)
		                  tid = omp_get_thread_num() + 1

		                  do lvl = 1, max_level
!$omp do schedule(dynamic,16)
		                    do pos = level_offset(lvl), level_offset(lvl+1) - 1
		                      ij = level_edges(pos)
		                      call process_pair (ij, retry, tiny, biglim, sumb_local, nrej_local, &
		                           & iused_ws(:,tid), latoms_ws(:,tid), storei_ws(:,tid), storej_ws(:,tid))
		                    end do
!$omp end do
		                  end do
!$omp end parallel
		                  elapsed = Max (1.d-12, omp_get_wtime() - t_start)
		                end if

		                ! Eq. (1) reused for level-sweep path.
		                edge_cost = elapsed / Dble (Max (1, nactive))
		                if (tune_level_samples(tune_slot) == 0) then
		                  tune_level_cost(tune_slot) = edge_cost
		                else
		                  ! Eq. (2) EWMA update for level mode.
		                  tune_level_cost(tune_slot) = (1.d0-tune_alpha)*tune_level_cost(tune_slot) + tune_alpha*edge_cost
		                end if
		                tune_level_samples(tune_slot) = tune_level_samples(tune_slot) + 1

		                sumb = sumb_local
		                nrej = nrej_local
		              end if
			            else
			              use_parallel_diagg2 = .false.
			            end if

		            if (allocated(level_edges)) deallocate (level_edges)
		            if (allocated(level_pos)) deallocate (level_pos)
		            if (allocated(level_offset)) deallocate (level_offset)
		            if (allocated(level_count)) deallocate (level_count)
	          end if
	        end if

	        deallocate (edge_level, active_edges, rem_edges, next_edges, last_occ, last_vir)
	      end if
#endif

    if (.not. use_parallel_diagg2) then
      outer_loop: do ij = 1, nij
        i = ifmo(1, ij)
        j = ifmo(2, ij)
	        if (Abs (fmo(ij)) >= tiny) then
	          c = fmo(ij) * const
	          d = eigs(j) - eigv(i) - shift
	          if (biglim >= 0.d0) then
	            if (d == 0.d0) then
	              if (Abs (c) <= 0.d0) cycle outer_loop
	            else if (Abs (c/d) < biglim) then
	              cycle outer_loop
	            end if
	          end if
	            ncfj = ncf(j)
	            ncei = nce(i)
	        !
	        !  STORE LMOS FOR POSSIBLE REJECTION, IF LMOS EXPAND TOO MUCH.
        !
            jlr = ncocc(j) + 1
            if (j /= nocc) then
              jur = ncocc(j+1)
              jncf = nncf(j+1)
            else
              jur = cocc_dim
              jncf = icocc_dim
            end if
            jur = Min (jlr+norbs-1, jur)
            ilr = ncvir(i) + 1
            if (i /= nvir) then
              iur = ncvir(i+1)
              incv = nnce(i+1)
	            else
	              iur = cvir_dim
	              incv = icvir_dim
	            end if
	            iur = Min (ilr+norbs-1, iur)
	            !
	            ! Eq. (E2): rollback backups need only the currently populated
	            ! coefficient spans (not the full reserved windows).
	            !   n_j = Sum_{a in supp(occ_j)} n_AO(a)
	            !   n_i = Sum_{a in supp(vir_i)} n_AO(a)
	            !
	            njlen = 0
	            do lf = nncf(j) + 1, nncf(j) + ncf(j)
	              njlen = njlen + iorbs(icocc(lf))
	            end do
	            nilen = 0
	            do le = nnce(i) + 1, nnce(i) + nce(i)
	              nilen = nilen + iorbs(icvir(le))
	            end do
	            storej(1:njlen) = cocc(jlr:jlr+njlen-1)
	            storei(1:nilen) = cvir(ilr:ilr+nilen-1)
            !
            !   STORAGE DONE.
            !
            lij = lij + 1
            e = Sign (Sqrt(4.d0*c*c+d*d), d)
            alpha = Sqrt (0.5d0*(1.d0+d/e))
            do
              beta = -Sign (Sqrt(1.d0-alpha*alpha), c)
              sumb = sumb + Abs (beta)
              !
              ! IDENTIFY THE ATOMS IN THE OCCUPIED LMO.  ATOMS NOT USED ARE
              ! FLAGGED BY '-1' IN IUSED.
              !
              mlf = 0
              !
              do lf = nncf(j) + 1, nncf(j) + ncf(j)
                ii = icocc(lf)
                iused(ii) = mlf
                mlf = mlf + iorbs(ii)
              end do
              loopi = ncvir(i)
              loopj = ncocc(j)
              mle = 0
           !
           !      ROTATION OF PSEUDO-EIGENVECTORS
           !
              do le = nnce(i) + 1, nnce(i) + nce(i)
                mie = icvir(le)

                latoms(mie) = .true.
                mlff = iused(mie) + loopj
                if (iused(mie) >= 0) then
                  !
                  !  TWO BY TWO ROTATION OF ATOMS WHICH ARE COMMON
                  !  TO OCCUPIED LMO J AND VIRTUAL LMO I
                  !
                  do mlee = mle + 1 + loopi, mle + iorbs(mie) + loopi
                    mlff = mlff + 1
                    a = cocc(mlff)
                    b = cvir(mlee)
                    cocc(mlff) = alpha * a + beta * b
                    cvir(mlee) = alpha * b - beta * a
                  end do
                else
                  !
                  !   FILLED  LMO ATOM 'MIE' DOES NOT EXIST.
                  !   CHECK IF IT SHOULD EXIST
                  !
                  sum = 0.d0
                  do mlee = mle + 1 + loopi, mle + iorbs(mie) + loopi
                    sum = sum + (beta*cvir(mlee)) ** 2
                  end do
                  if (sum > thresh) then
                    !
                    if (nncf(j)+ncf(j) >= jncf) go to 1000
                    if (mlf+iorbs(mie)+loopj > jur) go to 1000
                    !
                    !  YES, OCCUPIED LMO ATOM 'MIE' SHOULD EXIST
                    !
                    ncf(j) = ncf(j) + 1
                    icocc(nncf(j)+ncf(j)) = mie
                    !
                    iused(mie) = mlf
                    mlf = mlf + iorbs(mie)
                    !
                    !   PUT INTENSITY INTO OCCUPIED LMO ATOM 'MIE'
                    !
                    mlff = iused(mie) + loopj
                    do mlee = mle + 1 + loopi, mle + iorbs(mie) + loopi
                      mlff = mlff + 1
                      cocc(mlff) = beta * cvir(mlee)
                      cvir(mlee) = alpha * cvir(mlee)
                    end do
                  end if
                end if
                mle = mle + iorbs(mie)
              end do
              !
              !  NOW CHECK ALL ATOMS WHICH WERE IN THE OCCUPIED LMO
              !  WHICH ARE NOT IN THE VIRTUAL LMO, TO SEE IF THEY
              !  SHOULD BE IN THE VIRTUAL LMO.
              !
              do lf = nncf(j) + 1, nncf(j) + ncf(j)
                ii = icocc(lf)

                if ( .not. latoms(ii)) then
                  sum = 0.d0
                  do mlff = iused(ii) + loopj + 1, iused(ii) + loopj + &
                       & iorbs(ii)
                    sum = sum + (beta*cocc(mlff)) ** 2
                  end do
                  if (sum > thresh) then
                    if (nnce(i)+nce(i) >= incv) go to 1000
                    if (mle+iorbs(ii)+loopi > iur) go to 1000
                    !
                    !  YES, VIRTUAL  LMO ATOM 'II' SHOULD EXIST
                    !
                    nce(i) = nce(i) + 1
                    icvir(nnce(i)+nce(i)) = ii
                    latoms(ii) = .true.
                    !
                    !   PUT INTENSITY INTO VIRTUAL  LMO ATOM 'II'
                    !
                    mlff = iused(ii) + loopj
                    do mlee = mle + 1 + loopi, mle + iorbs(ii) + loopi
                      mlff = mlff + 1
                       !
                      cvir(mlee) = -beta * cocc(mlff)
                      cocc(mlff) = alpha * cocc(mlff)
                    end do
                    mle = mle + iorbs(ii)
                  end if
                end if
              end do
            exit
        1000  continue
              nrej = nrej + 1
                !
                !   THE ARRAY BOUNDS WERE GOING TO BE EXCEEDED.
                !   TO PREVENT THIS, RESET THE LMOS.
                !
	              ! Previous-commit improvement kept here: contiguous slice rollback replaces scalar loops.
	              cocc(jlr:jlr+njlen-1) = storej(1:njlen)
	              cvir(ilr:ilr+nilen-1) = storei(1:nilen)
	              do le = nnce(i) + 1, nnce(i) + nce(i)
	                latoms(icvir(le)) = .false.
	              end do
	              do lf = nncf(j) + 1, nncf(j) + ncf(j)
	                iused(icocc(lf)) = -1
	              end do
	              ncf(j) = ncfj
	              nce(i) = ncei
              if (retry) then
                  !
                  !   HALF THE ROTATION ANGLE.  WILL THIS PREVENT THE
                  !   ARRAY BOUND FROM BEING EXCEEDED?
                  !
                alpha = 0.5d0 * (alpha+1.d0)
              else
                cycle outer_loop
              end if
            end do
          !
          !  RESET COUNTERS WHICH HAVE BEEN SET.
          !
            do le = nnce(i) + 1, nnce(i) + nce(i)
              mie = icvir(le)
              latoms(mie) = .false.
            end do
          !
	            do lf = nncf(j) + 1, nncf(j) + ncf(j)
	              iused(icocc(lf)) = -1
	            end do
	        end if
	      end do outer_loop
	    end if
    nrejct(2) = nrejct(1)
    nrejct(1) = nrej
    if (times) then
      call timer (" AFTER DIAGG2 IN ITER")
    end if
  contains
	  subroutine process_occ_block (j_occ, nedges, edges, retry, tiny, biglim, sumb_acc, nrej_acc, &
	       & iused_t, latoms_t, storei_t, storej_t, cwork, cwork2, atomwork)
    !***********************************************************************
    !
	    ! Occupied-centered multi-virtual DIAGG2 update for one occupied LMO j_occ and
	    ! up to k virtual neighbors (nedges <= diagg2_block_cap) inside one DIAGG2 level.
	    !
	    ! Motivation:
	    ! - High-degree occupied LMOs can force many levels when each level allows
	    !   only one (i,j) per occupied. Processing (j + up to k virtuals) in one unit
	    !   reduces barrier depth for such cases.
	    !
	    ! Math (same per-pair Jacobi rotation as process_pair):
	    ! - For each pair (j_occ, i_p), apply:
	    !     a' = alpha_p*a + beta_p*b
	    !     b' = alpha_p*b - beta_p*a
	    !   where a = occ(j_occ), b = vir(i_p). (Eq. matches process_pair.)
	    !
	    ! Pruning/expansion semantics (MOZYME):
	    ! - Only intersection atoms are rotated.
	    ! - For atoms present in only one LMO, apply the rotation (and insert the atom)
	    !   iff:
	    !     sum = Sum_AO (beta_p*coeff_AO)^2 = beta_p^2 * ||coeff_atom||_2^2 > thresh.   (Eq. P1)
	    !   Otherwise leave that atom's coefficients untouched (no alpha scaling), matching process_pair.
	    !
	    ! Safety:
	    ! - All updates are computed in scratch first. If reserved storage would be exceeded,
	    !   fall back to calling process_pair for each edge (original behavior).
	    !***********************************************************************
    integer, intent (in) :: j_occ, nedges
    integer, dimension (:), intent (in) :: edges
    logical, intent (in) :: retry
    double precision, intent (in) :: tiny, biglim
    double precision, intent (inout) :: sumb_acc
	    integer, intent (inout) :: nrej_acc
	    integer, dimension (numat), intent (inout) :: iused_t
	    logical, dimension (numat), intent (inout) :: latoms_t
	    double precision, dimension (:), intent (inout) :: storei_t, storej_t
		    double precision, dimension (:,:), intent (inout) :: cwork, cwork2
		    integer, dimension (:), intent (inout) :: atomwork

		    integer, parameter :: mask_base = 32
		    integer, parameter :: enc_base = mask_base*mask_base
		    integer :: p, q
		    integer :: ij, i_vir, col, ncols
		    integer :: lf, le, atom, len, row
		    integer :: base, pos
		    integer :: nunion, m
		    integer :: jlr, jur, jncf
		    integer :: ilr, iur, incv
		    integer :: add_atoms, add_len
		    integer :: ncf_new, nce_new
		    integer :: attempt, max_attempt
		    integer :: ao
		    integer :: enc, off, mask_orig, mask_now
		    integer :: ncf_old
		    integer :: nce_old(diagg2_block_cap)
		    logical :: ok
		    logical :: occ_has, vir_has
		    double precision :: a, b, c, d, e, alpha, beta, beta2
		    double precision :: norm2
		    double precision :: alpha_vec(diagg2_block_cap), beta_vec(diagg2_block_cap), c_vec(diagg2_block_cap)
		    integer :: i_list(diagg2_block_cap)

	    if (nedges <= 0) return
	    if (nedges == 1) then
	      call process_pair (edges(1), retry, tiny, biglim, sumb_acc, nrej_acc, iused_t, latoms_t, storei_t, storej_t)
	      return
	    end if

	    do p = 1, nedges
	      ij = edges(p)
	      i_list(p) = ifmo(1, ij)
	      if (ifmo(2, ij) /= j_occ) then
	        ! Defensive fallback if schedule bookkeeping is inconsistent.
	        goto 9000
	      end if
	    end do

		    ! Build union atom list and per-atom presence masks.
		    ! iused_t(atom) temporarily stores a bitmask of which columns contain the atom:
		    !   bit 0: occupied (j_occ), bit p: virtual i_list(p).
		    nunion = 0
		    do lf = nncf(j_occ) + 1, nncf(j_occ) + ncf(j_occ)
		      atom = icocc(lf)
		      if (iused_t(atom) < 0) then
		        nunion = nunion + 1
		        atomwork(nunion) = atom
		        iused_t(atom) = 1
		      else
		        iused_t(atom) = Ior (iused_t(atom), 1)
		      end if
		    end do
		    do p = 1, nedges
		      i_vir = i_list(p)
		      do le = nnce(i_vir) + 1, nnce(i_vir) + nce(i_vir)
		        atom = icvir(le)
		        if (iused_t(atom) < 0) then
		          nunion = nunion + 1
		          atomwork(nunion) = atom
		          iused_t(atom) = Ishft (1, p)
		        else
		          iused_t(atom) = Ior (iused_t(atom), Ishft (1, p))
		        end if
		      end do
		    end do

		    ! Encode AO row offsets and masks into iused_t(atom):
		    !   iused_t(atom) = off*enc_base + mask_orig*mask_base + mask_now,
		    ! where off is 0-based AO offset, mask_now starts as mask_orig and can change on insertion.
		    !
		    m = 0
		    do q = 1, nunion
		      atom = atomwork(q)
		      len = iorbs(atom)
		      mask_orig = iused_t(atom)
		      iused_t(atom) = m*enc_base + mask_orig*mask_base + mask_orig
		      m = m + len
		    end do

		    ncols = 1 + nedges
		    if (m <= 0) goto 8000
		    if (m > size(cwork,1) .or. ncols > size(cwork,2)) then
	      goto 9000
	    end if

	    cwork(1:m, 1:ncols) = 0.d0

	    ! Pack occupied column.
	    base = ncocc(j_occ)
	    pos = 0
		    do lf = nncf(j_occ) + 1, nncf(j_occ) + ncf(j_occ)
		      atom = icocc(lf)
		      len = iorbs(atom)
		      row = (iused_t(atom) / enc_base) + 1
		      cwork(row:row+len-1, 1) = cocc(base+pos+1:base+pos+len)
		      pos = pos + len
		    end do

	    ! Pack each virtual column.
	    do p = 1, nedges
	      i_vir = i_list(p)
	      col = p + 1
	      base = ncvir(i_vir)
	      pos = 0
		      do le = nnce(i_vir) + 1, nnce(i_vir) + nce(i_vir)
		        atom = icvir(le)
		        len = iorbs(atom)
		        row = (iused_t(atom) / enc_base) + 1
		        cwork(row:row+len-1, col) = cvir(base+pos+1:base+pos+len)
		        pos = pos + len
		      end do
		    end do

		    ! Backup packed coefficients for retry damping (attempt 2).
		    cwork2(1:m, 1:ncols) = cwork(1:m, 1:ncols)

	    ! Compute initial Jacobi rotation parameters for each (j_occ, i_list(p)).
	    do p = 1, nedges
	      ij = edges(p)
	      i_vir = i_list(p)
	      if (tiny >= 0.d0) then
	        if (Abs (fmo(ij)) < tiny) then
	          alpha_vec(p) = 1.d0
	          beta_vec(p) = 0.d0
	          c_vec(p) = 0.d0
	          cycle
	        end if
	      end if
	      c = fmo(ij) * const
	      d = eigs(j_occ) - eigv(i_vir) - shift
	      if (biglim >= 0.d0) then
	        if (d == 0.d0) then
	          if (Abs (c) <= 0.d0) then
	            alpha_vec(p) = 1.d0
	            beta_vec(p) = 0.d0
	            c_vec(p) = c
	            cycle
	          end if
	        else if (Abs (c/d) < biglim) then
	          alpha_vec(p) = 1.d0
	          beta_vec(p) = 0.d0
	          c_vec(p) = c
	          cycle
	        end if
	      end if
	      e = Sign (Sqrt (Max (0.d0, 4.d0*c*c + d*d)), d)
	      if (e == 0.d0) then
	        alpha = 1.d0
	        beta = 0.d0
	      else
	        alpha = Sqrt (0.5d0*(1.d0 + d/e))
	        beta = -Sign (Sqrt (Max (0.d0, 1.d0 - alpha*alpha)), c)
	      end if
	      alpha_vec(p) = alpha
	      beta_vec(p) = beta
	      c_vec(p) = c
		    end do

		    ncf_old = ncf(j_occ)
		    do p = 1, nedges
		      nce_old(p) = nce(i_list(p))
		    end do

		    max_attempt = 1
		    if (retry) max_attempt = 2
		    ok = .false.
		    do attempt = 1, max_attempt
		      ! Reset scratch state for this attempt.
		      cwork(1:m, 1:ncols) = cwork2(1:m, 1:ncols)
		      do q = 1, nunion
		        atom = atomwork(q)
		        enc = iused_t(atom)
		        off = enc / enc_base
		        mask_orig = Mod (enc / mask_base, mask_base)
		        iused_t(atom) = off*enc_base + mask_orig*mask_base + mask_orig
		      end do

		      !
		      ! Apply the same sequence of Jacobi rotations as process_pair, but with
		      ! MOZYME pruning semantics enforced per atom (Eq. P1 insertion rule).
		      !
		      do q = 1, nunion
		        atom = atomwork(q)
		        enc = iused_t(atom)
		        off = enc / enc_base
		        mask_orig = Mod (enc / mask_base, mask_base)
		        mask_now = Mod (enc, mask_base)
		        row = off + 1
		        len = iorbs(atom)

		        do p = 1, nedges
		          alpha = alpha_vec(p)
		          beta = beta_vec(p)
		          if (beta == 0.d0) cycle
		          col = p + 1
		          occ_has = Btest (mask_now, 0)
		          vir_has = Btest (mask_now, p)

		          if (occ_has .and. vir_has) then
		            do ao = 0, len - 1
		              a = cwork(row+ao, 1)
		              b = cwork(row+ao, col)
		              cwork(row+ao, 1) = alpha*a + beta*b
		              cwork(row+ao, col) = alpha*b - beta*a
		            end do
		          else if ((.not. occ_has) .and. vir_has) then
		            beta2 = beta*beta
		            norm2 = 0.d0
		            do ao = 0, len - 1
		              b = cwork(row+ao, col)
		              norm2 = norm2 + b*b
		            end do
		            if (beta2*norm2 > thresh) then
		              mask_now = Ibset (mask_now, 0)
		              do ao = 0, len - 1
		                b = cwork(row+ao, col)
		                cwork(row+ao, 1) = beta*b
		                cwork(row+ao, col) = alpha*b
		              end do
		            end if
		          else if (occ_has .and. (.not. vir_has)) then
		            beta2 = beta*beta
		            norm2 = 0.d0
		            do ao = 0, len - 1
		              a = cwork(row+ao, 1)
		              norm2 = norm2 + a*a
		            end do
		            if (beta2*norm2 > thresh) then
		              mask_now = Ibset (mask_now, p)
		              do ao = 0, len - 1
		                a = cwork(row+ao, 1)
		                cwork(row+ao, col) = -beta*a
		                cwork(row+ao, 1) = alpha*a
		              end do
		            end if
		          end if
		        end do
		        iused_t(atom) = off*enc_base + mask_orig*mask_base + mask_now
		      end do

		      !
		      ! Capacity check for occupied + virtuals after applying insertion decisions.
		      !
		      ok = .true.

		      ! Occupied bounds (same logic as process_pair).
	      jlr = ncocc(j_occ) + 1
	      if (j_occ /= nocc) then
	        jur = ncocc(j_occ+1)
	        jncf = nncf(j_occ+1)
	      else
	        jur = cocc_dim
	        jncf = icocc_dim
	      end if
		      jur = Min (jlr+norbs-1, jur)

		      add_atoms = 0
		      add_len = 0
		      do lf = nncf(j_occ) + 1, nncf(j_occ) + ncf_old
		        latoms_t(icocc(lf)) = .true.
		      end do
		      do q = 1, nunion
		        atom = atomwork(q)
		        if (.not. latoms_t(atom)) then
		          mask_now = Mod (iused_t(atom), mask_base)
		          if (Btest(mask_now, 0)) then
		            add_atoms = add_atoms + 1
		            add_len = add_len + iorbs(atom)
		          end if
		        end if
		      end do
		      do lf = nncf(j_occ) + 1, nncf(j_occ) + ncf_old
		        latoms_t(icocc(lf)) = .false.
		      end do
			      ncf_new = ncf_old + add_atoms
			      if (nncf(j_occ) + ncf_new > jncf) ok = .false.

			      ! Coefficient capacity: current occupied used length is sum iorbs over its atom list.
			      pos = 0
			      do lf = nncf(j_occ) + 1, nncf(j_occ) + ncf_old
		        pos = pos + iorbs(icocc(lf))
		      end do
		      if (ncocc(j_occ) + pos + add_len > jur) ok = .false.

		      ! Virtual bounds.
		      do p = 1, nedges
		        i_vir = i_list(p)
	        ilr = ncvir(i_vir) + 1
	        if (i_vir /= nvir) then
	          iur = ncvir(i_vir+1)
	          incv = nnce(i_vir+1)
	        else
	          iur = cvir_dim
	          incv = icvir_dim
	        end if
	        iur = Min (ilr+norbs-1, iur)

		        add_atoms = 0
		        add_len = 0
		        do le = nnce(i_vir) + 1, nnce(i_vir) + nce_old(p)
		          latoms_t(icvir(le)) = .true.
		        end do
		        do q = 1, nunion
		          atom = atomwork(q)
		          if (.not. latoms_t(atom)) then
		            mask_now = Mod (iused_t(atom), mask_base)
		            if (Btest(mask_now, p)) then
		              add_atoms = add_atoms + 1
		              add_len = add_len + iorbs(atom)
		            end if
		          end if
		        end do
		        do le = nnce(i_vir) + 1, nnce(i_vir) + nce_old(p)
		          latoms_t(icvir(le)) = .false.
		        end do
		        nce_new = nce_old(p) + add_atoms
		        if (nnce(i_vir) + nce_new > incv) ok = .false.
		        pos = 0
		        do le = nnce(i_vir) + 1, nnce(i_vir) + nce_old(p)
		          pos = pos + iorbs(icvir(le))
		        end do
		        if (ncvir(i_vir) + pos + add_len > iur) ok = .false.
		      end do

		      if (ok) exit

		      if (attempt < max_attempt) then
		        ! Retry by damping each rotation angle (same rule as process_pair retry path).
		        do p = 1, nedges
		          alpha_vec(p) = 0.5d0 * (alpha_vec(p) + 1.d0)
		          beta_vec(p) = -Sign (Sqrt (Max (0.d0, 1.d0 - alpha_vec(p)*alpha_vec(p))), c_vec(p))
		        end do
		      end if
		    end do

	    if (.not. ok) then
	      goto 9000
	    end if

	    ! Accumulate the same rotation measure used in process_pair.
	    do p = 1, nedges
	      sumb_acc = sumb_acc + Abs (beta_vec(p))
	    end do

		    !
		    ! Write back occupied coefficients, preserving atom order and appending new atoms by union order.
		    !
		    base = ncocc(j_occ)
		    pos = 0
		    do lf = nncf(j_occ) + 1, nncf(j_occ) + ncf_old
		      atom = icocc(lf)
		      len = iorbs(atom)
		      row = (iused_t(atom) / enc_base) + 1
		      cocc(base+pos+1:base+pos+len) = cwork(row:row+len-1, 1)
		      pos = pos + len
		      latoms_t(atom) = .true.
		    end do
		    do q = 1, nunion
		      atom = atomwork(q)
		      if (.not. latoms_t(atom)) then
		        mask_now = Mod (iused_t(atom), mask_base)
		        if (Btest(mask_now, 0)) then
		          ncf(j_occ) = ncf(j_occ) + 1
		          icocc(nncf(j_occ)+ncf(j_occ)) = atom
		          len = iorbs(atom)
		          row = (iused_t(atom) / enc_base) + 1
		          cocc(base+pos+1:base+pos+len) = cwork(row:row+len-1, 1)
		          pos = pos + len
		          latoms_t(atom) = .true.
		        end if
		      end if
		    end do
		    do q = 1, nunion
		      latoms_t(atomwork(q)) = .false.
		    end do

	    !
	    ! Write back each virtual orbital.
	    !
		    do p = 1, nedges
		      i_vir = i_list(p)
		      col = p + 1
		      base = ncvir(i_vir)
		      pos = 0
		      do le = nnce(i_vir) + 1, nnce(i_vir) + nce_old(p)
		        atom = icvir(le)
		        len = iorbs(atom)
		        row = (iused_t(atom) / enc_base) + 1
		        cvir(base+pos+1:base+pos+len) = cwork(row:row+len-1, col)
		        pos = pos + len
		        latoms_t(atom) = .true.
		      end do
		      do q = 1, nunion
		        atom = atomwork(q)
		        if (.not. latoms_t(atom)) then
		          mask_now = Mod (iused_t(atom), mask_base)
		          if (Btest(mask_now, p)) then
		            nce(i_vir) = nce(i_vir) + 1
		            icvir(nnce(i_vir)+nce(i_vir)) = atom
		            len = iorbs(atom)
		            row = (iused_t(atom) / enc_base) + 1
		            cvir(base+pos+1:base+pos+len) = cwork(row:row+len-1, col)
		            pos = pos + len
		            latoms_t(atom) = .true.
		          end if
		        end if
		      end do
		      do q = 1, nunion
		        latoms_t(atomwork(q)) = .false.
	      end do
	    end do

8000  continue
	    ! Reset union markers (required before returning to scheduler).
	    do q = 1, nunion
	      iused_t(atomwork(q)) = -1
	      latoms_t(atomwork(q)) = .false.
	    end do
	    return

9000  continue
	    ! Fall back to original per-pair kernel (includes retry/rollback behavior).
	    do q = 1, nunion
	      iused_t(atomwork(q)) = -1
	      latoms_t(atomwork(q)) = .false.
	    end do
	    do p = 1, nedges
	      call process_pair (edges(p), retry, tiny, biglim, sumb_acc, nrej_acc, iused_t, latoms_t, storei_t, storej_t)
	    end do
	    return

	  end subroutine process_occ_block

	  subroutine process_pair (ij, retry, tiny, biglim, sumb_acc, nrej_acc, iused_t, latoms_t, storei_t, storej_t)
    integer, intent (in) :: ij
    logical, intent (in) :: retry
    double precision, intent (in) :: tiny, biglim
    double precision, intent (inout) :: sumb_acc
	    integer, intent (inout) :: nrej_acc
	    integer, dimension (numat), intent (inout) :: iused_t
	    logical, dimension (numat), intent (inout) :: latoms_t
	    double precision, dimension (:), intent (inout) :: storei_t, storej_t

    integer :: i, j, ii, jur, l, ilr, incv, iur, jlr, jncf, k, le, lf, loopi, loopj
    integer :: mie, mle, mlee, mlf, mlff, ncei, ncfj
    integer :: njlen, nilen
    double precision :: a, alpha, b, beta, c, d, e, sum

	    i = ifmo(1, ij)
	    j = ifmo(2, ij)
	    if (tiny >= 0.d0) then
	      if (Abs (fmo(ij)) < tiny) return
	    end if

	    c = fmo(ij) * const
	    d = eigs(j) - eigv(i) - shift
	    !
	    ! Eq. (E1): pair is active when |c/d| >= biglim (or |c|>0 when d=0).
	    ! Guarding d=0 keeps scheduler and kernel predicates identical.
	    !
	    if (biglim >= 0.d0) then
	      if (d == 0.d0) then
	        if (Abs (c) <= 0.d0) return
	      else if (Abs (c/d) < biglim) then
	        return
	      end if
	    end if

    ncfj = ncf(j)
    ncei = nce(i)
    !
    !  STORE LMOS FOR POSSIBLE REJECTION, IF LMOS EXPAND TOO MUCH.
    !
    jlr = ncocc(j) + 1
    if (j /= nocc) then
      jur = ncocc(j+1)
      jncf = nncf(j+1)
    else
      jur = cocc_dim
      jncf = icocc_dim
    end if
    jur = Min (jlr+norbs-1, jur)
    ilr = ncvir(i) + 1
    if (i /= nvir) then
      iur = ncvir(i+1)
      incv = nnce(i+1)
    else
      iur = cvir_dim
      incv = icvir_dim
    end if
	    iur = Min (ilr+norbs-1, iur)
	    !
	    ! Eq. (E2): backup sizes are active LMO support lengths:
	    !   n_j = Sum_{a in supp(occ_j)} n_AO(a), n_i = Sum_{a in supp(vir_i)} n_AO(a)
	    ! This lowers rollback copy traffic while keeping the same pre-rotation
	    ! state for all currently populated coefficients.
	    !
	    njlen = 0
	    do lf = nncf(j) + 1, nncf(j) + ncf(j)
	      njlen = njlen + iorbs(icocc(lf))
	    end do
	    nilen = 0
	    do le = nnce(i) + 1, nnce(i) + nce(i)
	      nilen = nilen + iorbs(icvir(le))
	    end do
	    !
	    ! Safety guard: required backup lengths must satisfy
	    !   njlen <= size(storej_t), nilen <= size(storei_t)
	    ! to avoid out-of-bounds in rollback buffers.
	    !
	    if (njlen > size(storej_t) .or. nilen > size(storei_t)) then
	      nrej_acc = nrej_acc + 1
	      return
	    end if
	    ! Previous-commit improvement kept here: contiguous backup copy for rollback path.
	    storej_t(1:njlen) = cocc(jlr:jlr+njlen-1)
	    storei_t(1:nilen) = cvir(ilr:ilr+nilen-1)
    !
    !   STORAGE DONE.
    !
    e = Sign (Sqrt(4.d0*c*c+d*d), d)
    alpha = Sqrt (0.5d0*(1.d0+d/e))

    do
      beta = -Sign (Sqrt(1.d0-alpha*alpha), c)
      sumb_acc = sumb_acc + Abs (beta)
      !
      ! IDENTIFY THE ATOMS IN THE OCCUPIED LMO.  ATOMS NOT USED ARE
      ! FLAGGED BY '-1' IN IUSED.
      !
      mlf = 0
      do lf = nncf(j) + 1, nncf(j) + ncf(j)
        ii = icocc(lf)
        iused_t(ii) = mlf
        mlf = mlf + iorbs(ii)
      end do
      loopi = ncvir(i)
      loopj = ncocc(j)
      mle = 0
      !
      !      ROTATION OF PSEUDO-EIGENVECTORS
      !
      do le = nnce(i) + 1, nnce(i) + nce(i)
        mie = icvir(le)
        latoms_t(mie) = .true.
        mlff = iused_t(mie) + loopj
        if (iused_t(mie) >= 0) then
          !
          !  TWO BY TWO ROTATION OF ATOMS WHICH ARE COMMON
          !  TO OCCUPIED LMO J AND VIRTUAL LMO I
          !
          do mlee = mle + 1 + loopi, mle + iorbs(mie) + loopi
            mlff = mlff + 1
            a = cocc(mlff)
            b = cvir(mlee)
            cocc(mlff) = alpha * a + beta * b
            cvir(mlee) = alpha * b - beta * a
          end do
        else
          !
          !   FILLED  LMO ATOM 'MIE' DOES NOT EXIST.
          !   CHECK IF IT SHOULD EXIST
          !
          sum = 0.d0
          do mlee = mle + 1 + loopi, mle + iorbs(mie) + loopi
            sum = sum + (beta*cvir(mlee)) ** 2
          end do
          if (sum > thresh) then
            if (nncf(j)+ncf(j) >= jncf) exit
            if (mlf+iorbs(mie)+loopj > jur) exit
            !
            !  YES, OCCUPIED LMO ATOM 'MIE' SHOULD EXIST
            !
            ncf(j) = ncf(j) + 1
            icocc(nncf(j)+ncf(j)) = mie
            !
            iused_t(mie) = mlf
            mlf = mlf + iorbs(mie)
            !
            !   PUT INTENSITY INTO OCCUPIED LMO ATOM 'MIE'
            !
            mlff = iused_t(mie) + loopj
            do mlee = mle + 1 + loopi, mle + iorbs(mie) + loopi
              mlff = mlff + 1
              cocc(mlff) = beta * cvir(mlee)
              cvir(mlee) = alpha * cvir(mlee)
            end do
          end if
        end if
        mle = mle + iorbs(mie)
      end do

      if (le <= nnce(i) + nce(i)) then
        ! Rejected due to array bounds in occupied expansion.
        nrej_acc = nrej_acc + 1
        ! Contiguous restore (same values as pre-rotation state).
	        cocc(jlr:jlr+njlen-1) = storej_t(1:njlen)
	        cvir(ilr:ilr+nilen-1) = storei_t(1:nilen)
	        do le = nnce(i) + 1, nnce(i) + nce(i)
	          latoms_t(icvir(le)) = .false.
	        end do
	        do lf = nncf(j) + 1, nncf(j) + ncf(j)
	          iused_t(icocc(lf)) = -1
	        end do
	        ncf(j) = ncfj
	        nce(i) = ncei
        if (retry) then
          alpha = 0.5d0 * (alpha+1.d0)
          cycle
        end if
        return
      end if

      !
      !  NOW CHECK ALL ATOMS WHICH WERE IN THE OCCUPIED LMO
      !  WHICH ARE NOT IN THE VIRTUAL LMO, TO SEE IF THEY
      !  SHOULD BE IN THE VIRTUAL LMO.
      !
      do lf = nncf(j) + 1, nncf(j) + ncf(j)
        ii = icocc(lf)
        if ( .not. latoms_t(ii)) then
          sum = 0.d0
          do mlff = iused_t(ii) + loopj + 1, iused_t(ii) + loopj + iorbs(ii)
            sum = sum + (beta*cocc(mlff)) ** 2
          end do
          if (sum > thresh) then
            if (nnce(i)+nce(i) >= incv) exit
            if (mle+iorbs(ii)+loopi > iur) exit
            !
            !  YES, VIRTUAL  LMO ATOM 'II' SHOULD EXIST
            !
            nce(i) = nce(i) + 1
            icvir(nnce(i)+nce(i)) = ii
            latoms_t(ii) = .true.
            !
            !   PUT INTENSITY INTO VIRTUAL  LMO ATOM 'II'
            !
            mlff = iused_t(ii) + loopj
            do mlee = mle + 1 + loopi, mle + iorbs(ii) + loopi
              mlff = mlff + 1
              cvir(mlee) = -beta * cocc(mlff)
              cocc(mlff) = alpha * cocc(mlff)
            end do
            mle = mle + iorbs(ii)
          end if
        end if
      end do

      if (lf <= nncf(j) + ncf(j)) then
        ! Rejected due to array bounds in virtual expansion.
        nrej_acc = nrej_acc + 1
        ! Contiguous restore (same values as pre-rotation state).
	        cocc(jlr:jlr+njlen-1) = storej_t(1:njlen)
	        cvir(ilr:ilr+nilen-1) = storei_t(1:nilen)
	        do le = nnce(i) + 1, nnce(i) + nce(i)
	          latoms_t(icvir(le)) = .false.
	        end do
	        do lf = nncf(j) + 1, nncf(j) + ncf(j)
	          iused_t(icocc(lf)) = -1
	        end do
	        ncf(j) = ncfj
	        nce(i) = ncei
        if (retry) then
          alpha = 0.5d0 * (alpha+1.d0)
          cycle
        end if
        return
      end if

      exit
    end do

    ! Reset counters which have been set.
    do le = nnce(i) + 1, nnce(i) + nce(i)
      mie = icvir(le)
      latoms_t(mie) = .false.
    end do
    do lf = nncf(j) + 1, nncf(j) + ncf(j)
      iused_t(icocc(lf)) = -1
    end do

  end subroutine process_pair
  end subroutine diagg2
