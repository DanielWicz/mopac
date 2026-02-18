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
#ifdef _OPENMP
    integer, parameter :: mode_level = 1, mode_task = 2
    double precision, parameter :: tune_alpha = 0.25d0, tune_switch_margin = 1.02d0
    integer, parameter :: tune_explore_period = 24
    integer :: alloc_stat, max_threads, tid, diagg2_threads, width_hint
    integer :: lvl, max_level, nactive, pos
    integer :: batch_start, batch_end, batch_target, batch_edges
    integer :: batch_min, batch_max
    integer :: i_task, j_task
    integer :: mode_selected, mode_used, tune_slot
    integer :: tune_calls_slot
    logical :: task_capable, force_explore
    double precision :: t_start, elapsed, edge_cost
    integer, allocatable :: edge_level(:), last_occ(:), last_vir(:)
    integer, allocatable :: level_count(:), level_offset(:), level_pos(:), level_edges(:)
    integer, allocatable, save :: iused_ws(:,:)
    logical, allocatable, save :: latoms_ws(:,:)
    double precision, allocatable, save :: storei_ws(:,:), storej_ws(:,:)
    double precision, allocatable :: sumb_thread(:)
    integer, allocatable :: nrej_thread(:)
    integer, save :: ws_numat = 0, ws_norbs = 0, ws_nthreads = 0
    logical :: need_ws_resize
    logical :: active
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
      allocate (edge_level(nij), last_occ(nocc), last_vir(nvir), stat=alloc_stat)
      if (alloc_stat /= 0) then
        call mopend("Insufficient memory to run DIAGG2 (parallel scheduling)")
        use_parallel_diagg2 = .false.
      else
        edge_level(:) = 0
        last_occ(:) = 0
        last_vir(:) = 0
        max_level = 0
        nactive = 0

        do ij = 1, nij
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
            lvl = Max (last_occ(j), last_vir(i)) + 1
            edge_level(ij) = lvl
            last_occ(j) = lvl
            last_vir(i) = lvl
            if (lvl > max_level) max_level = lvl
            nactive = nactive + 1
          end if
        end do

        if (nactive > 0) then
          allocate (level_count(max_level), level_offset(max_level+1), level_pos(max_level), &
               & level_edges(nactive), stat=alloc_stat)
          if (alloc_stat /= 0) then
            call mopend("Insufficient memory to run DIAGG2 (parallel scheduling)")
            use_parallel_diagg2 = .false.
          else
            level_count(:) = 0
            do ij = 1, nij
              lvl = edge_level(ij)
              if (lvl > 0) level_count(lvl) = level_count(lvl) + 1
            end do

            level_offset(1) = 1
            do lvl = 1, max_level
              level_offset(lvl+1) = level_offset(lvl) + level_count(lvl)
            end do

            level_pos(:) = level_offset(1:max_level)
	            do ij = 1, nij
	              lvl = edge_level(ij)
	              if (lvl > 0) then
	                pos = level_pos(lvl)
	                level_edges(pos) = ij
	                level_pos(lvl) = pos + 1
	              end if
	            end do

	            ! Previous-commit improvement kept here: adaptive team sizing from level width.
	            ! Eq. (0): width_hint = ceil(nactive / max_level)
	            !         diagg2_threads = min(max_threads, 2*width_hint, nactive)
	            width_hint = Max (1, (nactive + max_level - 1) / max_level)
	            diagg2_threads = Min (max_threads, Max (1, 2 * width_hint))
	            diagg2_threads = Min (diagg2_threads, nactive)

	            if (diagg2_threads > 1) then
	              need_ws_resize = .true.
	              ! Previous-commit improvement kept here: persistent per-thread scratch reuse.
	              ! Allocation is skipped when dimensions already satisfy:
	              !   ws_rows >= required_rows and ws_cols >= diagg2_threads
	              if (allocated(iused_ws) .and. allocated(latoms_ws) .and. allocated(storei_ws) .and. allocated(storej_ws)) then
	                if (size(iused_ws,1) >= numat .and. size(iused_ws,2) >= diagg2_threads .and. &
	                   & size(latoms_ws,1) >= numat .and. size(latoms_ws,2) >= diagg2_threads .and. &
	                   & size(storei_ws,1) >= norbs .and. size(storei_ws,2) >= diagg2_threads .and. &
	                   & size(storej_ws,1) >= norbs .and. size(storej_ws,2) >= diagg2_threads) then
	                  need_ws_resize = .false.
	                end if
	              end if

	              if (need_ws_resize) then
	                if (allocated(iused_ws)) deallocate (iused_ws)
	                if (allocated(latoms_ws)) deallocate (latoms_ws)
	                if (allocated(storei_ws)) deallocate (storei_ws)
	                if (allocated(storej_ws)) deallocate (storej_ws)

	                ! Scratch footprint scales as O((numat + norbs)*threads), allocated once and reused.
	                allocate (iused_ws(numat, diagg2_threads), latoms_ws(numat, diagg2_threads), &
	                     & storei_ws(norbs, diagg2_threads), storej_ws(norbs, diagg2_threads), stat=alloc_stat)
	                if (alloc_stat /= 0) then
	                  call mopend("Insufficient memory to run DIAGG2 (parallel scratch workspace)")
	                end if
	                ws_numat = numat
	                ws_norbs = norbs
	                ws_nthreads = diagg2_threads
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
	              task_capable = (max_level > 1 .and. nactive >= 2*diagg2_threads)
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
!$omp parallel num_threads(diagg2_threads) default(shared) private(ij, lvl, pos, tid, i_task, j_task, batch_start, batch_end, batch_edges)
	                tid = omp_get_thread_num() + 1
	                iused_ws(:, tid) = -1
	                latoms_ws(:, tid) = .false.
!$omp single
	                batch_start = 1
	                do while (batch_start <= max_level)
	                  batch_end = batch_start
	                  batch_edges = 0
	                  do while (batch_end <= max_level .and. batch_edges < batch_target)
	                    batch_edges = batch_edges + level_count(batch_end)
	                    batch_end = batch_end + 1
	                  end do
	                  batch_end = Max (batch_start, batch_end - 1)
	                  do lvl = batch_start, batch_end
	                    do pos = level_offset(lvl), level_offset(lvl+1) - 1
	                      ij = level_edges(pos)
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
	                  end do
!$omp taskwait
	                  batch_start = batch_end + 1
	                end do
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
	                sumb_local = 0.d0
	                nrej_local = 0
	                t_start = omp_get_wtime()
!$omp parallel num_threads(diagg2_threads) default(shared) private(ij, lvl, pos, tid) &
!$omp& reduction(+:sumb_local, nrej_local)
	                tid = omp_get_thread_num() + 1
	                iused_ws(:, tid) = -1
	                latoms_ws(:, tid) = .false.

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

	            deallocate (level_edges, level_pos, level_offset, level_count)
	          end if
	        end if

        deallocate (edge_level, last_occ, last_vir)
      end if
    end if
#endif

    if (.not. use_parallel_diagg2) then
      outer_loop: do ij = 1, nij
        i = ifmo(1, ij)
        j = ifmo(2, ij)
        if (Abs (fmo(ij)) >= tiny) then
          c = fmo(ij) * const
          d = eigs(j) - eigv(i) - shift
          if (Abs (c/d) >= biglim) then
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
	            njlen = jur - jlr + 1
	            nilen = iur - ilr + 1
	            storej(1:njlen) = cocc(jlr:jur)
	            storei(1:nilen) = cvir(ilr:iur)
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
	              cocc(jlr:jur) = storej(1:njlen)
	              cvir(ilr:iur) = storei(1:nilen)
              ncf(j) = ncfj
              nce(i) = ncei
              do k = 1, numat
                iused(k) = -1
                latoms(k) = .false.
              end do
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
        end if
      end do outer_loop
    end if
    nrejct(2) = nrejct(1)
    nrejct(1) = nrej
    if (times) then
      call timer (" AFTER DIAGG2 IN ITER")
    end if
  contains
  subroutine process_pair (ij, retry, tiny, biglim, sumb_acc, nrej_acc, iused_t, latoms_t, storei_t, storej_t)
    integer, intent (in) :: ij
    logical, intent (in) :: retry
    double precision, intent (in) :: tiny, biglim
    double precision, intent (inout) :: sumb_acc
    integer, intent (inout) :: nrej_acc
    integer, dimension (numat), intent (inout) :: iused_t
    logical, dimension (numat), intent (inout) :: latoms_t
    double precision, dimension (norbs), intent (inout) :: storei_t, storej_t

    integer :: i, j, ii, jur, l, ilr, incv, iur, jlr, jncf, k, le, lf, loopi, loopj
    integer :: mie, mle, mlee, mlf, mlff, ncei, ncfj
    integer :: njlen, nilen
    double precision :: a, alpha, b, beta, c, d, e, sum

    i = ifmo(1, ij)
    j = ifmo(2, ij)
    if (Abs (fmo(ij)) < tiny) return

    c = fmo(ij) * const
    d = eigs(j) - eigv(i) - shift
    if (Abs (c/d) < biglim) return

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
    njlen = jur - jlr + 1
    nilen = iur - ilr + 1
    ! Previous-commit improvement kept here: contiguous backup copy for rollback path.
    storej_t(1:njlen) = cocc(jlr:jur)
    storei_t(1:nilen) = cvir(ilr:iur)
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
        cocc(jlr:jur) = storej_t(1:njlen)
        cvir(ilr:iur) = storei_t(1:nilen)
        ncf(j) = ncfj
        nce(i) = ncei
        iused_t(:) = -1
        latoms_t(:) = .false.
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
        cocc(jlr:jur) = storej_t(1:njlen)
        cvir(ilr:iur) = storei_t(1:nilen)
        ncf(j) = ncfj
        nce(i) = ncei
        iused_t(:) = -1
        latoms_t(:) = .false.
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
