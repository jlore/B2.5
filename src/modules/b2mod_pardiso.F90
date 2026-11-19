module b2mod_pardiso
  use b2mod_types, only: R8
#if defined(SOLPS_PARDISO) && !defined(LEGACYCOMP) && !defined(_OPENMP)
  use b2mod_types, only: I8
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite ! IGNORE
#endif
  implicit none
  private
  public :: pardiso_require, pardiso_solve, pardiso_finish
#if defined(SOLPS_PARDISO) && !defined(LEGACYCOMP) && !defined(_OPENMP)
  include 'mkl_pardiso.fi'
  public :: pardiso_backward_error

  type solver_context
    type(MKL_PARDISO_HANDLE) :: pt(64)
    integer :: iparm(64), n = 0
    integer :: solves = 0, analyses = 0
    integer :: refinements = 0, perturbed_pivots = 0
    real(kind=R8) :: seconds = 0.0_R8, max_error = 0.0_R8
    integer, allocatable :: ia(:), ja(:), perm(:)
    real(kind=R8), allocatable :: av(:), b(:), x(:)
    logical :: active = .false.
  end type solver_context

  type context_entry
    character(len=80) :: name
    integer :: stencil
    type(solver_context) :: solver
    type(context_entry), pointer :: next => null()
  end type context_entry

  type(context_entry), pointer, save :: first => null()
#endif

contains

  subroutine pardiso_require()
    ! The MKL interface uses default LP64 integers and Fortran 2003 IMPORT.
    ! Contexts are process-local; OpenMP builds retain the other solvers.
#if !defined(SOLPS_PARDISO)
    call xerrab('PARDISO requires a build with SOLPS_PARDISO and MKL')
#elif defined(LEGACYCOMP)
    call xerrab('PARDISO requires a Fortran 2003 compiler')
#elif defined(_OPENMP)
    call xerrab('Structured PARDISO requires a build without OpenMP')
#else
    call xertst(bit_size(0) == 32,'PARDISO requires LP64 integers')
#endif
  end subroutine pardiso_require

  subroutine pardiso_solve(stencil, name, n, nz, rows, cols, values, rhs, x)
    integer, intent(in) :: stencil, n, nz, rows(nz), cols(nz)
    character(len=*), intent(in) :: name
    real(kind=R8), intent(in) :: values(nz), rhs(n)
    real(kind=R8), intent(out) :: x(n)
#if defined(SOLPS_PARDISO) && !defined(LEGACYCOMP) && !defined(_OPENMP)
    type(context_entry), pointer :: entry
    integer, allocatable :: ia(:), ja(:)
    real(kind=R8), allocatable :: av(:)
    integer :: used, flag, release_flag
    real(kind=R8) :: started, finished

    call cpu_time(started)
    x = 0.0_R8
    call pardiso_require()
    if (n < 1 .or. nz < 1 .or. len_trim(name) > 80) then
      call release_all(release_flag)
      call xerrab('Invalid PARDISO matrix size or equation name')
      return
    end if
    ! Equation names include the species identifier supplied by the caller.
    ! Stencil identity also distinguishes alternative assemblies of an equation.
    entry => first
    do while (associated(entry))
      if (entry%stencil == stencil .and. entry%name == name) exit
      entry => entry%next
    end do
    if (.not.associated(entry)) then
      allocate(entry)
      entry%name = name
      entry%stencil = stencil
      entry%next => first
      first => entry
    end if
    allocate(ia(n+1),ja(nz),av(nz))
    call triplet_to_csr(n,rows,cols,values,ia,ja,av,used,flag)
    if (flag == 0) then
      call context_solve(entry%solver,n,ia,ja(1:used),av(1:used),rhs,x,flag)
    end if
    deallocate(ia,ja,av)
    call cpu_time(finished)
    entry%solver%seconds = entry%solver%seconds+finished-started
    if (flag /= 0) then
      x = 0.0_R8
      write(*,'(a,a,a,i6,a,i8)') 'PARDISO failure: ',trim(name), &
        ' stencil=',stencil,' status=',flag
      call release_all(release_flag)
      call xerrab('PARDISO failed; no correction accepted')
    end if
#else
    ! Unavailable backends fail explicitly rather than changing solver choice.
    x = 0.0_R8
    call pardiso_require()
#endif
  end subroutine pardiso_solve

  subroutine pardiso_finish()
#if defined(SOLPS_PARDISO) && !defined(LEGACYCOMP) && !defined(_OPENMP)
    integer :: flag
    call release_all(flag)
    call xertst(flag == 0,'PARDISO memory release failed')
#endif
  end subroutine pardiso_finish

#if defined(SOLPS_PARDISO) && !defined(LEGACYCOMP) && !defined(_OPENMP)
  subroutine release_all(flag)
    integer, intent(out) :: flag
    type(context_entry), pointer :: entry, next
    integer :: status
    flag = 0
    entry => first
    do while (associated(entry))
      next => entry%next
      write(*,'(a,a,a,i2,4(a,i8),2(a,es13.5))') &
        'PARDISO context: ',trim(entry%name),' stencil=',entry%stencil, &
        ' solves=',entry%solver%solves,' analyses=',entry%solver%analyses, &
        ' refinements=',entry%solver%refinements, &
        ' perturbed_pivots=',entry%solver%perturbed_pivots, &
        ' max_error=',entry%solver%max_error, &
        ' cpu_seconds=',entry%solver%seconds
      call release_context(entry%solver,status)
      if (status /= 0) flag = status
      deallocate(entry)
      entry => next
    end do
    nullify(first)
  end subroutine release_all

  subroutine release_context(context, flag)
    ! Phase -1 is required even after an unsuccessful library phase.
    ! Context work arrays isolate library release from the accepted correction.
    type(solver_context), intent(inout) :: context
    integer, intent(out) :: flag
    flag = 0
    if (context%active) then
      call pardiso(context%pt,1,1,11,-1,context%n,context%av,context%ia, &
        context%ja,context%perm,1,context%iparm,0,context%b,context%x,flag)
      context%active = .false.
      context%n = 0
      deallocate(context%ia,context%ja,context%av,context%perm,context%b,context%x)
      context%pt%dummy = 0_I8
    end if
  end subroutine release_context

  subroutine triplet_to_csr(n, rows, cols, values, ia, ja, av, used, flag)
    ! Preserve the assembled row/column orientation and one-based indices.
    ! Sort each short stencil row, sum duplicates and retain cancellation zeros.
    integer, intent(in) :: n, rows(:), cols(:)
    real(kind=R8), intent(in) :: values(:)
    integer, intent(out) :: ia(:), ja(:), used, flag
    real(kind=R8), intent(out) :: av(:)
    integer, allocatable :: counts(:), offsets(:), workcols(:)
    real(kind=R8), allocatable :: workvals(:)
    integer :: k, row, p, q, start, finish, col, nz
    real(kind=R8) :: value

    flag = -100
    used = 0
    nz = size(values)
    if (n < 1 .or. nz < 1) return
    if (size(rows) /= nz .or. size(cols) /= nz) return
    if (size(ia) /= n+1 .or. size(ja) < nz .or. size(av) < nz) return
    if (any(rows < 1) .or. any(rows > n)) return
    if (any(cols < 1) .or. any(cols > n)) return
    if (.not.all(ieee_is_finite(values))) return
    allocate(counts(n),offsets(n+1),workcols(nz),workvals(nz))
    counts = 0
    do k = 1, nz
      counts(rows(k)) = counts(rows(k))+1
    end do
    offsets(1) = 1
    do row = 1, n
      offsets(row+1) = offsets(row)+counts(row)
    end do
    counts = offsets(1:n)
    do k = 1, nz
      p = counts(rows(k))
      workcols(p) = cols(k)
      workvals(p) = values(k)
      counts(rows(k)) = p+1
    end do
    do row = 1, n
      start = offsets(row)
      finish = offsets(row+1)-1
      do p = start+1, finish
        col = workcols(p)
        value = workvals(p)
        q = p-1
        do while (q >= start)
          if (workcols(q) <= col) exit
          workcols(q+1) = workcols(q)
          workvals(q+1) = workvals(q)
          q = q-1
        end do
        workcols(q+1) = col
        workvals(q+1) = value
      end do
      ia(row) = used+1
      do p = start, finish
        if (used >= ia(row)) then
          if (ja(used) == workcols(p)) then
            av(used) = av(used)+workvals(p)
            cycle
          end if
        end if
        used = used+1
        ja(used) = workcols(p)
        av(used) = workvals(p)
      end do
    end do
    ia(n+1) = used+1
    flag = 0
    if (.not.all(ieee_is_finite(av(1:used)))) flag = -101
    deallocate(counts,offsets,workcols,workvals)
  end subroutine triplet_to_csr

  subroutine context_solve(context, n, ia, ja, av, rhs, x, flag)
    type(solver_context), intent(inout) :: context
    integer, intent(in) :: n, ia(:), ja(:)
    real(kind=R8), intent(in) :: av(:), rhs(:)
    real(kind=R8), intent(out) :: x(:)
    integer, intent(out) :: flag
    integer, parameter :: phases(3) = (/11,22,33/)
    integer :: k, phase, release_flag
    real(kind=R8) :: error

    x = 0.0_R8
    flag = -100
    error = huge(1.0_R8)
    if (.not.all(ieee_is_finite(rhs))) return
    call release_context(context,flag)
    if (flag /= 0) return
    allocate(context%ia(n+1),context%ja(size(ja)),context%av(size(av)), &
      context%perm(n),context%b(n),context%x(n))
    context%n = n
    context%ia = ia
    context%ja = ja
    context%av = av
    context%b = rhs
    context%x = 0.0_R8
    context%perm = 0
    context%pt%dummy = 0_I8
    ! Explicit parameters match the WG real unsymmetric solver. Entry 9
    ! controls refinement tolerance; entry 10 controls pivot perturbation.
    context%iparm = 0
    context%iparm(1) = 1
    context%iparm(2) = 2
    context%iparm(9) = 13
    context%iparm(10) = 13
    context%iparm(11) = 1
    context%iparm(13) = 1
    context%iparm(17) = -1
    context%iparm(18) = -1
    context%active = .true.
    ! Scaling and matching depend on current coefficients. Each system
    ! therefore receives a new analysis, factorization and solve.
    do k = 1, 3
      phase = phases(k)
      if (phase == 11) context%analyses = context%analyses+1
      call pardiso(context%pt,1,1,11,phase,n,context%av,context%ia, &
        context%ja,context%perm,1,context%iparm,0,context%b,context%x,flag)
      if (flag /= 0) exit
    end do
    if (flag == 0) then
      error = pardiso_backward_error(n,ia,ja,av,rhs,context%x)
      if (error > 1.0e-10_R8) flag = -102
    end if
    if (flag == 0) then
      x = context%x
      context%solves = context%solves+1
      context%refinements = context%refinements+context%iparm(7)
      context%perturbed_pivots = context%perturbed_pivots+context%iparm(14)
      context%max_error = max(context%max_error,error)
    end if
    call release_context(context,release_flag)
    if (flag == 0) flag = release_flag
    if (flag /= 0) x = 0.0_R8
  end subroutine context_solve

  function pardiso_backward_error(n, ia, ja, av, rhs, x) result(error)
    ! Componentwise backward error: max |b-Ax|/(|b|+|A||x|).
    ! Bounds admit ordinary products only when every nonzero product is
    ! normal and the sum of absolute terms cannot overflow. Extreme rows
    ! use fractions and a common radix exponent before forming products.
    integer, intent(in) :: n, ia(:), ja(:)
    real(kind=R8), intent(in) :: av(:), rhs(:), x(:)
    real(kind=R8) :: error, lower, upper, ax, denom, term, b, ratio
    integer :: row, k, rowexp, termexp
    logical :: ordinary

    error = huge(1.0_R8)
    if (.not.all(ieee_is_finite(x))) return
    if (.not.all(ieee_is_finite(av))) return
    if (.not.all(ieee_is_finite(rhs))) return
    lower = sqrt(tiny(1.0_R8))
    upper = sqrt(huge(1.0_R8)/(2.0_R8*(real(size(ja),R8)+1.0_R8)))
    ordinary = maxval(abs(av)) <= upper .and. maxval(abs(x)) <= upper
    ordinary = ordinary .and. minval(abs(av),mask=av /= 0.0_R8) >= lower
    ordinary = ordinary .and. minval(abs(x),mask=x /= 0.0_R8) >= lower
    ordinary = ordinary .and. maxval(abs(rhs)) <= huge(1.0_R8)/2.0_R8
    ordinary = ordinary .and. minval(abs(rhs),mask=rhs /= 0.0_R8) >= tiny(1.0_R8)
    error = 0.0_R8
    do row = 1, n
      ax = 0.0_R8
      if (ordinary) then
        b = rhs(row)
        denom = abs(b)
        do k = ia(row), ia(row+1)-1
          term = av(k)*x(ja(k))
          ax = ax+term
          denom = denom+abs(term)
        end do
      else
        rowexp = -huge(0)
        if (rhs(row) /= 0.0_R8) rowexp = exponent(rhs(row))
        do k = ia(row), ia(row+1)-1
          if (av(k) == 0.0_R8 .or. x(ja(k)) == 0.0_R8) cycle
          rowexp = max(rowexp,exponent(av(k))+exponent(x(ja(k))))
        end do
        if (rowexp == -huge(0)) cycle
        b = scale(rhs(row),-rowexp)
        denom = abs(b)
        do k = ia(row), ia(row+1)-1
          if (av(k) == 0.0_R8 .or. x(ja(k)) == 0.0_R8) cycle
          termexp = exponent(av(k))+exponent(x(ja(k)))-rowexp
          term = scale(fraction(av(k))*fraction(x(ja(k))),termexp)
          ax = ax+term
          denom = denom+abs(term)
        end do
      end if
      ratio = abs(b-ax)
      if (denom > 0.0_R8) ratio = ratio/denom
      error = max(error,ratio)
    end do
  end function pardiso_backward_error
#endif

end module b2mod_pardiso
