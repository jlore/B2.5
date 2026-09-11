module b2mod_input_fields
  use b2mod_types
  use b2mod_constants
  use b2us_plasma
  use b2mod_b2cmpa
  implicit none
  private
  public :: transport_input_fields
  external find_file, xerrab

contains

!-----------------------------------------------------------------------

  subroutine transport_input_fields(nCv, ns, pl, rt, dv, &
    partial_allowed, idna, idpa, ivla, ivma, ivsa, &
    ihci, ihvi, ihce, ihve, isig, ialf)
    implicit none
!***********************************************************************
!     Subroutine for reading 2D (nCv) transport coefficients fields
!     provided by an external source, such as GENE-X.
!     IMPORTANT: Partial assignment is allowed when switch=2.
!     Input file: b2.transport.fields
!     TODO: Add this file to the SOLPS case archive, conversion, and
!       SimDB file manifests once its grid-transfer policy is defined.
!     TODO: Generate differentiated versions of this reader.  Until
!       then, differentiated executables reject this input explicitly.
!     A species mapping has the form coefficient_ns(target)=source.
!     For example, hci_ns(2:6)=1 assigns the species-1 HCI field to
!     species 2 through 6.
!
!    List of implemented transport coefficients (from b2tqna.F):
!    - dna: particle transport coefficients for grad(na).
!    - dpa: particle transport coefficients for grad(pa).
!    - vla: anomalous velocity.
!    - vma: anomalous velocity in parallel momentum balance equation.
!    - vsa: viscosity.
!    - hci: species thermal conductivity.
!    - hvi: all-atom thermal strange velocity.
!    - hce: electron thermal conductivity.
!    - hve: electron thermal strange velocity.
!    - sig: electrical conductivity.
!    - alf: thermo-electric coefficient.
!
!    Particle transport.
!      fnax(is) =
!        -dna(is)*(d/dx)na(is)
!        -dpa(is)*(d/dx)pa(is)
!        +vla(0)*na(is)
!      fnay(is) =
!        -dna(is)*(d/dy)na(is)
!        -dpa(is)*(d/dy)pa(is)
!        +vla(1)*na(is)
!
!    Parallel momentum transport.
!      momx(is) = vma(0,is)*na(is)*up(is) - vsa(is)*(d/dx)up(is)
!      momy(is) = vma(1,is)*na(is)*up(is) - vsa(is)*(d/dy)up(is)
!
!    Ion and neutral heat flux.
!      fhix = -sum(hci)*(d/dx)ti + hvi(0)*ni*ti
!      fhiy = -sum(hci)*(d/dy)ti + hvi(1)*ni*ti
!
!    Electron heat flux.
!      fhex = alf*te*ehx - hce*(d/dx)te + hve(0)*ne*te
!      fhey = alf*te*ehy - hce*(d/dy)te + hve(1)*ne*te
!
!    Electric current.
!      fchx = sig*ehx - alf*(d/dx)te
!      fchy = sig*ehy - alf*(d/dy)te
!
!    Here, ehx and ehy are the components of the modified electric
!    field:
!    ehx = -(d/dx)po + (1/(qe*ne))*(d/dx)(ne*te),
!    ehy = -(d/dy)po + (1/(qe*ne))*(d/dy)(ne*te).
!
!***********************************************************************

!.....Input arguments (unchanged on exit)
    integer, intent(in) :: nCv, ns
    type (B2Plasma), intent (in) :: pl
    type (B2Derivatives), intent (in) :: dv
    type (B2Rates), intent (in) :: rt
    logical, intent(in) :: partial_allowed

!.....Input-output arguments (possibly changed on exit)
    real(kind=R8), intent(inout) :: idna(nCv, 0:ns-1)
    real(kind=R8), intent(inout) :: idpa(nCv, 0:ns-1)
    real(kind=R8), intent(inout) :: ivla(nCv, 0:1, 0:ns-1)
    real(kind=R8), intent(inout) :: ivma(nCv, 0:1, 0:ns-1)
    real(kind=R8), intent(inout) :: ivsa(nCv, 0:ns-1)
    real(kind=R8), intent(inout) :: ihci(nCv, 0:ns-1)
    real(kind=R8), intent(inout) :: ihvi(nCv, 0:1)
    real(kind=R8), intent(inout) :: ihce(nCv)
    real(kind=R8), intent(inout) :: ihve(nCv, 0:1)
    real(kind=R8), intent(inout) :: isig(nCv)
    real(kind=R8), intent(inout) :: ialf(nCv)

!.....Internal variables
    integer :: is, idir, iunit, ios, close_ios
    character(len=256) :: filename
    character(len=512) :: iomsg, close_iomsg
    logical :: file_ok
    logical, save :: input_fields_loaded = .false.
    logical, save :: cached_partial_allowed = .false.
    integer, save :: cached_nCv = -1, cached_ns = -1
    real(kind=R8), parameter :: nml_unset = -huge(1.0_R8)
    real(kind=R8), allocatable, save :: dna(:,:), dpa(:,:)
    real(kind=R8), allocatable, save :: vla(:,:,:), vma(:,:,:)
    real(kind=R8), allocatable, save :: vsa(:,:), hci(:,:)
    real(kind=R8), allocatable, save :: hvi(:,:), hce(:)
    real(kind=R8), allocatable, save :: hve(:,:), sig(:), alf(:)
    logical, allocatable, save :: new_dna(:,:), new_dpa(:,:)
    logical, allocatable, save :: new_vla(:,:,:), new_vma(:,:,:)
    logical, allocatable, save :: new_vsa(:,:), new_hci(:,:)
    logical, allocatable, save :: new_hvi(:,:), new_hce(:)
    logical, allocatable, save :: new_hve(:,:), new_sig(:), new_alf(:)
    real(kind=R8), allocatable :: species_field(:,:)
    real(kind=R8), allocatable :: vector_species_field(:,:,:)
    logical, allocatable :: new_species_field(:,:)
    logical, allocatable :: new_vector_species_field(:,:,:)
    integer, allocatable, save :: dna_ns(:), dpa_ns(:)
    integer, allocatable, save :: vla_ns(:,:), vma_ns(:,:)
    integer, allocatable, save :: vsa_ns(:), hci_ns(:)
#ifdef LEGACYCOMP
    integer newunit
    external newunit
#endif
    namelist /transport/ &
      dna, dpa, vla, vma, vsa, hci, hvi, hce, hve, sig, alf, &
      dna_ns, dpa_ns, vla_ns, vma_ns, vsa_ns, hci_ns

!.....Read and validate the input only once.  The raw effective
!     coefficients are cached, while their conversion using the plasma
!     state is repeated below on every call.
    if (input_fields_loaded) then
      if (nCv.ne.cached_nCv .or. ns.ne.cached_ns) &
        call xerrab('Transport field dimensions changed during the run')
      if (partial_allowed.neqv.cached_partial_allowed) &
        call xerrab('b2tqna_input_fields changed during the run')
    else

!.....Initialization
      filename = 'b2.transport.fields'
      allocate(dna(nCv,0:ns-1), dpa(nCv,0:ns-1), &
        vla(nCv,0:1,0:ns-1), vma(nCv,0:1,0:ns-1), &
        vsa(nCv,0:ns-1), hci(nCv,0:ns-1), hvi(nCv,0:1), &
        hce(nCv), hve(nCv,0:1), sig(nCv), alf(nCv), stat=ios)
      if (ios.ne.0) call xerrab( &
        'Could not allocate two-dimensional transport fields')
      allocate(new_dna(nCv,0:ns-1), new_dpa(nCv,0:ns-1), &
        new_vla(nCv,0:1,0:ns-1), new_vma(nCv,0:1,0:ns-1), &
        new_vsa(nCv,0:ns-1), new_hci(nCv,0:ns-1), &
        new_hvi(nCv,0:1), new_hce(nCv), new_hve(nCv,0:1), &
        new_sig(nCv), new_alf(nCv), stat=ios)
      if (ios.ne.0) call xerrab( &
        'Could not allocate two-dimensional transport masks')
      allocate(dna_ns(0:ns-1), dpa_ns(0:ns-1), &
        vla_ns(0:1,0:ns-1), vma_ns(0:1,0:ns-1), &
        vsa_ns(0:ns-1), hci_ns(0:ns-1), stat=ios)
      if (ios.ne.0) call xerrab( &
        'Could not allocate transport species mappings')
      dna = nml_unset
      dpa = nml_unset
      vla = nml_unset
      vma = nml_unset
      vsa = nml_unset
      hci = nml_unset
      hvi = nml_unset
      hce = nml_unset
      hve = nml_unset
      sig = nml_unset
      alf = nml_unset
      dna_ns = -1
      dpa_ns = -1
      vla_ns = -1
      vma_ns = -1
      vsa_ns = -1
      hci_ns = -1


!.....Find b2.transport.fields and gather the provided namelist values
!     In normal MPI operation only rank zero executes the B2 solver;
!     the other ranks remain in the EIRENE service loop.  Therefore no
!     MPI broadcast belongs here.
      call find_file(filename, file_ok)
      if (.not.file_ok) call xerrab('No '//trim(filename))

      iomsg = ''
#ifdef LEGACYCOMP
      iunit = newunit()
      open(unit=iunit, file=trim(filename), status='old', &
        action='read', form='formatted', iostat=ios)
#else
      open(newunit=iunit, file=trim(filename), status='old', &
        action='read', form='formatted', iostat=ios, iomsg=iomsg)
#endif
      if (ios.ne.0) then
#ifdef LEGACYCOMP
        call xerrab('Could not open '//trim(filename))
#else
        call xerrab('Could not open '//trim(filename)//': '// &
          trim(iomsg))
#endif
      endif

#ifdef LEGACYCOMP
      read(iunit, nml=transport, iostat=ios)
#else
      read(iunit, nml=transport, iostat=ios, iomsg=iomsg)
#endif
      if (ios.ne.0) then
        close(iunit, iostat=close_ios)
#ifdef LEGACYCOMP
        call xerrab('Could not read transport namelist from '// &
          trim(filename))
#else
        call xerrab('Could not read transport namelist from '// &
          trim(filename)//': '//trim(iomsg))
#endif
      endif

      close_iomsg = ''
#ifdef LEGACYCOMP
      close(iunit, iostat=close_ios)
#else
      close(iunit, iostat=close_ios, iomsg=close_iomsg)
#endif
      if (close_ios.ne.0) then
#ifdef LEGACYCOMP
        call xerrab('Could not close '//trim(filename))
#else
        call xerrab('Could not close '//trim(filename)//': '// &
          trim(close_iomsg))
#endif
      endif


!.....Find which values were set by the namelist
      new_dna = (dna .ne. nml_unset)
      new_dpa = (dpa .ne. nml_unset)
      new_vla = (vla .ne. nml_unset)
      new_vma = (vma .ne. nml_unset)
      new_vsa = (vsa .ne. nml_unset)
      new_hci = (hci .ne. nml_unset)
      new_hvi = (hvi .ne. nml_unset)
      new_hce = (hce .ne. nml_unset)
      new_hve = (hve .ne. nml_unset)
      new_sig = (sig .ne. nml_unset)
      new_alf = (alf .ne. nml_unset)

!.....Reject values that cannot be used safely by the transport model
      if (any(new_dna.and..not.(abs(dna).le.huge(1.0_R8)))) &
        call xerrab('Non-finite DNA value in '//trim(filename))
      if (any(new_dpa.and..not.(abs(dpa).le.huge(1.0_R8)))) &
        call xerrab('Non-finite DPA value in '//trim(filename))
      if (any(new_vla.and..not.(abs(vla).le.huge(1.0_R8)))) &
        call xerrab('Non-finite VLA value in '//trim(filename))
      if (any(new_vma.and..not.(abs(vma).le.huge(1.0_R8)))) &
        call xerrab('Non-finite VMA value in '//trim(filename))
      if (any(new_vsa.and..not.(abs(vsa).le.huge(1.0_R8)))) &
        call xerrab('Non-finite VSA value in '//trim(filename))
      if (any(new_hci.and..not.(abs(hci).le.huge(1.0_R8)))) &
        call xerrab('Non-finite HCI value in '//trim(filename))
      if (any(new_hvi.and..not.(abs(hvi).le.huge(1.0_R8)))) &
        call xerrab('Non-finite HVI value in '//trim(filename))
      if (any(new_hce.and..not.(abs(hce).le.huge(1.0_R8)))) &
        call xerrab('Non-finite HCE value in '//trim(filename))
      if (any(new_hve.and..not.(abs(hve).le.huge(1.0_R8)))) &
        call xerrab('Non-finite HVE value in '//trim(filename))
      if (any(new_sig.and..not.(abs(sig).le.huge(1.0_R8)))) &
        call xerrab('Non-finite SIG value in '//trim(filename))
      if (any(new_alf.and..not.(abs(alf).le.huge(1.0_R8)))) &
        call xerrab('Non-finite ALF value in '//trim(filename))

!.....Diffusivities, viscosity, and conductivity must be non-negative.
!     Velocity-like fields and the thermo-electric coefficient may be
!     signed.
      if (any(new_dna.and.(dna.lt.0.0_R8))) &
        call xerrab('Negative DNA value in '//trim(filename))
      if (any(new_dpa.and.(dpa.lt.0.0_R8))) &
        call xerrab('Negative DPA value in '//trim(filename))
      if (any(new_vsa.and.(vsa.lt.0.0_R8))) &
        call xerrab('Negative VSA value in '//trim(filename))
      if (any(new_hci.and.(hci.lt.0.0_R8))) &
        call xerrab('Negative HCI value in '//trim(filename))
      if (any(new_hce.and.(hce.lt.0.0_R8))) &
        call xerrab('Negative HCE value in '//trim(filename))
      if (any(new_sig.and.(sig.lt.0.0_R8))) &
        call xerrab('Negative SIG value in '//trim(filename))


!.....Reject incomplete nCv fields unless partial fields are allowed
      if (.not. partial_allowed) then
        do is=0,ns-1
          if (any(new_dna(:,is)).and..not.all(new_dna(:,is))) then
            write(*,*) 'Incomplete DNA field for species ', is
            call xerrab('Incomplete DNA field in '//trim(filename))
          endif
          if (any(new_dpa(:,is)).and..not.all(new_dpa(:,is))) then
            write(*,*) 'Incomplete DPA field for species ', is
            call xerrab('Incomplete DPA field in '//trim(filename))
          endif
          do idir=0,1
            if (any(new_vla(:,idir,is)).and. &
                .not.all(new_vla(:,idir,is))) then
              write(*,*) 'Incomplete VLA field for direction, species ', &
                idir, is
              call xerrab('Incomplete VLA field in '//trim(filename))
            endif
            if (any(new_vma(:,idir,is)).and. &
                .not.all(new_vma(:,idir,is))) then
              write(*,*) 'Incomplete VMA field for direction, species ', &
                idir, is
              call xerrab('Incomplete VMA field in '//trim(filename))
            endif
          enddo
          if (any(new_vsa(:,is)).and..not.all(new_vsa(:,is))) then
            write(*,*) 'Incomplete VSA field for species ', is
            call xerrab('Incomplete VSA field in '//trim(filename))
          endif
          if (any(new_hci(:,is)).and..not.all(new_hci(:,is))) then
            write(*,*) 'Incomplete HCI field for species ', is
            call xerrab('Incomplete HCI field in '//trim(filename))
          endif
        enddo
        do idir=0,1
          if (any(new_hvi(:,idir)).and. &
              .not.all(new_hvi(:,idir))) then
            write(*,*) 'Incomplete HVI field for direction ', idir
            call xerrab('Incomplete HVI field in '//trim(filename))
          endif
          if (any(new_hve(:,idir)).and. &
              .not.all(new_hve(:,idir))) then
            write(*,*) 'Incomplete HVE field for direction ', idir
            call xerrab('Incomplete HVE field in '//trim(filename))
          endif
        enddo
        if (any(new_hce).and..not.all(new_hce)) &
          call xerrab('Incomplete HCE field in '//trim(filename))
        if (any(new_sig).and..not.all(new_sig)) &
          call xerrab('Incomplete SIG field in '//trim(filename))
        if (any(new_alf).and..not.all(new_alf)) &
          call xerrab('Incomplete ALF field in '//trim(filename))
      endif


!.....Validate species assignments.  -1 means no assignment.
      if (any(dna_ns.lt.-1).or.any(dna_ns.ge.ns)) &
        call xerrab('Invalid species index in dna_ns')
      if (any(dpa_ns.lt.-1).or.any(dpa_ns.ge.ns)) &
        call xerrab('Invalid species index in dpa_ns')
      if (any(vla_ns.lt.-1).or.any(vla_ns.ge.ns)) &
        call xerrab('Invalid species index in vla_ns')
      if (any(vma_ns.lt.-1).or.any(vma_ns.ge.ns)) &
        call xerrab('Invalid species index in vma_ns')
      if (any(vsa_ns.lt.-1).or.any(vsa_ns.ge.ns)) &
        call xerrab('Invalid species index in vsa_ns')
      if (any(hci_ns.lt.-1).or.any(hci_ns.ge.ns)) &
        call xerrab('Invalid species index in hci_ns')

!.....Assign fields from source species to target species.  Keep a
!     snapshot so chained assignments do not depend on species order.
      if (any(dna_ns.ge.0).or.any(dpa_ns.ge.0).or. &
          any(vla_ns.ge.0).or.any(vma_ns.ge.0).or. &
          any(vsa_ns.ge.0).or.any(hci_ns.ge.0)) then
        allocate(species_field(nCv,0:ns-1), &
          new_species_field(nCv,0:ns-1), &
          vector_species_field(nCv,0:1,0:ns-1), &
          new_vector_species_field(nCv,0:1,0:ns-1), stat=ios)
        if (ios.ne.0) call xerrab( &
          'Could not allocate transport species-mapping workspace')

        species_field = dna
        new_species_field = new_dna
        do is=0,ns-1
          if (dna_ns(is).ge.0) then
            if (any(new_species_field(:,dna_ns(is)))) then
              where (new_species_field(:,dna_ns(is)))
                dna(:,is)=species_field(:,dna_ns(is))
              endwhere
              new_dna(:,is)=new_dna(:,is).or. &
                new_species_field(:,dna_ns(is))
            else
              write(*,'(A, I0, A, I0)') &
                'Cannot set DNA field of species ', is, &
                ' to the UNASSIGNED field of species ', dna_ns(is)
              call xerrab('Invalid DNA species mapping')
            endif
          endif
        enddo

        species_field = dpa
        new_species_field = new_dpa
        do is=0,ns-1
          if (dpa_ns(is).ge.0) then
            if (any(new_species_field(:,dpa_ns(is)))) then
              where (new_species_field(:,dpa_ns(is)))
                dpa(:,is)=species_field(:,dpa_ns(is))
              endwhere
              new_dpa(:,is)=new_dpa(:,is).or. &
                new_species_field(:,dpa_ns(is))
            else
              write(*,'(A, I0, A, I0)') &
                'Cannot set DPA field of species ', is, &
                ' to the UNASSIGNED field of species ', dpa_ns(is)
              call xerrab('Invalid DPA species mapping')
            endif
          endif
        enddo

        vector_species_field = vla
        new_vector_species_field = new_vla
        do is=0,ns-1
          do idir=0,1
            if (vla_ns(idir,is).ge.0) then
              if (any(new_vector_species_field(:,idir, &
                  vla_ns(idir,is)))) then
                where (new_vector_species_field(:,idir, &
                  vla_ns(idir,is)))
                  vla(:,idir,is)=vector_species_field(:,idir, &
                    vla_ns(idir,is))
                endwhere
                new_vla(:,idir,is)=new_vla(:,idir,is).or. &
                  new_vector_species_field(:,idir,vla_ns(idir,is))
              else
                write(*,'(A, I0, A, I0, A, I0)') &
                  'Cannot set VLA direction ', idir, &
                  ' field of species ', is, &
                  ' to the UNASSIGNED field of species ', vla_ns(idir,is)
                call xerrab('Invalid VLA species mapping')
              endif
            endif
          enddo
        enddo

        vector_species_field = vma
        new_vector_species_field = new_vma
        do is=0,ns-1
          do idir=0,1
            if (vma_ns(idir,is).ge.0) then
              if (any(new_vector_species_field(:,idir, &
                  vma_ns(idir,is)))) then
                where (new_vector_species_field(:,idir, &
                  vma_ns(idir,is)))
                  vma(:,idir,is)=vector_species_field(:,idir, &
                    vma_ns(idir,is))
                endwhere
                new_vma(:,idir,is)=new_vma(:,idir,is).or. &
                  new_vector_species_field(:,idir,vma_ns(idir,is))
              else
                write(*,'(A, I0, A, I0, A, I0)') &
                  'Cannot set VMA direction ', idir, &
                  ' field of species ', is, &
                  ' to the UNASSIGNED field of species ', vma_ns(idir,is)
                call xerrab('Invalid VMA species mapping')
              endif
            endif
          enddo
        enddo

        species_field = vsa
        new_species_field = new_vsa
        do is=0,ns-1
          if (vsa_ns(is).ge.0) then
            if (any(new_species_field(:,vsa_ns(is)))) then
              where (new_species_field(:,vsa_ns(is)))
                vsa(:,is)=species_field(:,vsa_ns(is))
              endwhere
              new_vsa(:,is)=new_vsa(:,is).or. &
                new_species_field(:,vsa_ns(is))
            else
              write(*,'(A, I0, A, I0)') &
                'Cannot set VSA field of species ', is, &
                ' to the UNASSIGNED field of species ', vsa_ns(is)
              call xerrab('Invalid VSA species mapping')
            endif
          endif
        enddo

        species_field = hci
        new_species_field = new_hci
        do is=0,ns-1
          if (hci_ns(is).ge.0) then
            if (any(new_species_field(:,hci_ns(is)))) then
              where (new_species_field(:,hci_ns(is)))
                hci(:,is)=species_field(:,hci_ns(is))
              endwhere
              new_hci(:,is)=new_hci(:,is).or. &
                new_species_field(:,hci_ns(is))
            else
              write(*,'(A, I0, A, I0)') &
                'Cannot set HCI field of species ', is, &
                ' to the UNASSIGNED field of species ', hci_ns(is)
              call xerrab('Invalid HCI species mapping')
            endif
          endif
        enddo
      endif

!.....The sentinel is needed only while reading the namelist.  Replace
!     omitted values with zero before caching so masked conversions
!     cannot perform arithmetic on the sentinel.
      where(.not.new_dna) dna = 0.0_R8
      where(.not.new_dpa) dpa = 0.0_R8
      where(.not.new_vla) vla = 0.0_R8
      where(.not.new_vma) vma = 0.0_R8
      where(.not.new_vsa) vsa = 0.0_R8
      where(.not.new_hci) hci = 0.0_R8
      where(.not.new_hvi) hvi = 0.0_R8
      where(.not.new_hce) hce = 0.0_R8
      where(.not.new_hve) hve = 0.0_R8
      where(.not.new_sig) sig = 0.0_R8
      where(.not.new_alf) alf = 0.0_R8

      cached_nCv = nCv
      cached_ns = ns
      cached_partial_allowed = partial_allowed
      input_fields_loaded = .true.
    endif

!.....Convert effective diffusivities to the internal coefficients
!     used by b2tqna.  DNA is already a diffusivity, while VLA, VMA,
!     HVI, and HVE are already velocities, so they remain unchanged.
    where(new_dna) idna = dna
    where(new_vla) ivla = vla
    where(new_vma) ivma = vma
    where(new_hvi) ihvi = hvi
    where(new_hve) ihve = hve
    do is=0,ns-1
      where (new_dpa(:,is))
        idpa(:,is) = dpa(:,is)/(rt%rza(:,is)*pl%te+pl%ti)
      endwhere
      where(new_hci(:,is)) ihci(:,is) = hci(:,is)*pl%na(:,is)
      where(new_vsa(:,is)) &
        ivsa(:,is) = vsa(:,is)*mp*am(is)*pl%na(:,is)
    enddo
    where(new_hce) ihce = hce*dv%ne
    where(new_sig) isig = sig*qe*dv%ne
    where(new_alf) ialf = alf*dv%ne*sqrt(qe/pl%te)

    return
  end subroutine transport_input_fields

end module b2mod_input_fields
