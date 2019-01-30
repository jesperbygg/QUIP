
  !*************************************************************************
  !*
  !*  Potential_Switch routines
  !*
  !*************************************************************************

  recursive subroutine Potential_Switch_Initialise(this, args_str, pot1, pot2, r1, r2, mpi, error)
    type(Potential_Switch), intent(inout) :: this
    character(len=*), intent(in) :: args_str
    type(Potential), intent(in), target :: pot1, pot2
    real(dp) :: r1, r2
    type(MPI_Context), intent(in), optional :: mpi
    integer, intent(out), optional :: error

    INIT_ERROR(error)

    call finalise(this)

    this%pot1 => pot1
    this%pot2 => pot2
    this%r1 = r1
    this%r2 = r2

    if (present(mpi)) this%mpi = mpi

  end subroutine Potential_Switch_Initialise

  recursive subroutine Potential_Switch_Finalise(this)
    type(Potential_Switch), intent(inout) :: this
    
    nullify(this%pot1)
    nullify(this%pot2)

  end subroutine Potential_Switch_Finalise

  recursive subroutine Potential_Switch_Print(this, file)
    type(Potential_Switch), intent(inout) :: this
    type(Inoutput), intent(inout), optional :: file

    call print('Potential_Switch:', file=file)
    call print('', file=file)
    call print('r1: '//this%r1//'', file=file)
    call print('r2: '//this%r2//'', file=file)
    call print('', file=file)
    if (associated(this%pot1)) then
       call print('Potential 1:', file=file)
       call print(this%pot1, file=file)
       call print('', file=file)
    else
       call print('Potential 1 not initialised', file=file)
       call print('', file=file)
    end if
    if (associated(this%pot2)) then
       call print('Potential 2:', file=file)
       call Print(this%pot2, file=file)
       call print('', file=file)
    else
       call print('Potential 2 not initialised', file=file)
       call print('', file=file)
    end if

  end subroutine Potential_Switch_Print

  recursive subroutine Potential_Switch_Calc(this, at, args_str, error)
    type(Potential_Switch), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    character(*), intent(in), optional :: args_str
    integer, intent(out), optional :: error

    real(dp) :: energy, virial(3,3)
    real(dp), pointer :: at_force_ptr(:,:), at_local_energy_ptr(:), at_local_virial_ptr(:,:)

    real(dp) :: my_e_1, my_e_2, energy_1, energy_2
    real(dp), allocatable :: my_local_e_1(:), my_local_e_2(:)
    real(dp), allocatable :: my_f_1(:,:), my_local_virial_1(:,:)
    real(dp), allocatable :: my_f_2(:,:), my_local_virial_2(:,:)
    real(dp) :: my_virial_1(3,3), my_virial_2(3,3)
    type(Dictionary) :: params
    character(STRING_LENGTH) :: calc_energy, calc_force, calc_local_energy, calc_virial, calc_local_virial, calc_args_pot1, calc_args_pot2, my_args_str
    logical :: store_contributions

    real(dp) :: rmin, r, dr(3), drmin(3), S, Sp, r1, r2
    integer :: i, j, ij

    INIT_ERROR(error)

    call initialise(params)
    call param_register(params,"energy", "", calc_energy, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params,"force", "", calc_force, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params,"virial", "", calc_virial, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params,"local_energy", "", calc_local_energy, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params,"local_virial", "", calc_local_virial, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params,"calc_args_pot1", "", calc_args_pot1, help_string="additional args_str to pass along to pot1")
    call param_register(params,"calc_args_pot2", "", calc_args_pot2, help_string="additional args_str to pass along to pot2")
    call param_register(params,"store_contributions", "F", store_contributions, help_string="if true, store contributions to sum with _pot1 and _pot2 suffixes")
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_Switch_calc args_str')) then
       RAISE_ERROR('Potential_Switch_calc failed to parse args_str="'//trim(args_str)//'"', error)
    endif
    call finalise(params)

    my_args_str = optional_default("", args_str)


    ! Simple switching scheme: using minimum interatomic distances in cutoff environment
    ! TODO make 'calc_local_energy' true, or whatever it should be when we need it (why 'local_energy' ??)

    ! check that we have read in r1, r2 and that they are reasonaby chosen (r1<r2)
    r1 = this%r1
    r2 = this%r2
    if (r1 < 0.0) then
      RAISE_ERROR('Potential_Switch_calc: bad or missing value for r1 in init_args (0<r1<r2), '//r1//'', error)
    end if
    if (r2 < 0.0) then
      RAISE_ERROR('Potential_Switch_calc: bad or missing value for r2 in init_args (0<r1<r2), '//r2//'', error)
    end if
    if (r1 > r2) then
      RAISE_ERROR('Potential_Switch_calc: r1 should not be larger than r2', error)
    end if

    ! calculate and save energies/forces/virials from both pots
    call calc(this%pot1, at, args_str=trim(my_args_str)//" "//calc_args_pot1, error=error)
    PASS_ERROR(error)
    ! TODO is this needed?
    call assign_property_pointer(at, trim(calc_local_energy), at_local_energy_ptr, error=error)
    PASS_ERROR(error)
    allocate(my_local_e_1(at%N))
    my_local_e_1 = at_local_energy_ptr
    call get_param_value(at, trim(calc_energy), energy)
    energy_1 = energy
    if (len_trim(calc_force) > 0) then
      ! TODO is this needed?
      call assign_property_pointer(at, trim(calc_force), at_force_ptr, error=error)
      PASS_ERROR(error)
      allocate(my_f_1(3, at%N))
      my_f_1 = at_force_ptr
    end if
    if (len_trim(calc_virial) > 0) then
      RAISE_ERROR('Potential_Switch_calc: virials not yet implemented..', error)
      allocate(my_local_virial_1(9, at%N))
    end if

    call calc(this%pot2, at, args_str=trim(my_args_str)//" "//calc_args_pot2, error=error)
    PASS_ERROR(error)
    ! TODO is this needed?
    call assign_property_pointer(at, trim(calc_local_energy), at_local_energy_ptr, error=error)
    PASS_ERROR(error)
    allocate(my_local_e_2(at%N))
    my_local_e_2 = at_local_energy_ptr
    call get_param_value(at, trim(calc_energy), energy)
    energy_2 = energy
    if (len_trim(calc_force) > 0) then
      call assign_property_pointer(at, trim(calc_force), at_force_ptr, error=error)
      PASS_ERROR(error)
      allocate(my_f_2(3, at%N))
      my_f_2 = at_force_ptr
    end if
    if (len_trim(calc_virial) > 0) then
      RAISE_ERROR('Potential_Switch_calc: virials not yet implemented..', error)
      allocate(my_local_virial_2(9, at%N))
    end if

    call print("Potential_switch my_e_1_orig " // energy_1, PRINT_VERBOSE)
    call print("Potential_switch my_e_2_orig " // energy_2, PRINT_VERBOSE)

    ! so the switching/scaling for each atom
    ! TODO parallelise?
    do i = 1, at%N
      ! find minimum r_ij
      rmin = 1000.0_dp
      do ij = 1, n_neighbours(at, i)
        j = neighbour(at, i, ij, distance=r, cosines=dr)
        if (r < rmin) then
          rmin = r
          drmin = dr
        endif
      end do

      ! energies
      S = Switching_Function(rmin, r1, r2)
      my_local_e_1(i) = (1.0_dp - S) * my_local_e_1(i)
      my_local_e_2(i) = S * my_local_e_2(i)
      at_local_energy_ptr(i) = my_local_e_1(i) +  my_local_e_2(i)

      ! forces
      if (len_trim(calc_force) > 0) then
         Sp = Switching_Function_Deriv(rmin, r1, r2)
         ! F_1_scaled = (1-S)*F_1 - r_vector / r * S' * V_1
         ! F_2_scaled = S*F_2 + r_vector / r * S' * V_2
         my_f_1(1, i) = (1.0_dp - S) * my_f_1(1, i) - drmin(1) * Sp * my_local_e_1(i) / rmin
         my_f_1(2, i) = (1.0_dp - S) * my_f_1(2, i) - drmin(2) * Sp * my_local_e_1(i) / rmin
         my_f_1(3, i) = (1.0_dp - S) * my_f_1(3, i) - drmin(3) * Sp * my_local_e_1(i) / rmin
         my_f_2(1, i) = S * my_f_2(1, i) + drmin(1) * Sp * my_local_e_2(i) / rmin
         my_f_2(2, i) = S * my_f_2(2, i) + drmin(2) * Sp * my_local_e_2(i) / rmin
         my_f_2(3, i) = S * my_f_2(3, i) + drmin(3) * Sp * my_local_e_2(i) / rmin
         at_force_ptr(1, i) = my_f_1(1, i) + my_f_2(1, i)
         at_force_ptr(2, i) = my_f_1(2, i) + my_f_2(2, i)
         at_force_ptr(3, i) = my_f_1(3, i) + my_f_2(3, i)
      end if

      ! TODO virials
      if (len_trim(calc_virial) > 0) then
        RAISE_ERROR('Potential_Switch_calc: virials not yet implemented..', error)
        call assign_property_pointer(at, trim(calc_local_virial), at_local_virial_ptr, error=error)
        PASS_ERROR(error)
        my_local_virial_1 = at_local_virial_ptr
      end if

    end do


    ! save total quantities and save what we want
    ! final modified quantities saved as at_local_energy_ptr, at_force_ptr
    my_e_1 = sum(my_local_e_1)
    my_e_2 = sum(my_local_e_2)
    energy = my_e_1 + my_e_2
    call print("Potential_switch my_e_1_switch " // my_e_1, PRINT_VERBOSE)
    call print("Potential_switch my_e_2_switch " // my_e_2, PRINT_VERBOSE)
    call set_param_value(at, trim(calc_energy), energy)
    if (len_trim(calc_energy) > 0) then
      ! TODO does this store the right numbers (and where?)
      if (store_contributions) call set_param_value(at, trim(calc_energy)//"_pot1", my_e_1)
      if (store_contributions) call set_param_value(at, trim(calc_energy)//"_pot2", my_e_2)
    end if
    if (len_trim(calc_local_energy) > 0) then
      ! TODO what does this do, and does it save the right energies? Add pot2
      if (store_contributions) call add_property(at, trim(calc_local_energy)//"_pot1", at_local_energy_ptr, overwrite=.true.)
    endif
    if (len_trim(calc_force) > 0) then
       ! TODO does this store the right numbers (and where?)
       if (store_contributions) call add_property(at, trim(calc_force)//"_pot1", at_force_ptr, overwrite=.true.)
    endif
    if (len_trim(calc_virial) > 0) then
       call get_param_value(at, trim(calc_virial), my_virial_1)
       if (store_contributions) call set_param_value(at, trim(calc_virial)//"_pot1", my_virial_1)
    endif
    if (len_trim(calc_local_virial) > 0) then
       if (store_contributions) call add_property(at, trim(calc_local_virial)//"_pot1", at_local_virial_ptr, overwrite=.true.)
    endif


    if (allocated(my_local_e_1)) deallocate(my_local_e_1)
    if (allocated(my_local_e_2)) deallocate(my_local_e_2)
    if (allocated(my_f_1)) deallocate(my_f_1)
    if (allocated(my_f_2)) deallocate(my_f_2)
    if (allocated(my_local_virial_1)) deallocate(my_local_virial_1)
    if (allocated(my_local_virial_2)) deallocate(my_local_virial_1)

  end subroutine Potential_Switch_Calc

  recursive function Potential_Switch_Cutoff(this)
    type(Potential_Switch), intent(in) :: this
    real(dp) :: potential_switch_cutoff

    if(associated(this%pot1) .and. associated(this%pot2)) then
       potential_switch_cutoff = max(cutoff(this%pot1), cutoff(this%pot2))
    else
       potential_switch_cutoff = 0.0_dp
    endif

  end function Potential_Switch_Cutoff

  ! Fermi switching
  !real(dp) function Switching_Function(r)
    !!implicit none
    !real(dp), intent(in) :: r
    !real(dp) :: bf, rf

    !bf = 10.0_dp
    !rf = 1.0_dp
    !Switching_Function = 1.0_dp / (1.0_dp + exp(-bf * (r - rf)))

  !end function Switching_Function

  !real(dp) function Switching_Function_Deriv(r)
    !!implicit none
    !real(dp), intent(in) :: r
    !real(dp) :: bf, rf

    !bf = 10.0_dp
    !rf = 1.0_dp
    !Switching_Function_Deriv = bf * exp(-bf * (r - rf)) / (1.0_dp + exp(-bf * (r - rf)))**2

  !end function Switching_Function_Deriv

  ! cosine 'cutoff-function' switching
  !real(dp) function Switching_Function(r)
    !!implicit none
    !real(dp), intent(in) :: r
    !real(dp) :: r1, r2, pi

    !pi = 4.0_dp * atan(1.0_dp)
    !! hard-coded switching interval
    !r1 = 0.5_dp
    !r2 = 1.0_dp
    !if (r <= r1) then
      !Switching_Function = 0.0_dp
    !else if (r >= r2) then
      !Switching_Function = 1.0_dp
    !else
      !Switching_Function = 0.5_dp * (cos(pi * (r - r1) / (r2 - r1) + pi) + 1.0_dp)
    !end if

  !end function Switching_Function

  !real(dp) function Switching_Function_Deriv(r)
    !!implicit none
    !real(dp), intent(in) :: r
    !real(dp) :: r1, r2, pi

    !pi = 4.0_dp * atan(1.0_dp)
    !! hard-coded switching interval
    !r1 = 0.5_dp
    !r2 = 1.0_dp
    !if (r <= r1) then
      !Switching_Function_Deriv = 0.0_dp
    !else if (r >= r2) then
      !Switching_Function_Deriv = 0.0_dp
    !else
      !Switching_Function_Deriv = 0.5_dp * pi * sin(pi * (r1 - r) / (r1 - r2)) / (r2 - r1)
    !end if

  !end function Switching_Function_Deriv

  ! Perriot polynomial cutoff function (inverted)
  real(dp) function Switching_Function(r, r1, r2)
    !implicit none
    real(dp), intent(in) :: r, r1, r2
    real(dp) :: chi

    chi = (r - r1) / (r2 - r1)
    if (r < r1) then
      Switching_Function = 0.0_dp
    else if (r > r2) then
      Switching_Function = 1.0_dp
    else
      Switching_Function = 6.0_dp * chi**5 - 15.0_dp * chi**4 + 10.0_dp * chi**3
    end if

  end function Switching_Function

  real(dp) function Switching_Function_Deriv(r, r1, r2)
    !implicit none
    real(dp), intent(in) :: r, r1, r2
    real(dp) :: chi

    chi = (r - r1) / (r2 - r1)
    if (r < r1) then
      Switching_Function_Deriv = 0.0_dp
    else if (r > r2) then
      Switching_Function_Deriv = 0.0_dp
    else
      Switching_Function_Deriv = 30.0_dp * chi**4 - 60.0_dp * chi**3 + 30.0_dp * chi**2
    end if

  end function Switching_Function_Deriv
