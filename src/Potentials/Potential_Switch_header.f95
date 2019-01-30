
  public :: Potential_Switch
  type Potential_Switch
     type(MPI_context) :: mpi
     type(Potential), pointer :: pot1 => null() 
     type(Potential), pointer :: pot2 => null() 
     real(dp) :: r1, r2
  end type Potential_Switch

  interface Initialise
     module procedure Potential_Switch_Initialise
  end interface

  interface Finalise
     module procedure Potential_Switch_Finalise
  end interface

  interface Print
     module procedure Potential_Switch_Print
  end interface

  interface Cutoff
     module procedure Potential_Switch_Cutoff
  end interface

  interface Calc
     module procedure Potential_Switch_Calc
  end interface
