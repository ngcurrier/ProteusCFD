module private_complex_module
  
  type private_complex
     private
     real :: real, imaginary
  end type private_complex

  interface operator(*)
     module procedure pc_mult
  end interface operator(*)

  contains
  
  subroutine pc_init(this,real,imaginary)
    ! initialize private_complex variable from reals
    type (private_complex), intent(out) :: this
    real, intent(in) :: real, imaginary    
    this%real = real
    this%imaginary = imaginary
  end subroutine pc_init
  
  subroutine pc_display(this,c)
    ! display value of private_complex variable with label
    type (private_complex), intent(in) :: this
    character*(*), intent(in) :: c
    print *, c, this%real, this%imaginary
  end subroutine pc_display

  type (private_complex) function pc_mult(a,b)
    type (private_complex), intent(in) :: a,b
    pc_mult%real = a%real*b%real - a%imaginary*b%imaginary
    pc_mult%imaginary = a%real*b%imaginary + a%imaginary*b%real
  end function pc_mult
    
end module private_complex_module
      
program main
  use private_complex_module
  type (private_complex) :: a, b, c, d
  ! initialize sample variables
  call pc_init(a,1.,-1.)
  call pc_init(b,-1.,2.)
  ! perform multiplication
  c = pc_mult(a,b)
  d = a*b
  ! display result
  call pc_display(c,'c=')
  call pc_display(d,'d=')
  stop
end program main
