MODULE Equations
    IMPLICIT NONE

    REAL(8), PARAMETER :: DT = 0.001
    INTEGER, PARAMETER :: nSteps = 100000
    INTEGER, PARAMETER :: NDIM = 3
    REAL, PARAMETER :: S0(1) = (/28.0/)
    REAL, PARAMETER :: X0(NDIM) = (/1.0,1.0,1.0/)

    CONTAINS

SUBROUTINE Step(x, s)
    IMPLICIT NONE

    REAL(8), INTENT(inout) :: x(3)
    REAL(8), INTENT(in) :: s(1)

    REAL(8), PARAMETER :: dt = 0.001

    REAL(8) :: dx(3)
    dx(1) = 10 * (x(2) - x(1))
    dx(2) = x(1) * (s(1) - x(3)) - x(2)
    dx(3) = x(1) * x(2) - 8. / 3 * x(3)
    x(:) = x(:) + dt * dx(:)
END SUBROUTINE

SUBROUTINE Adjoint(x, y, s)
    IMPLICIT NONE
    REAL(8), INTENT(in) :: x(3)
    REAL(8), INTENT(inout) :: y(3)
    REAL(8), INTENT(in) :: s(1)

    REAL(8) :: dy(3)
    dy(1) = -10 * y(1) + (s(1) - x(3)) * y(2) + x(2) * y(3)
    dy(2) = 10 * y(1) - y(2) + x(1) * y(3)
    dy(3) = -x(1) * y(2) - 8./3 * y(3)
    y = y + DT * dy
END SUBROUTINE

REAL(8) FUNCTION gradContribution(x, y, s)
    REAL(8), INTENT(in) :: x(3), y(3), s(1)
    gradContribution = x(1) * y(2) * DT
END FUNCTION

END MODULE
