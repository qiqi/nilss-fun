MODULE Equations
    IMPLICIT NONE

    REAL(8), PARAMETER :: DT = 0.001
    INTEGER, PARAMETER :: nSteps = 100000
    INTEGER, PARAMETER :: NDIM = 2
    REAL, PARAMETER :: S0(1) = (/-1.0/)
    REAL, PARAMETER :: X0(NDIM) = (/1.0,1.0/)

    CONTAINS

SUBROUTINE Step(x, s)
    IMPLICIT NONE

    REAL(8), INTENT(inout) :: x(2)
    REAL(8), INTENT(in) :: s(1)

    REAL(8), PARAMETER :: dt = 0.001

    REAL(8) :: spd, dx(2)
    spd = exp(s(1) * x(2))
    dx(1) = x(2)
    dx(2) = -x(1) + (1 - x(1) * x(1)) * x(2)
    x(:) = x(:) + dt * dx(:) * spd
END SUBROUTINE

SUBROUTINE Adjoint(x, y, s)
    IMPLICIT NONE
    REAL(8), INTENT(in) :: x(2)
    REAL(8), INTENT(inout) :: y(2)
    REAL(8), INTENT(in) :: s(1)

    REAL(8) :: spd, dx(2), dy(2)
    spd = exp(s(1) * x(2))
    dx(1) = x(2)
    dx(2) = -x(1) + (1 - x(1) * x(1)) * x(2)

    dy(2) = y(1) + (1 - x(1) * x(1)) * y(2)
    dy(1) = -y(2) - 2 * x(1) * x(2) * y(2)
    dy(2) = dy(2) + s(1) * dot_product(dx, y)
    y = y + DT * dy * spd
END SUBROUTINE

REAL(8) FUNCTION gradContribution(x, y, s)
    REAL(8), INTENT(in) :: x(2), y(2), s(1)
    REAL(8) :: spd, dx(2)
    spd = exp(s(1) * x(2))
    dx(1) = x(2)
    dx(2) = -x(1) + (1 - x(1) * x(1)) * x(2)
    gradContribution = x(2) * dot_product(dx, y) * spd * DT
END FUNCTION

END MODULE
