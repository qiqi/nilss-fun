PROGRAM Lorenz
    IMPLICIT NONE

    INTEGER, PARAMETER :: nSteps = 100000
    INTEGER :: iStep
    REAL(8) :: x(2), history(2,nSteps), s(1)

    s = (/0.0/)

    x(:) = (/1.0, 1.0/)
    DO iStep = 1, 1000
        CALL Step(x, s)
    END DO

    DO iStep = 1, nSteps
        history(:,iStep) = x(:)
        CALL Step(x, s)
    END DO

    Open(1, file="solution.bin", form="unformatted", access="stream")
    Write(1) history
    Close(1)
CONTAINS

SUBROUTINE Step(x, s)
    IMPLICIT NONE
    REAL(8), INTENT(inout) :: x(2)
    REAL(8), INTENT(in) :: s(1)

    REAL(8), PARAMETER :: dt = 0.001

    REAL(8) :: dx(2)
    dx(1) = x(2) * (1 + s(1) * x(1))
    dx(2) = (-x(1) + (1 - x(1) * x(1)) * x(2)) * (1 + s(1) * x(1))
    x(:) = x(:) + dt * dx(:)
END SUBROUTINE

END PROGRAM
