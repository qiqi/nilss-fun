PROGRAM Verification

    Use Equations

    IMPLICIT NONE

    REAL(8), PARAMETER :: EPS=1E-5
    INTEGER, PARAMETER :: Steps = 20
    INTEGER :: iStep
    REAL(8) :: x(NDIM), y(NDIM), history(NDIM,Steps+1), s(1)
    REAL(8) :: dJds_fd, dJds_adj

    s = S0
    x = X0
    DO iStep = 1, Steps
        history(:,iStep) = x(:)
        CALL Step(x, s)
    END DO
    history(:,Steps+1) = x(:)

    s = S0
    s(1) = s(1) + EPS
    x(:) = X0
    DO iStep = 1, Steps
        CALL Step(x, s)
    END DO

    dJds_fd = (x(NDIM) - history(NDIM,Steps+1)) / EPS

    y(:) = 0.0
    y(NDIM) = 1.0
    dJds_adj = 0.0
    DO iStep = Steps, 1, -1
        x = history(:,iStep)
        dJds_adj = dJds_adj + gradContribution(x, y, s)
        CALL Adjoint(x, y, s)
    END DO

    PRINT *, dJds_fd, dJds_adj

END PROGRAM
