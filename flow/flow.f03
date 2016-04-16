PROGRAM Flow

    USE Equations

    IMPLICIT NONE

    INTEGER :: iStep
    REAL(8) :: x(NDIM), history(NDIM,nSteps), s(1)

    s = S0
    x(:) = X0
    DO iStep = 1, 1000
        CALL Step(x, s)
    END DO

    DO iStep = 1, nSteps
        history(:,iStep) = x(:)
        CALL Step(x, s)
    END DO

    Open(1, file="objective.bin", form="unformatted", access="stream", &
         convert='big_endian')
    Write(1) history(NDIM,:)
    Close(1)

    Open(1, file="solution.bin", form="unformatted", access="stream", &
         convert='big_endian')
    Write(1) history
    Close(1)

END PROGRAM
