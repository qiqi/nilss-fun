PROGRAM Adj

    USE Equations

    IMPLICIT NONE

    INTEGER :: iStep, nStep0, nStep1, objective
    REAL(8) :: y(NDIM,2), history(NDIM,nSteps), s(1)
    REAL(8) :: dxdt(NDIM,2)
    REAL(8) :: adj_history(NDIM,nSteps)
    CHARACTER(len=128) :: arg
    REAL(8) :: zAvg, z, grad

    s = S0

    if (command_argument_count() .ne. 3) then
        print *, "Need 3 arguments"
        call exit(-1)
    end if

    CALL get_command_argument(1, arg)
    Read(arg, '(i10)') nStep1

    CALL get_command_argument(2, arg)
    Read(arg, '(i10)') nStep0

    CALL get_command_argument(3, arg)
    Read(arg, '(i10)') objective

    Open(1, file="../flow/solution.bin", form="unformatted", access="stream", &
            status="old")
    Read(1) history
    Close(1)

    Open(1, file="input.bin",form="unformatted",access="stream", status="old", &
         convert="big_endian")
    Read(1) y
    Close(1)

    ! Projection
    ! dxdt = (history(:,nStep1) - history(:,nStep1-1)) / DT
    ! dxdt_normalized = dxdt / sqrt(sum(dxdt * dxdt))
    ! y(:,1) = y(:,1) - sum(dxdt_normalized * y(:,1)) * dxdt_normalized
    ! y(:,2) = y(:,2) - sum(dxdt_normalized * y(:,2)) * dxdt_normalized

    ! Add time dilation contribtion
    ! if (objective .ne. 0) then
    !     zAvg = sum(history(3,:)) / nSteps
    !     z = history(3,nStep1)
    !     y = y + (zAvg - z) / 10 * dxdt / sum(dxdt * dxdt)
    !     print *, (zAvg - z) / 10 * dxdt / sum(dxdt * dxdt)
    ! end if

    ! Run adjoint backwards
    grad = 0.0
    DO iStep = nStep1-1, nStep0, -1
        grad = grad + gradContribution(history(1,iStep), y(:,1), s)
        grad = grad + gradContribution(history(1,iStep), y(:,2), s)
        CALL Adjoint2P(history(:, iStep), y(:,1), y(:,2), s)
        adj_history(:,iStep) = y(:,1)
        if (objective .ne. 0) then
            y(NDIM,1) = y(NDIM,1) + DT / 2
            y(NDIM,2) = y(NDIM,1) + DT / 2
        end if
    END DO

    dxdt = (history(:,nStep0+1:nStep0+2) - history(:,nStep0:nStep0+1)) / DT

    ! dxdt = (history(:,nStep0) - history(:,nStep0-1)) / DT
    ! dxdt_normalized = dxdt / sqrt(sum(dxdt * dxdt))
    ! y(:,1) = y(:,1) - sum(dxdt_normalized * y(:,1)) * dxdt_normalized
    ! y(:,2) = y(:,2) - sum(dxdt_normalized * y(:,2)) * dxdt_normalized
    
    ! Add time dilation contribution
    
    ! if (objective .ne. 0) then
    !     dxdt = (history(:,nStep0+1) - history(:,nStep0)) / DT
    !     z = history(3,nStep0)
    !     y = y - (zAvg - z) / 10 * dxdt / sum(dxdt * dxdt)
    !     print *, (zAvg - z) / 10 * dxdt / sum(dxdt * dxdt)
    !     print *, y
    ! end if

    Open(1, file="output.bin",form="unformatted",access="stream", &
         convert='big_endian')
    Write(1) y
    Close(1)

    Open(1, file="dxdt.bin",form="unformatted",access="stream", &
         convert='big_endian')
    Write(1) dxdt
    Close(1)

    Open(1, file="grad.txt",form="formatted",access="stream")
    Write(1, *) grad
    Close(1)

    ! Open(1, file="adj.txt",form="formatted",access="stream")
    ! Write(1, *) adj_history(:,nStep0:nStep1-1)
    ! Close(1)

CONTAINS
    SUBROUTINE Adjoint2P(x, y, yf, s)
        IMPLICIT NONE
        REAL(8), INTENT(in) :: x(NDIM)
        REAL(8), INTENT(inout) :: y(NDIM), yf(NDIM)
        REAL(8), INTENT(in) :: s(1)
        CALL Adjoint(x, y, s)
        CALL Adjoint(x, yf, s)
    END SUBROUTINE
END PROGRAM
