## ============================================ ##

function FGStark_Nstep_propagator( X0, acc, dtau, N, order, 
    cSun, caseSun, alphaSun, caseStark, divSafeguard, 
    accuracyFlag, Xf, accuracyOut, statusFlag ) 

    # dealing with acceleration  
    q = 10^(-3) ; 

    Xk   = X0 ; 
    Xkp1 = Xf ; 
    X_hist        = Xk ; 
    accuracy_hist = accuracyOut ; 
    status_hist   = statusFlag ; 
    acc_hist      = [  ] ; 

    for i = 1 : N 

        println("i = ", i)

        # rv to oe 
        oe_k = cart2kep( Xk[1:6], 1 ) ; 
        nu_k = oe_k[6] ; 
        r_norm = norm( Xk[1:3] ) ;

        # orient acc in same direction as velocity 
        v_vec = Xk[4:6] ; 
        acc   = q / r_norm * ( 1 + cos(nu_k) ) * v_vec ;  
        acc   = 0 * acc ; 

        # propagate one dtau step 
        Xkp1, accuracyOut, statusFlag = FGStark_oneStep_propagator( Xk, acc, dtau, order, 
            cSun, caseSun, alphaSun, caseStark, divSafeguard, 
            accuracyFlag, Xkp1, accuracyOut, statusFlag ) ;     

        # save hist 
        X_hist        = [ X_hist ; Xkp1 ] ; 
        accuracy_hist = [ accuracy_hist ; accuracyOut ] ; 
        status_hist   = [ status_hist ; statusFlag ] ; 
        acc_hist      = [ acc_hist ; norm(acc) ] ; 

        # ready for next iter 
        Xk = Xkp1 ; 

    end 

    return X_hist, accuracy_hist, status_hist, acc_hist

end 

export FGStark_Nstep_propagator 

## ============================================ ##

function FGStark_oneStep_propagator( 
    X0, acc, dtau, order, 
    cSun, caseSun, alphaSun, caseStark, divSafeguard, 
    accuracyFlag, Xf, accuracyOut, statusFlag )

    MAXORD = 30 ; 

    #         !-- Declaration of variables
    #         ! I/O
    #         real*8, intent(in)  :: X0[7]          ! X0: Initial Conditions. Dimension 7x1, containing r0 (3D), v0 (3D), t0 
    #         real*8, intent(in)  :: acc[3], dtau   ! acc: Stark acceleration (3D), dtau: delta tau for integration
    #         real*8, intent(in)  :: cSun           ! Sundman transformation constant (dt = cSun * r**alphaSun * dtau). Typically cSun = 1
    #         integer, intent(in) :: caseSun        ! Select a case for the Sundman transformation: 0: alpha = 0, tau = time
    #                                               !     1: alpha = 1, tau = eccentric anomaly ; This is the case resulting in best efficiency of the series
    #                                               !     2: alpha = 2, tau = true anomaly
    #                                               !     3: alpha = 3/2, tau = intermediate anomaly
    #                                               !     4: alpha = alphaSun, generic Sundman transformation: user provided alpha
    #         integer, intent(in) :: order          ! order: order of the desired Taylor Series integration
    #         real*8, intent(in)  :: alphaSun       ! Value of alpha when case 4 is selected. For other cases, alphaSun is not used
    #         logical, intent(in) :: caseStark      ! Indicates whether a Stark or Kepler propagation should be computed
    #         logical, intent(in) :: divSafeguard   ! Indicates whether divergence tests should be realized
    #         logical, intent(in) :: accuracyFlag[4]! Indicates which accuracy measurements are to be computed (1 = Hamiltonian, 2 = delR, 3 = delV, 4 = delT)
    # 
    #         real*8, intent(out) :: Xf[7]          ! Xf: Output vector, dimension 7x1
    #         real*8, intent(out) :: accuracyOut[4] ! Accuracy measurements: change in Hamiltonian, last term error estimates for r, v, p
    #         integer, intent(out) :: statusFlag    ! Indicates possible overflow/divergence of the series. If statusFlag = -1, the routine exited before finishing the
    #                                               ! computation. The step size has to be reduced.


    #         ! Internal variables
    #         integer, parameter  :: rng = RANGE(1.e0)
    #         real*8, parameter  :: frst = 10.e0**(rng-10)
    # 
    #         real*8  :: dtaun, dtaunprime, normR, v2
    #         real*8  :: hamil, hamilPlus1
    #         real*8  :: delRvec[3], delVvec[3]
    #         real*8  :: Ftot, Gtot, Htot, Ttot, dtaunprimej 
    #         real*8  :: Ftotprime, Gtotprime, Htotprime
    #         real*8  :: first[4], new[4], newprime[3], numerator
    #         integer :: j, i, k, set[4]
    #         real*8  :: FGHT(4*MAXORD), coeff[4]

    ## internal variables 

    rng  = 307 ; 

    frst = 10.e0^(rng-10) ; 

    dtaun = 0 ; dtaunprime = 0 ; normR = 0 ; v2 = 0 ; 

    hamil = 0 ; hamilPlus1 = 0 ; 

    delRvec = zeros(3,1) ; delVvec = zeros(3,1) ; 

    Ftot  = 0 ; Gtot = 0 ; Htot = 0 ; Ttot = 0 ; dtaunprimej = 0 ; 

    Ftotprime = 0 ; Gtotprime = 0 ; Htotprime = 0 ; 

    first = zeros(4,1) ; new = zeros(4,1) ; newprime = zeros(3,1) ; numerator = 0 ; 

    j = 0 ; i = 0 ; k = 0 ; set = zeros(4,1) ; 

    FGHT = zeros( 4 * MAXORD, 1 ) ; coeff = zeros(4,1) ; 


    ## Begin subroutine

    # Variables needed for overflow safeguard
    # first = frst ; 
    # set   = 0 ; 
    # statusFlag  = 0 ; 
    # accuracyOut = 0.e0 ; 

    # Initialize the sums
    Ftot = 0.e0 ; 
    Gtot = 0.e0 ; 
    Htot = 0.e0 ; 
    Ttot = 0.e0 ; 
    Ftotprime = 0.e0 ; 
    Gtotprime = 0.e0 ; 
    Htotprime = 0.e0 ; 

    dtaun = 1.e0 ; 
    dtaunprime = 1.e0 ; 

    # Compute the Hamiltonian
    if accuracyFlag[1] 
        normR = sqrt( X0[1]*X0[1] + X0[2]*X0[2] + X0[3]*X0[3] ) ; 
        v2    = X0[4]*X0[4] + X0[5]*X0[5] + X0[6]*X0[6] ; 
        hamil = 0.5e0 * v2 - 1.e0/normR - dot(X0[1:3], acc) ;     
    end 

    # Call the Maple generated file to get the values of
    # the F&G Stark series coefficients. Call depending on the
    # caseSun and caseStark 
    if caseSun == 0 
        FGHT = starkCoeffs_time(X0, acc, order, cSun, FGHT) ; 

    elseif caseSun == 1 # The independent variable is proportional to the eccentric anomaly
        FGHT = starkCoeffs_ecc(X0, acc, order, cSun, FGHT) ; 

    elseif caseSun == 2 # The independent variable is proportional to the true anomaly
        FGHT = starkCoeffs_true(X0, acc, order, cSun, FGHT) ; 

    elseif caseSun == 3 # The independent variable is proportional to the intermediate anomaly
        FGHT = starkCoeffs_inter(X0, acc, order, cSun, FGHT) ; 

    elseif caseSun == 4 # Generic Sundman transformation: dt = cSun * r**alphaSun * dtau
        FGHT = starkCoeffs_generic(X0, acc, order, cSun, alphaSun, FGHT) ; 

    else
        println("ERROR: Please select a Sundman case between 0 and 4") ; 

    end  
        

    # Form the Taylor Series from the coefficients
    k = 1 ; 
    for j = 1 : order

        dtaun = dtaun*dtau ; 

        # Compute position TS
        coeff[1] = FGHT[k] ; 
        coeff[2] = FGHT[k+1] ; 
        coeff[3] = FGHT[k+2] ; 
        coeff[4] = FGHT[k+3] ; 
        k = k + 4 ; 

        new[1] = coeff[1]*dtaun ; 
        new[2] = coeff[2]*dtaun ; 
        new[3] = coeff[3]*dtaun ; 
        new[4] = coeff[4]*dtaun ; 


        if (divSafeguard) 
            
            # Heuristically check convergence of the series. If the new
            # term in the series is bigger than 100* the very first term, we
            # exit the integrator
            for i = 1 : 4

            # Set the "first" value (first non zero value taken by each
            # series)
            if ( (set[i] == 0) && ( abs(new[i]) > 0.e0 ) ) 
                first[i] = 100.e0*(1.e0 + abs(new[i])) ; 
                set[i] = 1 ; 
            end 

            # Compare the new value to "first"
            if (abs(new[i]) > first[i]) 
                statusFlag = -1 ; 
                break 
            end 
            
            end 
        end 

        Ftot = Ftot + new[1] ; 
        Gtot = Gtot + new[2] ; 
        Htot = Htot + new[3] ; 
        Ttot = Ttot + new[4] ; 

        # Compute velocity TS: derivative of the position.
        dtaunprimej = dtaunprime*j ; 

        newprime[1] = coeff[1]*dtaunprimej ; 
        newprime[2] = coeff[2]*dtaunprimej ; 
        newprime[3] = coeff[3]*dtaunprimej ; 

        Ftotprime = Ftotprime + newprime[1] ; 
        Gtotprime = Gtotprime + newprime[2] ; 
        Htotprime = Htotprime + newprime[3] ; 

        dtaunprime = dtaunprime*dtau ; 

    end 

    # Finish the TS integration: multiply by the basis vectors

    # r = F*r0 + G*v0 + H*thrust
    Xf[1:3] = X0[1:3] + Ftot*X0[1:3] + Gtot*X0[4:6] + Htot * acc[:] ; 

    # v = Fdot*r0 + Gdot*v0 + Hdot*thrust
    Xf[4:6] = Ftotprime * X0[1:3] + Gtotprime * X0[4:6] + Htotprime * acc[:] ; 

    # Careful: to go from r to v, differentiation wrt 
    # time. However, the coefficients have been
    # differentiated with respect to tau. v = dr/dt =
    # dr/dtau*dtau/dt = dr/dtau * 1/(cSun * r**alphaSun)
    normR = sqrt( Xf[1]*Xf[1] + Xf[2]*Xf[2] + Xf[3]*Xf[3] ) ; 

    if caseSun == 0
        numerator = (cSun) ; 
    elseif caseSun == 1 
        numerator = (cSun*normR) ; 
    elseif caseSun == 2 
        numerator = (cSun*normR^2) ; 
    elseif caseSun == 3 
        numerator = (cSun*sqrt(normR^3)) ; 
    elseif caseSun == 4 
        numerator = (cSun*normR^alphaSun) ; 
    end 

    Xf[4:6] = Xf[4:6] / numerator ; 

    # time = time0 + dtime
    Xf[7] = X0[7] + Ttot ; 
        
    # Accuracy computations:
    if (accuracyFlag[1]) 
        # Compute new Hamiltonian
        v2 = Xf[4]*Xf[4] + Xf[5]*Xf[5] + Xf[6]*Xf[6] ; 
        hamilPlus1 = 0.5e0 * v2 - 1.e0 / normR -dot(Xf[1:3], acc[:]) ; 
        accuracyOut[1] = abs(hamilPlus1 - hamil) ; 
    end 
    if (accuracyFlag[2]) 
        delRvec = new[1] * X0[1:3] + new[2] * X0[4:6] + new[3] * acc[:] ; 
        accuracyOut[2] = sqrt( delRvec[1]*delRvec[1] + delRvec[2]*delRvec[2] + delRvec[3]*delRvec[3]) ; 
    end 
    if (accuracyFlag[3]) 
        delVvec = ( newprime[1] * X0[1:3] + newprime[2] * X0[4:6] + newprime[3] * acc[:] ) / numerator ; 
        accuracyOut[3] = sqrt( delVvec[1]*delVvec[1] + delVvec[2]*delVvec[2] + delVvec[3]*delVvec[3] ) ; 
    end 
    if (accuracyFlag[4]) 
        accuracyOut[4] = abs(new[4]) ; 
    end 

    return Xf, accuracyOut, statusFlag 

end 

export FGStark_oneStep_propagator 

## ============================================ ##

function taup_caseSun( caseSun, a, e, mu, cSun ) 

    # tau period 
    n = sqrt( mu / a^3 ) ; 

    # 0: alpha = 0, tau = time 
    if caseSun == 0 
        tau_p = 2 * pi / ( n * cSun ) ; 
        alphaSun = 0 ; 

    # 1: alpha = 1, tau = eccentric anomaly 
    elseif caseSun == 1 
        tau_p = 2 * pi / ( n * cSun * a ) ;  
        alphaSun = 1 ; 

    # 2: alpha = 2, tau = true anomaly 
    elseif caseSun == 2 
        tau_p = 2 * pi / ( n * cSun * sqrt( a * ( 1 - e^2 ) ) ) ; 
        alphaSun = 2 ; 

    # 3: alpha = 3/2, tau = intermediate anomaly 
    elseif caseSun == 3 
        M     = 2 * e / ( 1 + e ) ; 
        tau_p = 4 * ellipke( M ) / ( cSun * sqrt( mu * ( 1 + e ) ) ) ; 
        alphaSun = 3/2 ;     

    # 4: alpha = alphaSun, generic Sundman transformation: user provided alpha  
    else 
        tau_p = 2 * pi ; 

    end 

    return tau_p, alphaSun
end 

export taup_caseSun 
