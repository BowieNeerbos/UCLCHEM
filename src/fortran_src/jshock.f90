
!J-shock paramterization
!Based on James et al. 2019 A&A 634
!https://ui.adsabs.harvard.edu/abs/2020A%26A...634A..17J/abstract
MODULE jshock_mod
    USE physicscore
    USE network
    USE constants
    IMPLICIT NONE
 
    
    integer :: coflag

    REAL(dp) :: tstart,maxTemp,vMin,mfp,tCool,tShock,d,dMax,maxDens
    REAL(dp) :: t_lambda, n_lambda

    REAL(dp) :: z2,vs,v0,zn,vn,at,z3,tsat
    REAL(dp) :: ucm,z1,driftVel,vi,tempi,vn0,zn0,vA,dlength
    REAL(dp) :: grainRadius5,dens6,grainNumberDensity,dzv,start_vel
    REAL(dp), allocatable :: tn(:),ti(:),tgc(:),tgr(:),tg(:)
    !variables for the collisional and radiative heating of grains
    REAL(dp) :: mun,tgc0,Frs,tgr0,tgr1,tgr2,tau100,trs0,G0
    REAL(dp) :: coshinv1,coshinv2,zmax,a1,eta,eps,epso,sConst

    integer :: inrad,projectiles(6)
    REAL(dp), PARAMETER ::nu0=3.0d15,K_BOLTZ_CGS=1.38d-16,bm0=1.e-6,bt=6.
    REAL(dp), PARAMETER :: GAS_DUST_NUMBER_RATIO=1.14d-12,CODES_TEMP=130.0
    REAL(dp), PARAMETER :: grainRadius=1.0d-5
    !*******************************************************************

CONTAINS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Checks inputs make sense and then calculates a few constants and!
    ! sets up variables for the shock paramterization that follows    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE initializePhysics
        INTEGER :: iLoop

        !Reset variables for python wrap.
        coflag=0 !should reset sputtering
        
        cloudSize=(rout-rin)*pc

        if (freefall) THEN
            write(*,*) "Cannot have freefall on during jshock"
            Write(*,*) "setting freefall=0 and continuing"
            freefall=.False.
        ENDIF
        density=initialDens

        ! Determine the maximum temperature
        maxTemp = (5e3)*(vs/10)**2


        ! Determine minimum velocity
        vMin = ((-2.058e-07*(vs**4) + 3.844e-05*(vs**3) - 0.002478*(vs**2) + 0.06183*(vs) - 0.4254)**2)**0.5

        ! Determine the shock width (of the order of the mean free path)
        mfp = ((SQRT(2.0)*(1e3)*(pi*(2.4e-8)**2))**(-1))/1d4
        tShock = mfp/(vs*1d5)
        ! Determine shock width
        tCool = (1/initialDens)*1d6*(60*60*24*365)
        ! Determine the maximum density attained
        maxDens = vs*initialDens*(1d2)
        ! Determine the rate constants
        t_lambda = LOG(maxTemp/initialTemp)
        n_lambda = LOG(maxDens/initialDens)


        if (allocated(tn)) deallocate(tn,ti,tgc,tgr,tg)
        allocate(tn(points),ti(points),tgc(points),tgr(points),tg(points))
        mun=2*mh
        grainRadius5=grainRadius/4.e-5
        dens6=density(dstep)/1.e6
        currentTimeOld=0.0

        !Need to find the location of the sputtering projectiles in species arrays
        DO iLoop=1,SIZE(specName)
            IF (specName(iLoop).eq."H2") projectiles(1)=iLoop
            IF (specName(iLoop).eq."HE") projectiles(2)=iLoop
            IF (specName(iLoop).eq."C") projectiles(3)=iLoop
            IF (specName(iLoop).eq."O") projectiles(4)=iLoop
            IF (specName(iLoop).eq."SI") projectiles(5)=iLoop
            IF (specName(iLoop).eq."CO") projectiles(6)=iLoop
        END DO

        !tsat proportional to 1/pre-shock density. Fit to tsats from Jimenez-Serra 2008.
        tsat=(-15.38729*vs*vs*vs)+(2069.56962*vs*vs)-(90272.826991*vs)+1686858.54278
        tsat=tsat/initialDens        
    END SUBROUTINE initializePhysics

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Called every time loop in main.f90. Sets the timestep for the next output from   !
    !UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90 !
    !but the integrator itself chooses an integration timestep.                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updateTargetTime
        IF (timeInYears .gt. 1e6) THEN
            targetTime=(timeInYears+1e5)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 1.0d4) THEN
            targetTime=(timeInYears+1000)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 1.0d3) THEN
            targetTime=(timeInYears+100.)*SECONDS_PER_YEAR
        ELSE IF (timeInYears*SECONDS_PER_YEAR .lt. tShock) THEN
            targetTime=currentTime+0.05*tShock
        ELSE
            targetTime=1.1*currentTime
        END IF
        ! ELSE IF (timeInYears .gt. 1.) THEN
        !     targetTime=(timeInYears+1.)*SECONDS_PER_YEAR
        ! ELSE IF (timeInYears .gt. 0.001) THEN
        !     targetTime=(timeInYears+0.1)*SECONDS_PER_YEAR
        ! ELSE IF (timeInYears .gt. 0.00001) THEN
        !     targetTime=(timeInYears+0.0001)*SECONDS_PER_YEAR
        ! ELSE IF (timeInYears .gt. 0.000001) THEN
        !     targetTime=(timeInYears+0.0000001)*SECONDS_PER_YEAR
        ! ELSE IF  (timeInYears.gt.0.0) THEN
        !     targetTime=(timeInYears+0.0000000001)*SECONDS_PER_YEAR
        ! ELSE
        !     targetTime=3.16d-03
        ! ENDIF
    END SUBROUTINE updateTargetTime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate shock properties for current time and set density, temperature and Av  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updatePhysics
        !calculate column density. Remember dstep counts from edge of core in to centre
        !calculate column density. Remember dstep counts from edge of core in to centre
        IF (dstep .lt. points) THEN
            !column density of current point + column density of all points further out
            coldens(dstep)=(cloudSize/real(points))*density(dstep)
            coldens(dstep)=coldens(dstep)+sum(coldens(dstep:points))
        ELSE
            coldens(dstep)=cloudSize/real(points)*density(dstep)
        END IF
      
        !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        av(dstep)= baseAv +coldens(dstep)/1.6d21

        ! Determine the shock velocity at the current time
        v0 = vs*(exp(LOG(vMin/vs)*(currentTime/(finalTime*60*60*24*365))))
        IF (v0 .lt. vMin) THEN
            v0 = vMin
        END IF

        ! Determine whether shock is still increasing the temperature
        ! Or whether it is in the post-shock cooling phase
        ! Or whether the temperature is now constant
        IF (currentTime .le. tShock) THEN
            tn(dstep) = ((currentTime/tShock)**2)*(maxTemp) + initialTemp
            density = (((currentTime/tShock)**3)*(4*initialDens))
            WHERE (density .lt. initialDens) density = initialDens
        ELSE IF (currentTime .gt. tShock .AND. currentTime .le. tCool) THEN
            ! Otherwise we're in the cooling phase
            tn(dstep) = maxTemp*EXP(-t_lambda*(currentTime/(tCool)))
            density = (4*initialDens)*EXP(n_lambda*(currentTime/(tCool)))

            ! Ensure the gas does not cool below around 10 K
            IF (tn(dstep) .le. 10) THEN
                tn(dstep) = 10
            END IF

            where(density .gt. maxDens) density = maxDens
        ELSE
            tn(dstep) = 10
            density = maxDens
        END IF
        gasTemp(dstep)=tn(dstep)
        dustTemp(dstep)=gasTemp(dstep)
        IF (timeInYears .gt. 0) THEN
            write(92,1234) tn(dstep),density(dstep),timeInYears
            1234 format(3(e16.9))
        ENDIF
    END SUBROUTINE updatePhysics


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This subroutine must be in every physics module.                                !
    ! It receives the abundance array and performs any sublimation related activity   !
    ! In hot core that means following thermalEvaporation subroutine.                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE sublimation(abund)
        REAL(dp) :: abund(nspec+1,points)
        INTENT(INOUT) :: abund

        IF (coflag .ne. 1) THEN
            IF (gasTemp(dstep) .gt. CODES_TEMP) THEN
                coflag=1
                abund(gasIceList,dstep)=abund(gasIceList,dstep)+abund(iceList,dstep)
                abund(iceList,dstep)=1d-30
            ELSE
                IF ((sum(abund(iceList,dstep)) .gt. 1d-25) .AND. (driftVel .gt. 0)) CALL sputtering(abund)
            END IF
        END IF
        WHERE(abund<1.0d-30) abund=1.0d-30
    END SUBROUTINE sublimation


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine that will sputter the ices based in their entirety
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE sputtering(abund)
        REAL(dp) :: abund(nspec+1,points)
        INTENT(INOUT) :: abund
        abund(gasIcelist,dstep)=abund(gasIcelist,dstep)+abund(iceList,dstep)
        abund(iceList,dstep)=abund(iceList,dstep)-abund(iceList,dstep)
    END SUBROUTINE

END MODULE jshock_mod
