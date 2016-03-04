	FUNCTION deltat( tempr )

c    *******************************************************************
c    ** calculates width of polymer brush as a function of temp       **
c    ** it gives delta in terms of radius of the core                 **
c    *******************************************************************

        double precision deltat, tempr, tstar

        tstar = ( tempr - 5 ) / 80d0
        deltat = ( 4.5 + 87.36 * exp ( -exp ( 1.6 * tstar + 0.15 ) ) )

        end
c
c