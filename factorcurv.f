        double precision function Factorcurv (i)

        use mparameters
        use mgraft
        use mparameters_chain
        use mparameters_monomer

        implicit none
        integer i
        real*8 radio

        radio = posgraft(1)

        factorcurv = 1/(2*pi*(delta**3)*(dfloat(i) - 0.5d0))
        return
        end
 
