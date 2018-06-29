"""
Shocks is a module that consists of canonical hydrodynamical 
equations used for modeling astrophysical plasmas


Classes:

Shocks -- A class that emcompasses a multitude of hydrodynamical 
          shock equations relevant for plasma calculations.


Functions:

rh_density -- Returns the density jump relation derived from the
              Rankine-Hugoniot jump relations

rh_temp --   Returns the temperature jump relations derived from 
             the Rankine-Hugoniot jump relations

"""


class MHD:
    """
    Stores results of NEI simulation involving shock dynamics
    """

    def __init__(self):
        pass

    def rh_density(self, init_dens, gamma, mach):
        """
        Returns the density ratio according to the Rankine-Hugoniot 
        jump conditions

        Parameters
        ------
        init_dens: ~astropy.units.Quantity,
                   The initial density of the plasma pre-shock
        gamma: float,
               The specific heats ratios of the system
        mach: int,
              The mach number of the system

        Returns
        ------
        den_ratio: array-like
                   The density solution to the mass conservation equation as 
                   defined by the Rankine-Hugoniot relations.
        """

        dens_ratio = ((gamma+1)*mach**2)/(2+(gamma-1)*mach**2)

        final_dens = dens_ratio*init_dens

        return final_dens 

    def rh_temp(self, init_temp, gamma, mach):
        """
        Returns the temperature ratio according to the Rankine-Hugoniot 
        jump conditions

        Parameters
        ------
        init_temp: ~astropy.units.Quantity,
                   The initial temperature of the plasma pre-shock
        gamma: float,
               The specific heats ratios of the system
        mach: int,
              The mach number of the system

        Returns
        ------
        final_temp: array-like
                   The temperature solutions to the energy conservation equation as 
                   defined by the Rankine-Hugoniot relations.
        """


        temp_ratio = ( ( ( (gamma + 1)+ 2 * gamma * (mach**2 - 1) ) *
                        ( (gamma + 1) + (gamma-1)*(mach**2 - 1)) )  / 
                        ((gamma + 1)**2 * mach**2 )   )

        final_temp = temp_ratio * init_temp
        
        return final_temp