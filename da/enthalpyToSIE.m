function SIE = enthalpyToSIE(enthalpy)

SIE = -enthalpy(1,:) - 18; % only transform the first component!
