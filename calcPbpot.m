function [Pbpot units] = calcPbpot(Epur, Es, Pbs, units)


Pbpot = Pbs*(1-exp(-Epur/Es));      %mgC/mgChl/time  

units.Pbpot = units.Pbs;