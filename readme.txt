10.05.23

To do:

sw_reset_Photons_escape: set to .TRUE. and adjust main prog accordingly

Python: brems opacity; DONE (10.05.23)

Photon split; DONE (10.05.23)

Tambov change as photons escape

energy conservation (DONE)

even grid in mass rather than xi (DONE)

vary initial temperature

!!!###########################################################
!!!###### To write down:

Thermalization condition: how long is it valid. Informs simulation start time. Can one do analytics before that?

Temperature evolution: balance between free-free cooling and Compton heating?

Effective opacity/photosphere, resulting escaping spectrum

Bolometric light curve: looks like declining all the way, like t^-1 until t_diff


!!!####### Details to fix:

Zz, Zav, mu_mol, Gaunt factor

Increase mass somewhat (see eq. 13 in IM23)

!!!### For setting up realistic simulation:

What density profile to use? Simular to SNe (Narar, Sari 2010)? Probably not.

NS10: initial kinetic and thermal energies should be EQUAL! Simul 20 happens to be close!


Dep on initial radius, e.g. test_14 vs test_16: no strong effect since both get fully thermalized (see T(t) plot)

Dep in initial temp, e.g. 14c vs 16: X-ray LC-s very different since 16 hasn't thermalized early enough. Later, differences are smaller, but the high T simul doesn't completely thermalize (at least hasn't by 800 sec [pooleli]).








