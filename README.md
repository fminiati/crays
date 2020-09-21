# crays
Multithreaded Lagrangian cosmic-rays transport code.
It time-integrates the Boltzmann's equation for transport of charged cosmic-ray particles in 
momentum space.
The input file specifies:
1. number of threads
2. number of resolution elements and duration of the integration
3. initial conditions, e.g. a time-dependent source term defining the injection mechanism, 
   which is either instantaneous (delta-function) or contniuous over a time period. 
4. time dependent from a plasma model which includes a specification of the plasma 
  density, temperature and turbulence parameters, used to calculate the atransport coefficienct, 
  such as energy losses and advection and diffusion in momentum space, determined . 
5. how often to compute the synchrotron emission spectrum and its frequency range.

