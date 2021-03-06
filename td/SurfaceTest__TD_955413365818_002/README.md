SurfaceEnergyTestDriver
=======================

Surface energy test adapted for the OpenKIM ecosystem.  Calculates various and
low and high symmetry surfaces.  Using the energies at these cuts, we construct
a broken bond model.  This models essentially counts the step density around
each surface, creating cusps that capture all surface energies with very few
parameters.

This test outputs the typical properties as well as xyz files containing the 
reconstructed surfaces as well as a plot of the final broken bond fit.

By default, the test only calculates four surfaces, but can be extended to
any number of surfaces.  200 unique surfaces are included in IndexList.pkl
which is commented out in runner for the time being.
