This is WIP project.

Lots of NR codes operate on a (rectangular?) grid
and likewise they keep track of only the metric and coordinate information
with no information about the background coordinate system.

This is an attempt to reconstruct the surface of the coordinate space using only the numerical metric information.
I am integrating the basis vectors using the connection coefficients, and integrating the points in background coordinate space using the basis vectors.

I've tried 
* forward-Euler integration
* RK4 integration
* explicitly solving linear dynamic systems (assuming the connection is constant), which I only have available for certain classifications of matrices.

I am testing it against known metrics.
Right now just the polar metric.

Here's the difference in error between the analytical connection and the numerically computed connection:


![alt](conn%20numeric%20vs%20analytic.png)

So far I can get two rings of a polar surface to compute correctly, but the third one busts.
