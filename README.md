This is a simple toy 2D antiplane quasi-dynamic earthquake sequence simulator utilizing Runge-Kutta 4/5 to solve for the time history evolution of slip and stress.

The primary file is QDBIM2D which calls the ode file DieterichRuinaRegAging which implements the state evolution of slip, shear stress, slip rate and the state variable from the regularized Dieterich-Ruina rate-and-state friction law.

The additional ode files related to hacked versions of ode23 and ode45 to allow for the periodic output of the solution vector to disk, as otherwise the vector exceeds machine memory limits for large spatial domains and long time periods.
