This is a simple quasi-dynamic boundary element code for performing 2D antiplane simulations of sequences of earthquakes and aseismic slip (SEAS) utilizing Runge-Kutta 4/5 adaptive time-stepping to solve for the time history evolution of slip and stress.

The main file is QDBIM2D which calls the ode function  DieterichRuinaRegAging that describes the state evolution of slip, shear stress, slip rate and the rate-and-state state variable following the regularized Dieterich-Ruina rate-and-state friction formulation using the aging law.

The physical problem setup in QDBIM2D is consistent with benchmark problem BP1-QD from the Southern California Earthquake Center Working Group for Advancing Simulations of Earthquakes and Aseismic Slip: https://strike.scec.org/cvws/seas/benchmark_descriptions.html

The include directory contains additional functions related to modified versions of ode23 and ode45 to allow for the periodic output of the solution vector to disk, as otherwise the state vector history can easily exceed machine memory limits for large spatial domains and long time periods.
