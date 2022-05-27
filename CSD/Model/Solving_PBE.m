function [t, L, C, T] = Solving_PBE(iter_L, iter_Z, size_L, size_Z, V_in, U, n_0, C_0, C_sat, MASS_DENS, K_V, K_g, K_b, g, b)
% Solver PBE: 
% Input:
%       - iter_L: Amount of discretization in L-plane  
%       - size_L: Biggest size that crystal can get [m_crystal]
%       - iter_z: Amount of discretizations over Z
%       - U: Average speed of fluid (Plug flow at this moment) [m(spatial)/s]
%       - n_0: Initial condition of the smallest fraction of crystals [(# of crystals)/(m_crystal * m^3_space)]
%       - C_0: Initial concentration of the solution [kg / m^3_space] 
%       - C_sat: Saturated concentration of a solution at the temperature
%                at which the crystallization happens. [kg/ m^3_space]
%       - MASS_DENS: Mass density  of paracetamol [kg/ m^3_crystal]
%       - K_V: Form factor [-]
%       - K_g: Growth constant [variable]
%       - K_b: Birth constant [variable]
%       - g: Order of growth [-]
%       - b: Order of birth [-]
%
%
% Output:
%       Nog niet bepaald

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the test section. All values in this section are just used for
% testing and should not be standing here if the test is over. It should be
% written as more readable code.

T_0 = 40 + 273;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Making a grid with equal sized bins
delta_L = size_L/iter_L;

% Making a grid with equal sized volumes in tube
delta_Z = size_Z/iter_Z;

% All Z values at which the ODE-Solve will evaluate the differential equations 
z_values = [0: delta_Z: size_Z];


% Initial conditions for all values that will be evaluated.
% The first value can be chosen, but at this moment in time all other
% values start as 0. This simply says that at time 0 (or axial space 0)
% there are no crystals bigger than the smallest fraction.
InitCond = [n_0, zeros(1, iter_L), C_0, T_0];


% The solver of the system of equations 
[t,y] = ode15s(@(t,y) Stelsel(t, y, U, delta_L, iter_L, size_L, V_in, delta_Z, C_0, C_sat, MASS_DENS, K_V, K_g, K_b, g, b), z_values, InitCond);


% Takes the values for sizes out of the output of the solver
L = y(:, 1:end-2);


% Takes the values for the concentrations out of the ouput from the solver
C = y(:, end-1);

% Takes the values of the temperature out of the output from the solver
T = y(:, end);

