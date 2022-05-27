function [z, L, C, T] = CSD_model(K_g, g, K_b, b)
% CSD_model is the same function as the main.m in model, but only in
% function format. This helps in automating the comparison between model
% and experiments. Only the print functions are ommited, because this would
% only slow down calculations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hereafter code is copy pasted

%{
Assumptions:
    - Plug Flow
    - Steady State
    - Isothermic solution
    - No agglomeration
    - No breakage
    - Simple birth and growth terms
    - No seeding and no solids at the entrance
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L plane

% Standard value of 100. DO NOT CHANGE THIS VALUE ATM.
iter_L = 100;         % Amount of discretization in the L plane 
                      % (excl. the 0th point, because that is the boundary
                      % condition).
                      % Meaning that iter_L of 100 will lead to 101 bins.

size_L = 500*10^(-6); % This is the maximum size of the crystal. It is used
                      % together with iter_L to make intervals. In this
                      % case all intervals are equidistant 
                      % Unit: m_crystal

delta_L = size_L/ iter_L; % Discretization size of L-plane. This will be 
                          % the width of the bins in the L-plane. [m_crystal] 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Z plane

iter_Z = 10;          % Amount of discretization in the Z plane. 

size_Z = 0.10;        % The lenght of the tube in meters_space.

delta_Z = size_Z/iter_Z; % Discretization size of Z-plane. [m_space]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants

% Flow rate [m^3_space/s] 
Q = 1*10^-6 /60;          % Experimental value: Q = +- 1 ml/min

% Inner diameter of tubing [m_space]
D_in = 1.6*10^-3;

% Area of inner tubing [m^2_space]
A_in = pi*D_in/4;

% Volume of discretized piece of tubing [m^3_space]
V_in = A_in * delta_Z ;

% Average linear speed of fluid [m_space /s]
U = Q/A_in;

% Starting conditions for Growth [m_crystal/ s]  and Birth [(# of crystals)/ (m^3_space * s)]
G_0 = 1;
B_0 = 1;

% Starting condition for smallest fraction of crystals. Used as boundary
% condition in the ODE solver. Best set at 0. Otherwise G_0/B_0
n_0 = 0;

% Initial concentration of the solution [kg/ m^3_space]. Solution used
% supersaturates at 40 degrees Celcius. This uses 2.4464 g/100 ml solution.
C_0 = 24.464;

% Saturation concentration of the solution [kg/ m^3_space]. The correlating
% temperature is 20 degrees Celcius. Normally a saturated solution should
% have 1.182146 g/100 ml
C_sat = 12.66;

% Mass density of paracetamol [kg / m^3_crystal]. Information from KU
% Leuven database.
MASS_DENS = 1.293;

% Form factor [-]. Best to leave at 1. Unless another shape of crystals is
% assumed
K_V = 1;

%{
% Order of growth [-]. 
% Ranges between 1 and 2 (Rasmuson, Handbook of industrial crystallization p173)
% Similar values in Fujiwara paper.
g = 1;

% Order of birth [-]. 
% Primary nucleation between 5 and 15 (Rasmuson, Handbook of industrial crystallization p173)
% Secundary nucleation between 2 and 4 (Rasmuson, Handbook of industrial crystallization p173)
% Fujiwara: 5 - 7.5
b = 2; % Van Gerven uses 2

% Birth constant [variable]
% Fujiwara: between 8.298323E17 and 7.43139E21

k_sec1 = 1.15 * 10^11; %  (For Van Gerven use 1.15 * 10^11)
k_sec2 = NaN;
K_b = k_sec1 * C_sat^b; % Concentration should have power 'b'

% Growth constant [variable]
% When referencing Fujiwara: the used paper is Determination of the Kinetic Parameters for the Crystallization of Paracetamol
% from Water Using Metastable Zone Width Experiments.
% Fujiwara: between 4.79587E-3 and 5.613476E-2
%K_g = 5.613476E-10;

k_g1 = 3.21 * 10^(-5); % (For Van Gerven use 3.21 * 10^(-4))
k_g2 = 2.58 * 10^4;
R = 8.314;
T = 293; % K
K_g = k_g1 * exp(-k_g2 /R /T);

% K_g = 8.0693E-10; % Just takes the above value directly.
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the PBM

[z, L, C, T] = Solving_PBE(iter_L, iter_Z, size_L, size_Z, V_in, U, n_0, C_0, C_sat, MASS_DENS, K_V, K_g, K_b, g, b);

end