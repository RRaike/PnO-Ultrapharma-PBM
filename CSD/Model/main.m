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

clc
clf
clear

%{
b = 2;
g = 1;

% k_sec1 = 1.15 * 10^11;
% k_sec2 = NaN;
% K_b = k_sec1 * C_sat^b;
K_b = 1E18;

K_g = 10^-6*10^5;    % Just takes the above value directly.

[z, L, C, T] = CSD_model(K_g, g, K_b, b);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boolean to print figures or not.

% If Show_PBM_Figure is true, then the numerical solution of the PBM will
% be shown. All distributions are on one plot to be able to compare.
Show_PBM_Figure = true;

% If Gauss_Figure is True (or 1) then all Gaussian curves will be printed
% and saves separately.
Show_Gauss_Figure = false;

% The figures of Rsquare, Means and the Standard Deviations (STDs) show all
% of these values in 1 plot of their own. Unlike the Gauss_Figure it
% agglomerates information of multple Gaussian figures.
Show_Rsquare_Figure = false;
Show_Means_Figure = false;
Show_STDs_Figure = false;


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

iter_Z = 4;          % Amount of discretization in the Z plane. 

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the PBM

[z, L, C, T] = Solving_PBE(iter_L, iter_Z, size_L, size_Z, V_in, U, n_0, C_0, C_sat, MASS_DENS, K_V, K_g, K_b, g, b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting of the data reveived from the solver.

if Show_PBM_Figure
   Fig_PBM = figure();
    for i = 2: length(z)
        delta_L = size_L/(iter_L);
        plot(0:delta_L:size_L, L(i, 1:end),'DisplayName', num2str(i))
        axis([0, 1.5*10^-4, 0, 4.5*10^17])
    
        h = legend( {'25% through reactor', '50% through reactor', '75% through reactor', 'end of reactor'},...
            'location','northeast');
        xlabel("Length of crystal [m_{crystal}]")
        ylabel(["Number probability density function"; ...
                "[# crystals/ (m^3_{space} *m_{crystal})]"])
        %title(["Number probability density as function of"," crystal length at different points in the reactor"])

        stringInAnnotationHeaders = {"Order of birth:", 
                                  "Order of growth:",
                                  "Birth constant:",
                                  "Growth constant:"
                                  };
        annotation('textbox', [.56, .42, .3, .3], 'String', stringInAnnotationHeaders, 'FitBoxToText', 'on', 'LineStyle', 'none')

        stringInAnnotationValues = {string(b),
                                       string(g),
                                       string(sprintf('%.3e',K_b)),
                                       string(sprintf('%.3e' ,K_g))
                                       };
        annotation('textbox', [.76, .42, .3, .3], 'String', stringInAnnotationValues, 'FitBoxToText', 'on', 'LineStyle', 'none')
        hold on
    end


    fig_Concentration = figure();
    plot(z, C)

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit a Gaussian curve trough the data of the solver and get its
% characteristics.


% Initiation of characteristics. This optimizes the program.
Rsquare = zeros(iter_Z,1);
Means = zeros(iter_Z,1);
STDs = zeros(iter_Z,1);


% Iterate over every non-boundary condition. The boundary condition is the
% first value in the matrix and after which iter_Z iterations are done,
% thus iterate until iter_Z + 1.
for k = 2:(iter_Z+1)

L_bins = [0:delta_L:size_L].';

% Get the to-fit data as the n function gotten from the L-matrix. These are
% the values of n at a certain Z as function of L. 
Data = L(k,:).';


% Fit function. Built-in Matlab function.
[bell, goodness] = fit(L_bins, Data, 'gauss1');


% Save the values in the right variable.
Rsquare(k-1,1) = goodness.rsquare;
Means(k-1,1) = bell.b1;
STDs(k-1,1) = bell.c1/sqrt(2);


if Show_Gauss_Figure 
    Gauss_fig = figure;
    hold on 
    plot(L_bins, Data);
    plot(bell);
    %saveas(Gauss_fig,sprintf('Gauss %d.png',k));
    %hold off
    %clf
end

end


if Show_Rsquare_Figure
% Plotting R_squared values as scatterplot. The X-axis are the Z-values
% over which is iterated. The Y-axis displays the fractional R^2 value.
    Rsquare_fig = figure;
        scatter([delta_Z: delta_Z: size_Z], Rsquare)
        xlabel("Z-value; Spatial coordinate in tube [m_{space}]")
        ylabel("R^2 value corresponding with n(Z, L)")
        title("R^square value of Gaussian fit on n(Z, L)")
end

if Show_Means_Figure
    Means_fig = figure;
        scatter([delta_Z: delta_Z: size_Z], Means)
        xlabel("Z-value; Spatial coordinate in tube [m_{space}]")
        ylabel("Mean value of the Gaussian distribution [m_{crystal}]")
        title("Mean value of Gaussian fit on n(Z, L)")
end

if Show_STDs_Figure
    STDs_fig = figure;
        scatter([delta_Z: delta_Z: size_Z], STDs)
        xlabel("Z-value; Spatial coordinate in tube [m_{space}]")
        ylabel("Standard deviation of the Gaussian distribution [m_{crystal}]")
        title("Standard deviation of Gaussian fit on n(Z, L)")
end



