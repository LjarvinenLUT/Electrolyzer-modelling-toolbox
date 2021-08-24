function [F,R,n_e] = getConstants()
% GETCONSTANTS Provides the accurate values of Faraday's constant,
% universal gas constant and number of electrons transferred in a single
% electrochemical water splitting reaction.

N_A = 6.02214076e23; % 1/mol, Avogadro constant
k_B = 1.380649e-23; % J/K, Boltzmann constant
Q_e = 1.602176634e-19; % C, Elemental charge

F = Q_e*N_A; % C/mol, Faraday constant
R = k_B*N_A; % J/K/mol, Universal gas constant
n_e = 2; % Number of electrons transferred in one reaction

end