%% This script takes the necessary input data for running the 
%% FVDAM code

clear all;

% file_input_title=input('Read the title of the input datafile with the extension (i.e. file type viz. .fgm/.dat/.inp etc..)\n','s');
% LOP=input('\nRead the loading option.\n');
% theta=input('Read the angle between the global-x direction and principle-1 direction (in degrees).\n');
% load_initial=input('\nRead the starting load for the given LOP\n');
% load_max=input('\nRead the maximum load for the given LOP\n');
% incr_NUM=input('\nRead the number of increments to be used in reaching from initial load to maximum load\n');
% Input the loading increment numbers after which data needs to plotted/recorded.\n');
% t=input('Read the number of integration points to be used in the Legendre-Gauss Quadrature\n');
% para_M=input('Read the number of terms to be used in Legendre polynomial expansion\n');


file_input_title='generated_ruc_46_bal_insitu.fgm';      % Name of input data file with extension
LOP=4;                                        % Loading option in the principle material coordinate system
theta=0;                                      % Angle between global-x and principle-1
load_initial=0;                               % Initial macroscopic strain (corresponds to stress-free state)
load_max=0.005;                               % Final macroscopic strain
incr_NUM=25;                                  % Total number of steps to reach the final load from initial load.
incr_register=[10 25];                        % Register stresses, displacements, effective, hydrostatic and plastic strains at following increments;
t=5;                                          % integration points grid in each subcell 5 for 5x5, 3 for 3x3, etc..   
para_M=3;                                     % Number of Legendre ploynomials
iter_limit=25;                                % Maximum iterations allowed at each load step
error_tolerance=2;                            % Maximum (%) error allowed in the incremental plastic strain between final two successive iterations

%% DAMAGE MODELING

% % crack_num=input('Read the number of cracks\n');
% crack_num=0;
% crack_dir=[];
% crack_subcell_beta=[];
% crack_subcell_gama1=[];
% crack_subcell_gama2=[];
% crack_subcell_gama=[];
% crack_subcell_beta1=[];
% crack_subcell_beta2=[];
% for i=1:crack_num
%     crack_dir(i)=input(sprintf('Read the direction (2 for X_2 and 3 for X_3) of the %ith crack\n',i));
%     if crack_dir(i) == 2
%         crack_subcell_gama(i)=input(sprintf('Read the index (gamma) of the subcell on top of which the %ith crack lies\n',i));
%         crack_subcell_beta1(i)=input(sprintf('Read the beginning index (beta) fot the %ith crack\n',i));
%         crack_subcell_beta2(i)=input(sprintf('Read the ending index (beta) fot the %ith crack\n',i));
%     end
%     if crack_dir(i) == 3
%         crack_subcell_beta(i)=input(sprintf('Read the index (beta) of the subcell to the right of which the %ith crack lies\n',i));
%         crack_subcell_gama1(i)=input(sprintf('Read the beginning index (gama) fot the %ith crack\n',i));
%         crack_subcell_gama2(i)=input(sprintf('Read the ending index (gama) fot the %ith crack\n',i));
%     end
% end

crack_yes=0;                                             % Input 1 if crack/debondng present or 0 if absent
crack_num=0;                                             % Input total number of cracks if present else put 0
crack_dir=[];                                            % input 2 for horizontal crack and 3 for vertical crack, An array OF 2 and 3 for multiple cracks, empty if no crack
crack_subcell_beta=[];                                   % Subcell index (gamma) on top of which the i^th crack lies, empty if no crack
crack_subcell_gama1=[];                                  % Beginning index (beta) fot the i^th crack, empty if no crack
crack_subcell_gama2=[];                                  % Ending index (beta) fot the i^th crack, empty if no crack
crack_subcell_gama=[];                                   % Subcell index (beta) to the right of which the i^th crack lies, empty if no crack
crack_subcell_beta1=[];                                  % Beginning index (gama) fot the i^th crack, empty if no crack
crack_subcell_beta2=[];                                  % Ending index (gama) fot the i^th crack, empty if no crack


if length(findstr('.fgm',file_input_title)) > 0
    string_title=regexprep(file_input_title, '.fgm', '');
elseif length(findstr('.dat',file_input_title)) > 0
    string_title=regexprep(file_input_title, '.dat', '');
elseif length(findstr('.inp',file_input_title)) > 0
    string_title=regexprep(file_input_title, '.inp', '');
elseif length(findstr('.fvdam',file_input_title)) > 0
    string_title=regexprep(file_input_title, '.fvdam', '');
end

file_output_title=strcat(string_title, '.out');

fvdam_global;                                            % EXECUTES THE FVDAM CODE
