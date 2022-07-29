function my_function()
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
    LOP=2;                                        % Loading option in the principle material coordinate system
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
    
    tic; % time starts ticking

fid=fopen(file_output_title,'w');   % CREATING THE OUTPUT FILE CONTAINING THE LOGO, RELEVANT INPUT AND OUTPUT DATA.
% %__________________________________________________________________________________________________________________________________
fprintf(fid,'\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf(fid,'\t|\t\tFINITE VOLUME DIRECT AVERAGING MICROMECHANICS (FVDAM) MODEL\t\t|\n');
fprintf(fid,'\t|\t\t\t\t\t\t\t\t\t\t\t\t\t|\n');
fprintf(fid,'\t| DEVELOPED BY:\t\t\t\t\t\t\t\t\t\t\t|\n');
fprintf(fid,'\t|\t\t\tYOGESH BANSAL (GRADUATE RESEARCH ASSISTANT)\t\t\t|\n');
fprintf(fid,'\t|\t\t\t\t\tunder the supervision of\t\t\t\t|\n');
fprintf(fid,'\t|\t\t\tMAREK-JERZY PINDERA (PROFESSOR OF APPLIED MECHANICS)\t\t|\n');
fprintf(fid,'\t|\t\t\tDEPARTMENT OF CIVIL ENGINEERING, UNIVERSITY OF VIRGINIA\t|\n');
fprintf(fid,'\t|\t\t\t\t\t(SEPTEMBER 1st, 2004)\t\t\t\t\t|\n');
fprintf(fid,'\t|\t\t\t\t\t\t\t\t\t\t\t\t\t|\n');
fprintf(fid,'\t|\t\t\t(Developed under NASA Grant NAG3-2524 with Dr. S.M. Arnold\t|\n');
fprintf(fid,'\t|\t\t\t as the Technical Monitor)\t\t\t\t\t\t|\n');
fprintf(fid,'\t|\t\t\t\t\t\t\t\t\t\t\t\t\t|\n');
fprintf(fid,'\t| FUNCTIONALITY: Determination of macroscopic response and micro-level\t\t|\n');
fprintf(fid,'\t|\t\t\tdisplacement, stress & strain fields of periodic multiphase\t|\n');
fprintf(fid,'\t|\t\t\tmaterials with phase alignment along a common direction.\t|\n');
fprintf(fid,'\t|\t\t\t\t\t\t\t\t\t\t\t\t\t|\n');
fprintf(fid,'\t| REFERENCES:Yogesh Bansal and M-J. Pindera: Testing the Predictive Capability|\n');
fprintf(fid,'\t|\t\t of High-Fidelity Generalized Method of Cells Using an Efficient\t|\n');
fprintf(fid,'\t|\t\t Reformulation, NASA Contractor Report 213043, April 2004.\t\t|\n');
fprintf(fid,'\t|\t\t\t\t\t\t\t\t\t\t\t\t\t|\n');
fprintf(fid,'\t|\t\tYogesh Bansal and M-J. Pindera: Efficient Reformulation of the\t|\n');
fprintf(fid,'\t|\t\tThermoelastic Higher-Order Theory for FGMs, J. Thermal Stresses,\t|\n');
fprintf(fid,'\t|\t\tVol. 26(11/12), 2003, pp. 1055-1092.\t\t\t\t\t|\n');
fprintf(fid,'\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n');
% %__________________________________________________________________________________________________________________________________
fclose(fid);


fid_1=fopen(file_input_title,'r');    % OPENING AND READING DATA FROM THE FORMATTED INPUT FILE USED IN ORIGINAL HFGMC2D
% %__________________________________________________________________________________________________________________________________
% Reading number of materials specified (NMAT)
signal=0;
tline=fgetl(fid_1);
while signal ~= 1
    if length(findstr('Nmat',tline)) > 0
        data = sscanf(tline,'%f',[1,1]);
        NMAT=data;
        signal =1;
    else
        tline=fgetl(fid_1);
        if feof(fid_1) == 1
            error('Error in reading total number of materials, could not find the string "Nmat" in the input datafile')
        end
    end
end
clear data;

% Reading material properties for all the materials specified.
for k=1:NMAT

    % Material Number and Name (MAT.name)
    signal=0;
    tline=fgetl(fid_1);
    while signal ~= 1
        if length(findstr('''',tline)) > 1
            MAT(k).name = sscanf(tline,'%s');
            signal =1;
        else
            tline=fgetl(fid_1);
            if feof(fid_1) == 1
                error('Error in reading material number and name, could not find quotation marks "''" in the input datafile')
            end
        end
    end

    % Material Model (MM)
    signal=0;
    tline=fgetl(fid_1);
    while signal ~= 1
        if length(findstr('plasticity',tline)) > 0
            MM(k) = sscanf(tline,'%f',[1,1]);
            signal =1;
        else
            tline=fgetl(fid_1);
            if feof(fid_1) == 1
                error('Error in reading material model, could not find the string "plasticity" in the input datafile')
            end
        end
    end
    clear data

    % Direction of Anisotropy (ID)
    signal=0;
    tline=fgetl(fid_1);
    while signal ~= 1
        if length(findstr('Id',tline)) > 0
            ID(k) = sscanf(tline,'%f',[1,1]);
            signal =1;
        else
            tline=fgetl(fid_1);
            if feof(fid_1) == 1
                error('Error in reading direction of anisotropy for the current material, could not find the string "Id" in the input datafile')
            end
        end
    end
    clear data

    tline=fgetl(fid_1);                         % Reading the line containing "Young's Modulus in the Axial Direction", "Young's Modulus in the Transverse Direction", and "Axial Shear Modulus".
    data = sscanf(tline,'%f',[3,1]);
    [EA(k)]=data(1);                            % Young's Modulus in the Axial Direction
    [ET(k)]=data(2);                            % Young's Modulus in the Transverse Direction
    [GA(k)]=data(3);                            % Axial Shear Modulus

    tline=fgetl(fid_1);                           % Reading the line containing "Longitudinal Poisson's Ratio", and "Transverse Poisson's Ratio".
    data = sscanf(tline,'%f',[2,1]);
    [NUA(k)]=data(1);                           % Longitudinal Poisson's Ratio
    [NUT(k)]=data(2);                           % Transverse Poisson's Ratio

    tline=fgetl(fid_1);                         % Reading the line containing "CTE in the Axial Direction", and "CTE in the Transverse Direction".
    data = sscanf(tline,'%f',[2,1]);
    [ALFA(k)]=data(1);                          % Coefficient of Thermal Expansion in the Axial Direction
    [ALFT(k)]=data(2);                          % Coefficient of Thermal Expansion in the Transverse Direction

    tline=fgetl(fid_1);                         % Reading the line containing "Yield Strength", and "Hardening Parameter".
    data = sscanf(tline,'%f',[2,1]);
    [Y(k)]=data(1);                             % Yield Strength
    [HP(k)]=data(2);                            % Hardening Parameter
end

% Reading Number of subcells (m) in X_2 direction (Aboudi's X_3 direction), and (n) in X_3 direction (Aboudi's X_2 direction)
signal=0;
tline=fgetl(fid_1);
while signal ~= 1
    if length(findstr('Nbeta',tline)) > 0
        data = sscanf(tline,'%f',[2,1]);
        n=data(1);
        m=data(2);
        signal =1;
    else
        tline=fgetl(fid_1);
        if feof(fid_1) == 1
            error('Error in reading the number of subcells in X_2 and X_3 direction, could not find the string "Nbeta" in the input datafile')
        end
    end
end

% Reading the type of Geometry, Ind, (user supplied or vf specified), user supplied by default
signal=0;
tline=fgetl(fid_1);
while signal ~= 1
    if length(findstr('Ind',tline)) > 0
        Ind = sscanf(tline,'%f',[1,1]);
        signal =1;
    else
        tline=fgetl(fid_1);
        if feof(fid_1) == 1
            error('Error in reading the type of geometry, could not find the string "Ind" in the input datafile')
        end
    end
end


% Reading the height (L) of each subcell in X_3 direction (Aboudi's X_2 direction)
for i=1:n
    [num]=fscanf(fid_1,'%d',[1,1]);
    [L(num)]=fscanf(fid_1,'%f',[1,1]);
end

% Reading the length (H) of each subcell in X_2 direction (Aboudi's X_3 direction)
[H]=fscanf(fid_1,'%f',[m,1]);

% Reading the material assigned to each subcell
for i=1:n
    [index]=fscanf(fid_1,'%d',[1,1]);
    [k]=fscanf(fid_1,'%f',[m,1]);
    K(index,:)=k';
end
fclose(fid_1);
% %__________________________________________________________________________________________________________________________________
% DATA READ FROM THE INPUT FILE


fid_2=fopen(file_output_title,'a');   % CREATING THE OUTPUT FILE CONTAINING RELEVANT INPUT DATA AND RESULTING EFFECTIVE PROPERTIES.
% %__________________________________________________________________________________________________________________________________
fprintf(fid_2,'TOTAL NUMBER OF MATERIALS SPECIFIED\t\t%d\n', NMAT);

for i=1:NMAT
    fprintf(fid_2,'\n%s\n', MAT(i).name);
    fprintf(fid_2,'DIRECTION OF ANISOTROPY\t\t\t%d\n', ID(i));
    fprintf(fid_2,'Axial Young''s Modulus (EA)\t\t%d\n', EA(i));
    fprintf(fid_2,'Transverse Young''s Modulus (ET)\t%d\n', ET(i));
    fprintf(fid_2,'Axial Shear Modulus (GA)\t\t%d\n', GA(i));
    fprintf(fid_2,'Axial Poisson''s ratio (NUA)\t\t%d\n', NUA(i));
    fprintf(fid_2,'Transverse Poisson''s ratio (NUT)\t%d\n', NUT(i));
    fprintf(fid_2,'Axial CTE\t\t\t\t\t%d\n', ALFA(i));
    fprintf(fid_2,'Transverse CTE\t\t\t\t%d\n', ALFT(i));
    fprintf(fid_2,'Yield Stress\t\t\t\t%d\n', Y(i));
    fprintf(fid_2,'Hardening Parameter\t\t\t%d\n', HP(i));
end

fprintf(fid_2,'\n\nTOTAL NUMBER OF SUBCELLS IN X_2 (--->) DIRECTION\t%d\n\n', m);
fprintf(fid_2,'                                  ^ \n');
fprintf(fid_2,'TOTAL NUMBER OF SUBCELLS IN X_3 ( | ) DIRECTION\t\t%d\n', n);

fprintf(fid_2,'\n\nLENGTH OF EACH SUBCELL IN X_2 (--->) DIRECTION\n');
for num=1:m
    fprintf(fid_2,'%G\n',H(num));
end
fprintf(fid_2,'\n\nHEIGHT OF EACH SUBCELL IN X_3 DIRECTION\n');
for num=1:n
    fprintf(fid_2,'%-5d\t%G\n',num,L(num));
end

fprintf(fid_2,'\n\nMICROSTRUCTURE\n\n');
for i=n:-1:1
    fprintf(fid_2,'%-5d',i);
    for j=1:m
        fprintf(fid_2,'%2d',K(i,j));
    end
    fprintf(fid_2,'\n');
end
fclose(fid_2);
% %__________________________________________________________________________________________________________________________________

% % CREATING x2.dat WHICH CONTAINS HORIZONTAL SUBCELL COORDINATES
% %__________________________________________________________________________
fid_3=fopen('x2.dat','w');
x(1)=0;
fprintf(fid_3,'%G\n',x(1));
for i=1:m
    x(i+1)=x(i)+H(i);
    fprintf(fid_3,'%G\n',x(i+1));
end
fclose(fid_3);
% %__________________________________________________________________________

% % CREATING x3.dat WHICH CONTAINS VERTICAL SUBCELL COORDINATES
% %__________________________________________________________________________
fid_4=fopen('x3.dat','w');
y(1)=0;
fprintf(fid_4,'%G\n',y(1));
for i=1:n
    y(i+1)=y(i)+L(i);
    fprintf(fid_4,'%G\n',y(i+1));
end
fclose(fid_4);
% %__________________________________________________________________________

% Reassignment of material type, length, and height to all the subcells
clear k;
i=1:n;
temp_num=(i-1)*m;
for j=1:m
    k(temp_num+j)=K(i,j);
    h(temp_num+j)=H(j);
    l(temp_num+j)=L;
end
% clear K H L;

% % Plotting the microstructure for the given problem
% %__________________________________________________________________________
% load x2.dat
% load x3.dat
% ms=figure;
% % set(ms,'PaperPosition',[0.0 0.0 8.5 11.0],'units','inches');
% % subplot(2,1,1);
% hold on;
% % axis([x2(1) x2(m+1) x3(1) x3(n+1)]);
% % Plot
% for j=1:n
%     for i=1:m
%         % Define rectangle using corner points stored in xnew
%         rectx = [x2(i) x2(i+1) x2(i+1) x2(i)];
%         recty = [x3(j) x3(j) x3(j+1) x3(j+1)];
%         if k((j-1)*m+i) == 1
%             fill(rectx,recty,[0.9 0.9 0.9]);
%         elseif k((j-1)*m+i) == 2
%             fill(rectx,recty,[0.5 0.5 0.5]);
%         elseif k((j-1)*m+i) == 3
%             fill(rectx,recty,[0.3 0.3 0.3]);
%         end
%     end
% end
% for i=1:crack_num
%     if crack_dir(i) == 2
%         line([x2(crack_subcell_beta1(i)) x2(crack_subcell_beta2(i)+1)],[x3(crack_subcell_gama(i)+1) x3(crack_subcell_gama(i)+1)],'color', [0.99 0.99 0.99], 'linewidth', 2);
%     end
%     if crack_dir(i) == 3
%         line([x2(crack_subcell_beta(i)+1) x2(crack_subcell_beta(i)+1)],[x3(crack_subcell_gama1(i)) x3(crack_subcell_gama2(i)+1)], 'color', [0.99 0.99 0.99], 'linewidth', 2);
%     end
% end
% axis([0.5 1 0.4 0.8]);
% % axis image;
% axis off;
% % eval(['title(''\bf{Mesh Discretization and Microstructure} (' num2str(m) ' x ' num2str(n) ')'', ''fontsize'',14);']);
% % print(strcat(string_title,'_ms.tif'), '-dtiff', '-r300');
% % close;
% % %__________________________________________________________________________
% break;


% Assigning Engineering Constants in X_1, X_2, and X_3 directions based on the Direction of Anisotropy
for j=1:NMAT
    if ID(j)==1
        E11_mat(j)=EA(j);
        E22_mat(j)=ET(j);
        E33_mat(j)=ET(j);
        NU12_mat(j)=NUA(j);
        NU21_mat(j)=E22_mat(j)*NU12_mat(j)/E11_mat(j);
        NU13_mat(j)=NUA(j);
        NU31_mat(j)=E33_mat(j)*NU13_mat(j)/E11_mat(j);
        NU23_mat(j)=NUT(j);
        NU32_mat(j)=E33_mat(j)*NU23_mat(j)/E22_mat(j);
        G12_mat(j)=GA(j);
        G13_mat(j)=GA(j);
        G23_mat(j)=E22_mat(j)/(2*(1+NU23_mat(j)));
        ALF11_mat(j)=ALFA(j);
        ALF22_mat(j)=ALFT(j);
        ALF33_mat(j)=ALFT(j);
    elseif ID(j) == 3
        E22_mat(j)=EA(j);
        E11_mat(j)=ET(j);
        E33_mat(j)=ET(j);
        NU21_mat(j)=NUA(j);
        NU12_mat(j)=E11_mat(j)*NU21_mat(j)/E22_mat(j);
        NU23_mat(j)=NUA(j);
        NU32_mat(j)=E33_mat(j)*NU23_mat(j)/E22_mat(j);
        NU13_mat(j)=NUT(j);
        NU31_mat(j)=E33_mat(j)*NU13_mat(j)/E11_mat(j);
        G12_mat(j)=GA(j);
        G23_mat(j)=GA(j);
        G13_mat(j)=E11_mat(j)/(2*(1+NU13_mat(j)));
        ALF22_mat(j)=ALFA(j);
        ALF11_mat(j)=ALFT(j);
        ALF33_mat(j)=ALFT(j);
    elseif ID(j) == 2
        E33_mat(j)=EA(j);
        E11_mat(j)=ET(j);
        E22_mat(j)=ET(j);
        NU13_mat(j)=NUA(j);
        NU31_mat(j)=E33_mat(j)*NU13_mat(j)/E11_mat(j);
        NU23_mat(j)=NUA(j);
        NU32_mat(j)=E33_mat(j)*NU23_mat(j)/E22_mat(j);
        NU12_mat(j)=NUT(j);
        NU21_mat(j)=E22_mat(j)*NU12_mat(j)/E11_mat(j);
        G13_mat(j)=GA(j);
        G23_mat(j)=GA(j);
        G12_mat(j)=E11_mat(j)/(2*(1+NU12_mat(j)));
        ALF33_mat(j)=ALFA(j);
        ALF11_mat(j)=ALFT(j);
        ALF22_mat(j)=ALFT(j);
    end
end


% Reassignment of material properties to all the real subcells
E11=E11_mat(k);
E22=E22_mat(k);
E33=E33_mat(k);
v12=NU12_mat(k);
v21=NU21_mat(k);
v13=NU13_mat(k);
v31=NU31_mat(k);
v23=NU23_mat(k);
v32=NU32_mat(k);
G12=G12_mat(k);
G13=G13_mat(k);
G23=G23_mat(k);
alpha11=ALF11_mat(k);
alpha22=ALF22_mat(k);
alpha33=ALF33_mat(k);
Hp=HP(k);
Y_Sigma=Y(k);

delta=(1-v12.*v21-v23.*v32-v13.*v31-2*v21.*v32.*v13)./(E11.*E22.*E33);

C11=(1-v32.*v23)./(E22.*E33.*delta);
C22=(1-v31.*v13)./(E11.*E33.*delta);
C33=(1-v12.*v21)./(E11.*E22.*delta);

C12=(v21+v31.*v23)./(E33.*E22.*delta);
C13=(v31+v21.*v32)./(E33.*E22.*delta);
C23=(v32+v31.*v12)./(E33.*E11.*delta);

C44=G23;
C55=G13;
C66=G12;

C22bar=C22+C44.*h.^2./l.^2;
C23bar=C33.*h.^2./l.^2+C44;
C32bar=C22.*l.^2./h.^2+C44;
C33bar=C33+C44.*l.^2./h.^2;
C55bar=C55+C66.*l.^2./h.^2;
C66bar=C66+C55.*h.^2./l.^2;

[crack_row, crack_col, crack_faces_ver, crack_faces_hor, ul, ur, ub, ut] = crack(m, n, crack_num, crack_dir, crack_subcell_beta, crack_subcell_gama1,....
    crack_subcell_gama2, crack_subcell_gama, crack_subcell_beta1, crack_subcell_beta2);


% "_i" represents the inplane quantities
% "_o" represents the out of plane quantities.

% Evaluating the inplane global stiffness matrix (K_i) and out-of-plane global stiffness matrix (K_o). They can be assembled seperately as the inplane and out of plane quantities are uncoupled
[K_i, K_o] = stiffness_cracks(m, n, h, l, C22, C33, C12, C13, C23, C44, C55, C66, C22bar, C23bar, C32bar, C33bar, C55bar, C66bar, crack_row, crack_col, crack_faces_ver, crack_faces_hor, ul, ur, ub, ut);

K_i_r=K_i;          % K_i and K_o are preserved and reassigned to K_i_r and K_o_r which will be reduced as the corner subcell is fixed.
K_o_r=K_o;

NS=m*n;             % Total number of subcells

% If fixing one corner subcell
p=ur(NS)-n;
q=ut(NS)-m;
temp_i=[ul(1) ul(m*(n-1)+1)-(n-1) p+ul(1) p+ul(m*(n-1)+1)-(n-1)  2*p+ub(1) 2*p+ub(m)-(m-1)  2*p+q+ub(1) 2*p+q+ub(m)-(m-1)];
temp_o=[ul(1) ul(m*(n-1)+1)-(n-1) p+ub(1) p+ub(m)-(m-1)];

K_i_r(temp_i,:)=[];
K_i_r(:,temp_i)=[];
K_o_r(temp_o,:)=[];
K_o_r(:,temp_o)=[];


% spparms('spumoni',2);
% p_i=colamd(K_i_r);
% [L_i,U_i,P_i] = lu(K_i_r(:,p_i)); % lu factorization of the inplane equations
% spparms('spumoni',0);
% p_o=colamd(K_o_r);
% [L_o,U_o,P_o] = lu(K_o_r(:,p_o)); % lu factorization of the out-of-plane equations

[L_i,U_i,P_i,Q_i] = lu(K_i_r); % lu factorization of the inplane equations
[L_o,U_o,P_o,Q_o] = lu(K_o_r); % lu factorization of the out-of-plane equations

dim=2*(ur(NS)-n)+2*(ut(NS)-m);        % Dimension of the global stiffness matrix

% Initializing the inplane and out-of plane surface-averaged displacements and tractions for the
Avu_ini=zeros(dim,1);
Avu_ino=zeros(dim/2,1);

Sigma_i=sparse(dim,1);
Sigma_o=sparse(dim/2,1);


%%_________________________________________________________________________
Ave11=0.01;
Ave22=0;
Ave33=0;
Ave23=0;
Ave13=0;
Ave12=0;

clear R1 R2 R3 R4 R5 R6
R1=(Ave11*C12+Ave22*C22+Ave33*C23);
R2=2*Ave23*C44;
R3=2*Ave23*C44;
R4=(Ave11*C13+Ave22*C23+Ave33*C33);
R5=2*Ave13*C55;
R6=2*Ave12*C66;

[A_i, A_o] = loading(m, n, R1, R2, R3, R4, R5, R6, crack_row, crack_col, crack_faces_ver, crack_faces_hor, ul, ur, ub, ut);

Sigma_i_r=Sigma_i+A_i;
Sigma_i_r(temp_i,:)=[];
Sigma_o_r=Sigma_o+A_o;
Sigma_o_r(temp_o,:)=[];

% [Avu_i]=inverse_next(L_i, U_i, P_i, p_i, Sigma_i_r,Avu_ini,temp_i);
% [Avu_o]=inverse_next(L_o, U_o, P_o, p_o, Sigma_o_r,Avu_ino,temp_o);
[Avu_i]=inverse_next(L_i, U_i, P_i, Q_i, Sigma_i_r,Avu_ini,temp_i);
[Avu_o]=inverse_next(L_o, U_o, P_o, Q_o, Sigma_o_r,Avu_ino,temp_o);

[W2010, W3010, W2001, W3001, W1010, W1001] = microvariables(m, n, h, l, Avu_i, Avu_o, ul, ur, ub, ut);

e11(1:NS)=Ave11;
e22=Ave22+W2010;
e33=Ave33+W3001;
e23=Ave23+(W2001+W3010)/2;
e13=Ave13+W1001/2;
e12=Ave12+W1010/2;

A11=e11/Ave11;
A21=e22/Ave11;
A31=e33/Ave11;
A41=2*e23/Ave11;
A51=2*e13/Ave11;
A61=2*e12/Ave11;
%%_________________________________________________________________________


%%_________________________________________________________________________
Ave11=0;
Ave22=0.01;
Ave33=0;
Ave23=0;
Ave13=0;
Ave12=0;

clear R1 R2 R3 R4 R5 R6
R1=(Ave11*C12+Ave22*C22+Ave33*C23);
R2=2*Ave23*C44;
R3=2*Ave23*C44;
R4=(Ave11*C13+Ave22*C23+Ave33*C33);
R5=2*Ave13*C55;
R6=2*Ave12*C66;

[A_i, A_o] = loading(m, n, R1, R2, R3, R4, R5, R6, crack_row, crack_col, crack_faces_ver, crack_faces_hor, ul, ur, ub, ut);

Sigma_i_r=Sigma_i+A_i;
Sigma_i_r(temp_i,:)=[];
Sigma_o_r=Sigma_o+A_o;
Sigma_o_r(temp_o,:)=[];

% [Avu_i]=inverse_next(L_i, U_i, P_i, Q_i, p_i, Sigma_i_r,Avu_ini,temp_i);
% [Avu_o]=inverse_next(L_o, U_o, P_o, Q_o, p_o, Sigma_o_r,Avu_ino,temp_o);
[Avu_i]=inverse_next(L_i, U_i, P_i, Q_i, Sigma_i_r,Avu_ini,temp_i);
[Avu_o]=inverse_next(L_o, U_o, P_o, Q_o, Sigma_o_r,Avu_ino,temp_o);

[W2010, W3010, W2001, W3001, W1010, W1001] = microvariables(m, n, h, l, Avu_i, Avu_o, ul, ur, ub, ut);

e11(1:NS)=Ave11;
e22=Ave22+W2010;
e33=Ave33+W3001;
e23=Ave23+(W2001+W3010)/2;
e13=Ave13+W1001/2;
e12=Ave12+W1010/2;

A12=e11/Ave22;
A22=e22/Ave22;
A32=e33/Ave22;
A42=2*e23/Ave22;
A52=2*e13/Ave22;
A62=2*e12/Ave22;
%%_________________________________________________________________________


%%_________________________________________________________________________
Ave11=0;
Ave22=0;
Ave33=0.01;
Ave23=0;
Ave13=0;
Ave12=0;

clear R1 R2 R3 R4 R5 R6
R1=(Ave11*C12+Ave22*C22+Ave33*C23);
R2=2*Ave23*C44;
R3=2*Ave23*C44;
R4=(Ave11*C13+Ave22*C23+Ave33*C33);
R5=2*Ave13*C55;
R6=2*Ave12*C66;

[A_i, A_o] = loading(m, n, R1, R2, R3, R4, R5, R6, crack_row, crack_col, crack_faces_ver, crack_faces_hor, ul, ur, ub, ut);

Sigma_i_r=Sigma_i+A_i;
Sigma_i_r(temp_i,:)=[];
Sigma_o_r=Sigma_o+A_o;
Sigma_o_r(temp_o,:)=[];

% [Avu_i]=inverse_next(L_i, U_i, P_i, Q_i, p_i, Sigma_i_r,Avu_ini,temp_i);
% [Avu_o]=inverse_next(L_o, U_o, P_o, Q_o, p_o, Sigma_o_r,Avu_ino,temp_o);
[Avu_i]=inverse_next(L_i, U_i, P_i, Q_i, Sigma_i_r,Avu_ini,temp_i);
[Avu_o]=inverse_next(L_o, U_o, P_o, Q_o, Sigma_o_r,Avu_ino,temp_o);

[W2010, W3010, W2001, W3001, W1010, W1001] = microvariables(m, n, h, l, Avu_i, Avu_o, ul, ur, ub, ut);

e11(1:NS)=Ave11;
e22=Ave22+W2010;
e33=Ave33+W3001;
e23=Ave23+(W2001+W3010)/2;
e13=Ave13+W1001/2;
e12=Ave12+W1010/2;

A13=e11/Ave33;
A23=e22/Ave33;
A33=e33/Ave33;
A43=2*e23/Ave33;
A53=2*e13/Ave33;
A63=2*e12/Ave33;
%%_________________________________________________________________________


%%_________________________________________________________________________
Ave11=0;
Ave22=0;
Ave33=0;
Ave23=0.01;
Ave13=0;
Ave12=0;

clear R1 R2 R3 R4 R5 R6
R1=(Ave11*C12+Ave22*C22+Ave33*C23);
R2=2*Ave23*C44;
R3=2*Ave23*C44;
R4=(Ave11*C13+Ave22*C23+Ave33*C33);
R5=2*Ave13*C55;
R6=2*Ave12*C66;

[A_i, A_o] = loading(m, n, R1, R2, R3, R4, R5, R6, crack_row, crack_col, crack_faces_ver, crack_faces_hor, ul, ur, ub, ut);

Sigma_i_r=Sigma_i+A_i;
Sigma_i_r(temp_i,:)=[];
Sigma_o_r=Sigma_o+A_o;
Sigma_o_r(temp_o,:)=[];

% [Avu_i]=inverse_next(L_i, U_i, P_i, Q_i, p_i, Sigma_i_r,Avu_ini,temp_i);
% [Avu_o]=inverse_next(L_o, U_o, P_o, Q_o, p_o, Sigma_o_r,Avu_ino,temp_o);
[Avu_i]=inverse_next(L_i, U_i, P_i, Q_i, Sigma_i_r,Avu_ini,temp_i);
[Avu_o]=inverse_next(L_o, U_o, P_o, Q_o, Sigma_o_r,Avu_ino,temp_o);

[W2010, W3010, W2001, W3001, W1010, W1001] = microvariables(m, n, h, l, Avu_i, Avu_o, ul, ur, ub, ut);

e11(1:NS)=Ave11;
e22=Ave22+W2010;
e33=Ave33+W3001;
e23=Ave23+(W2001+W3010)/2;
e13=Ave13+W1001/2;
e12=Ave12+W1010/2;

A14=e11/(2*Ave23);
A24=e22/(2*Ave23);
A34=e33/(2*Ave23);
A44=e23/Ave23;
A54=e13/Ave23;
A64=e12/Ave23;
%%_________________________________________________________________________


%%_________________________________________________________________________
Ave11=0;
Ave22=0;
Ave33=0;
Ave23=0;
Ave13=0.01;
Ave12=0;

clear R1 R2 R3 R4 R5 R6
R1=(Ave11*C12+Ave22*C22+Ave33*C23);
R2=2*Ave23*C44;
R3=2*Ave23*C44;
R4=(Ave11*C13+Ave22*C23+Ave33*C33);
R5=2*Ave13*C55;
R6=2*Ave12*C66;

[A_i, A_o] = loading(m, n, R1, R2, R3, R4, R5, R6, crack_row, crack_col, crack_faces_ver, crack_faces_hor, ul, ur, ub, ut);

Sigma_i_r=Sigma_i+A_i;
Sigma_i_r(temp_i,:)=[];
Sigma_o_r=Sigma_o+A_o;
Sigma_o_r(temp_o,:)=[];

% [Avu_i]=inverse_next(L_i, U_i, P_i, Q_i, p_i, Sigma_i_r,Avu_ini,temp_i);
% [Avu_o]=inverse_next(L_o, U_o, P_o, Q_o, p_o, Sigma_o_r,Avu_ino,temp_o);
[Avu_i]=inverse_next(L_i, U_i, P_i, Q_i, Sigma_i_r,Avu_ini,temp_i);
[Avu_o]=inverse_next(L_o, U_o, P_o, Q_o, Sigma_o_r,Avu_ino,temp_o);

[W2010, W3010, W2001, W3001, W1010, W1001] = microvariables(m, n, h, l, Avu_i, Avu_o, ul, ur, ub, ut);

e11(1:NS)=Ave11;
e22=Ave22+W2010;
e33=Ave33+W3001;
e23=Ave23+(W2001+W3010)/2;
e13=Ave13+W1001/2;
e12=Ave12+W1010/2;

A15=e11/(2*Ave13);
A25=e22/(2*Ave13);
A35=e33/(2*Ave13);
A45=e23/Ave13;
A55=e13/Ave13;
A65=e12/Ave13;
%%_________________________________________________________________________


%%_________________________________________________________________________
Ave11=0;
Ave22=0;
Ave33=0;
Ave23=0;
Ave13=0;
Ave12=0.01;

clear R1 R2 R3 R4 R5 R6
R1=(Ave11*C12+Ave22*C22+Ave33*C23);
R2=2*Ave23*C44;
R3=2*Ave23*C44;
R4=(Ave11*C13+Ave22*C23+Ave33*C33);
R5=2*Ave13*C55;
R6=2*Ave12*C66;

[A_i, A_o] = loading(m, n, R1, R2, R3, R4, R5, R6, crack_row, crack_col, crack_faces_ver, crack_faces_hor, ul, ur, ub, ut);

Sigma_i_r=Sigma_i+A_i;
Sigma_i_r(temp_i,:)=[];
Sigma_o_r=Sigma_o+A_o;
Sigma_o_r(temp_o,:)=[];

% [Avu_i]=inverse_next(L_i, U_i, P_i, Q_i, p_i, Sigma_i_r,Avu_ini,temp_i);
% [Avu_o]=inverse_next(L_o, U_o, P_o, Q_o, p_o, Sigma_o_r,Avu_ino,temp_o);
[Avu_i]=inverse_next(L_i, U_i, P_i, Q_i, Sigma_i_r,Avu_ini,temp_i);
[Avu_o]=inverse_next(L_o, U_o, P_o, Q_o, Sigma_o_r,Avu_ino,temp_o);

[W2010, W3010, W2001, W3001, W1010, W1001] = microvariables(m, n, h, l, Avu_i, Avu_o, ul, ur, ub, ut);

e11(1:NS)=Ave11;
e22=Ave22+W2010;
e33=Ave33+W3001;
e23=Ave23+(W2001+W3010)/2;
e13=Ave13+W1001/2;
e12=Ave12+W1010/2;

A16=e11/(2*Ave12);
A26=e22/(2*Ave12);
A36=e33/(2*Ave12);
A46=e23/Ave12;
A56=e13/Ave12;
A66=e12/Ave12;
%%_________________________________________________________________________



% clear A;
% for i=1:NS
%     A(:,:,i)=[A11(i) A12(i) A13(i) A14(i) A15(i) A16(i); A21(i) A22(i) A23(i) A24(i) A25(i) A26(i); A31(i) A32(i) A33(i) A34(i) A35(i) A36(i); A41(i) A42(i) A43(i) A44(i) A45(i) A46(i); A51(i) A52(i) A53(i) A54(i) A55(i) A56(i); A61(i) A62(i) A63(i) A64(i) A65(i) A66(i);];
%     Cnext(:,:,i)=[C11(i) C12(i) C13(i) 0 0 0; C12(i) C22(i) C23(i) 0 0 0; C13(i) C23(i) C33(i) 0 0 0; 0 0 0 C44(i) 0 0; 0 0 0 0 C55(i) 0; 0 0 0 0 0 C66(i)];
%     Gamma(:,:,i)=[C11(i)*alpha11(i)+C12(i)*alpha22(i)+C13(i)*alpha33(i);C12(i)*alpha11(i)+C22(i)*alpha22(i)+C23(i)*alpha33(i);C13(i)*alpha11(i)+C23(i)*alpha22(i)+C33(i)*alpha33(i); 0; 0; 0];
% end
%
% Cstar=zeros(6);
% for i=1:NS
%     Cstar=Cstar+h(i)*l(i)*Cnext(:,:,i)*A(:,:,i);
% end
% Gammastar=zeros(6,1);
% for i=1:NS
%     Gammastar=Gammastar+h(i)*l(i)*A(:,:,i)'*Gamma(:,:,i);
% end

H=0;
for i=1:m
    H=H+h(i);
end
L=0;
for i=1:m:(n-1)*m+1
    L=L+l(i);
end

effective_stiffness;
Cstar=Cstar/(H*L);
Gammastar=Gammastar/(H*L);
alphastar=Cstar\Gammastar;

Sstar=inv(Cstar);
E11=1/Sstar(1,1);
E22=1/Sstar(2,2);
E33=1/Sstar(3,3);
G23=1/Sstar(4,4);
G13=1/Sstar(5,5);
G12=1/Sstar(6,6);
v12=-Sstar(2,1)*E11;
v13=-Sstar(3,1)*E11;
v23=-Sstar(3,2)*E22;

fid=fopen(file_output_title,'a');
fprintf(fid,'\nEFFECTIVE STIFFNESS\n');
fprintf(fid,'%E\t%E\t%E\t%E\t%E\t%E\n',Cstar);

fprintf(fid,'\n\nEFFECTIVE MODULI\n');
fprintf(fid,'E11 = %14.4E\n',E11);
fprintf(fid,'E22 = %14.4E\n',E22);
fprintf(fid,'E33 = %14.4E\n',E33);
fprintf(fid,'G23 = %14.4E\n',G23);
fprintf(fid,'G13 = %14.4E\n',G13);
fprintf(fid,'G12 = %14.4E\n',G12);
fprintf(fid,'NU23 = %14.4E\n',v23);
fprintf(fid,'NU13 = %14.4E\n',v13);
fprintf(fid,'NU12 = %14.4E\n',v12);
fclose(fid);


%%_________________________________________________________________________
% Integrating the plastic strain contribution numerically within a subcell
% using t x t Gauss-Quadrature (Gauss-Legendre) method.

% % For demonstration, if t = 3, then the representation used for the
% corresponding points, their weights and coordinates are
% %   3,1    3,2   3,3
% %   2,2    2,2   2,3
% %   1,1    1,2   1,3
%
% %   13    23    33
% %   12    22    32
% %   11    21    31
%
% %   7    8    9
% %   4    5    6
% %   1    2    3
%
% % Coordinates of those points
% array_x2=[-0.774596669241483 0 0.774596669241483];
% array_x3=[-0.774596669241483 0 0.774596669241483];
%
% % weights associated with each point
% w1=0.55555555556;
% w2=0.88888888889;
% w3=0.55555555556;
% weight=[w1*w1 w1*w2 w1*w3;....
%     w2*w1 w2*w2 w2*w3;....
%     w3*w1 w3*w2 w3*w3];

% % t=input('Read the number of integration points to be used in the Legendre-Gauss Quadrature\n');
[array_x,array_w]=lgwt(t,-1,1);
for i=1:t
    for j = 1:t
        weight(i,j) = array_w(i)*array_w(j);
    end
end
array_x=flipdim(array_x,1);
array_x=array_x';
array_w=array_w';

% % Number of terms in the Legendre polynomial expansion
% % para_M=input('Read the number of terms to be used in Legendre polynomial expansion\n');
para_N=para_M;

% % Legendre polynomials LP_2 as a function of x2 and LP_3 as a function of
% % x3, local subcell coordinates.
for para_m=0:para_M
    LP2=legendre(para_m,array_x);
    LP_2(para_m+1,:)=LP2(1,:);
end
for para_n=0:para_N
    LP3=legendre(para_n,array_x);
    LP_3(para_n+1,:)=LP3(1,:);
end

load_current=load_initial;
load_array=load_initial;
av_stress_macro=zeros(6,1);
av_strain_macro=zeros(6,1);
av_stress_macro_bar=zeros(6,1);
av_strain_macro_bar=zeros(6,1);
incr_register_length=length(incr_register);

fid=fopen(file_output_title,'a');
fprintf(fid,'\nOFF-AXIS ANGLE BETWEEN GLOBAL X AND PRINCIPLE 1 (FIBER) DIRECTION\t\t%d DEGREES\n', theta);
fprintf(fid,'\n\n\nELASTIC-PLASTIC ANALYSIS\n');
fprintf(fid,'____________________________________________________________________\n');
fprintf(fid,'\nLoading option\t%d\n',LOP);
if LOP == 1
    fprintf(fid,'EPS11 Applied,\tSIG22=SIG33=SIG23=SIG13=SIG12=0,\tDELTA_T=0\n');
elseif LOP == 2
    fprintf(fid,'EPS22 Applied,\tSIG11=SIG33=SIG23=SIG13=SIG12=0,\tDELTA_T=0\n');
elseif LOP == 3
    fprintf(fid,'EPS33 Applied,\tSIG11=SIG22=SIG23=SIG13=SIG12=0,\tDELTA_T=0\n');
elseif LOP == 4
    fprintf(fid,'EPS23 Applied,\tSIG11=SIG22=SIG33=SIG13=SIG12=0,\tDELTA_T=0\n');
elseif LOP == 5
    fprintf(fid,'EPS13 Applied,\tSIG11=SIG22=SIG33=SIG23=SIG12=0,\tDELTA_T=0\n');
elseif LOP == 6
    fprintf(fid,'EPS12 Applied,\tSIG11=SIG22=SIG33=SIG23=SIG13=0,\tDELTA_T=0\n');
elseif LOP == 7
    fprintf(fid,'EPS11=0, EPS22=-EPS33 Applied,\tSIG23=SIG13=SIG12=0,\tDELTA_T=0\n');
elseif LOP == 8
    fprintf(fid,'DELTA_T Applied,\tSIG11=SIG22=SIG33=SIG23=SIG13=SIG12=0\n');
elseif LOP == 9
    fprintf(fid,'EPS22 Applied, EPS11=0,\tSIG33=SIG23=SIG13=SIG12=0,\tDELTA_T=0\n');
elseif LOP == 10
    fprintf(fid,'EPS11=0, EPS22=EPS33 Applied,\tSIG23=SIG13=SIG12=0,\tDELTA_T=0\n');
elseif LOP == 11
    fprintf(fid,'EPS22=EPS33 Applied,\tSIG11=SIG23=SIG13=SIG12=0,\tDELTA_T=0\n');
end


if LOP == 8
    fprintf(fid,'\nThermal Loading\n');
    fprintf(fid,'Initial Load\t%d degrees\n', load_initial);
    fprintf(fid,'Final Load\t\t%d degrees\n', load_max);
else
    fprintf(fid,'\nMechanical Loading\n');
    fprintf(fid,'\nInitial Load\t%d\n', load_initial);
    fprintf(fid,'Final Load\t\t%d\tper-cent\n', load_max*100);
end
fprintf(fid,'Total Number of Load Increments\t%d\n', incr_NUM);
fprintf(fid,'Maximum Number of Iterations allowed at each Load Increment\t\t%d\n', iter_limit);
fprintf(fid,'Data recorded after increment\t');
if incr_register_length == 0
    fprintf(fid,'NONE\n');
else
    for i=1:incr_register_length-1
        fprintf(fid,'%d, ', incr_register(i));
    end
    fprintf(fid,'%d\n', incr_register(incr_register_length));
end
fclose(fid);

%% Conversion from row vector to column vector
%%_________________________________________________________________________
mu=C44';
h=h';
l=l';
C11=C11';
C22=C22';
C33=C33';
C44=C44';
C55=C55';
C66=C66';
C22bar=C22bar';
C33bar=C33bar';
C55bar=C55bar';
C66bar=C66bar';
C12=C12';
C13=C13';
C23=C23';
alpha11=alpha11';
alpha22=alpha22';
alpha33=alpha33';
Y_Sigma=Y_Sigma';
Hp=Hp';
%%_________________________________________________________________________

%% Initializing total plastic strains and plastic strain increments to zero
for j=1:t
    for i=1:t
        eps_pl_pre(i,j).eff = sparse(NS,1);
        eps_pl_pre(i,j).dir11=sparse(NS,1);
        eps_pl_pre(i,j).dir22=sparse(NS,1);
        eps_pl_pre(i,j).dir33=sparse(NS,1);
        eps_pl_pre(i,j).dir23=sparse(NS,1);
        eps_pl_pre(i,j).dir13=sparse(NS,1);
        eps_pl_pre(i,j).dir12=sparse(NS,1);

        eps_pl(i,j).dir11=sparse(NS,1);
        eps_pl(i,j).dir22=sparse(NS,1);
        eps_pl(i,j).dir33=sparse(NS,1);
        eps_pl(i,j).dir23=sparse(NS,1);
        eps_pl(i,j).dir13=sparse(NS,1);
        eps_pl(i,j).dir12=sparse(NS,1);
    end
end


em=cos(pi/180*theta);
en=sin(pi/180*theta);
T1=[em^2 en^2 0 0 0 2*em*en; en^2 em^2 0 0 0 -2*em*en; 0 0 1 0 0 0; 0 0 0 em -en 0; 0 0 0 en em 0; -em*en em*en 0 0 0 em^2-en^2];
T2=[em^2 en^2 0 0 0 em*en; en^2 em^2 0 0 0 -em*en; 0 0 1 0 0 0; 0 0 0 em -en 0; 0 0 0 en em 0; -2*em*en 2*em*en 0 0 0 em^2-en^2];
T1_inv=inv(T1);
T2_inv=inv(T2);
Cstar_bar=T1_inv*Cstar*T2;

av_stress_macro_at_each_iter=av_stress_macro;
load_array_at_each_iter=load_array;
strain_thermal=sparse(6,1);
Rbar=sparse(NS,6);
D=sparse(NS,6);
field_at_incr=[];
delta_T=0;
toggle_matrix=sparse(3,3);
z=1;
for incr_num=1:incr_NUM
    incr_num
    load_current = load_current + (load_max-load_initial)/incr_NUM

    if LOP == 8
        delta_T=load_current;
    end

    for j=1:t
        for i=1:t
            deps_pl_pre(i,j).eff=ones(NS,1);
            deps_pl(i,j).dir11=sparse(NS,1);
            deps_pl(i,j).dir22=sparse(NS,1);
            deps_pl(i,j).dir33=sparse(NS,1);
            deps_pl(i,j).dir23=sparse(NS,1);
            deps_pl(i,j).dir13=sparse(NS,1);
            deps_pl(i,j).dir12=sparse(NS,1);
        end
    end

    for iter=1:iter_limit

        %                 sig_inel=sparse(6,1);
        %                 for i=1:NS
        %                     sig_inel=sig_inel + h(i)*l(i)*(Cnext(:,:,i)*D(:,i) - Rbar(:,i) - Gamma(:,:,i)*delta_T);
        %                 end
        effective_sig_inel;
        sig_inel = -sig_inel/(H*L);
        sig_inel_bar=T1_inv*sig_inel;

        strain_thermal_bar=T2_inv*strain_thermal;
        [av_macro_strains_bar, av_macro_stresses_bar] = traction_loading_plastic(LOP, load_current, Cstar_bar, sig_inel_bar, strain_thermal_bar);
        av_macro_strains=T2*av_macro_strains_bar;

        Ave11=av_macro_strains(1);
        Ave22=av_macro_strains(2);
        Ave33=av_macro_strains(3);
        Ave23=av_macro_strains(4)/2;
        Ave13=av_macro_strains(5)/2;
        Ave12=av_macro_strains(6)/2;

        R1=(Ave11*C12+Ave22*C22+Ave33*C23) - (alpha11.*C12+alpha22.*C22+alpha33.*C23)*delta_T;
        R2=2*Ave23*C44;
        R3=2*Ave23*C44;
        R4=(Ave11*C13+Ave22*C23+Ave33*C33)- (alpha11.*C13+alpha22.*C23+alpha33.*C33)*delta_T;
        R5=2*Ave13*C55;
        R6=2*Ave12*C66;

        [tau] = stress_coeff(para_M, para_N, mu, weight, LP_2, LP_3, eps_pl, t);
        G1_1=sparse(NS,1);
        G2_1=sparse(NS,1);
        G3_1=sparse(NS,1);
        for para_m=1:2:para_M
            G1_1 = G1_1 + sqrt(1+2*para_m)*tau(para_m+1,1).dir22;
            G2_1 = G2_1 + sqrt(1+2*para_m)*tau(para_m+1,1).dir23;
            G3_1 = G3_1 + sqrt(1+2*para_m)*tau(para_m+1,1).dir12;
        end
        G1_2=sparse(NS,1);
        G2_2=sparse(NS,1);
        G3_2=sparse(NS,1);
        for para_n=1:2:para_N
            G1_2 = G1_2 + sqrt(1+2*para_n)*tau(1,para_n+1).dir23;
            G2_2 = G2_2 + sqrt(1+2*para_n)*tau(1,para_n+1).dir33;
            G3_2 = G3_2 + sqrt(1+2*para_n)*tau(1,para_n+1).dir13;
        end

        G1_3p=sparse(NS,1);
        G2_3p=sparse(NS,1);
        for para_m=0:para_M
            G1_3p = G1_3p + sqrt(1+2*para_m)*tau(para_m+1,1).dir22;
            G2_3p = G2_3p + sqrt(1+2*para_m)*tau(para_m+1,1).dir23;
        end
        G1_3m=sparse(NS,1);
        G2_3m=sparse(NS,1);
        for para_m=0:para_M
            G1_3m = G1_3m + (-1)^(para_m+1)*sqrt(1+2*para_m)*tau(para_m+1,1).dir22;
            G2_3m = G2_3m + (-1)^(para_m+1)*sqrt(1+2*para_m)*tau(para_m+1,1).dir23;
        end
        G3_3p=sparse(NS,1);
        G4_3p=sparse(NS,1);
        for para_n=0:para_N
            G3_3p = G3_3p + sqrt(1+2*para_n)*tau(1,para_n+1).dir23;
            G4_3p = G4_3p + sqrt(1+2*para_n)*tau(1,para_n+1).dir33;
        end
        G3_3m=sparse(NS,1);
        G4_3m=sparse(NS,1);
        for para_n=0:para_N
            G3_3m = G3_3m + (-1)^(para_n+1)*sqrt(1+2*para_n)*tau(1,para_n+1).dir23;
            G4_3m = G4_3m + (-1)^(para_n+1)*sqrt(1+2*para_n)*tau(1,para_n+1).dir33;
        end
        G5_3p=sparse(NS,1);
        G5_3m=sparse(NS,1);
        for para_n=0:para_N
            G5_3p = G5_3p + sqrt(1+2*para_n)*tau(1,para_n+1).dir13;
            G5_3m = G5_3m + (-1)^(para_n+1)*sqrt(1+2*para_n)*tau(1,para_n+1).dir13;
        end
        G6_3p=sparse(NS,1);
        G6_3m=sparse(NS,1);
        for para_m=0:para_M
            G6_3p = G6_3p + sqrt(1+2*para_m)*tau(para_m+1,1).dir12;
            G6_3m = G6_3m + (-1)^(para_m+1)*sqrt(1+2*para_m)*tau(para_m+1,1).dir12;
        end

        G1p = (C22./C22bar).*(G1_1+(h./l).*G1_2) - G1_3p;
        G1m = (C22./C22bar).*(G1_1+(h./l).*G1_2) - G1_3m;

        G2p = (l./h).*(C44./C33bar).*((l./h).*G2_1+G2_2) - G2_3p;
        G2m = (l./h).*(C44./C33bar).*((l./h).*G2_1+G2_2) - G2_3m;

        G3p = (h./l).*(C44./C22bar).*(G1_1+(h./l).*G1_2) - G3_3p;
        G3m = (h./l).*(C44./C22bar).*(G1_1+(h./l).*G1_2) - G3_3m;

        G4p = (C33./C33bar).*((l./h).*G2_1+G2_2) - G4_3p;
        G4m = (C33./C33bar).*((l./h).*G2_1+G2_2) - G4_3m;

        G5p = (h./l).*(C55./C66bar).*(G3_1+(h./l).*G3_2) - G5_3p;
        G5m = (h./l).*(C55./C66bar).*(G3_1+(h./l).*G3_2) - G5_3m;

        G6p = (C66./C66bar).*(G3_1+(h./l).*G3_2) - G6_3p;
        G6m = (C66./C66bar).*(G3_1+(h./l).*G3_2) - G6_3m;

        [A_i, A_o] = loading_plastic(m, n, R1, R2, R3, R4, R5, R6, G1p, G1m, G2p, G2m, G3p, G3m, G4p, G4m, G5p, G5m, G6p, G6m, crack_row, crack_col, crack_faces_ver, crack_faces_hor, ul, ur, ub, ut);
        Sigma_i_r=Sigma_i+A_i;
        Sigma_i_r(temp_i,:)=[];
        Sigma_o_r=Sigma_o+A_o;
        Sigma_o_r(temp_o,:)=[];

        [Avu_i]=inverse_next(L_i, U_i, P_i, Q_i, Sigma_i_r,Avu_ini,temp_i);
        [Avu_o]=inverse_next(L_o, U_o, P_o, Q_o, Sigma_o_r,Avu_ino,temp_o);

        [eps, Av_e] = total_strains(m, n, h, l, G1_1, G2_1, G3_1, G1_2, G2_2, G3_2, C22, C33, C44, C55, C66, C22bar, C33bar, C55bar, C66bar, Avu_i, Avu_o, Ave11, Ave22, Ave33, Ave23, Ave13, Ave12, ul, ur, ub, ut, array_x);

        %         for i=1:NS
        %             D(:,i) = Av_e(:,i) - A(:,:,i)*av_macro_strains;
        %         end
        %         Rbar=[tau(1,1).dir11 tau(1,1).dir22 tau(1,1).dir33 tau(1,1).dir23 tau(1,1).dir13 tau(1,1).dir12];
        %         Rbar=Rbar';    % transverse
        D_11= Av_e(:,1) - (A11*av_macro_strains(1)+A12*av_macro_strains(2)+A13*av_macro_strains(3)+A14*av_macro_strains(4)+A15*av_macro_strains(5)+A16*av_macro_strains(6))';
        D_21= Av_e(:,2) - (A21*av_macro_strains(1)+A22*av_macro_strains(2)+A23*av_macro_strains(3)+A24*av_macro_strains(4)+A25*av_macro_strains(5)+A26*av_macro_strains(6))';
        D_31= Av_e(:,3) - (A31*av_macro_strains(1)+A32*av_macro_strains(2)+A33*av_macro_strains(3)+A34*av_macro_strains(4)+A35*av_macro_strains(5)+A36*av_macro_strains(6))';
        D_41= Av_e(:,4) - (A41*av_macro_strains(1)+A42*av_macro_strains(2)+A43*av_macro_strains(3)+A44*av_macro_strains(4)+A45*av_macro_strains(5)+A46*av_macro_strains(6))';
        D_51= Av_e(:,5) - (A51*av_macro_strains(1)+A52*av_macro_strains(2)+A53*av_macro_strains(3)+A54*av_macro_strains(4)+A55*av_macro_strains(5)+A56*av_macro_strains(6))';
        D_61= Av_e(:,6) - (A61*av_macro_strains(1)+A62*av_macro_strains(2)+A63*av_macro_strains(3)+A64*av_macro_strains(4)+A65*av_macro_strains(5)+A66*av_macro_strains(6))';

        D=[D_11 D_21 D_31 D_41 D_51 D_61];
        Rbar=[tau(1,1).dir11 tau(1,1).dir22 tau(1,1).dir33 tau(1,1).dir23 tau(1,1).dir13 tau(1,1).dir12];

        %         sig_inel=sparse(6,1);
        %         for i=1:NS
        %             sig_inel=sig_inel + h(i)*l(i)*(Cnext(:,:,i)*D(:,i) - Rbar(:,i) - Gamma(:,:,i)*delta_T);
        %         end
        %         sig_inel = -sig_inel/(H*L)


        for i=1:t
            for j=1:t
                e(i,j).dir11=eps(i,j).dir11-(eps(i,j).dir11+eps(i,j).dir22+eps(i,j).dir33)/3-eps_pl_pre(i,j).dir11;
                e(i,j).dir22=eps(i,j).dir22-(eps(i,j).dir11+eps(i,j).dir22+eps(i,j).dir33)/3-eps_pl_pre(i,j).dir22;
                e(i,j).dir33=eps(i,j).dir33-(eps(i,j).dir11+eps(i,j).dir22+eps(i,j).dir33)/3-eps_pl_pre(i,j).dir33;
                e(i,j).dir23=eps(i,j).dir23-eps_pl_pre(i,j).dir23;
                e(i,j).dir13=eps(i,j).dir13-eps_pl_pre(i,j).dir13;
                e(i,j).dir12=eps(i,j).dir12-eps_pl_pre(i,j).dir12;

                ebar(i,j).eff = sqrt( 2*(e(i,j).dir11.^2+e(i,j).dir22.^2+e(i,j).dir33.^2)/3 + 4*(e(i,j).dir12.^2+e(i,j).dir13.^2+e(i,j).dir23.^2)/3 );
                deps_pl(i,j).eff = sqrt( 2*(deps_pl(i,j).dir11.^2+deps_pl(i,j).dir22.^2+deps_pl(i,j).dir33.^2)/3 + 4*(deps_pl(i,j).dir12.^2+deps_pl(i,j).dir13.^2+deps_pl(i,j).dir23.^2)/3 );

                % % deps_pl_eff(i,j).diff is the difference in the incremental effective plastic
                % % strains (deps_pl, at each point within each subcell) obtained from two
                % % successive iterations for a particular load increment.
                deps_pl_eff(i,j).diff = abs(deps_pl_pre(i,j).eff - deps_pl(i,j).eff);
                deps_pl_pre(i,j).eff = deps_pl(i,j).eff;

                eps_pl(i,j).eff = eps_pl_pre(i,j).eff + deps_pl(i,j).eff;
                sigma(i,j).bar=Y_Sigma+Hp.*eps_pl(i,j).eff;
                dlamda(i,j).expr=1 - sigma(i,j).bar./(3*mu.*ebar(i,j).eff);

                toggle(i,j).temp = dlamda(i,j).expr > 0;
                toggle_matrix(i,j) = any(toggle(i,j).temp);

                deps_pl(i,j).dir11=dlamda(i,j).expr.*e(i,j).dir11.*toggle(i,j).temp;
                deps_pl(i,j).dir22=dlamda(i,j).expr.*e(i,j).dir22.*toggle(i,j).temp;
                deps_pl(i,j).dir33=dlamda(i,j).expr.*e(i,j).dir33.*toggle(i,j).temp;
                deps_pl(i,j).dir23=dlamda(i,j).expr.*e(i,j).dir23.*toggle(i,j).temp;
                deps_pl(i,j).dir13=dlamda(i,j).expr.*e(i,j).dir13.*toggle(i,j).temp;
                deps_pl(i,j).dir12=dlamda(i,j).expr.*e(i,j).dir12.*toggle(i,j).temp;

                eps_pl(i,j).dir11=eps_pl_pre(i,j).dir11+deps_pl(i,j).dir11;
                eps_pl(i,j).dir22=eps_pl_pre(i,j).dir22+deps_pl(i,j).dir22;
                eps_pl(i,j).dir33=eps_pl_pre(i,j).dir33+deps_pl(i,j).dir33;
                eps_pl(i,j).dir23=eps_pl_pre(i,j).dir23+deps_pl(i,j).dir23;
                eps_pl(i,j).dir13=eps_pl_pre(i,j).dir13+deps_pl(i,j).dir13;
                eps_pl(i,j).dir12=eps_pl_pre(i,j).dir12+deps_pl(i,j).dir12;
            end
        end

        % % Check for beginning of yielding. Breaks out of the loop if the
        % % yielding has not begun at any point.
        % % ___________________________
        %         if toggle(3,3).temp == 0
        if toggle_matrix == 0
            break;
        end
        % % ___________________________


        % % Check for convergence (if yielding has begun)
        % % _________________________________________________________________
        % % max_diff is the maximum of deps_pl_eff(i,j).diff (see defn. of deps_pl_eff(i,j).diff above) within each subcell
        max_diff=sparse(NS,1);
        for i=1:t
            for j=1:t
                max_diff=max(max_diff, deps_pl_eff(i,j).diff);
            end
        end

        % % ave_deps_pl_eff is the average value of incremental effective
        % plastic strain in each subcell at each iteration.
        ave_deps_pl_eff=sparse(NS,1);
        for i=1:t
            for j=1:t
                ave_deps_pl_eff = ave_deps_pl_eff + deps_pl_pre(i,j).eff;
            end
        end
        ave_deps_pl_eff=ave_deps_pl_eff/t^2;


        if max_diff <= 0.01*error_tolerance*ave_deps_pl_eff
            break;
        end

        load_array_at_each_iter=[load_array_at_each_iter load_current];
        av_stress_macro_at_each_iter=[av_stress_macro_at_each_iter av_macro_stresses_bar];

    end

    iterations(incr_num)=iter;
    iter

    for i=1:t
        for j=1:t
            eps_pl_pre(i,j).eff=eps_pl_pre(i,j).eff+deps_pl(i,j).eff;
            eps_pl_pre(i,j).dir11=eps_pl_pre(i,j).dir11+deps_pl(i,j).dir11;
            eps_pl_pre(i,j).dir22=eps_pl_pre(i,j).dir22+deps_pl(i,j).dir22;
            eps_pl_pre(i,j).dir33=eps_pl_pre(i,j).dir33+deps_pl(i,j).dir33;
            eps_pl_pre(i,j).dir23=eps_pl_pre(i,j).dir23+deps_pl(i,j).dir23;
            eps_pl_pre(i,j).dir13=eps_pl_pre(i,j).dir13+deps_pl(i,j).dir13;
            eps_pl_pre(i,j).dir12=eps_pl_pre(i,j).dir12+deps_pl(i,j).dir12;
        end
    end

    g=incr_num;
    if z <= incr_register_length
        if  g == incr_register(z)
            %             [field_at_incr(g).Sigma22, field_at_incr(g).Sigma33, field_at_incr(g).Sigma23, field_at_incr(g).U2, field_at_incr(g).U3] = field_plastic(m, n, h, l, G1_1, G2_1, G3_1, G1_2, G2_2, G3_2, C11, C22, C33, C12, C13, C23, C44, C55, C66, C22bar, C33bar, C55bar, C66bar, alpha11, alpha22, alpha33, Avu_i, Avu_o, Ave11, Ave22, Ave33, Ave23, Ave13, Ave12, ul, ur, ub, ut, eps_pl_pre, delta_T, array_x);
            %             [field_at_incr(g).Sigma11, field_at_incr(g).Sigma22, field_at_incr(g).Sigma33, field_at_incr(g).Sigma23, field_at_incr(g).Sigma13, field_at_incr(g).Sigma12, field_at_incr(g).Sigma_HS, field_at_incr(g).Sigma_eff, field_at_incr(g).Strain_p_eff] = field_plastic(m, n, h, l, G1_1, G2_1, G3_1, G1_2, G2_2, G3_2, C11, C22, C33, C12, C13, C23, C44, C55, C66, C22bar, C33bar, C55bar, C66bar, alpha11, alpha22, alpha33, Avu_i, Avu_o, Ave11, Ave22, Ave33, Ave23, Ave13, Ave12, ul, ur, ub, ut, eps_pl_pre, delta_T, array_x);
            %             [field_at_incr(g).Sigma_HS, field_at_incr(g).Sigma_eff, field_at_incr(g).Strain_p_eff] = field_plastic(m, n, h, l, G1_1, G2_1, G3_1, G1_2, G2_2, G3_2, C11, C22, C33, C12, C13, C23, C44, C55, C66, C22bar, C33bar, C55bar, C66bar, alpha11, alpha22, alpha33, Avu_i, Avu_o, Ave11, Ave22, Ave33, Ave23, Ave13, Ave12, ul, ur, ub, ut, eps_pl_pre, delta_T, array_x);
                        [field_at_incr(g).Sigma11, field_at_incr(g).Sigma22, field_at_incr(g).Sigma33, field_at_incr(g).Sigma23, field_at_incr(g).Sigma13, field_at_incr(g).Sigma12, field_at_incr(g).Sigma_HS, field_at_incr(g).Sigma_eff, field_at_incr(g).Strain_p_eff, field_at_incr(g).U1, field_at_incr(g).U2, field_at_incr(g).U3] = field_plastic(m, n, h, l, G1_1, G2_1, G3_1, G1_2, G2_2, G3_2, C11, C22, C33, C12, C13, C23, C44, C55, C66, C22bar, C33bar, C55bar, C66bar, alpha11, alpha22, alpha33, Avu_i, Avu_o, Ave11, Ave22, Ave33, Ave23, Ave13, Ave12, ul, ur, ub, ut, eps_pl_pre, delta_T, array_x);
            %             [field_at_incr(g).Strain11, field_at_incr(g).Strain22, field_at_incr(g).Strain33, field_at_incr(g).Strain23, field_at_incr(g).Strain13, field_at_incr(g).Strain12] = struct2double(m, n, eps, t);
            z=z+1;
        end
    end

    av_macro_stresses=T1*av_macro_stresses_bar;
    load_array = [load_array load_current];
    av_stress_macro=[av_stress_macro av_macro_stresses];
    av_strain_macro=[av_strain_macro av_macro_strains];
    av_stress_macro_bar=[av_stress_macro_bar av_macro_stresses_bar];
    av_strain_macro_bar=[av_strain_macro_bar av_macro_strains_bar];
end

x2=[];
length_total=0;
for j=1:m
    x2=cat(2,x2,array_x*h(j)/2+h(j)/2+length_total);
    length_total=length_total+h(j);
end

x3=[];
height=0;
for t=1:n
    j=m*(t-1)+1;
    x3=cat(2,x3,array_x*l(j)/2+l(j)/2+height);
    height=height+l(j);
end
% %%_________________________________________________________________________

time_elapsed=toc % time stops ticking

fid=fopen(file_output_title,'a');
fprintf(fid,'\nIncrement Number\t\tLoad\t\tNumber of Iterations\tSIGMA11\t\tSIGMA22\t\tSIGMA33\t\tSIGMA23\t\tSIGMA13\t\tSIGMA12 \n');
for i=1:incr_NUM
    fprintf(fid,'\t%d\t\t%12.4E\t\t%d\t\t\t%12.4E\t%12.4E\t%12.4E\t%12.4E\t%12.4E\t%12.4E\n', i, load_array(i+1), iterations(i), av_stress_macro(1,i+1), av_stress_macro(2,i+1), av_stress_macro(3,i+1), av_stress_macro(4,i+1), av_stress_macro(5,i+1), av_stress_macro(6,i+1));
end
fprintf(fid,'\nTotal time taken to run the current analysis\t\t%d seconds, i.e. %d minutes\n', time_elapsed, time_elapsed/60);
fclose(fid);

if LOP == 8
    string_title = strcat(string_title,'_thermal');
else
    string_title = strcat(string_title,'_mechanical');
end
save(string_title, 'x2', 'x3', 'field_at_incr', 'load_array', 'av_stress_macro', 'av_strain_macro', 'av_stress_macro_bar', 'av_strain_macro_bar','load_array_at_each_iter','av_stress_macro_at_each_iter');
% if length(strfind(version,'R14')) > 0
%     save(string_title, 'x2', 'x3', 'field_at_incr', 'load_array', 'av_stress_macro', 'av_strain_macro', 'av_stress_macro_bar', 'av_strain_macro_bar','load_array_at_each_iter','av_stress_macro_at_each_iter','-v6');
% else
%     save(string_title, 'x2', 'x3', 'field_at_incr', 'load_array', 'av_stress_macro', 'av_strain_macro', 'av_stress_macro_bar', 'av_strain_macro_bar','load_array_at_each_iter','av_stress_macro_at_each_iter');
% end

% [u11, u22, u33, S11, S22, S33, S23, S13, S12, X2, X3, u1, u2, u3, Sigma11, Sigma22, Sigma33, Sigma23, Sigma13, Sigma12, x2, x3] = field(m, n, h, l, C11, C22, C33, C12, C13, C23, C44, C55, C66, C22bar, C23bar, C32bar, C33bar, C55bar, C66bar, Avu_i, Avu_o, Ave11, Ave22, Ave33, Ave23, Ave13, Ave12, ul, ur, ub, ut);
% % [u2, u3, x2, x3] = field_special(m, n, h, l, C11, C22, C33, C12, C13, C23, C44, C55, C66, C22bar, C23bar, C32bar, C33bar, C55bar, C66bar, Avu_i, Avu_o, Ave11, Ave22, Ave33, Ave23, Ave13, Ave12, ul, ur, ub, ut);
% save(strcat(string_title,'_deformation'), 'u11', 'u22', 'u33', 'S11', 'S22', 'S33', 'S23', 'S13', 'S12', 'X2', 'X3', 'u1', 'u2', 'u3', 'Sigma11', 'Sigma22', 'Sigma33', 'Sigma23', 'Sigma13', 'Sigma12', 'x2', 'x3', 'Cstar');


% % m is the number of subcells in the horizontal direction
% % n is the number of subcells in the vertical direction
% % size(x2) = t*m (integration points)
% % size(x3) = t*n (integration points)
                                          % EXECUTES THE FVDAM CODE
end