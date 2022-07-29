function [crack_row, crack_col, crack_faces_ver, crack_faces_hor, ul, ur, ub, ut] = crack(m, n, crack_num, crack_dir, crack_subcell_beta, crack_subcell_gama1,....
    crack_subcell_gama2, crack_subcell_gama, crack_subcell_beta1, crack_subcell_beta2)

% Creates an array containing the subcell interfaces which constitute a crack

% m: Total number of subcells in the X_2 direction in the entire grid 
% n: Total number of subcells in the X_3 direction in the entire grid 
% crack_num: Total number of cracks
% crack_dir: Direction of each crack

% CRACK IS BEING MODELED CONSIDERING IT LIES BETWEEN THE SUBCELLS NUMBERED CRACK_SUBCELL_1 AND CRACK_SUBCELL_2, 
% ON THE TOP IF IT IS IN THE X_2 DIRECTION AND ON THE RIGHT IF IT IS IN THE X_3 DIRECTION  

crack_faces_ver=[];   % vertical subcell interfaces constituting the crack
crack_faces_hor=[];   % horizontal subcell interfaces constituting the crack

for k=1:crack_num    
    if crack_dir(k) == 2
        crack_faces_hor=[crack_faces_hor (crack_subcell_beta1(k)-1)*(n+1)+crack_subcell_gama(k)+1:n+1:(crack_subcell_beta2(k)-1)*(n+1)+crack_subcell_gama(k)+1];        
    elseif crack_dir(k) == 3
        crack_faces_ver=[crack_faces_ver (crack_subcell_gama1(k)-1)*(m+1)+crack_subcell_beta(k)+1:m+1:(crack_subcell_gama2(k)-1)*(m+1)+crack_subcell_beta(k)+1];        
    end
end

crack_row=[];   % rows containing cracks
crack_col=[];   % columns containing cracks
for k=1:crack_num
    if crack_dir(k) == 2
        crack_col=[crack_col crack_subcell_beta1(k):crack_subcell_beta2(k)];        
    elseif crack_dir(k) == 3
        crack_row=[crack_row crack_subcell_gama1(k):crack_subcell_gama2(k)];        
    end
end
crack_row=unique(crack_row);
crack_col=unique(crack_col);



% p=n*(m+1);
% q=m*(n+1);
% Ttemp_crack=[crack_faces_ver p+crack_faces_hor];
% temp_crack=[crack_faces_ver p+crack_faces_ver 2*p+crack_faces_hor 2*p+q+crack_faces_hor];


t=1;
for i=0:n-1
    for j=1:m
        k=i*m+j;
        x=i*(m+1)+j;
        if ismember(x,crack_faces_ver) == 0
            ul(k) = t;
            ur(k) = t+1;
            t=t+1;
        else
            t=t+1;
            ul(k) = t;
            ur(k) = t+1;
            t=t+1;
        end
    end
    t=t+1;
end

t=1;
for i=0:m-1
    for j=1:n
        k=i+1+m*(j-1);
        x=i*(n+1)+j;
        if ismember(x,crack_faces_hor) == 0
            ub(k) = t;
            ut(k) = t+1;
            t=t+1;
        else
            t=t+1;
            ub(k) = t;
            ut(k) = t+1;
            t=t+1;
        end
    end
    t=t+1;
end
