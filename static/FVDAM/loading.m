function [A_i, A_o] = loading(m, n, R1, R2, R3, R4, R5, R6, crack_row, crack_col, crack_int_ver, crack_int_hor, ul, ur, ub, ut);

% Generates the loading due to applied macroscopic strains


%% In Plane
%__________________________________________________________________________
% 1
A_i_1=zeros(ur(m*n)-n,1);
for i=0:(n-1)
    t=i*m;
    if ismember(i+1,crack_row) == 0           
        A_i_1(ul(t+1)-i)=R1(t+1)-R1(t+m);
        for k=2:m
            A_i_1(ul(t+k)-i)=-R1(t+k-1)+R1(t+k);    
        end
    else
        A_i_1(ul(t+1)-i)=R1(t+1)-R1(t+m);
        for k=2:m
            if ismember(i*(m+1)+k, crack_int_ver) == 0
                A_i_1(ul(t+k)-i)=-R1(t+k-1)+R1(t+k);    
            else
                A_i_1(ur(t+k-1)-i)=-R1(t+k-1);    
                A_i_1(ul(t+k)-i)=R1(t+k);    
            end
        end
    end
end

% 2
A_i_2=zeros(ur(m*n)-n,1);
for i=0:(n-1)
    t=i*m;
    if ismember(i+1,crack_row) == 0           
        A_i_2(ul(t+1)-i)=R2(t+1)-R2(t+m);
        for k=2:m
            A_i_2(ul(t+k)-i)=-R2(t+k-1)+R2(t+k);    
        end
    else
        A_i_2(ul(t+1)-i)=R2(t+1)-R2(t+m);
        for k=2:m
            if ismember(i*(m+1)+k, crack_int_ver) == 0
                A_i_2(ul(t+k)-i)=-R2(t+k-1)+R2(t+k);    
            else
                A_i_2(ur(t+k-1)-i)=-R2(t+k-1);    
                A_i_2(ul(t+k)-i)=R2(t+k);    
            end
        end
    end
end

% 3
A_i_3=zeros(ut(m*n)-m,1);
for i=0:(m-1)
    if ismember(i+1,crack_col) == 0
        A_i_3(ub(i+1)-i)=R3(i+1)-R3(i+1+m*(n-1));
        for k=2:n
            A_i_3(ub(i+1+m*(k-1))-i)=-R3(i+1+m*(k-2))+R3(i+1+m*(k-1));
        end
    else
        A_i_3(ub(i+1)-i)=R3(i+1)-R3(i+1+m*(n-1));
        for k=2:n
            if ismember(i*(n+1)+k, crack_int_hor) == 0
                A_i_3(ub(i+1+m*(k-1))-i)=-R3(i+1+m*(k-2))+R3(i+1+m*(k-1));
            else
                A_i_3(ut(i+1+m*(k-2))-i)=-R3(i+1+m*(k-2));
                A_i_3(ub(i+1+m*(k-1))-i)=R3(i+1+m*(k-1));
            end
        end
    end
end

% 4
A_i_4=zeros(ut(m*n)-m,1);
for i=0:(m-1)
    if ismember(i+1,crack_col) == 0
        A_i_4(ub(i+1)-i)=R4(i+1)-R4(i+1+m*(n-1));
        for k=2:n
            A_i_4(ub(i+1+m*(k-1))-i)=-R4(i+1+m*(k-2))+R4(i+1+m*(k-1));
        end
    else
        A_i_4(ub(i+1)-i)=R4(i+1)-R4(i+1+m*(n-1));
        for k=2:n
            if ismember(i*(n+1)+k, crack_int_hor) == 0
                A_i_4(ub(i+1+m*(k-1))-i)=-R4(i+1+m*(k-2))+R4(i+1+m*(k-1));
            else
                A_i_4(ut(i+1+m*(k-2))-i)=-R4(i+1+m*(k-2));
                A_i_4(ub(i+1+m*(k-1))-i)=R4(i+1+m*(k-1));
            end
        end
    end
end
%__________________________________________________________________________


%% Out of Plane
%__________________________________________________________________________
% 1
A_o_1=zeros(ur(m*n)-n,1);
for i=0:(n-1)
    t=i*m;
    if ismember(i+1,crack_row) == 0           
        A_o_1(ul(t+1)-i)=R6(t+1)-R6(t+m);
        for k=2:m
            A_o_1(ul(t+k)-i)=-R6(t+k-1)+R6(t+k);    
        end
    else
        A_o_1(ul(t+1)-i)=R6(t+1)-R6(t+m);
        for k=2:m
            if ismember(i*(m+1)+k, crack_int_ver) == 0
                A_o_1(ul(t+k)-i)=-R6(t+k-1)+R6(t+k);    
            else
                A_o_1(ur(t+k-1)-i)=-R6(t+k-1);    
                A_o_1(ul(t+k)-i)=R6(t+k);    
            end
        end
    end
end

% 2
A_o_2=zeros(ut(m*n)-m,1);
for i=0:(m-1)
    if ismember(i+1,crack_col) == 0
        A_o_2(ub(i+1)-i)=R5(i+1)-R5(i+1+m*(n-1));
        for k=2:n
            A_o_2(ub(i+1+m*(k-1))-i)=-R5(i+1+m*(k-2))+R5(i+1+m*(k-1));
        end
    else
        A_o_2(ub(i+1)-i)=R5(i+1)-R5(i+1+m*(n-1));
        for k=2:n
            if ismember(i*(n+1)+k, crack_int_hor) == 0
                A_o_2(ub(i+1+m*(k-1))-i)=-R5(i+1+m*(k-2))+R5(i+1+m*(k-1));
            else
                A_o_2(ut(i+1+m*(k-2))-i)=-R5(i+1+m*(k-2));
                A_o_2(ub(i+1+m*(k-1))-i)=R5(i+1+m*(k-1));
            end
        end
    end
end
%__________________________________________________________________________


A_i=[A_i_1; A_i_2; A_i_3; A_i_4];
A_o=[A_o_1; A_o_2];
