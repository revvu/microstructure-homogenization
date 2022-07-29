function [K, Knext] = stiffness_cracks(m, n, h, l, C22, C33, C12, C13, C23, C44, C55, C66, C22bar, C23bar, C32bar, C33bar, C55bar, C66bar, crack_row, crack_col, crack_int_ver, crack_int_hor, ul, ur, ub, ut)

% Generates the global stiffness matrix

% This function is strictly for a rectangular grid. It however can contain
% any number of cracks. The assembly of global stiffness matrix is
% cumbersome and may be difficult to follow. However, the global stiffness matrix can
% be assembled in several different and more efficient ways.

% Calculating the elements of the local stiffness matrix for each subcell
% %-----------------------------------------------------------------------------------------------------------------
L11=1./h.*(4-3*C22./C22bar).*C22;
L12=1./h.*(2-3*C22./C22bar).*C22;
L13=L12;
L14=L11;    

N11=-3*h.*C22.*C44./(C22bar.*l.^2);
N12=N11;
N13=N11;
N14=N11;    

O11=C23./l;
O12=-O11;
O13=-O11;
O14=O11;    

M21=1./h.*(4-3*C44./C23bar).*C44;
M22=1./h.*(2-3*C44./C23bar).*C44;
M23=M22;
M24=M21;    

N21=C44./l;
N22=-N21; 
N23=-N21;
N24=N21;    

O21=-3*h.*C33.*C44./(C23bar.*l.^2);
O22=O21;
O23=O21;
O24=O21;    

L31=-3*l.*C22.*C44./(C32bar.*h.^2);
L32=L31;
L33=L31;
L34=L31;    

M31=C44./h;
M32=-M31;
M33=-M31;
M34=M31;

N31=1./l.*(4-3*C44./C32bar).*C44;
N32=1./l.*(2-3*C44./C32bar).*C44;
N33=N32;
N34=N31;    

L41=C23./h;
L42=-L41;
L43=-L41;
L44=L41;    

M41=-3*l.*C33.*C44./(C33bar.*h.^2);
M42=M41;
M43=M41;
M44=M41;    

O41=1./l.*(4-3*C33./C33bar).*C33;
O42=1./l.*(2-3*C33./C33bar).*C33;
O43=O42;
O44=O41;    

U11=1./h.*(4-3*C66./C66bar).*C66;
U12=1./h.*(2-3*C66./C66bar).*C66;
U13=U12;
U14=U11;

V11=-3*C55.*C66./(C55bar.*h);
V12=V11;
V13=V11;
V14=V11;

U21=-3*C55.*C66./(C66bar.*l);
U22=U21;
U23=U21;
U24=U21;

V21=1./l.*(4-3*C55./C55bar).*C55;
V22=1./l.*(2-3*C55./C55bar).*C55;
V23=V22;
V24=V21;
% %-----------------------------------------------------------------------------------------------------------------


% Generating the individual submatrices of the global stiffness matrix
% %_________________________________________________________________________________________________________________
% [1,1]
K_11=[];
for i=0:(n-1)
    t=i*m;
    K11=sparse(ur(t+m)-ul(t+1),ur(t+m)-ul(t+1));
    if ismember(i+1,crack_row) == 0
        K11(1,1)=L14(t+1)+L11(t+m);
        K11(1,2)=L13(t+1);
        K11(1,m)=L12(t+m);
        for k=2:m-1
            K11(k, k-1)=L12(t+k-1);
            K11(k, k)=L11(t+k-1)+L14(t+k);
            K11(k, k+1)=L13(t+k);
        end
        K11(m,m-1)=L12(t+m-1);
        K11(m,m)=L11(t+m-1)+L14(t+m);
        K11(m,1)=L13(t+m);
    else
        K11(1,1)=L14(t+1)+L11(t+m);
        K11(1,2)=L13(t+1);
        K11(1,ur(t+m)-ul(t+1))=L12(t+m);
        k=2;
        for kp=2:m-1
            if ismember(i*(m+1)+kp, crack_int_ver) == 0
                K11(k, k-1)=L12(t+kp-1);
                K11(k, k)=L11(t+kp-1)+L14(t+kp);
                K11(k, k+1)=L13(t+kp);
                k=k+1;
            else
                K11(k, k-1)=L12(t+kp-1);
                K11(k, k)=L11(t+kp-1);
                k=k+1;
                K11(k, k)=L14(t+kp);
                K11(k, k+1)=L13(t+kp);
                k=k+1;
            end
        end
        if ismember(i*(m+1)+m, crack_int_ver) == 0
            K11(k,k-1)=L12(t+m-1);
            K11(k,k)=L11(t+m-1)+L14(t+m);
            K11(k,1)=L13(t+m);
        else
            K11(k,k-1)=L12(t+m-1);
            K11(k,k)=L11(t+m-1);
            k=k+1;
            K11(k,k)=L14(t+m);
            K11(k,1)=L13(t+m);
        end
    end
    [x,y]=size(K11);
    [p,q]=size(K_11);
    K_11=[K_11 sparse(p,y); sparse(x,q) K11];
end

%[1,2]
[p,q]=size(K_11);
K_12=sparse(p,q);

%[1,3]
K_13=[];
for i=0:(n-2)
    t=i*m;
    K13=sparse(ur(m*(i+1))-ul(t+1),ut(m*n)-m);
    if ismember(i+1,crack_row) == 0   
        K13(1,ub(t+1))=N14(t+1);
        K13(1,ut(t+1))=N13(t+1);
        K13(1,ub(t+m)-(m-1))=N12(t+m);
        K13(1,ut(t+m)-(m-1))=N11(t+m);        
        for k=2:m
            K13(k, ub(t+k-1)-(k-2))=N12(t+k-1);
            K13(k, ut(t+k-1)-(k-2))=N11(t+k-1);
            K13(k, ub(t+k)-(k-1))=N14(t+k);
            K13(k, ut(t+k)-(k-1))=N13(t+k);
        end
    else
        K13(1,ub(t+1))=N14(t+1);
        K13(1,ut(t+1))=N13(t+1);
        K13(1,ub(t+m)-(m-1))=N12(t+m);
        K13(1,ut(t+m)-(m-1))=N11(t+m);
        k=2;
        for kp=2:m
            if ismember(i*(m+1)+kp, crack_int_ver) == 0
                K13(k, ub(t+kp-1)-(kp-2))=N12(t+kp-1);
                K13(k, ut(t+kp-1)-(kp-2))=N11(t+kp-1);
                K13(k, ub(t+kp)-(kp-1))=N14(t+kp);
                K13(k, ut(t+kp)-(kp-1))=N13(t+kp);
                k=k+1;
            else
                K13(k, ub(t+kp-1)-(kp-2))=N12(t+kp-1);
                K13(k, ut(t+kp-1)-(kp-2))=N11(t+kp-1);
                k=k+1;
                K13(k, ub(t+kp)-(kp-1))=N14(t+kp);
                K13(k, ut(t+kp)-(kp-1))=N13(t+kp);
                k=k+1;
            end
        end
    end
    K_13=[K_13; K13];
end
i=n-1;
t=i*m;
K13=sparse(ur(m*n)-ul(t+1),ut(m*n)-m);
if ismember(n,crack_row) == 0   
    K13(1,ub(t+1))=N14(t+1);
    K13(1,1)=N13(t+1);
    K13(1,ub(t+m)-(m-1))=N12(t+m);
    K13(1,ub(m)-(m-1))=N11(t+m);        
    for k=2:m
        K13(k, ub(t+k-1)-(k-2))=N12(t+k-1);
        K13(k, ub(k-1)-(k-2))=N11(t+k-1);
        K13(k, ub(t+k)-(k-1))=N14(t+k);
        K13(k, ub(k)-(k-1))=N13(t+k);
    end
else
    K13(1,ub(t+1))=N14(t+1);
    K13(1,1)=N13(t+1);
    K13(1,ub(t+m)-(m-1))=N12(t+m);
    K13(1,ub(m)-(m-1))=N11(t+m);        
    k=2;
    for kp=2:m
        if ismember(i*(m+1)+kp, crack_int_ver) == 0
            K13(k, ub(t+kp-1)-(kp-2))=N12(t+kp-1);
            K13(k, ub(kp-1)-(kp-2))=N11(t+kp-1);
            K13(k, ub(t+kp)-(kp-1))=N14(t+kp);
            K13(k, ub(kp)-(kp-1))=N13(t+kp);
            k=k+1;
        else
            K13(k, ub(t+kp-1)-(kp-2))=N12(t+kp-1);
            K13(k, ub(kp-1)-(kp-2))=N11(t+kp-1);
            k=k+1;
            K13(k, ub(t+kp)-(kp-1))=N14(t+kp);
            K13(k, ub(kp)-(kp-1))=N13(t+kp);
            k=k+1;
        end
    end
end
K_13=[K_13; K13];

%[1,4]
K_14=[];
for i=0:(n-2)
    t=i*m;
    K14=sparse(ur(m*(i+1))-ul(t+1),ut(m*n)-m);
    if ismember(i+1,crack_row) == 0   
        K14(1,ub(t+1))=O14(t+1);
        K14(1,ut(t+1))=O13(t+1);
        K14(1,ub(t+m)-(m-1))=O12(t+m);
        K14(1,ut(t+m)-(m-1))=O11(t+m);        
        for k=2:m
            K14(k, ub(t+k-1)-(k-2))=O12(t+k-1);
            K14(k, ut(t+k-1)-(k-2))=O11(t+k-1);
            K14(k, ub(t+k)-(k-1))=O14(t+k);
            K14(k, ut(t+k)-(k-1))=O13(t+k);
        end
    else
        K14(1,ub(t+1))=O14(t+1);
        K14(1,ut(t+1))=O13(t+1);
        K14(1,ub(t+m)-(m-1))=O12(t+m);
        K14(1,ut(t+m)-(m-1))=O11(t+m);
        k=2;
        for kp=2:m
            if ismember(i*(m+1)+kp, crack_int_ver) == 0
                K14(k, ub(t+kp-1)-(kp-2))=O12(t+kp-1);
                K14(k, ut(t+kp-1)-(kp-2))=O11(t+kp-1);
                K14(k, ub(t+kp)-(kp-1))=O14(t+kp);
                K14(k, ut(t+kp)-(kp-1))=O13(t+kp);
                k=k+1;
            else
                K14(k, ub(t+kp-1)-(kp-2))=O12(t+kp-1);
                K14(k, ut(t+kp-1)-(kp-2))=O11(t+kp-1);
                k=k+1;
                K14(k, ub(t+kp)-(kp-1))=O14(t+kp);
                K14(k, ut(t+kp)-(kp-1))=O13(t+kp);
                k=k+1;
            end
        end
    end
    K_14=[K_14; K14];
end
i=n-1;
t=i*m;
K14=sparse(ur(m*n)-ul(t+1),ut(m*n)-m);
if ismember(n,crack_row) == 0   
    K14(1,ub(t+1))=O14(t+1);
    K14(1,1)=O13(t+1);
    K14(1,ub(t+m)-(m-1))=O12(t+m);
    K14(1,ub(m)-(m-1))=O11(t+m);        
    for k=2:m
        K14(k, ub(t+k-1)-(k-2))=O12(t+k-1);
        K14(k, ub(k-1)-(k-2))=O11(t+k-1);
        K14(k, ub(t+k)-(k-1))=O14(t+k);
        K14(k, ub(k)-(k-1))=O13(t+k);
    end
else
    K14(1,ub(t+1))=O14(t+1);
    K14(1,1)=O13(t+1);
    K14(1,ub(t+m)-(m-1))=O12(t+m);
    K14(1,ub(m)-(m-1))=O11(t+m);        
    k=2;
    for kp=2:m
        if ismember(i*(m+1)+kp, crack_int_ver) == 0
            K14(k, ub(t+kp-1)-(kp-2))=O12(t+kp-1);
            K14(k, ub(kp-1)-(kp-2))=O11(t+kp-1);
            K14(k, ub(t+kp)-(kp-1))=O14(t+kp);
            K14(k, ub(kp)-(kp-1))=O13(t+kp);
            k=k+1;
        else
            K14(k, ub(t+kp-1)-(kp-2))=O12(t+kp-1);
            K14(k, ub(kp-1)-(kp-2))=O11(t+kp-1);
            k=k+1;
            K14(k, ub(t+kp)-(kp-1))=O14(t+kp);
            K14(k, ub(kp)-(kp-1))=O13(t+kp);
            k=k+1;
        end
    end
end
K_14=[K_14; K14];


%[2,1]
[p,q]=size(K_12);
K_21=sparse(p,q);

%[2,2]
K_22=[];
for i=0:(n-1)
    t=i*m;
    K22=sparse(ur(t+m)-ul(t+1),ur(t+m)-ul(t+1));
    if ismember(i+1,crack_row) == 0
        K22(1,1)=M24(t+1)+M21(t+m);
        K22(1,2)=M23(t+1);
        K22(1,m)=M22(t+m);
        for k=2:m-1
            K22(k, k-1)=M22(t+k-1);
            K22(k, k)=M21(t+k-1)+M24(t+k);
            K22(k, k+1)=M23(t+k);
        end
        K22(m,m-1)=M22(t+m-1);
        K22(m,m)=M21(t+m-1)+M24(t+m);
        K22(m,1)=M23(t+m);
    else
        K22(1,1)=M24(t+1)+M21(t+m);
        K22(1,2)=M23(t+1);
        K22(1,ur(t+m)-ul(t+1))=M22(t+m);
        k=2;
        for kp=2:m-1
            if ismember(i*(m+1)+kp, crack_int_ver) == 0
                K22(k, k-1)=M22(t+kp-1);
                K22(k, k)=M21(t+kp-1)+M24(t+kp);
                K22(k, k+1)=M23(t+kp);
                k=k+1;
            else
                K22(k, k-1)=M22(t+kp-1);
                K22(k, k)=M21(t+kp-1);
                k=k+1;
                K22(k, k)=M24(t+kp);
                K22(k, k+1)=M23(t+kp);
                k=k+1;
            end
        end
        if ismember(i*(m+1)+m, crack_int_ver) == 0
            K22(k,k-1)=M22(t+m-1);
            K22(k,k)=M21(t+m-1)+M24(t+m);
            K22(k,1)=M23(t+m);
        else
            K22(k,k-1)=M22(t+m-1);
            K22(k,k)=M21(t+m-1);
            k=k+1;
            K22(k,k)=M24(t+m);
            K22(k,1)=M23(t+m);
        end
    end
    [x,y]=size(K22);
    [p,q]=size(K_22);
    K_22=[K_22 sparse(p,y); sparse(x,q) K22];
end

%[2,3]
K_23=[];
for i=0:(n-2)
    t=i*m;
    K23=sparse(ur(m*(i+1))-ul(t+1),ut(m*n)-m);
    if ismember(i+1,crack_row) == 0   
        K23(1,ub(t+1))=N24(t+1);
        K23(1,ut(t+1))=N23(t+1);
        K23(1,ub(t+m)-(m-1))=N22(t+m);
        K23(1,ut(t+m)-(m-1))=N21(t+m);        
        for k=2:m
            K23(k, ub(t+k-1)-(k-2))=N22(t+k-1);
            K23(k, ut(t+k-1)-(k-2))=N21(t+k-1);
            K23(k, ub(t+k)-(k-1))=N24(t+k);
            K23(k, ut(t+k)-(k-1))=N23(t+k);
        end
    else
        K23(1,ub(t+1))=N24(t+1);
        K23(1,ut(t+1))=N23(t+1);
        K23(1,ub(t+m)-(m-1))=N22(t+m);
        K23(1,ut(t+m)-(m-1))=N21(t+m);
        k=2;
        for kp=2:m
            if ismember(i*(m+1)+kp, crack_int_ver) == 0
                K23(k, ub(t+kp-1)-(kp-2))=N22(t+kp-1);
                K23(k, ut(t+kp-1)-(kp-2))=N21(t+kp-1);
                K23(k, ub(t+kp)-(kp-1))=N24(t+kp);
                K23(k, ut(t+kp)-(kp-1))=N23(t+kp);
                k=k+1;
            else
                K23(k, ub(t+kp-1)-(kp-2))=N22(t+kp-1);
                K23(k, ut(t+kp-1)-(kp-2))=N21(t+kp-1);
                k=k+1;
                K23(k, ub(t+kp)-(kp-1))=N24(t+kp);
                K23(k, ut(t+kp)-(kp-1))=N23(t+kp);
                k=k+1;
            end
        end
    end
    K_23=[K_23; K23];
end
i=n-1;
t=i*m;
K23=sparse(ur(m*n)-ul(t+1),ut(m*n)-m);
if ismember(n,crack_row) == 0   
    K23(1,ub(t+1))=N24(t+1);
    K23(1,1)=N23(t+1);
    K23(1,ub(t+m)-(m-1))=N22(t+m);
    K23(1,ub(m)-(m-1))=N21(t+m);        
    for k=2:m
        K23(k, ub(t+k-1)-(k-2))=N22(t+k-1);
        K23(k, ub(k-1)-(k-2))=N21(t+k-1);
        K23(k, ub(t+k)-(k-1))=N24(t+k);
        K23(k, ub(k)-(k-1))=N23(t+k);
    end
else
    K23(1,ub(t+1))=N24(t+1);
    K23(1,1)=N23(t+1);
    K23(1,ub(t+m)-(m-1))=N22(t+m);
    K23(1,ub(m)-(m-1))=N21(t+m);        
    k=2;
    for kp=2:m
        if ismember(i*(m+1)+kp, crack_int_ver) == 0
            K23(k, ub(t+kp-1)-(kp-2))=N22(t+kp-1);
            K23(k, ub(kp-1)-(kp-2))=N21(t+kp-1);
            K23(k, ub(t+kp)-(kp-1))=N24(t+kp);
            K23(k, ub(kp)-(kp-1))=N23(t+kp);
            k=k+1;
        else
            K23(k, ub(t+kp-1)-(kp-2))=N22(t+kp-1);
            K23(k, ub(kp-1)-(kp-2))=N21(t+kp-1);
            k=k+1;
            K23(k, ub(t+kp)-(kp-1))=N24(t+kp);
            K23(k, ub(kp)-(kp-1))=N23(t+kp);
            k=k+1;
        end
    end
end
K_23=[K_23; K23];

%[2,4]
K_24=[];
for i=0:(n-2)
    t=i*m;
    K24=sparse(ur(m*(i+1))-ul(t+1),ut(m*n)-m);
    if ismember(i+1,crack_row) == 0   
        K24(1,ub(t+1))=O24(t+1);
        K24(1,ut(t+1))=O23(t+1);
        K24(1,ub(t+m)-(m-1))=O22(t+m);
        K24(1,ut(t+m)-(m-1))=O21(t+m);        
        for k=2:m
            K24(k, ub(t+k-1)-(k-2))=O22(t+k-1);
            K24(k, ut(t+k-1)-(k-2))=O21(t+k-1);
            K24(k, ub(t+k)-(k-1))=O24(t+k);
            K24(k, ut(t+k)-(k-1))=O23(t+k);
        end
    else
        K24(1,ub(t+1))=O24(t+1);
        K24(1,ut(t+1))=O23(t+1);
        K24(1,ub(t+m)-(m-1))=O22(t+m);
        K24(1,ut(t+m)-(m-1))=O21(t+m);
        k=2;
        for kp=2:m
            if ismember(i*(m+1)+kp, crack_int_ver) == 0
                K24(k, ub(t+kp-1)-(kp-2))=O22(t+kp-1);
                K24(k, ut(t+kp-1)-(kp-2))=O21(t+kp-1);
                K24(k, ub(t+kp)-(kp-1))=O24(t+kp);
                K24(k, ut(t+kp)-(kp-1))=O23(t+kp);
                k=k+1;
            else
                K24(k, ub(t+kp-1)-(kp-2))=O22(t+kp-1);
                K24(k, ut(t+kp-1)-(kp-2))=O21(t+kp-1);
                k=k+1;
                K24(k, ub(t+kp)-(kp-1))=O24(t+kp);
                K24(k, ut(t+kp)-(kp-1))=O23(t+kp);
                k=k+1;
            end
        end
    end
    K_24=[K_24; K24];
end
i=n-1;
t=i*m;
K24=sparse(ur(m*n)-ul(t+1),ut(m*n)-m);
if ismember(n,crack_row) == 0   
    K24(1,ub(t+1))=O24(t+1);
    K24(1,1)=O23(t+1);
    K24(1,ub(t+m)-(m-1))=O22(t+m);
    K24(1,ub(m)-(m-1))=O21(t+m);        
    for k=2:m
        K24(k, ub(t+k-1)-(k-2))=O22(t+k-1);
        K24(k, ub(k-1)-(k-2))=O21(t+k-1);
        K24(k, ub(t+k)-(k-1))=O24(t+k);
        K24(k, ub(k)-(k-1))=O23(t+k);
    end
else
    K24(1,ub(t+1))=O24(t+1);
    K24(1,1)=O23(t+1);
    K24(1,ub(t+m)-(m-1))=O22(t+m);
    K24(1,ub(m)-(m-1))=O21(t+m);        
    k=2;
    for kp=2:m
        if ismember(i*(m+1)+kp, crack_int_ver) == 0
            K24(k, ub(t+kp-1)-(kp-2))=O22(t+kp-1);
            K24(k, ub(kp-1)-(kp-2))=O21(t+kp-1);
            K24(k, ub(t+kp)-(kp-1))=O24(t+kp);
            K24(k, ub(kp)-(kp-1))=O23(t+kp);
            k=k+1;
        else
            K24(k, ub(t+kp-1)-(kp-2))=O22(t+kp-1);
            K24(k, ub(kp-1)-(kp-2))=O21(t+kp-1);
            k=k+1;
            K24(k, ub(t+kp)-(kp-1))=O24(t+kp);
            K24(k, ub(kp)-(kp-1))=O23(t+kp);
            k=k+1;
        end
    end
end
K_24=[K_24; K24];


%[3,1]
K_31=[];
for i=0:(m-2)
    K31=sparse(ut(i+1+m*(n-1))-ub(i+1),ur(m*n)-n);
    if ismember(i+1,crack_col) == 0   
        K31(1,ul(i+1))=L34(i+1);
        K31(1,ur(i+1))=L33(i+1);
        K31(1,ul(i+1+m*(n-1))-(n-1))=L32(i+1+m*(n-1));
        K31(1,ur(i+1+m*(n-1))-(n-1))=L31(i+1+m*(n-1));
        for k=2:n
            K31(k, ul(i+1+m*(k-2))-(k-2))=L32(i+1+m*(k-2));
            K31(k, ur(i+1+m*(k-2))-(k-2))=L31(i+1+m*(k-2));
            K31(k, ul(i+1+m*(k-1))-(k-1))=L34(i+1+m*(k-1));
            K31(k, ur(i+1+m*(k-1))-(k-1))=L33(i+1+m*(k-1));
        end
    else
        K31(1,ul(i+1))=L34(i+1);
        K31(1,ur(i+1))=L33(i+1);
        K31(1,ul(i+1+m*(n-1))-(n-1))=L32(i+1+m*(n-1));
        K31(1,ur(i+1+m*(n-1))-(n-1))=L31(i+1+m*(n-1));
        k=2;
        for kp=2:n
            if ismember(i*(n+1)+kp, crack_int_hor) == 0
                K31(k, ul(i+1+m*(kp-2))-(kp-2))=L32(i+1+m*(kp-2));
                K31(k, ur(i+1+m*(kp-2))-(kp-2))=L31(i+1+m*(kp-2));
                K31(k, ul(i+1+m*(kp-1))-(kp-1))=L34(i+1+m*(kp-1));
                K31(k, ur(i+1+m*(kp-1))-(kp-1))=L33(i+1+m*(kp-1));
                k=k+1;
            else
                K31(k, ul(i+1+m*(kp-2))-(kp-2))=L32(i+1+m*(kp-2));
                K31(k, ur(i+1+m*(kp-2))-(kp-2))=L31(i+1+m*(kp-2));
                k=k+1;
                K31(k, ul(i+1+m*(kp-1))-(kp-1))=L34(i+1+m*(kp-1));
                K31(k, ur(i+1+m*(kp-1))-(kp-1))=L33(i+1+m*(kp-1));
                k=k+1;
            end
        end
    end
    K_31=[K_31; K31];
end
i=m-1;
K31=sparse(ut(i+1+m*(n-1))-ub(i+1),ur(m*n)-n);
if ismember(m,crack_col) == 0   
    K31(1,ul(i+1))=L34(i+1);
    K31(1,ul(1))=L33(i+1);
    K31(1,ul(i+1+m*(n-1))-(n-1))=L32(i+1+m*(n-1));
    K31(1,ul(m*(n-1)+1)-(n-1))=L31(i+1+m*(n-1));
    for k=2:n
        K31(k, ul(i+1+m*(k-2))-(k-2))=L32(i+1+m*(k-2));
        K31(k, ul(1+m*(k-2))-(k-2))=L31(i+1+m*(k-2));
        K31(k, ul(i+1+m*(k-1))-(k-1))=L34(i+1+m*(k-1));
        K31(k, ul(1+m*(k-1))-(k-1))=L33(i+1+m*(k-1));
    end
else
    K31(1,ul(i+1))=L34(i+1);
    K31(1,ul(1))=L33(i+1);
    K31(1,ul(i+1+m*(n-1))-(n-1))=L32(i+1+m*(n-1));
    K31(1,ul(1+m*(n-1))-(n-1))=L31(i+1+m*(n-1));
    k=2;
    for kp=2:n
        if ismember(i*(n+1)+kp, crack_int_hor) == 0
            K31(k, ul(i+1+m*(kp-2))-(kp-2))=L32(i+1+m*(kp-2));
            K31(k, ul(1+m*(kp-2))-(kp-2))=L31(i+1+m*(kp-2));
            K31(k, ul(i+1+m*(kp-1))-(kp-1))=L34(i+1+m*(kp-1));
            K31(k, ul(1+m*(kp-1))-(kp-1))=L33(i+1+m*(kp-1));
            k=k+1;
        else
            K31(k, ul(i+1+m*(kp-2))-(kp-2))=L32(i+1+m*(kp-2));
            K31(k, ul(1+m*(kp-2))-(kp-2))=L31(i+1+m*(kp-2));
            k=k+1;
            K31(k, ul(i+1+m*(kp-1))-(kp-1))=L34(i+1+m*(kp-1));
            K31(k, ul(1+m*(kp-1))-(kp-1))=L33(i+1+m*(kp-1));
            k=k+1;
        end
    end
end
K_31=[K_31; K31];

%[3,2]
K_32=[];
for i=0:(m-2)
    K32=sparse(ut(i+1+m*(n-1))-ub(i+1),ur(m*n)-n);
    if ismember(i+1,crack_col) == 0   
        K32(1,ul(i+1))=M34(i+1);
        K32(1,ur(i+1))=M33(i+1);
        K32(1,ul(i+1+m*(n-1))-(n-1))=M32(i+1+m*(n-1));
        K32(1,ur(i+1+m*(n-1))-(n-1))=M31(i+1+m*(n-1));
        for k=2:n
            K32(k, ul(i+1+m*(k-2))-(k-2))=M32(i+1+m*(k-2));
            K32(k, ur(i+1+m*(k-2))-(k-2))=M31(i+1+m*(k-2));
            K32(k, ul(i+1+m*(k-1))-(k-1))=M34(i+1+m*(k-1));
            K32(k, ur(i+1+m*(k-1))-(k-1))=M33(i+1+m*(k-1));
        end
    else
        K32(1,ul(i+1))=M34(i+1);
        K32(1,ur(i+1))=M33(i+1);
        K32(1,ul(i+1+m*(n-1))-(n-1))=M32(i+1+m*(n-1));
        K32(1,ur(i+1+m*(n-1))-(n-1))=M31(i+1+m*(n-1));
        k=2;
        for kp=2:n
            if ismember(i*(n+1)+kp, crack_int_hor) == 0
                K32(k, ul(i+1+m*(kp-2))-(kp-2))=M32(i+1+m*(kp-2));
                K32(k, ur(i+1+m*(kp-2))-(kp-2))=M31(i+1+m*(kp-2));
                K32(k, ul(i+1+m*(kp-1))-(kp-1))=M34(i+1+m*(kp-1));
                K32(k, ur(i+1+m*(kp-1))-(kp-1))=M33(i+1+m*(kp-1));
                k=k+1;
            else
                K32(k, ul(i+1+m*(kp-2))-(kp-2))=M32(i+1+m*(kp-2));
                K32(k, ur(i+1+m*(kp-2))-(kp-2))=M31(i+1+m*(kp-2));
                k=k+1;
                K32(k, ul(i+1+m*(kp-1))-(kp-1))=M34(i+1+m*(kp-1));
                K32(k, ur(i+1+m*(kp-1))-(kp-1))=M33(i+1+m*(kp-1));
                k=k+1;
            end
        end
    end
    K_32=[K_32; K32];
end
i=m-1;
K32=sparse(ut(i+1+m*(n-1))-ub(i+1),ur(m*n)-n);
if ismember(m,crack_col) == 0   
    K32(1,ul(i+1))=M34(i+1);
    K32(1,ul(1))=M33(i+1);
    K32(1,ul(i+1+m*(n-1))-(n-1))=M32(i+1+m*(n-1));
    K32(1,ul(m*(n-1)+1)-(n-1))=M31(i+1+m*(n-1));
    for k=2:n
        K32(k, ul(i+1+m*(k-2))-(k-2))=M32(i+1+m*(k-2));
        K32(k, ul(1+m*(k-2))-(k-2))=M31(i+1+m*(k-2));
        K32(k, ul(i+1+m*(k-1))-(k-1))=M34(i+1+m*(k-1));
        K32(k, ul(1+m*(k-1))-(k-1))=M33(i+1+m*(k-1));
    end
else
    K32(1,ul(i+1))=M34(i+1);
    K32(1,ul(1))=M33(i+1);
    K32(1,ul(i+1+m*(n-1))-(n-1))=M32(i+1+m*(n-1));
    K32(1,ul(1+m*(n-1))-(n-1))=M31(i+1+m*(n-1));
    k=2;
    for kp=2:n
        if ismember(i*(n+1)+kp, crack_int_hor) == 0
            K32(k, ul(i+1+m*(kp-2))-(kp-2))=M32(i+1+m*(kp-2));
            K32(k, ul(1+m*(kp-2))-(kp-2))=M31(i+1+m*(kp-2));
            K32(k, ul(i+1+m*(kp-1))-(kp-1))=M34(i+1+m*(kp-1));
            K32(k, ul(1+m*(kp-1))-(kp-1))=M33(i+1+m*(kp-1));
            k=k+1;
        else
            K32(k, ul(i+1+m*(kp-2))-(kp-2))=M32(i+1+m*(kp-2));
            K32(k, ul(1+m*(kp-2))-(kp-2))=M31(i+1+m*(kp-2));
            k=k+1;
            K32(k, ul(i+1+m*(kp-1))-(kp-1))=M34(i+1+m*(kp-1));
            K32(k, ul(1+m*(kp-1))-(kp-1))=M33(i+1+m*(kp-1));
            k=k+1;
        end
    end
end
K_32=[K_32; K32];

%[3,3]
K_33=[];
for i=0:(m-1)
    K33=sparse(ut(i+1+m*(n-1))-ub(i+1),ut(i+1+m*(n-1))-ub(i+1));
    if ismember(i+1,crack_col) == 0      
        K33(1,1)=N34(i+1)+N31(i+1+m*(n-1));
        K33(1,2)=N33(i+1);
        K33(1,n)=N32(i+1+m*(n-1));
        for k=2:n-1
            K33(k, k-1)=N32(i+1+m*(k-2));
            K33(k, k)=N31(i+1+m*(k-2))+N34(i+1+m*(k-1));
            K33(k, k+1)=N33(i+1+m*(k-1));
        end
        K33(n, n-1)=N32(i+1+m*(n-2));
        K33(n, n)=N31(i+1+m*(n-2))+N34(i+1+m*(n-1));
        K33(n, 1)=N33(i+1+m*(n-1));
    else
        K33(1,1)=N34(i+1)+N31(i+1+m*(n-1));
        K33(1,2)=N33(i+1);
        K33(1,ut(i+1+m*(n-1))-ub(i+1))=N32(i+1+m*(n-1));
        k=2;
        for kp=2:n-1
            if ismember(i*(n+1)+kp, crack_int_hor) == 0
                K33(k, k-1)=N32(i+1+m*(kp-2));
                K33(k, k)=N31(i+1+m*(kp-2))+N34(i+1+m*(kp-1));
                K33(k, k+1)=N33(i+1+m*(kp-1));
                k=k+1;
            else
                K33(k, k-1)=N32(i+1+m*(kp-2));
                K33(k, k)=N31(i+1+m*(kp-2));
                k=k+1;
                K33(k, k)=N34(i+1+m*(kp-1));
                K33(k, k+1)=N33(i+1+m*(kp-1));
                k=k+1;
            end
        end
        if ismember(i*(n+1)+n, crack_int_hor) == 0
            K33(k, k-1)=N32(i+1+m*(n-2));
            K33(k, k)=N31(i+1+m*(n-2))+N34(i+1+m*(n-1));
            K33(k, 1)=N33(i+1+m*(n-1));
        else
            K33(k, k-1)=N32(i+1+m*(n-2));
            K33(k, k)=N31(i+1+m*(n-2));
            k=k+1;
            K33(k, k)=N34(i+1+m*(n-1));
            K33(k, 1)=N33(i+1+m*(n-1));
        end
    end
    [x,y]=size(K33);
    [p,q]=size(K_33);
    K_33=[K_33 sparse(p,y); sparse(x,q) K33];
end

%[3,4]
[p,q]=size(K_33);
K_34=sparse(p,q);


%[4,1]
K_41=[];
for i=0:(m-2)
    K41=sparse(ut(i+1+m*(n-1))-ub(i+1),ur(m*n)-n);
    if ismember(i+1,crack_col) == 0   
        K41(1,ul(i+1))=L44(i+1);
        K41(1,ur(i+1))=L43(i+1);
        K41(1,ul(i+1+m*(n-1))-(n-1))=L42(i+1+m*(n-1));
        K41(1,ur(i+1+m*(n-1))-(n-1))=L41(i+1+m*(n-1));
        for k=2:n
            K41(k, ul(i+1+m*(k-2))-(k-2))=L42(i+1+m*(k-2));
            K41(k, ur(i+1+m*(k-2))-(k-2))=L41(i+1+m*(k-2));
            K41(k, ul(i+1+m*(k-1))-(k-1))=L44(i+1+m*(k-1));
            K41(k, ur(i+1+m*(k-1))-(k-1))=L43(i+1+m*(k-1));
        end
    else
        K41(1,ul(i+1))=L44(i+1);
        K41(1,ur(i+1))=L43(i+1);
        K41(1,ul(i+1+m*(n-1))-(n-1))=L42(i+1+m*(n-1));
        K41(1,ur(i+1+m*(n-1))-(n-1))=L41(i+1+m*(n-1));
        k=2;
        for kp=2:n
            if ismember(i*(n+1)+kp, crack_int_hor) == 0
                K41(k, ul(i+1+m*(kp-2))-(kp-2))=L42(i+1+m*(kp-2));
                K41(k, ur(i+1+m*(kp-2))-(kp-2))=L41(i+1+m*(kp-2));
                K41(k, ul(i+1+m*(kp-1))-(kp-1))=L44(i+1+m*(kp-1));
                K41(k, ur(i+1+m*(kp-1))-(kp-1))=L43(i+1+m*(kp-1));
                k=k+1;
            else
                K41(k, ul(i+1+m*(kp-2))-(kp-2))=L42(i+1+m*(kp-2));
                K41(k, ur(i+1+m*(kp-2))-(kp-2))=L41(i+1+m*(kp-2));
                k=k+1;
                K41(k, ul(i+1+m*(kp-1))-(kp-1))=L44(i+1+m*(kp-1));
                K41(k, ur(i+1+m*(kp-1))-(kp-1))=L43(i+1+m*(kp-1));
                k=k+1;
            end
        end
    end
    K_41=[K_41; K41];
end
i=m-1;
K41=sparse(ut(i+1+m*(n-1))-ub(i+1),ur(m*n)-n);
if ismember(m,crack_col) == 0   
    K41(1,ul(i+1))=L44(i+1);
    K41(1,ul(1))=L43(i+1);
    K41(1,ul(i+1+m*(n-1))-(n-1))=L42(i+1+m*(n-1));
    K41(1,ul(m*(n-1)+1)-(n-1))=L41(i+1+m*(n-1));
    for k=2:n
        K41(k, ul(i+1+m*(k-2))-(k-2))=L42(i+1+m*(k-2));
        K41(k, ul(1+m*(k-2))-(k-2))=L41(i+1+m*(k-2));
        K41(k, ul(i+1+m*(k-1))-(k-1))=L44(i+1+m*(k-1));
        K41(k, ul(1+m*(k-1))-(k-1))=L43(i+1+m*(k-1));
    end
else
    K41(1,ul(i+1))=L44(i+1);
    K41(1,ul(1))=L43(i+1);
    K41(1,ul(i+1+m*(n-1))-(n-1))=L42(i+1+m*(n-1));
    K41(1,ul(1+m*(n-1))-(n-1))=L41(i+1+m*(n-1));
    k=2;
    for kp=2:n
        if ismember(i*(n+1)+kp, crack_int_hor) == 0
            K41(k, ul(i+1+m*(kp-2))-(kp-2))=L42(i+1+m*(kp-2));
            K41(k, ul(1+m*(kp-2))-(kp-2))=L41(i+1+m*(kp-2));
            K41(k, ul(i+1+m*(kp-1))-(kp-1))=L44(i+1+m*(kp-1));
            K41(k, ul(1+m*(kp-1))-(kp-1))=L43(i+1+m*(kp-1));
            k=k+1;
        else
            K41(k, ul(i+1+m*(kp-2))-(kp-2))=L42(i+1+m*(kp-2));
            K41(k, ul(1+m*(kp-2))-(kp-2))=L41(i+1+m*(kp-2));
            k=k+1;
            K41(k, ul(i+1+m*(kp-1))-(kp-1))=L44(i+1+m*(kp-1));
            K41(k, ul(1+m*(kp-1))-(kp-1))=L43(i+1+m*(kp-1));
            k=k+1;
        end
    end
end
K_41=[K_41; K41];

%[4,2]
K_42=[];
for i=0:(m-2)
    K42=sparse(ut(i+1+m*(n-1))-ub(i+1),ur(m*n)-n);
    if ismember(i+1,crack_col) == 0   
        K42(1,ul(i+1))=M44(i+1);
        K42(1,ur(i+1))=M43(i+1);
        K42(1,ul(i+1+m*(n-1))-(n-1))=M42(i+1+m*(n-1));
        K42(1,ur(i+1+m*(n-1))-(n-1))=M41(i+1+m*(n-1));
        for k=2:n
            K42(k, ul(i+1+m*(k-2))-(k-2))=M42(i+1+m*(k-2));
            K42(k, ur(i+1+m*(k-2))-(k-2))=M41(i+1+m*(k-2));
            K42(k, ul(i+1+m*(k-1))-(k-1))=M44(i+1+m*(k-1));
            K42(k, ur(i+1+m*(k-1))-(k-1))=M43(i+1+m*(k-1));
        end
    else
        K42(1,ul(i+1))=M44(i+1);
        K42(1,ur(i+1))=M43(i+1);
        K42(1,ul(i+1+m*(n-1))-(n-1))=M42(i+1+m*(n-1));
        K42(1,ur(i+1+m*(n-1))-(n-1))=M41(i+1+m*(n-1));
        k=2;
        for kp=2:n
            if ismember(i*(n+1)+kp, crack_int_hor) == 0
                K42(k, ul(i+1+m*(kp-2))-(kp-2))=M42(i+1+m*(kp-2));
                K42(k, ur(i+1+m*(kp-2))-(kp-2))=M41(i+1+m*(kp-2));
                K42(k, ul(i+1+m*(kp-1))-(kp-1))=M44(i+1+m*(kp-1));
                K42(k, ur(i+1+m*(kp-1))-(kp-1))=M43(i+1+m*(kp-1));
                k=k+1;
            else
                K42(k, ul(i+1+m*(kp-2))-(kp-2))=M42(i+1+m*(kp-2));
                K42(k, ur(i+1+m*(kp-2))-(kp-2))=M41(i+1+m*(kp-2));
                k=k+1;
                K42(k, ul(i+1+m*(kp-1))-(kp-1))=M44(i+1+m*(kp-1));
                K42(k, ur(i+1+m*(kp-1))-(kp-1))=M43(i+1+m*(kp-1));
                k=k+1;
            end
        end
    end
    K_42=[K_42; K42];
end
i=m-1;
K42=sparse(ut(i+1+m*(n-1))-ub(i+1),ur(m*n)-n);
if ismember(m,crack_col) == 0   
    K42(1,ul(i+1))=M44(i+1);
    K42(1,ul(1))=M43(i+1);
    K42(1,ul(i+1+m*(n-1))-(n-1))=M42(i+1+m*(n-1));
    K42(1,ul(m*(n-1)+1)-(n-1))=M41(i+1+m*(n-1));
    for k=2:n
        K42(k, ul(i+1+m*(k-2))-(k-2))=M42(i+1+m*(k-2));
        K42(k, ul(1+m*(k-2))-(k-2))=M41(i+1+m*(k-2));
        K42(k, ul(i+1+m*(k-1))-(k-1))=M44(i+1+m*(k-1));
        K42(k, ul(1+m*(k-1))-(k-1))=M43(i+1+m*(k-1));
    end
else
    K42(1,ul(i+1))=M44(i+1);
    K42(1,ul(1))=M43(i+1);
    K42(1,ul(i+1+m*(n-1))-(n-1))=M42(i+1+m*(n-1));
    K42(1,ul(1+m*(n-1))-(n-1))=M41(i+1+m*(n-1));
    k=2;
    for kp=2:n
        if ismember(i*(n+1)+kp, crack_int_hor) == 0
            K42(k, ul(i+1+m*(kp-2))-(kp-2))=M42(i+1+m*(kp-2));
            K42(k, ul(1+m*(kp-2))-(kp-2))=M41(i+1+m*(kp-2));
            K42(k, ul(i+1+m*(kp-1))-(kp-1))=M44(i+1+m*(kp-1));
            K42(k, ul(1+m*(kp-1))-(kp-1))=M43(i+1+m*(kp-1));
            k=k+1;
        else
            K42(k, ul(i+1+m*(kp-2))-(kp-2))=M42(i+1+m*(kp-2));
            K42(k, ul(1+m*(kp-2))-(kp-2))=M41(i+1+m*(kp-2));
            k=k+1;
            K42(k, ul(i+1+m*(kp-1))-(kp-1))=M44(i+1+m*(kp-1));
            K42(k, ul(1+m*(kp-1))-(kp-1))=M43(i+1+m*(kp-1));
            k=k+1;
        end
    end
end
K_42=[K_42; K42];

%[4,3]
[p,q]=size(K_34);
K_43=sparse(p,q);

%[4,4]
K_44=[];
for i=0:(m-1)
    K44=sparse(ut(i+1+m*(n-1))-ub(i+1),ut(i+1+m*(n-1))-ub(i+1));
    if ismember(i+1,crack_col) == 0      
        K44(1,1)=O44(i+1)+O41(i+1+m*(n-1));
        K44(1,2)=O43(i+1);
        K44(1,n)=O42(i+1+m*(n-1));
        for k=2:n-1
            K44(k, k-1)=O42(i+1+m*(k-2));
            K44(k, k)=O41(i+1+m*(k-2))+O44(i+1+m*(k-1));
            K44(k, k+1)=O43(i+1+m*(k-1));
        end
        K44(n, n-1)=O42(i+1+m*(n-2));
        K44(n, n)=O41(i+1+m*(n-2))+O44(i+1+m*(n-1));
        K44(n, 1)=O43(i+1+m*(n-1));
    else
        K44(1,1)=O44(i+1)+O41(i+1+m*(n-1));
        K44(1,2)=O43(i+1);
        K44(1,ut(i+1+m*(n-1))-ub(i+1))=O42(i+1+m*(n-1));
        k=2;
        for kp=2:n-1
            if ismember(i*(n+1)+kp, crack_int_hor) == 0
                K44(k, k-1)=O42(i+1+m*(kp-2));
                K44(k, k)=O41(i+1+m*(kp-2))+O44(i+1+m*(kp-1));
                K44(k, k+1)=O43(i+1+m*(kp-1));
                k=k+1;
            else
                K44(k, k-1)=O42(i+1+m*(kp-2));
                K44(k, k)=O41(i+1+m*(kp-2));
                k=k+1;
                K44(k, k)=O44(i+1+m*(kp-1));
                K44(k, k+1)=O43(i+1+m*(kp-1));
                k=k+1;
            end
        end
        if ismember(i*(n+1)+n, crack_int_hor) == 0
            K44(k, k-1)=O42(i+1+m*(n-2));
            K44(k, k)=O41(i+1+m*(n-2))+O44(i+1+m*(n-1));
            K44(k, 1)=O43(i+1+m*(n-1));
        else
            K44(k, k-1)=O42(i+1+m*(n-2));
            K44(k, k)=O41(i+1+m*(n-2));
            k=k+1;
            K44(k, k)=O44(i+1+m*(n-1));
            K44(k, 1)=O43(i+1+m*(n-1));
        end    
    end
    [x,y]=size(K44);
    [p,q]=size(K_44);
    K_44=[K_44 sparse(p,y); sparse(x,q) K44];
end
% %_________________________________________________________________________________________________________________


% Out-pf-Plane Shear Stiffness Matrix
% %_________________________________________________________________________________________________________________
%[1,1]
Knext_11=[];
for i=0:(n-1)
    t=i*m;
    Knext11=sparse(ur(t+m)-ul(t+1),ur(t+m)-ul(t+1));
    if ismember(i+1,crack_row) == 0
        Knext11(1,1)=U14(t+1)+U11(t+m);
        Knext11(1,2)=U13(t+1);
        Knext11(1,m)=U12(t+m);
        for k=2:m-1
            Knext11(k, k-1)=U12(t+k-1);
            Knext11(k, k)=U11(t+k-1)+U14(t+k);
            Knext11(k, k+1)=U13(t+k);
        end
        Knext11(m,m-1)=U12(t+m-1);
        Knext11(m,m)=U11(t+m-1)+U14(t+m);
        Knext11(m,1)=U13(t+m);
    else
        Knext11(1,1)=U14(t+1)+U11(t+m);
        Knext11(1,2)=U13(t+1);
        Knext11(1,ur(t+m)-ul(t+1))=U12(t+m);
        k=2;
        for kp=2:m-1
            if ismember(i*(m+1)+kp, crack_int_ver) == 0
                Knext11(k, k-1)=U12(t+kp-1);
                Knext11(k, k)=U11(t+kp-1)+U14(t+kp);
                Knext11(k, k+1)=U13(t+kp);
                k=k+1;
            else
                Knext11(k, k-1)=U12(t+kp-1);
                Knext11(k, k)=U11(t+kp-1);
                k=k+1;
                Knext11(k, k)=U14(t+kp);
                Knext11(k, k+1)=U13(t+kp);
                k=k+1;
            end
        end
        if ismember(i*(m+1)+m, crack_int_ver) == 0
            Knext11(k,k-1)=U12(t+m-1);
            Knext11(k,k)=U11(t+m-1)+U14(t+m);
            Knext11(k,1)=U13(t+m);
        else
            Knext11(k,k-1)=U12(t+m-1);
            Knext11(k,k)=U11(t+m-1);
            k=k+1;
            Knext11(k,k)=U14(t+m);
            Knext11(k,1)=U13(t+m);
        end        
    end
    [x,y]=size(Knext11);
    [p,q]=size(Knext_11);
    Knext_11=[Knext_11 sparse(p,y); sparse(x,q) Knext11];
end

%[1,2]
Knext_12=[];
for i=0:(n-2)
    t=i*m;
    Knext12=sparse(ur(m*(i+1))-ul(t+1),ut(m*n)-m);
    if ismember(i+1,crack_row) == 0   
        Knext12(1,ub(t+1))=V14(t+1);
        Knext12(1,ut(t+1))=V13(t+1);
        Knext12(1,ub(t+m)-(m-1))=V12(t+m);
        Knext12(1,ut(t+m)-(m-1))=V11(t+m);        
        for k=2:m
            Knext12(k, ub(t+k-1)-(k-2))=V12(t+k-1);
            Knext12(k, ut(t+k-1)-(k-2))=V11(t+k-1);
            Knext12(k, ub(t+k)-(k-1))=V14(t+k);
            Knext12(k, ut(t+k)-(k-1))=V13(t+k);
        end
    else
        Knext12(1,ub(t+1))=V14(t+1);
        Knext12(1,ut(t+1))=V13(t+1);
        Knext12(1,ub(t+m)-(m-1))=V12(t+m);
        Knext12(1,ut(t+m)-(m-1))=V11(t+m);
        k=2;
        for kp=2:m
            if ismember(i*(m+1)+kp, crack_int_ver) == 0
                Knext12(k, ub(t+kp-1)-(kp-2))=V12(t+kp-1);
                Knext12(k, ut(t+kp-1)-(kp-2))=V11(t+kp-1);
                Knext12(k, ub(t+kp)-(kp-1))=V14(t+kp);
                Knext12(k, ut(t+kp)-(kp-1))=V13(t+kp);
                k=k+1;
            else
                Knext12(k, ub(t+kp-1)-(kp-2))=V12(t+kp-1);
                Knext12(k, ut(t+kp-1)-(kp-2))=V11(t+kp-1);
                k=k+1;
                Knext12(k, ub(t+kp)-(kp-1))=V14(t+kp);
                Knext12(k, ut(t+kp)-(kp-1))=V13(t+kp);
                k=k+1;
            end
        end
    end
    Knext_12=[Knext_12; Knext12];
end
i=n-1;
t=i*m;
Knext12=sparse(ur(m*n)-ul(t+1),ut(m*n)-m);
if ismember(n,crack_row) == 0   
    Knext12(1,ub(t+1))=V14(t+1);
    Knext12(1,1)=V13(t+1);
    Knext12(1,ub(t+m)-(m-1))=V12(t+m);
    Knext12(1,ub(m)-(m-1))=V11(t+m);        
    for k=2:m
        Knext12(k, ub(t+k-1)-(k-2))=V12(t+k-1);
        Knext12(k, ub(k-1)-(k-2))=V11(t+k-1);
        Knext12(k, ub(t+k)-(k-1))=V14(t+k);
        Knext12(k, ub(k)-(k-1))=V13(t+k);
    end
else
    Knext12(1,ub(t+1))=V14(t+1);
    Knext12(1,1)=V13(t+1);
    Knext12(1,ub(t+m)-(m-1))=V12(t+m);
    Knext12(1,ub(m)-(m-1))=V11(t+m);        
    k=2;
    for kp=2:m
        if ismember(i*(m+1)+kp, crack_int_ver) == 0
            Knext12(k, ub(t+kp-1)-(kp-2))=V12(t+kp-1);
            Knext12(k, ub(kp-1)-(kp-2))=V11(t+kp-1);
            Knext12(k, ub(t+kp)-(kp-1))=V14(t+kp);
            Knext12(k, ub(kp)-(kp-1))=V13(t+kp);
            k=k+1;
        else
            Knext12(k, ub(t+kp-1)-(kp-2))=V12(t+kp-1);
            Knext12(k, ub(kp-1)-(kp-2))=V11(t+kp-1);
            k=k+1;
            Knext12(k, ub(t+kp)-(kp-1))=V14(t+kp);
            Knext12(k, ub(kp)-(kp-1))=V13(t+kp);
            k=k+1;
        end
    end
end
Knext_12=[Knext_12; Knext12];


%[2,1]
Knext_21=[];
for i=0:(m-2)
    Knext21=sparse(ut(i+1+m*(n-1))-ub(i+1),ur(m*n)-n);
    if ismember(i+1,crack_col) == 0   
        Knext21(1,ul(i+1))=U24(i+1);
        Knext21(1,ur(i+1))=U23(i+1);
        Knext21(1,ul(i+1+m*(n-1))-(n-1))=U22(i+1+m*(n-1));
        Knext21(1,ur(i+1+m*(n-1))-(n-1))=U21(i+1+m*(n-1));
        for k=2:n
            Knext21(k, ul(i+1+m*(k-2))-(k-2))=U22(i+1+m*(k-2));
            Knext21(k, ur(i+1+m*(k-2))-(k-2))=U21(i+1+m*(k-2));
            Knext21(k, ul(i+1+m*(k-1))-(k-1))=U24(i+1+m*(k-1));
            Knext21(k, ur(i+1+m*(k-1))-(k-1))=U23(i+1+m*(k-1));
        end
    else
        Knext21(1,ul(i+1))=U24(i+1);
        Knext21(1,ur(i+1))=U23(i+1);
        Knext21(1,ul(i+1+m*(n-1))-(n-1))=U22(i+1+m*(n-1));
        Knext21(1,ur(i+1+m*(n-1))-(n-1))=U21(i+1+m*(n-1));
        k=2;
        for kp=2:n
            if ismember(i*(n+1)+kp, crack_int_hor) == 0
                Knext21(k, ul(i+1+m*(kp-2))-(kp-2))=U22(i+1+m*(kp-2));
                Knext21(k, ur(i+1+m*(kp-2))-(kp-2))=U21(i+1+m*(kp-2));
                Knext21(k, ul(i+1+m*(kp-1))-(kp-1))=U24(i+1+m*(kp-1));
                Knext21(k, ur(i+1+m*(kp-1))-(kp-1))=U23(i+1+m*(kp-1));
                k=k+1;
            else
                Knext21(k, ul(i+1+m*(kp-2))-(kp-2))=U22(i+1+m*(kp-2));
                Knext21(k, ur(i+1+m*(kp-2))-(kp-2))=U21(i+1+m*(kp-2));
                k=k+1;
                Knext21(k, ul(i+1+m*(kp-1))-(kp-1))=U24(i+1+m*(kp-1));
                Knext21(k, ur(i+1+m*(kp-1))-(kp-1))=U23(i+1+m*(kp-1));
                k=k+1;
            end
        end
    end
    Knext_21=[Knext_21; Knext21];
end
i=m-1;
Knext21=sparse(ut(i+1+m*(n-1))-ub(i+1),ur(m*n)-n);
if ismember(m,crack_col) == 0   
    Knext21(1,ul(i+1))=U24(i+1);
    Knext21(1,ul(1))=U23(i+1);
    Knext21(1,ul(i+1+m*(n-1))-(n-1))=U22(i+1+m*(n-1));
    Knext21(1,ul(m*(n-1)+1)-(n-1))=U21(i+1+m*(n-1));
    for k=2:n
        Knext21(k, ul(i+1+m*(k-2))-(k-2))=U22(i+1+m*(k-2));
        Knext21(k, ul(1+m*(k-2))-(k-2))=U21(i+1+m*(k-2));
        Knext21(k, ul(i+1+m*(k-1))-(k-1))=U24(i+1+m*(k-1));
        Knext21(k, ul(1+m*(k-1))-(k-1))=U23(i+1+m*(k-1));
    end
else
    Knext21(1,ul(i+1))=U24(i+1);
    Knext21(1,ul(1))=U23(i+1);
    Knext21(1,ul(i+1+m*(n-1))-(n-1))=U22(i+1+m*(n-1));
    Knext21(1,ul(1+m*(n-1))-(n-1))=U21(i+1+m*(n-1));
    k=2;
    for kp=2:n
        if ismember(i*(n+1)+kp, crack_int_hor) == 0
            Knext21(k, ul(i+1+m*(kp-2))-(kp-2))=U22(i+1+m*(kp-2));
            Knext21(k, ul(1+m*(kp-2))-(kp-2))=U21(i+1+m*(kp-2));
            Knext21(k, ul(i+1+m*(kp-1))-(kp-1))=U24(i+1+m*(kp-1));
            Knext21(k, ul(1+m*(kp-1))-(kp-1))=U23(i+1+m*(kp-1));
            k=k+1;
        else
            Knext21(k, ul(i+1+m*(kp-2))-(kp-2))=U22(i+1+m*(kp-2));
            Knext21(k, ul(1+m*(kp-2))-(kp-2))=U21(i+1+m*(kp-2));
            k=k+1;
            Knext21(k, ul(i+1+m*(kp-1))-(kp-1))=U24(i+1+m*(kp-1));
            Knext21(k, ul(1+m*(kp-1))-(kp-1))=U23(i+1+m*(kp-1));
            k=k+1;
        end
    end
end
Knext_21=[Knext_21; Knext21];

%[2,2]
Knext_22=[];
for i=0:(m-1)
    Knext22=sparse(ut(i+1+m*(n-1))-ub(i+1),ut(i+1+m*(n-1))-ub(i+1));
    if ismember(i+1,crack_col) == 0      
        Knext22(1,1)=V24(i+1)+V21(i+1+m*(n-1));
        Knext22(1,2)=V23(i+1);
        Knext22(1,n)=V22(i+1+m*(n-1));
        for k=2:n-1
            Knext22(k, k-1)=V22(i+1+m*(k-2));
            Knext22(k, k)=V21(i+1+m*(k-2))+V24(i+1+m*(k-1));
            Knext22(k, k+1)=V23(i+1+m*(k-1));
        end
        Knext22(n, n-1)=V22(i+1+m*(n-2));
        Knext22(n, n)=V21(i+1+m*(n-2))+V24(i+1+m*(n-1));
        Knext22(n, 1)=V23(i+1+m*(n-1));
    else
        Knext22(1,1)=V24(i+1)+V21(i+1+m*(n-1));
        Knext22(1,2)=V23(i+1);
        Knext22(1,ut(i+1+m*(n-1))-ub(i+1))=V22(i+1+m*(n-1));
        k=2;
        for kp=2:n-1
            if ismember(i*(n+1)+kp, crack_int_hor) == 0
                Knext22(k, k-1)=V22(i+1+m*(kp-2));
                Knext22(k, k)=V21(i+1+m*(kp-2))+V24(i+1+m*(kp-1));
                Knext22(k, k+1)=V23(i+1+m*(kp-1));
                k=k+1;
            else
                Knext22(k, k-1)=V22(i+1+m*(kp-2));
                Knext22(k, k)=V21(i+1+m*(kp-2));
                k=k+1;
                Knext22(k, k)=V24(i+1+m*(kp-1));
                Knext22(k, k+1)=V23(i+1+m*(kp-1));
                k=k+1;
            end
        end
        if ismember(i*(n+1)+n, crack_int_hor) == 0
            Knext22(k, k-1)=V22(i+1+m*(n-2));
            Knext22(k, k)=V21(i+1+m*(n-2))+V24(i+1+m*(n-1));
            Knext22(k, 1)=V23(i+1+m*(n-1));
        else
            Knext22(k, k-1)=V22(i+1+m*(n-2));
            Knext22(k, k)=V21(i+1+m*(n-2));
            k=k+1;
            Knext22(k, k)=V24(i+1+m*(n-1));
            Knext22(k, 1)=V23(i+1+m*(n-1));
        end
    end
    [x,y]=size(Knext22);
    [p,q]=size(Knext_22);
    Knext_22=[Knext_22 sparse(p,y); sparse(x,q) Knext22];
end


K = [K_11, K_12, K_13, K_14; K_21, K_22, K_23, K_24; K_31, K_32, K_33, K_34; K_41, K_42, K_43, K_44];
Knext=[Knext_11,Knext_12; Knext_21,Knext_22];
% %_________________________________________________________________________________________________________________
