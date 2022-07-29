% Takes the inverse of the reduced global stiffness matrix and calculates the surface averaged displacements
% function [Avu]=inverse_next(L, U, P, perm, Sigmat, Avu, temp)
function [Avu]=inverse_next(L, U, P, Q, Sigmat, Avu, temp)

% spparms('spumoni',2);
% z = P*Sigmat;
% y = L\z; % lower triangular solve
% x = U\y; % upper triangular solve
% B(perm)=x;
% spparms('spumoni',0);

y = L\(P*Sigmat);    % lower triangular solve
x = U\y;       % upper triangular solve
% B(perm)=x;
B=Q*x;

dim=length(Avu);
dimnext=length(temp);

j=1;
k=1;
for i=1:dim
    if i ~= temp(j)
        Avu(i)=B(k);
        k=k+1;
    elseif j ~= dimnext
        j=j+1;
    end
end
