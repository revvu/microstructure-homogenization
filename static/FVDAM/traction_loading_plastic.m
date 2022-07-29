function [av_macro_strains, av_macro_stresses] = traction_loading_plastic(LOP, load_current, Cstar, sig_inel, strain_thermal);

% strain=zeros(6,1);
stress=zeros(6,1);
strain=strain_thermal;

stress_r=stress;
Cstar_r=Cstar;
sig_inel_r=sig_inel;


if ismember(LOP, [1:6]) == 1

    strain(LOP)=strain(LOP)+load_current;

    Cstar_r(LOP,:)=[];   % Row reduction of global stiffness matrix corresponding to the known displacement boundary conditions
    stress_r(LOP,:)=[];  % Removing the unknown tractions at the boundaries, direct consequence of the above
    sig_inel_r(LOP,:)=[];
    stress_r=stress_r+sig_inel_r-(Cstar_r(:,LOP)*strain(LOP));
    Cstar_r(:,LOP)=[];   % Column reduction of global stiffness matrix corresponding to the known displacement boundary conditions
    strain_r=Cstar_r\stress_r;
    k=1;
    for j=1:6
        if j ~= LOP
            strain(j)=strain_r(k);
            k=k+1;
        end
    end


elseif LOP == 7

    strain(1)=0;
    strain(2)=load_current;
    strain(3)=-load_current;

    temp=[1 2 3];
    dim=length(temp);

    Cstar_r(temp,:)=[];
    stress_r(temp,:)=[];
    sig_inel_r(temp,:)=[];
    stress_r=stress_r+sig_inel_r;
    for i=1:dim
        stress_r=stress_r-(Cstar_r(:,temp(i))*strain(temp(i)));
    end
    Cstar_r(:,temp)=[];
    strain_r=Cstar_r\stress_r;

    j=1;
    k=1;
    for i=1:6
        if i ~= temp(j)
            strain(i)=strain_r(k);
            k=k+1;
        elseif j ~= dim
            j=j+1;
        end
    end


elseif LOP == 8
    strain=Cstar\(stress+sig_inel);


elseif LOP == 9

    strain(1)=0;
    strain(2)=load_current;

    temp=[1 2];
    dim=length(temp);

    Cstar_r(temp,:)=[];
    stress_r(temp,:)=[];
    sig_inel_r(temp,:)=[];
    stress_r=stress_r+sig_inel_r;
    for i=1:dim
        stress_r=stress_r-(Cstar_r(:,temp(i))*strain(temp(i)));
    end
    Cstar_r(:,temp)=[];
    strain_r=Cstar_r\stress_r;

    j=1;
    k=1;
    for i=1:6
        if i ~= temp(j)
            strain(i)=strain_r(k);
            k=k+1;
        elseif j ~= dim
            j=j+1;
        end
    end


elseif LOP == 10

    strain(1)=0;
    strain(2)=load_current;
    strain(3)=load_current;

    temp=[1 2 3];
    dim=length(temp);

    Cstar_r(temp,:)=[];
    stress_r(temp,:)=[];
    sig_inel_r(temp,:)=[];
    stress_r=stress_r+sig_inel_r;
    for i=1:dim
        stress_r=stress_r-(Cstar_r(:,temp(i))*strain(temp(i)));
    end
    Cstar_r(:,temp)=[];
    strain_r=Cstar_r\stress_r;

    j=1;
    k=1;
    for i=1:6
        if i ~= temp(j)
            strain(i)=strain_r(k);
            k=k+1;
        elseif j ~= dim
            j=j+1;
        end
    end


elseif LOP == 11

    strain(2)=load_current;
    strain(3)=load_current;

    temp=[2 3];
    dim=length(temp);

    Cstar_r(temp,:)=[];
    stress_r(temp,:)=[];
    sig_inel_r(temp,:)=[];
    stress_r=stress_r+sig_inel_r;
    for i=1:dim
        stress_r=stress_r-(Cstar_r(:,temp(i))*strain(temp(i)));
    end
    Cstar_r(:,temp)=[];
    strain_r=Cstar_r\stress_r;

    j=1;
    k=1;
    for i=1:6
        if i ~= temp(j)
            strain(i)=strain_r(k);
            k=k+1;
        elseif j ~= dim
            j=j+1;
        end
    end

elseif LOP == 12

    strain(1)=0;
    strain(2)=0;
    strain(3)=load_current;

    temp=[1 2 3];
    dim=length(temp);

    Cstar_r(temp,:)=[];
    stress_r(temp,:)=[];
    sig_inel_r(temp,:)=[];
    stress_r=stress_r+sig_inel_r;
    for i=1:dim
        stress_r=stress_r-(Cstar_r(:,temp(i))*strain(temp(i)));
    end
    Cstar_r(:,temp)=[];
    strain_r=Cstar_r\stress_r;

    j=1;
    k=1;
    for i=1:6
        if i ~= temp(j)
            strain(i)=strain_r(k);
            k=k+1;
        elseif j ~= dim
            j=j+1;
        end
    end

elseif LOP == 13

    strain(1)=0;
    strain(3)=load_current;

    temp=[1 3];
    dim=length(temp);

    Cstar_r(temp,:)=[];
    stress_r(temp,:)=[];
    sig_inel_r(temp,:)=[];
    stress_r=stress_r+sig_inel_r;
    for i=1:dim
        stress_r=stress_r-(Cstar_r(:,temp(i))*strain(temp(i)));
    end
    Cstar_r(:,temp)=[];
    strain_r=Cstar_r\stress_r;

    j=1;
    k=1;
    for i=1:6
        if i ~= temp(j)
            strain(i)=strain_r(k);
            k=k+1;
        elseif j ~= dim
            j=j+1;
        end
    end

end

av_macro_strains=strain;
av_macro_stresses=Cstar*strain - sig_inel;
