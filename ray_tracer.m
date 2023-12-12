for li=1:N_ray% N_ray is the total number of rays
    if AngD==1
        tet(li)=asind(rand()); %Initial polar angle of the ray with a Lamberian distribution
        fi=360*rand(); %Initial azimuthal angle
    elseif AngD==2
        tet(li)=rand(); %Initial polar angle of the ray with a uniform distribution
        fi=360*rand(); %Initial azimuthal angle
    else
        if li==1 && ff==1 && gg==1
            prompt={'Polar angle "teta" in degree:','Azimuthal angle "phi" in degree:'};
            def={'0','0'};
            d_title='Direction of radiation';
            n_lines=1;
            N=inputdlg(prompt,d_title,n_lines,def);
        end
        tet(li)=str2double(N(1)); % The polar angle of ray's direction
        fi=str2double(N(2)); % The azimuthal angle of ray's direction
    end
    jj=Start_plate; % Determining the source 
    col=[.1,.1,1]; %blue 
    n1=n01; n2=n02;
    if jj==10 
        n0=-cosd(tet(li)); % 3rd direction cosine for 'apex-to-bas' direction
    else 
        n0=cosd(tet(li)); % 3rd direction cosine for 'base-to-apex' direction
    end
    l0=sind(tet(li))*cosd(fi); %first direction cosine
    m0=sind(tet(li))*sind(fi);% 2nd direction cosine
    di=[l0;m0;n0]; % Unit vector of beam direction
    r00(1:3,li)=[(2*rand()-1);(2*rand()-1);+1/L(jj,3)]; %gives a random point on the source surface
    r0=r00(1:3,li);      r=r0;    tt=0;
    while r(3)<=1/L(10,3) && r(3)>=1/L(11,3)% condition for remaining between two surfaces
        r0=r;% sets the last position of the beam as new start point
        
        % Determining the distances from the current posiotion to each surface
        ti=([(L(1,4)-L(1,1)*r(1)-L(1,2)*r(2)-L(1,3)*r(3))/(di(1)*L(1,1)+di(2)*L(1,2)+di(3)*L(1,3)),
       (L(2,4)-L(2,1)*r(1)-L(2,2)*r(2)-L(2,3)*r(3))/(di(1)*L(2,1)+di(2)*L(2,2)+di(3)*L(2,3)),
       (L(3,4)-L(3,1)*r(1)-L(3,2)*r(2)-L(3,3)*r(3))/(di(1)*L(3,1)+di(2)*L(3,2)+di(3)*L(3,3)),
       (L(4,4)-L(4,1)*r(1)-L(4,2)*r(2)-L(4,3)*r(3))/(di(1)*L(4,1)+di(2)*L(4,2)+di(3)*L(4,3)),
       (L(5,4)-L(5,1)*r(1)-L(5,2)*r(2)-L(5,3)*r(3))/(di(1)*L(5,1)+di(2)*L(5,2)+di(3)*L(5,3)),
       (L(6,4)-L(6,1)*r(1)-L(6,2)*r(2)-L(6,3)*r(3))/(di(1)*L(6,1)+di(2)*L(6,2)+di(3)*L(6,3)),
       (L(7,4)-L(7,1)*r(1)-L(7,2)*r(2)-L(7,3)*r(3))/(di(1)*L(7,1)+di(2)*L(7,2)+di(3)*L(7,3)),
       (L(8,4)-L(8,1)*r(1)-L(8,2)*r(2)-L(8,3)*r(3))/(di(1)*L(8,1)+di(2)*L(8,2)+di(3)*L(8,3)),
       (L(9,4)-L(9,1)*r(1)-L(9,2)*r(2)-L(9,3)*r(3))/(di(1)*L(9,1)+di(2)*L(9,2)+di(3)*L(9,3)),
       (L(10,4)-L(10,1)*r(1)-L(10,2)*r(2)-L(10,3)*r(3))/(di(1)*L(10,1)+di(2)*L(10,2)+di(3)*L(10,3)),
       (L(11,4)-L(11,1)*r(1)-L(11,2)*r(2)-L(11,3)*r(3))/(di(1)*L(11,1)+di(2)*L(11,2)+di(3)*L(11,3)),]);
        
        ti(jj)=NaN;
        tt=tt+1; % tt is the number of surfaces that one ray encounters during its trip.
        for j=1:11% Selecting the correct sufaces
            r=r0+di*ti(j);
            if (sum((r-r0).*di)<=0 || norm(di*ti(j))<=1*eps || ti(j)==-inf)
                ti(j)=NaN;
            end
        end
        [m,jj]=(min(ti(1:N_surfaces)));%Finding the nearest surface
        r=r0+di*ti(jj);% Incident point on new surfsce
        Ns=sqrt(sum(L(jj,1:3).*L(jj,1:3)));
        % Is the incident point inside the piramid?
          A1=sum(L(1,1:3).*r')<=1+eps;
          A2=sum(L(2,1:3).*r')<=1+eps;
          A3=sum(L(3,1:3).*r')<=1+eps;
          A4=sum(L(4,1:3).*r')<=1+eps;
          A5=sum(L(5,1:3).*r')<=1+eps;
          A6=sum(L(6,1:3).*r')<=1+eps;
          A7=sum(L(7,1:3).*r')<=1+eps;
          A8=sum(L(8,1:3).*r')<=1+eps;
          A9=sum(L(9,1:3).*r')<=1+eps;
          A_piramid=[A1,A2,A3,A4,A5,A6,A7,A8,A9];

        if simulation_key==1 && rem(li,1000)==0.0 %Drawing the ray path
            line([r0(1),r(1)],[r0(2),r(2)],[r0(3),r(3)],'color',col); pause(0.05) %
        end
            switch jj%
                case {6}% Enclosure walls on x axis
                   %di(1)=-di(1);%
                   r(1)=-r(1); jj=7;     %   col=[.1,.1,.1]; %black 
                case {7}% Enclosure walls on x axis
                   %di(1)=-di(1);%
                   r(1)=-r(1); jj=6;                   
                case {8}%Enclosure walls on y axis
                   %di(2)=-di(2);%    
                   r(2)=-r(2); jj=9;%%    col=[.1,.1,.1]; %black 
                case {9}%Enclosure walls on y axis
                   %di(2)=-di(2);%    
                   r(2)=-r(2); jj=8;%di(2)=-di(2);%    %    col=[.1,.1,.1]; %black 
                case {10} %Atop source
                    if rand<=reflect
                        di(3)=-cosd(asind(rand));
                    else
                        r10=r10+1;
                        break;
                    end
                case {11} %Beneath source
                    if rand<=reflect
                        di(3)=cosd(asind(rand));
                    else
                        r11=r11+1;
                        break;
                    end
                case {1,2,3,4,5}% lateral surfaces of piramid
                    A_piramid(jj)=1;
                    if all(A_piramid)% If incident point is on piramid's surfaces
                        CI1=sum(di'.*L(jj,1:3))./Ns;%Cosine(teta) of the incident angle
                        SI1=sqrt(1-CI1^2);% Sine(teta)
                        SI2=n1/n2*SI1;% Sine(teta') of the refracted ray (Snell's law)
                        CI2=sign(CI1)*sqrt(1-SI2^2);% Cosine(teta')
                        if SI2<eps% If the incident is almost normal
                            n=n2/n1;
                            ro_p=((abs(CI1)-sqrt(n^2-SI1^2))/(abs(CI1)+sqrt(n^2-SI1^2)))^2;
                            ro_n=((n^2*abs(CI1)-sqrt(n^2-SI1.^2))/(n^2*abs(CI1)+sqrt(n^2-SI1^2)))^2;
                            ro=0.5*(ro_n+ro_p); %Fresnel's reflection index
                            k=rand;
                            if k<1*ro % Checking the Fresnel's reflection condition
                               col=[.9,.5,.4]; %red
                               di=-di;% Reflection
                            else
                                nn=n1;    n1=n2;  n2=nn;%Entering to other medium
                            end
                        elseif abs(SI2)>=1%The condition for total internal reflection
                            C=[-CI1; 2*SI1^2-1; 0];%Reflection angle vector
                            T=tra(L(jj,1:3)/Ns,di);%Coeficient matrix
                            if det(T)~=0
                                di=T\C; %Direction of the reflected ray
                            end
                            col=[.5,.8,.1]; %Green
                        else
                            n=n2/n1;%
                            ro_p=((abs(CI1)-sqrt(n^2-SI1^2))/(abs(CI1)+sqrt(n^2-SI1^2)))^2;
                            ro_n=((n^2*abs(CI1)-sqrt(n^2-SI1.^2))/(n^2*abs(CI1)+sqrt(n^2-SI1^2)))^2;
                            ro=0.5*(ro_n+ro_p);
                            k=rand;
                            if k<1*ro %Checking the Fresnel's reflection condition
                                C=[-CI1; 2*SI1^2-1; 0];%Reflection angle vector
                                col=[.9,.5,.4]; %red
                            else
                                C=[CI2; (CI1*CI2+SI1*SI2); 0];% Refraction angle vector
                                nn=n1;    n1=n2;  n2=nn;% Changing the mediums
                            end
                            T=tra(L(jj,1:3)/Ns,di);%Coeficient matrix
                            if det(T)~=0
                                di=T\C;% New direction
                            end
                        end
                    end
            end
        if tt>100 %Restricting the ray path length
            er=er+1;    
            break;  % Reterning to the main program for starting a new ray
        end
    end
end