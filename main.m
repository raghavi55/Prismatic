% This code calculates the transmittances for both "apex-to-base" and
% "base-to-apex" directions for a prismatic sheet.
clear; 
prompt={'Number of rays:','Refraction index:','Simulation, on=1, off=0','reflection=[0-1]','Apex_0','Angular distribution; 1=Lambertian, 2= Uniform, 3= directional'};
def={'10000','1.53','1','0.05','60','1'};
d_title='Input';
n_lines=1;
N=inputdlg(prompt,d_title,n_lines,def);
N_ray=str2double(N(1));     % Number of rays     
nr=str2double(N(2));    % n is the refraction index of prism
simulation_key=str2double(N(3));    % 1= showing the ray's trace, 0= don't show
reflect=str2double(N(4)); %Reflection coeficient of sources surfsces. A number between 0 and 1   
A0=str2double(N(5)); %The smalest apex angle
As=5; %Step for apex angles
st=(180-A0)/As+1;  
AngD=str2double(N(6)); %The angular distribution of initial radiation
for ff=1:st % A loop over different prism apex angles
    Alf(ff)=A0-0.001+As*(ff-1)
    A=1/tand(Alf(ff)/2); 
    
    %Definition of the geometrical structure 
    L=[0    0     -2       1;   % Base of piramid         1
       1    0    1/A       1;   % Lateral surfaces        2
      -1    0    1/A       1;   % Lateral surfaces        3
       0    1    1/A       1;   % Lateral surfaces        4
       0   -1    1/A       1;   % Lateral surfaces        5
       1    0      0       1;   % Lateral mirror          6
      -1    0      0       1;   % Lateral mirror          7
       0    1      0       1;   % Lateral mirror          8
       0   -1      0       1;   % Lateral mirror          9
       0    0   1/(A+0.25) 1;   % Atop source             10
       0    0  -1/(A+0.75) 1];  % Beneath source          11
    
    clf; 
    unit_cell; view(2,5);       % Drawing the 
    N_surfaces=11;    
    for gg=1:2 
        Start_plate=9+gg;    
        n01=1;      n02=nr;
        er=0;       r11=0;       r10=0;  rep=0;
        ray_tracer;
        Tcb(gg,ff)=r11/(r11+r10); %Transmission for apex-to-base direction
        Tbc(gg,ff)=r10/(r11+r10); %Transmission for base-to-apex direction
        EEE(gg,ff)=er/(r11+r10+er); %Errors
    end
end
res=[Alf;Tcb(1,:);Tbc(2,:)];
plot(res(1,:),res(2,:),'-',res(1,:),res(3,:),'--',res(1,:),EEE(2,:),'.',res(1,:),abs(res(3,:)-res(2,:)),'-.');
