%% Aero 302 Lab 3
% Section 7 Group 1

clear; clc;

% span = 26.453125 in
% +/- 1/16 in
% chord = 4.680 in
% +/- .001 in

span = 0.671909375;     %m
spansigma = 0.0015875;  %m
chord = 0.118872;       %m
chordsigma = 0.0000254; %m

s = span*chord;

lbf2N = 4.44822;

q20 = 1/2*1.225*20^2;

% import all of the data

fds = fileDatastore('C:\Users\joels\Documents\MATLAB\Lab 3 Data', 'ReadFcn', @importdata);

fullFileNames = fds.Files;

numFiles = length(fullFileNames);

% Loop over all files reading them in and plotting them.

for k = 1 : numFiles
    data(k) = load(fullFileNames{k},"-mat");
end

% Pull alpha info from filenames
for k = 1 : numFiles - 1
    alpha(k) = str2double(fullFileNames{k}(52));
    if fullFileNames{k}(53) ~= 'R'
        alpha(k) = (alpha(k) * 10) + str2double(fullFileNames{k}(53));
    end   
    if fullFileNames{k}(51) == 'N'
        alpha(k) = -alpha(k);
    end
end

for k = 1 : numFiles
    data(k).F = data(k).F * lbf2N;
end

% Calc lift force, drag force, cl, cd
for k = 1 : numFiles - 1
    L(:,k) = (-data(k).F(:,3));
    Cl(:,k) = L(:,k)/((mean(data(k).P(:,1))-mean(data(k).P(:,2)))*s);
    
    index = floor((alpha(k)+30)*73.46);
    tare = ((abs(mean(data(67).F(index,1)))/q20)*(mean(data(k).P(:,1))-mean(data(k).P(:,2))));
    data(k).F(:,1) = data(k).F(:,1) - tare;

    D(:,k) = data(k).F(:,1);
    Cd(:,k) = D(:,k)/((mean(data(k).P(:,1))-mean(data(k).P(:,2)))*s);
end

[alpha,I] = sort(alpha);
L  =  L(:,I);
D  =  D(:,I);
Cl = Cl(:,I);
Cd = Cd(:,I);

%% Plots

y = L./D;

figure(1)
hold on; grid on;
errorbar(alpha(1:3:60),mean(Cl(:,1:3:60)),std(Cl(:,1:3:60)),'r')
errorbar(alpha(2:3:60),mean(Cl(:,2:3:60)),std(Cl(:,2:3:60)),'g')
errorbar(alpha(3:3:60),mean(Cl(:,3:3:60)),std(Cl(:,3:3:60)),'b')
title("Coeffecient of Lift v AoA")
legend("Re1","Re2","Re3")
xlabel("Angle of Attack (Degrees)")
ylabel("Cl")

figure(2)
hold on; grid on;
errorbar(alpha(1:3:60),mean(Cd(:,1:3:60)),std(Cd(:,1:3:60)),'r')
errorbar(alpha(2:3:60),mean(Cd(:,2:3:60)),std(Cd(:,2:3:60)),'g')
errorbar(alpha(3:3:60),mean(Cd(:,3:3:60)),std(Cd(:,3:3:60)),'b')
title("Coeffecient of Drag v AoA")
legend("Re1","Re2","Re3")
xlabel("Angle of Attack (Degrees)")
ylabel("Cd")

figure(3)
hold on; grid on;
errorbar(alpha(1:3:60),mean(y(:,1:3:60)),std(y(:,1:3:60)),'r')
errorbar(alpha(2:3:60),mean(y(:,2:3:60)),std(y(:,2:3:60)),'g')
errorbar(alpha(3:3:60),mean(y(:,3:3:60)),std(y(:,3:3:60)),'b')
title("L/D vs alpha")
legend("Re1","Re2","Re3")
xlabel("Angle of Attack (Degrees)")
ylabel("L/D")

%% Derivative Function

function[sigma_Cl, sigma_Cd, sigma_ClCd] = sigmafunc(p0, ps, fx, fz, Cl, Cd)

  

% fx, fz, Cl, Cd must be row vectors of length 66 (1,66)

 

    qinf = p0 - ps;

    qinfp = mean(qinf);

    std_qinf = std(qinf);

    b = 67.1909375;

    c = 11.8872;

    S = b*c;

 

    partialL = 1/((qinfp)*(S));

    partialq = -fz./((S)*(qinfp)^2);

    partialS = -fz./((qinfp)*S^2);

    std_fz = std(fz);

    std_fx = std(fx);

 

    partialb = c;

    partialc = b;

    sigmab = 0.15875;

    sigmac = 0.00254;

 

    sigma_S = sqrt( (partialb*sigmab)^2 + (partialc*sigmac)^2 );

 

    sigma_Cl = zeros(1,66);

    for i = 1:length(Cl)

        sigma_Cl(1,i) = sqrt( (partialL*std_fz)^2 + (partialq(1,i)*std_qinf)^2 + (partialS(1,i)*sigma_S)^2 );

    end

 

    partialD = 1/((qinfp)*(S));

    partialqD = -fx./((S)*(qinfp)^2);

    partialSD = -fx./((qinfp)*S^2);

 

    sigma_Cd = zeros(1,66);

    for i = 1:length(Cd)

        sigma_Cd(1,i) = sqrt( (partialD*std_fx)^2 + (partialqD(1,i)*std_qinf)^2 + (partialSD(1,i)*sigma_S)^2 );

    end

 

    partialCl = 1./Cd;

    partialCd = -Cl./(Cd.^2);

 

    sigma_ClCd = zeros(1,66);

    for i = 1:length(sigma_ClCd)

        sigma_ClCd(1,i) = sqrt( (partialCl(1,i)*sigma_Cl(1,i))^2 + (partialCd(1,i)*sigma_Cd(1,i))^2 );

    end

 

end







