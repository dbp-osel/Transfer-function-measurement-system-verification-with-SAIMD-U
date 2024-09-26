%% MATLAB code prepared for Regulatory Science Tool - Tabletop E-field generator for in vitro transfer fucntion model validation for device test in 1.5T and 3T MRI
prompt = "Enter the frequency to verify (64 or 128) MHz?   "; % choose the frequency to calculate (either 64 or 128)
x=input(prompt);
B0=sprintf('%.1f',1.5*x/64); % Relevant main magnetic field strength to the input frequency
sprintf("%.1f"+B0+" T testing selected")

% Skipping extra range of data acquired outside SAIMD-U
sk_b=2; % skip two pts acquired before the tip
sk_e=2; % skip one point acquired after SAIMD-U and one point acquired when a piX excitor placed on top of the generic IPG can (excitor-can copuling)

% Read JSON data from a file
cd(sprintf('%sT_SAIMD-U',B0))  % Change directory to the folder containing data either 1.5T or 3.0T
jsonFileName=sprintf('SAIMD-U_%dMHz_PiXscan1.json',x)  % Time-domain E-field sensor (TDS) measurement at 10 mm distance from the SAIMD-U tip
jsonData = jsondecode(fileread(jsonFileName)); % Access JSON data and convert it to MATLAB variables 
Meas=jsonData.Ch1R(1+sk_b:end-sk_e)+ i*jsonData.Ch1I(1+sk_b:end-sk_e); % Convert data in complex value, removing data range outside SAIMD-U lead

jsonFileName=sprintf('BG_SAIMD-U_%dMHz_Incident BG_Incident1.json',x) % background measurement with TDS without SAIMD-U (background response)
jsonData = jsondecode(fileread(jsonFileName)); % Access JSON data and convert it to MATLAB variables 
Bg=jsonData.Ch1R(1+sk_b:end-sk_e)+ i*jsonData.Ch1I(1+sk_b:end-sk_e); % Convert data in complex value, removing data range outside SAIMD-U lead

% Transfer function model calucaltion (Background reponse removal and scaling)
TF=Meas-Bg; % Subtract the background response from the measurement (Additional process required for the piecewise excitation method)
nTF=TF./max(abs(TF)); % Scale transfer function model in magnitude that maximum equals one 
ph=unwrap(angle(nTF))-angle(nTF(1));  % unwrap the phase of transfer function model and phase offset that phase at the SAIMD-U tip eqauls zero

% Normalized root mean squared error (NRMSE) calculation (Compared to FDA OSEL laboratory meausrement data) 
fid = fopen(sprintf('Tip%dHPM.txt',x),'rt'); % read the target SAIMD-U transfer function model value measured at FDA laboratory
data = textscan(fid, '%f %f %f', 'HeaderLines',0);
fclose(fid);  
abs_Error=sqrt(sum((abs(nTF)-data{2}).^2)/(sum(data{2}.^2))) % NRMSE error in Magnitude 
pha_Error=sqrt(sum((ph-data{3}).^2)/(sum(data{3}.^2))) % NRMSE error in Phase

% (additional metric) Calculate correlation between measurement and target value
Mag_corr=corr(abs(nTF),data{2}) 
Pha_corr=corr(ph,data{3})

if abs_Error < 0.2 && pha_Error < 0.2  % Confirm whether the TF measurement system is verified
    "Transfer function measurement system is verified. Ready to measure transfer function model of your device"
else
    "Transfer function measurement system is not verified. Please recheck the measurement system including phanotm properties, target frequency, sample positions and measure the SAIMD-U model again."
end

figure; ah1(1)=subplot(2,1,1); fts=14; % Display subplot of transfer function model in magnitude with font size 14
plot(abs(nTF(:)),'-or');hold on; % Plot the scaled transfer function model measurement in magnitude
plot(abs(data{2}),'-k'); hold on; % Plot the measurement target value in magnitude
title ("Magnitude of SAIMD-U transfer function model at "+ B0 + " T (normalized max magnitude = 1)")
xlabel('Distance from the tip(cm)');ylabel ('Magnitude(1/m)'); fontsize(gca, fts,'points'); 
xl=length(nTF); axis([1 xl 0 max(abs(nTF)+0.1)]); grid on; 
legend('Measured SAIMD-U (Magnitude)','Target value in magnitude','Location','NorthEastOutside'); legend boxoff
hold on; ah2(2)=subplot(2,1,2);  % Subplot of transfer function model in phase 
plot(ph,'-ob'); hold on  % Plot the scaled SAIMD-U transfer function model in phase (phase offset: phase at the SAIMD-U tip = 0)
plot(data{3},'-k');  % Plot the SAIMD-U measurement target value in phase
axis([1 xl min(ph) 0.2]); grid on; 
title ("Phase of SAIMD-U transfer function model at " + B0 + " T (Phase offset: Phase at the SAIMD-U tip = 0)");
xlabel('Distance from the tip(cm)');ylabel ('Phase (rad)'); fontsize(gca, fts,'points'); 
legend('Measured SAIMD-U (Phase)','Target value in Phase','Location','NorthEastOutside'); legend boxoff; set(gcf,'WindowState','maximized')

cd ..  % Change to the initial directory

%% (Optional scaling step before Tier 3 assessment with human body models) Normalization that integral of transfer function model in magnitude equales one 
dl=0.01;% Distance between piX data measurement point, 10 mm = 0.01 m 
nTF=TF./trapz(dl,TF) % Scale the transfer function by the numerical integration using a Trapezoidal rule with unit spacing dl
if abs(trapz(dl,nTF))-1 < 1e-6  % Check whether integral of scaled transfer function in magnitude equaled 1
    "Transfer function model is scaled that integral of TF magnitude equaled 1"
else
    "Transfer function model is not scaled that integral of transfer function magnitude is not equaled 1"
end
figure; ah1(1)=subplot(2,1,1); % Display subplot of scaled transfer function in mangnitude
plot(abs(nTF(:)),'-or');grid on; hold on  % scale the y-axis unit from 1/cm (i.e.,10mm per point) into 1/m by multiplying 100
title ("Magnitude of SAIMD-U transfer function model at "+ B0 + " T (normalized integral of transfer function in magnitude = 1)")
xl=length(nTF); xlabel('Distance from the tip(cm)');ylabel ('Magnitude(1/m)'); fontsize(gca, fts,'points'); 
legend('Measured SAIMD-U (Magnitude)','Location','NorthEastOutside');legend boxoff; 
axis([1 xl 0 max(abs(nTF))+1]); grid on; 
hold on; ah2(2)=subplot(2,1,2); plot(ph,'-ob'); % Display subplot of scaled transfer function in phase
axis([1 xl min(ph) 0.2]); grid on; 
title ("Phase of SAIMD-U transfer function model at " + B0 + " T (Phase offset: Phase at the SAIMD-U tip = 0)");
xlabel('Distance from the tip(cm)');ylabel ('Phase(rad)'); fontsize(gca, fts,'points'); grid on;
legend('Measured SAIMD-U (Phase)','Location','NorthEastOutside'); legend boxoff; set(gcf,'WindowState','maximized')