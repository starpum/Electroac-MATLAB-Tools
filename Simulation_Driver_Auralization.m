% Simulation Dirver Auralization : a simple GUI-less matlab software to
% simulate a driver in different settings (closed box, vented box) ; in
% linear and non-linear mode

% The calculations are based on the state-space model for linear and
% non-linear systems
% References :
%               - Kleinder, Electroacoustics
%               - Agervist, State-space models simulations
%               - Panzer, AkAbaK manual

% Created by Kolya Corno for the The IlaB - all rights reserved. 
% Date created : Aug. 2019
% Last updated : Nov. 2019

% Inputs :  - Sampling frequency (Hz) [default 44100]
%           - Number of FFT points [default 4000]
%           - LEM driver parameters, in SI units [default = Dayton ND91-4 specs]
%           - Driver configuration [closed box ; vented box]
%           - Box volume (L) [default 1]
%           - Vent specifications (mm)  [default 0.01 length and 0.001 diameter]
%           - PR specifications [default based on Dayton audio ND65-PR]
%           - Type of study (linear/non-linear) and non-linear parameters
%           if needed. 

% Outputs : - Exrcusion / freqnecy plot in the driver's linear range
%           - Membrane displacement / frequency in the driver's linear AND
%           non-linear mode 
%           - SPL @1W@1m on-axis with acoustic center / frequency for both
%           linear and non-linear models
%           - Auralized simulated output using selected parameters (if user
%           requests it)

% TODO : - implement solutions for PR and other systems. 
%        - optimize the nested loops computation by switching to a full
%        matrical-based calculation
%        - optimise the auralization and wawwriting 
%        - Create a GUI (most likely in Python + open-source access) 


close all
clc

%% User defined constants / study

prompt = {'Sampling frequency (Hz) : ','Number of points : '};
dlgtitle = 'Study parameters (if you''re unsure, keep the default) : ';
definput = ({'44100','4000'});
dims = [1 30];
answer = inputdlg(prompt,dlgtitle,dims,definput);
Fsamp = str2double(answer{1});
N     = str2double(answer{2});

u = zeros(1,N); 
u(1) = 1;

rho = 1.2;
c = 344;

Tsamp = 1/Fsamp;
t     = 0:Tsamp:(N-1)*Tsamp;

f   = 0:Fsamp/N:Fsamp/2-Fsamp/N; % frequency vector
w   = 2*pi*f;

% Analysis_type_menu = menu('What study should be completed?','linear (faster)',...
%     'non-linear (slower, more accurate for high excursions, requieres LSI parameters');
% switch Analysis_type_menu
%     case 1
%         analysis_type = 1;
%     case 2
%         analysis_type = 2;
% end

analysis_type = 1; %Only linear solution implemented now

prompt = {'Distance of radiation (m) : '};
dlgtitle = 'Select the distance of radiation (in m) ';
definput = ({'1'});
dims = [1 15];
answer = inputdlg(prompt,dlgtitle,dims,definput);
rLS = str2double(answer{1});

%% Speaker definition (user input)

prompt_LEM = {'Re (Ohms)';'Le (mH)';'Bl (T/m)';'Mms (g)';'Fs (Hz)';...
    'Cms (mm/N)';'Qms';'Qes';'Qts';'Sd (cm2)'};
dlgtitle_LEM = 'Enter the TnS parameters for your driver : ';
definput_LEM = ({'4';'0.1';'2.5';'1';'100';'0.9';'1.2';'1.2';'1.2';'30'});
dims = [1 100];
answer_LEM = inputdlg(prompt_LEM,dlgtitle_LEM,dims,definput_LEM);

Re   = str2double(answer_LEM{1});
Le   = str2double(answer_LEM{2})*10^-3;  % Convertion to H
Bl   = str2double(answer_LEM{3});
Mms  = str2double(answer_LEM{4})*10^-3;  % Convertion to kg
Fs   = str2double(answer_LEM{5});
ws   = 2*pi*Fs;
Cms  = str2double(answer_LEM{6})*10^-3;  % Convertion to m/N
Kms  = 1/Cms;
Qms  = str2double(answer_LEM{7});
Qes  = str2double(answer_LEM{8});
Qts  = str2double(answer_LEM{9}); 
Rms  = (Bl^2*Qes)/(Re*Qms);
Sd   = str2double(answer_LEM{10})*10^-4; % Convertion to m2
Cas  = Sd^2*Cms;
Vas  = rho*c^2*Cas;

%% ------------------------------------------------------------------
%       Box parameters (see page 176 akabak manual + chapter 7 Kleiner, Electroacoustic for details) 

% box_type_inquiry = menu('Select the system type','infinite baffle','closed box','vented box','passive radiator');
% switch box_type_inquiry
%     case 1
%         box_type = 1;   %inf. baffle
%     case 2
%         box_type = 2;   %vented box 
%     case 3
%         box_type = 3;   %closed box
%     case 4
%         box_type = 4;   %PR
% end

box_type_inquiry = menu('Select the system type','closed box','vented box');
switch box_type_inquiry
    case 1
        box_type = 2;   %vented box 
    case 2
        box_type = 3;   %closed box
end

if box_type == 2 %closed box
    prompt_box = {'Box volume (L)'};
    dlgtitle_box = 'Enter the box specifications : ';
    definput_box = ({'1'});
    dims = [1 30];
    answer_box = inputdlg(prompt_box,dlgtitle_box,dims,definput_box);
    
    Vb = str2double(answer_box{1})*1000;    % Convertion to m3
    Cb = Vb/(rho*c^2);
    
elseif box_type == 3  %vented box
    prompt_box = {'Box volume (L)','Vent diameter (mm)','Vent length (mm)'};
    dlgtitle_box = 'Enter the box and vent specifications : ';
    definput_box = ({'1', '0.01', '0.001'});
    dims = [1 30];
    answer_box = inputdlg(prompt_box,dlgtitle_box,dims,definput_box); 
    
    Vb = str2double(answer_box{1})*1000;    % Convertion to m3
    Dv = str2double(answer_box{2})*10^3;    % Convertion to m
    Lv = str2double(answer_box{2})*10^3;    % Convertion to m
    
    Cb = Vb/(rho*c^2);
    Mv = rho*(Lv/(pi*Dv^2));
    Rv = rho*c/(pi*Dv^2);
    
elseif box_type == 4  %PR box
    prompt_box = {'Box volume (L)','Membrane surface Sd_pr (cm2)','PR mass Mms_pr (g)',...
        'PR compliance Cms_pr (mm/N)','PR resistance Rms_pr (kg/s)'};
    dlgtitle_box = 'Enter the box and PR specifications : ';
    definput_box = ({'1', '15','3.2', '1.45','0.56'});
    dims = [1 30];
    answer_box = inputdlg(prompt_box,dlgtitle_box,dims,definput_box); 
    
    Vb = str2double(answer_box{1})*1000;    % Convertion to m3
    Sd_pr = str2double(answer_box{2})*10^4;    % Convertion to m2
    Mms_pr = str2double(answer_box{3})*10^3;    % Convertion to kg
    Cms_pr = str2double(answer_box{4})*10^3;    % Convertion to m/N
    Rms_pr = str2double(answer_box{5});    
    
    Cb = Vb/(rho*c^2);
end

%% ------------------------------------------------------------------
%       Linear state space model  (F Agervist, Kleiner chapter 16)

if box_type == 1    % Inf baffle
        % TODO
        
elseif box_type == 2    % Closed box
    
    A=[-Re/Le,          0,  -Bl/Le; 
            0,          0,       1;     
       Bl/Mms, -1/(Mms*Cms), - Rms/Mms];
    
    B=[ 1/Le; 0; 0];
    
elseif box_type == 3    % Vented box
    A=[ -Re/Le,            0,    -Bl/Le,       0,      0; 
             0,            0,         1,       0,      0; 
        Bl/Mms, -1/(Mms*Cms), - Rms/Mms, -Sd/Mms,      0;
             0,            0,     Sd/Cb,       0,  -1/Cb;
             0,            0,         0,    1/Mv, -Rv/Mv];

    B=[ 1/Le; 0; 0; 0; 0];
    
elseif box_type == 4    % PR
        % TODO
end    
    

x = zeros(length(A),N);

% Matrix discretion (bilinear method)

Ad = inv(eye(size(A))-A*Tsamp/2)*(eye(size(A))+A*Tsamp/2);
Bd = inv(eye(size(A))-A*Tsamp/2)*B*Tsamp; 

% Difeomorphism / Algorythmic computation 
for k=2:N
    x(:,k)=Ad*x(:,k-1)+Bd*u(k-1);
end

% Extraction of datas from result matrix
current  = x(1,:);
disp     = x(2,:);
velocity = x(3,:);

QLS      = reshape(Sd*velocity,1,length(velocity));
pressure = x(length(A),:);
pressure = reshape(pressure,1,length(pressure));

% Transition to frequency domain 
ffti   = fft(current,N);
fftQLS = fft(QLS,N);

%       Radiated pressure 

% Derivative
a = 1;
b = [1 -1]*Fsamp;

dQLS = filter(b,a,QLS); % derivative dQLS/dt
tauLS = rLS/c;

NtauLS = round(Fsamp*tauLS);

propaLS = zeros(1,NtauLS);

propaLS(NtauLS) = rho/(4*pi*rLS);


prad = filter(propaLS,1,dQLS);

fftprad = fft(prad,N);
Lp = 20*log10(abs(fftprad(1:N/2)/2e-5));

%% Non-linear parameters (measured using Klippel LIS (no thermal effects))

    % Non-linear Bl

% pos = disp;    % Position of the voice coil in the gap (mm)
% Bl2 = -0.061744;
% 
% BlNL = Bl.*ones(size(pos))+Bl2.*pos.^2;
% 
%     % Non-linear Cms 
% 
% Cm1 = 0.0021830;
% Cm2 = -0.0026017e-3;
% 
% CmNL = Cms*(ones(size(pos))+Cm1*pos+Cm2*pos.^2);
% 
%     % Non-linear L 
% 
% L1 = -0.074168e-3;
% L2 = -0.004e-3;
% 
% LeNL = Le*(ones(size(pos))+L1*pos+L2*pos.^2);
% 
% %% ------------------------------------------------------------------
% %       Nonlinear state space controler (Bright, Pedersen, Rubak)
% 
% u_controled = (Re./BlNL).*((1./CmNL).*QLS + (Rms + BlNL.^2 ./ Re).*disp + Mms*u);
% 
% %% ------------------------------------------------------------------
% %       Nonlinear state space model  
% 
% xNL = zeros(3,N);
% 
% for i = 1:length(pos) 
%     ANL=[-Re/LeNL(i),               0,  -BlNL(i)/LeNL(i); 
%                  0,               0,                 1; 
%         BlNL(i)/Mms, -1/(Mms*CmNL(i)),           - Rms/Mms];
% 
%     BNL=[ 1/LeNL(i); 0; 0];
% 
%     % Matrix discretion - bilinear method
%     AdNL = inv(eye(size(ANL))-ANL*Tsamp/2)*(eye(size(ANL))+ANL*Tsamp/2);
%     BdNL = inv(eye(size(ANL))-ANL*Tsamp/2)*BNL*Tsamp; 
%    
% end
% 
%     % Algorythmic computation 
% for k=2:N
%     xNL(:,k)=AdNL*xNL(:,k-1)+BdNL*u(k-1);
% end   
% 
% % Extraction of datas from result matrix
% currentNL  = xNL(1,:);
% dispNL     = xNL(2,:);
% velocityNL = xNL(3,:);
% 
% QLSNL      = reshape(Sd*velocityNL,1,length(velocityNL));
% 
% % Transition to frequency domain 
% fftdispNL = fft(dispNL, N);
% fftiNL    = fft(currentNL,N);
% fftQLSNL  = fft(QLSNL,N);
% 
% %       Radiated pressure 
% 
% % Derivation 
% a = 1;
% b = [1 -1]*Fsamp;
% 
% dQLS = filter(b,a,QLSNL); % derivative dQLS/dt
% 
% rLS     = 0.05;
% tauLS   = rLS/c;
% NtauLS  = round(Fsamp*tauLS);
% propaLS = zeros(1,NtauLS);
% propaLS(NtauLS) = rho/(4*pi*rLS);
% 
% pradNL = filter(propaLS,1,dQLS);
% 
% fftpradNL = fft(pradNL,N);
% LpNL      = 20*log10(fftpradNL(1:N/2)/2e-5);

%% ------------------------------------------------------------------
%       Figures 

figure(1)
    subplot(211)
semilogx(f,20*log10(abs(fftQLS(1:N/2))))
grid on
xlabel('Frequency (Hz)'); ylabel('Amplitude (dB)');
xlim([50 1500])
title('Volume displacement : magnitude'); 
    subplot(212)
semilogx(f,angle(fftQLS(1:N/2)))
grid on
xlabel('Frequency (Hz)'); ylabel('Phase (rad)');
xlim([50 20000])
title('Volume displacement : phase');

figure(2)
    subplot(211)
semilogx(f,20*log10(abs(ffti(1:N/2))))
grid on
xlabel('Frequency (Hz)'); ylabel('Amplitude (dB)');
xlim([50 1500])
title('Intensity : magnitude'); 
    subplot(212)
semilogx(f,angle(ffti(1:N/2)))
grid on
xlabel('Frequency (Hz)'); ylabel('Phase (rad)');
xlim([50 20000])
title('Intensity : phase');

figure(3)
semilogx(f,Lp)
grid on
xlabel ('Frequency (Hz)'); ylabel('Amplitude (dB)');
xlim([50 20000])
title('Radiated pressure at selected distance')


% figure(1)
% semilogx(f,fftdispNL(1:N/2).*10^3)
% grid on
% xlabel('Frequency (Hz)'); ylabel('Displacement (mm)')
% title('Excursion');
% xlim([50 1500]) 

% figure(2)
%     subplot(211)
% semilogx(f,20*log10(abs(fftQLS(1:N/2))),f,20*log10(abs(fftQLSNL(1:N/2))))
% grid on
% xlabel('Frequency (Hz)'); ylabel('Amplitude (dB)');
% xlim([50 1500])
% title('TRF Magnitude');
% legend('Lineaar loudspaker Displacement','Nonlineaar loudspaker Displacement'); 
%     subplot(212)
% semilogx(f,angle(fftQLS(1:N/2)),f,angle(fftQLSNL(1:N/2)))
% grid on
% xlabel('Frequency (Hz)'); ylabel('Amplitude (dB)');
% xlim([50 20000])
% title('TRF Phase');
% legend('Lineaar loudspaker Displacement','Nonlineaar loudspaker Displacement'); 
% 
% figure(3)
% semilogx(f,Lp,f,LpNL)
% grid on
% xlabel ('Frequency (Hz)'); ylabel('Amplitude (dB)');
% xlim([50 20000])
% title('Radiated pressure')
% legend ('Linear radiated pressure','Nonlinear radiated pressure');