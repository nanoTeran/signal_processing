clear all;

disp('El programa se ha iniciado')
pause(1)
disp('Por favor ingrese el nombre del archivo de entrada, sin ninguna extension.')
pause(1.5)
disp('Este archivo debe ser del tipo [.mat].')
pause(1.5)
disp('Si desea usar el nombre de archivo predeterminado, [SenalPrueba], solo presione [Enter].')
pause(2)
SignalName = input('','s');
if isempty(SignalName)
    SignalName = 'SenalPrueba';        
end
SignalName = [SignalName,'.mat'];
Signal = load(SignalName);
SignalName = fieldnames(Signal);
SignalName = cell2mat(SignalName(1,1));
Signal = Signal.(SignalName);
[SignalSize(1),SignalSize(2)]= size(Signal);
if SignalSize(1)>SignalSize(2)
   Signal = transpose(Signal); 
end

fs=8000;
Ws = round(2*pi*fs);
t=linspace(0,length(Signal)/fs,length(Signal));
characters=blanks(length(Signal)/(fs*0.2));   %CREATING the Output with number tones of the signal

W=[697,770,852,941,1209,1336,1477];         %In Hz.
W=round((2*pi).*W);                         %In Rad/s.
W1=zeros(1,6);                             %Initialized Vector for storing the different W1.
Rp=3;
As=20;
N_array=zeros(1,7);                         %Initialized Vector for storing the optimal Order for each filter. 
Wc=zeros(2,7);                              %Initialized Array for storing the optimal Cut Frequency (?) for each filter.   

for i=1:6
   W1(i)=round(sqrt(W(i)*W(i+1)));         %W1 as the geometric mean of consecutive frequencies.
end

for i=2:6                                   %Calculate the order and Cut frequency for the bandpass filters and store them.
    [N_array(i),Wc_n]=buttord([W1(i-1),W1(i)],[W(i-1),W(i+1)],Rp,As,'s');
    Wc(1,i)=Wc_n(1);                        
    Wc(2,i)=Wc_n(2);
end
N=max(N_array);                             %Choose the maximum order of the desired bandpass filters.
[N_array(1),Wc(1)]=buttord(W1(1),W(2),Rp,As,'s');       %Calculate the optimal order and Cut frequency for Low pass
[N_array(7),Wc(13)]=buttord(W1(6),W(6),Rp,As,'s');      %And highpass filters. The orders are not used as they're too high.   


%%---- Analogic Filters ----%%
b = cell(1,7);                           %Cell array of b vectors (analogic)
a = cell(1,7);                           %Cell array of a vectors (analogic)

figure(1)
title('Analogic Filters 1');
hold on
[b{1},a{1}] = butter(N,Wc(1,1),'s');              %Bode plot of the first 4 filters in the first figure.
bode(b{1},a{1})

for i=2:6
    if i<5
        figure(1)
        hold on
        [b{i},a{i}] = butter(N,[Wc(1,i),Wc(2,i)],'bandpass','s');
        bode(b{i},a{i})
        hold off
    else
        figure(2)
        hold on
        title('Analogic Filters 2');
        [b{i},a{i}] = butter(N,[Wc(1,i),Wc(2,i)],'bandpass','s');
        bode(b{i},a{i})   
        hold off
    end 
end

figure(2)
hold on
[b{7},a{7}] = butter(N,Wc(1,7),'high','s');
bode(b{7},a{7})
hold off


%%---- IIR impinvar Method ----%%
bz =cell(1,7);                           %Cell array of b vectors (IIR impinvar digital)
az =cell(1,7);                           %Cell array of b vectors (IIR impinvar digital)

[bz{1},az{1}] = impinvar(b{1},a{1},fs);
figure(3)
title('IIF impinvar 1');
freqz (bz{1},az{1})

for i=2:6
    if i<5
        figure(3)
        title('IIF impinvar 1');
        hold on
        [bz{i},az{i}] = impinvar(b{i},a{i},fs);
        freqz (bz{i},az{i})
        hold off
    else
        figure(4)
        hold on
        title('IIF impinvar 2');
        [bz{i},az{i}] = impinvar(b{i},a{i},fs);
        freqz (bz{i},az{i})
        hold off
    end 
end

[bz{7},az{7}] = impinvar(b{7},a{7},fs);
figure(4)
hold on
freqz(bz{7},az{7})


%%----IIR Bilinear method -----%%
bzb = cell(1,7);                              %Cell array of b vectors (IIR bilinear method digital)
azb = cell(1,7);                              %Cell array of b vectors (IIR bilinear method digital)

Tbilinear = atan(pi*Wc./Ws).*(2./Wc);         %Using the tangent transformation of BILINEAR (previous knowing that
%                                             Wd = Wa*pi /(Ws/2)
fbilinear = 1./Tbilinear;


%%% Same process than other filters 
[bzb{1},azb{1}] = bilinear(b{1},a{1},fbilinear(1,1));
figure(7)
title('IIF bilinear 1');
freqz (bzb{1},azb{1})

for i=2:6
    if i<5
        figure(7)
        hold on
        title('IIF bilinear 1');
        [bzb{i},azb{i}] = bilinear(b{i},a{i},fbilinear(1,i));
        freqz (bzb{i},azb{i})
        hold off
    else
        figure(8)
        hold on
        title('IIF bilinear 2');
        [bzb{i},azb{i}] = bilinear(b{i},a{i},fbilinear(1,i));
        freqz (bzb{i},azb{i})
        hold off
    end 
end

[bzb{7},azb{7}] = bilinear(b{7},a{7},fbilinear(1,7));
figure(8)
hold on
freqz(bzb{7},azb{7})


%%-- FIR FILTERS --%%        
bfir =cell(1,7);                         %Cell array of b vectors (IIR digital) (numerators)

%Find N with Harris method
nf_array(1) = round((2*pi*fs*As)/(22*(W1(6)-W(6))));     %low-high Pass transition band wide
nf_array(2) = round(fs*As/(22*(W(7)-W1(6))));       %band     Pass transition band wide 
nf=max(nf_array);                                   %Select higher nf (number of coefficients)
if rem(nf,2)~=0                                     %Select even nf
    nf=nf+1;
end
%Windows Design
WR = ones(nf+1,1);                                 %Rectangular
WHn = hann(nf+1);                                  %hanning
WB1 = blackman(nf+1);                              %blackman
WBa = bartlett(nf+1);                              %Bartlett

%%% Same process than other filters 
%Designed Filters
figure(5)
hold on
title('FIR 1');
bfir{1} = fir1(nf,Wc(1,1)/(Ws/2));                           %Filtering LP               
freqz(bfir{1},1)                                             %FIL LP Graphic
for i=2:6
    if i<5
        bfir{i} = fir1(2*nf,[Wc(1,i) Wc(2,i)]/(Ws/2));       %Find FIR BP with 2*nf (double of LP nf)  
        hold on
        freqz(bfir{i},1)                                     %BP FIR graphic
    else
        hold off
        figure(6)                                            %Active figure(6) of the last (5-6-7) filters (last frequencies)
        hold on
        title('FIR 2');
        bfir{i} = fir1(2*nf,[Wc(1,i) Wc(2,i)]/(Ws/2));       %Find FIR BP with 2*nf (double of LP nf)  
        hold on
        freqz(bfir{i},1)
    end        
end
bfir{7} = fir1(nf,Wc(1,7)/(Ws/2),'high');                    %HP Filter                    
hold on                    
freqz(bfir{1},1)                                             %HP graphic
hold off

% % b = fir1(N,W_1(6)/(Ws/2),WR,'high');        %Rectangular
% % freqz(b,1)
% % b = fir1(N,W_1(6)/(Ws/2),WHn,'high');       %Hanning
% % freqz(b,1)
% % b = fir1(N,W_1(6)/(Ws/2),WB1,'high');       %Blackman
% % freqz(b,1)
% % b = fir1(N,W_1(6)/(Ws/2),WBa,'high');       %bartlett
% % freqz(b,1)


%%--- FILTERING SIGNAL ---%%     
y = cell(3,7);     % Store Filtered Signal outputs
%                     {1,:} Analogic Filtering 
%                     {2,:} IIR Filtering
%                     {3,:} FIR Filtering

%%% Analogic Filter
%figure 9
figure ('name','Filtered Signal with analogics' )

for i=1:8                             % i=1:7 it is the seven filters, i=8 it is for plot the ORIGINAL signal to compare 
    if i<8
        subplot(3,3,i)                %It is a sub"graphic" in figure 11 for each filter         
        title(i)
        SYS = tf(b{i},a{i});          % transform NUM/DEM ->SYS
        y{1,i} = transpose(lsim(SYS,Signal,t));  % y{1,:}= output of Analogic filtered signal
        plot(t,y{1,i})                %it graphs the t vs filtered signal 
    else
        subplot(3,3,9)                %Subplot for ORIGINAL Signal to compare
        title('Original Signal')               
        plot(t,Signal)                % Plot the input signal to compare in time
    end
    
end

%%% IIR  Filter
%figure 10
figure ('name','Filtered Signal with IIR ')
for i=1:8
    if i<7
        subplot(3,3,i)                         %It is a sub"graphic" in figure 12 for each filter
        title(i)
        y{2,i} = filtfilt(bz{i},az{i},Signal); % y{2,:} = output of IIR impinvar filtered signal
        plot(t,y{2,i}) %it graphs the t vs filtered signal
    elseif i==7
        subplot(3,3,i)                         %It is a sub"graphic" in figure 12 for each filter
        title(i)
        y{2,i} = filtfilt(bzb{i},azb{i},Signal); % y{2,:} = output of IIR bilinear filtered signal
        plot(t,y{2,i}) %it graphs the t vs filtered signal
    else 
        subplot(3,3,9)                         %Subplot for ORIGINAL Signal to compare
        title('Signal')
        plot(t,Signal)                         % Plot the input signal to compare in time
    end
    
end


%%% FIR Filter
    %figure (11)
figure ('name','Filtered Signal with FIR ')
N = 6;

for i=1:8

    if i<8
       
        hold on
        
        subplot(3,3,i)                            %It is a sub"graphic" in figure 13 for each filter
        title(i)
        y{3,i} = fftfilt(bfir{i},Signal,N);       % y{3,:} = output of FIR filtered signal
        plot(t,y{3,i})                            %it graphs the t vs filtered signal
        
    else
        hold on
        subplot(3,3,9)                            %Subplot for ORIGINAL Signal to compare
        title('Original Signal')
        plot(t,Signal)                            % Plot the input signal to compare in time   
    end
    
end

%------ Decoder ------%
% Now, we will find the power of each filtered Signal to find which is the most powerful tone
% it means, the highest frequencies of each part to detect which number it is sounding

%example with FIR filter : y{3,:}
% 7 : 7 filters 
% 2 : 2 filters to decide what number is the tone 
% tones: the number of tones in the signal 

tones =length(characters);                         %For "SenalPrueba"it is 10 because the input signal has 10 "sounds" (0123456789)
P = zeros(tones,7);                                %P: Store the Power of each filtered Signal tone
pos = zeros(tones,2);                              %pos:it stores the number of the filters with the
%                                                        highest power responses, therefore, which tone sounds

%next, we start to analize each tone of the signal with each filter to know what
%frequencies are sounding. "i" will go around each tone, while "j" will go
%around each filter for every tone

for i=1:tones                                         %i: it's each tone(10 tones in SenalPrueba)
    
    %post is variable to check the values of a range of the signal, it will analyze the
    %power in a range around the middle of each different tone(number) sounds.
    %example: if i=1; post will be the postion in the Signal array of the value in t=0.1; 
    %it will analyse neighboring values around that time, when it is sounding the '0' (941Hz + 1336Hz)
                     
    post=round(length(Signal)*(i-0.5)/(tones));    %Takes time to analize each tone
    for j=1:7                                      %It analizes the seven filters by tone
        
%         This cicle it is divided in 2 parts(2 if): first determines
%         which of first 4 filters has the highest power. Second determines
%         which of last 3 fliters has the highest power. Then it is save in
%         pos which filters are passing the sound, then, it is possible to know which "number" its sounding
        
        P(i,j)=(max([y{1,j}(1,post-3) y{1,j}(1,post-2) y{1,j}(1,post-1) y{1,j}(1,post) y{1,j}(1,post+1) y{1,j}(1,post+2) y{1,j}(1,post+3)]))^2; 
        %Save power of the higuest valued point of each section of the filter.
        if j<5 && (P(i,j) == max(P(i,1:4)))        %choose the max P of the first 4 filters
            pos(i,1)=j;
        end                                        %choose the max P of the last 3 filters
        if j>=5 && P(i,j) == max(P(i,5:7))
            pos(i,2)=j;
        end
    end
end

% The output of the filtered signal is saved in 'characters'. Example: with
% Signal "SenalPrueba.mat" it will store "0123456789". 
% if the Signal is "SenalPrueba2.mat" it will store "*#*#*#"
%Now, it's selected which characters are stored in the input signal
%(0123456789*#)

for i=1:tones                                       
    if pos(i,1) == 1                        %if the first filter was lowpass (697Hz)
        if pos(i,2) == 5                    %if the second filter was bandpass (1209Hz)
            characters(1,i) ='1';               %Asign '1' to the output 
        elseif pos(i,2) == 6                %Repeat same logic to all the other characters
            characters(1,i) ='2'; 
        else
            characters(1,i) ='3';
        end
    elseif pos(i,1) == 2 
        if pos(i,2) == 5
            characters(1,i) ='4';
        elseif pos(i,2) == 6
            characters(1,i) ='5';
        else
            characters(1,i) ='6';
        end
    elseif pos (i,1) ==3
        if pos(i,2) == 5
            characters(1,i) ='7';
        elseif pos(i,2) == 6
            characters(1,i) ='8';
        else
            characters(1,i) ='9';
        end
    elseif pos (i,1) ==4
        if pos(i,2) == 5
            characters(1,i) ='*';
        elseif pos(i,2) == 6
            characters(1,i) ='0';
        else
            characters(1,i) ='#';
        end
    end
end

disp([' La senal es: ' characters ])

