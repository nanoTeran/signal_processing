amp=3; 
fs=8000;  % sampling frequency
duration=0.2;
freq=[941 697 697 697 770 770 770 852 852 852 941 941; 1336 1209 1336 1477 1209 1336 1477 1209 1336 1477 1209 1477];
time=linspace(0,duration, duration*fs);

values = [1 5 9 7 5 3 10 0 11 0];  %*=10; #=11

Tones=[];
for i=1:length(values)
    Tones=[Tones (amp*sin(2*pi* freq(1,values(i)+1)*time)+amp*sin(2*pi* freq(2,values(i)+1)*time))];
end