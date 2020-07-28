% Define arduino
a = arduino();

% assign and call photodiode pins
configurePin(a,'A0','AnalogInput');
inPin0=configurePin(a,'A0');
configurePin(a,'A5','AnalogInput');
inPin5=configurePin(a,'A5');
configurePin(a,'A2','AnalogInput');
inPin2=configurePin(a,'A2');
configurePin(a,'D8','DigitalOutput');
outPin8=configurePin(a,'D8');
configurePin(a,'D11','DigitalOutput');
outPin11=configurePin(a,'D11');
configurePin(a,'D12','DigitalOutput');
outPin12=configurePin(a,'D12');

%set up vectors for recording and graphing values
GFP_vec = NaN(288,1);
Laser_vec = NaN(288,1);
t_vec = NaN(288,1);

%run for 24 hrs (1440 minutes). take measurements every 5 min and put
%variables into vectors
t_max = 1440;
t = 0;
j=1;
i=1;
while(t<t_max)
        t_vec(j) = t;
        writeDigitalPin(a,'D8',0);
        writeDigitalPin(a,'D11',0);
        writeDigitalPin(a,'D12',0);
        Laser = 0;
        GFP = 0;
        
        %Laser readings...
        
        writeDigitalPin(a,'D8',1);
        pause(2);   
        %take 1000 measurements
        for i=i:1000
            Laser = Laser+readVoltage(a,'A5');
            pause(0.01);
        end
        i=1;
        writeDigitalPin(a,'D8',0);
        %find avg value
        Laser = Laser/1000;
        %find OD from avg value
        Laser = (log10(4.218/Laser))/2.54;
        Laser_vec(j) = Laser;
        pause(2);
        
        %GFP readings..
        
        %take 1000 readings
        for i=i:1000
            writeDigitalPin(a,'D12',1);
            GFP = GFP+readVoltage(a,'A0');
            pause(0.01);
        end
        i=1;
        writeDigitalPin(a,'D12',0);
        %take avg value
        GFP = GFP/1000;
        GFP_vec(j) = GFP;
        j=j+1;
        pause(300);
        t=t+((300+10+10+1+1)/60);
end

%remove unused array elements
GFP_vec = rmmissing(GFP_vec);
Laser_vec = rmmissing(Laser_vec);
t_vec = rmmissing(t_vec);

%plot the results
plot(t_vec,GFP_vec,'g-')
hold on;
plot(t_vec,Laser_vec,'r-')
title('OD and GFP Readings of Bacterial Sample')
xlabel('Time (minutes)') 
ylabel('Voltage (V) [GFP] or AU [OD]')
legend({'GFP','OD'},'Location','southwest')
