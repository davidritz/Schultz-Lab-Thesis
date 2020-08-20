import time
import piplates.DAQC2plate as DAQC2
from datetime import datetime
import sys
import csv
import pandas as pd
import numpy as np
import threading
import queue


# write LOG
def write_log(name,t,OD,concentrations,dilutions):
    with open(name,"a") as file:
        file.write("\n%.4f" % t)
        for item in OD:
            file.write(", %.4f" % item)
        for item in concentrations:
            file.write(", %.4f" % item)
        for item in dilutions:
            file.write(", %.4f" % item)


# read OD
def read_OD(): 
    OD_temp = np.zeros(6)
    DAQC2.setDOUTbit(1,7)
    time.sleep(0.2)
    for i in cultures:
        tempVolt = []
        for n in range(100):
            try:
                tempVolt.append(DAQC2.getADC(0,i))
            except:
                print('error in acquiring voltage, n=%s, culture %s, time %.4f' % (str(n+1), str(i+1), float(time.time()-time_init)/60))
        if(len(tempVolt)>0):
            AvgVolt = sum(tempVolt) / len(tempVolt)
            try:
                od = np.log10(blanks[i]/AvgVolt)/2.54
            except:
                print('error in calculationg OD, culture %s, time %.4f' % (str(i+1), float(time.time()-time_init)/60))
                od = 0
            OD_temp[i]=od
    DAQC2.clrDOUTbit(1,7)

    return OD_temp
    

# calculate fraction of volume in container from medium 2 and dilution rate in volumes/hour
def calculate(conc_ini, t_pump1, t_pump2):
    #conc = conc_ini * np.exp( - rate * t_pump1)
    #conc = 1 + (conc - 1) * np.exp( - rate * t_pump2)
    conc = rate * t_pump2
    dil = rate * (t_pump1 + t_pump2)
    return conc, dil


# dispense
def dispense(i, t_pump1, t_pump2):
    
    time.sleep(0.2)
    
    if(t_pump1>0.1):
        DAQC2.setDOUTbit(0,i)
        time.sleep(0.1)
        DAQC2.setDOUTbit(0,6)
        time.sleep(t_pump1)
        DAQC2.clrDOUTbit(0,6)
        time.sleep(0.1)
        DAQC2.clrDOUTbit(0,i)
        time.sleep(0.2)
    
    if(t_pump2>0.1):
        DAQC2.setDOUTbit(1,i)
        time.sleep(0.1)
        DAQC2.setDOUTbit(1,6)
        time.sleep(t_pump2)
        DAQC2.clrDOUTbit(1,6)
        time.sleep(0.1)
        DAQC2.clrDOUTbit(1,i)
        time.sleep(0.2)
    
# check for pause
def check_pause():
    #print('\nenter pause / resume\n')
    global checkpause
    while True:
        while checkpause==False:
            var=input('\n\nenter pause: ')
            if(var=='pause'):
                checkpause=True

# controller
def controller(i):

    global concentrations, time_init, times_last, OD_last, time_last, peakdose, delay, flagpeak, time_next, conc

    if(mode=='growthcurve'):
        t_pump2 = 0
        t_pump1 = 0
    
    elif(mode=='chemostat'):
        t_pump2=0
        t_pump1 = cycle * dil_fixed / rate
    
    elif(mode=='turbidostat'):
        t_pump2=0
        t_pump1 = Kp * (OD[i] - OD_target) + Ki * np.sum( (times[1:]-times[0:-1]) * (OD_past[i,1:] - OD_target) )
        if(t_pump1<0):
            t_pump1=0
    
    elif(mode=='morbidostat'):
        n = np.argmax( times > times[-1]-600 )
        t_pump = cycle * dil_fixed / rate
        if OD[i] > OD_target and OD_past[i,-1] > OD_past[i,n] and times[-1] > 600:
            t_pump2 = drugdose / rate
            t_pump1 = max( 0 , t_pump - t_pump2)
        else:
            t_pump2=0
            t_pump1=t_pump

    elif (mode == 'mixmode'):
        if i % 2 == 1 :
            if time.time() > time_start + 10800: 
                t_pump2=0
                t_pump1 = Kp * (OD[i] - OD_target) + Ki * np.sum( (times[1:]-times[0:-1]) * (OD_past[i,1:] - OD_target) )
                if(t_pump1<0):
                    t_pump1=0
            else:                    
                t_pump2=0
                t_pump1 = Kp * (OD[i] - OD_target)
                if(t_pump1<0):
                    t_pump1=0

        else:  
            n = np.argmax( times > times[-1]-600 )
            t_pump = cycle * dil_fixed / rate
            if OD[i] > OD_target and OD_past[i,-1] > OD_past[i,n] and times[-1] > 600:
                t_pump2 = drugdose / rate
                t_pump1 = max( 0 , t_pump - t_pump2)
            else:
                t_pump2=0
                t_pump1=t_pump

    elif(mode=='morbgradient'):
        t_pump = cycle * dilution_rate_gradient[i] / rate
        n1 = np.where( times > max(0,times[-1]-3600) )
        n2 = np.where( OD_past[i,:] > 0 )
        n = np.intersect1d(n1,n2)
        if len(n)>1:
            times_last = np.squeeze(np.take(times, n))/3600
            OD_last = np.squeeze(np.log2(np.take(OD_past[i,:],n)))
            coefs = np.polyfit(times_last, OD_last, 1)
            slope = coefs[0]
        else:
            slope = 0

        if OD[i] > od_gradient[i] and slope > 0 :
            t_pump2 = drugdose / rate
            t_pump1 = t_pump - t_pump2
        else :
            t_pump2=0
            t_pump1=t_pump

    elif (mode == 'shocksteady'):
        if i % 2 == 1 :
            if time.time() > time_next[i]:
                if flagpeak[i] and conc[i] > 0.01:
                    t_pump2 = 0
                else:
                    peakdose[i] = peakdose[i] * max( 1 , np.sqrt(delay0/delay[i]) ) 
                    t_pump2 = peakdose[i] / rate
                t_pump1 = max( 0, 1/rate - t_pump2 )
                conc[i] = conc[i] * np.exp( - rate * t_pump1)
                conc[i] = 1 + (conc[i] - 1) * np.exp( - rate * t_pump2)
                delay[i] = time.time()
                flagpeak[i]=True
                time_next[i] = time_next[i] + shock_cycle * 3600
                print(delay)
                print(peakdose)
            else:                    
                t_pump2=0
                t_pump1 = max( 0 , Kp * (OD[i] - OD_target) )
                conc[i] = conc[i] * np.exp( - rate * t_pump1)
                if OD[i] > OD_target and flagpeak[i] :
                    delay[i] = time.time() - delay[i]
                    print(delay)
                    flagpeak[i] = False

        else:  
            n = np.argmax( times > times[-1]-1800 )
            t_pump = cycle * dil_fixed / rate
            if OD[i] > OD_target and OD_past[i,-1] > OD_past[i,n] and times[-1] > 600:
                t_pump2 = drugdose / rate
                t_pump1 = max( 0 , t_pump - t_pump2)
            else:
                t_pump2=0
                t_pump1=t_pump
    
    # put limits on t_pump
    t_tot=(t_pump1+t_pump2)
    if(t_tot>210):
        t_pump1=t_pump1*(210/t_tot)
        t_pump2=t_pump2*(210/t_tot)
    if(OD[i]>1):
        t_pump2 = 0
        t_pump1 = 0
    
    return t_pump1, t_pump2

# direct control
def direct():
    print("\n\ncommand: <component> <number> <on/off> (or exit)")
    print("\ncomponent: pinch, pump, laser, photo")
    print("\nnumber: optional, if empty returns all\n")

    condition = 1
    while(condition):
        cmd = input("\n: ")
        cmd = cmd.split(' ')
        try: 

            if(cmd[0] == 'exit'):
                condition = 0
                DAQC2.setDOUTall(0,0)
                DAQC2.setDOUTall(1,0)
                print('\n\nexiting direct control\n\n')
                break
            
            elif(cmd[0] == 'laser'): 
                if(cmd[1]=='on'):
                    DAQC2.setDOUTbit(1,7)
                elif(cmd[1]=='off'):
                    DAQC2.clrDOUTbit(1,7)
            
            elif(cmd[0] == 'pump'):
                if(cmd[1]=='on'):
                    for i in [0,1]:
                        DAQC2.setDOUTbit(i,6)           
                elif(cmd[1]=='off'):
                    for i in [0,1]:
                        DAQC2.clrDOUTbit(i,6)
                else:
                    j=int(cmd[1])-1
                    if(cmd[2]=='on'):
                        DAQC2.setDOUTbit(j,6)
                    elif(cmd[2]=='off'):
                        DAQC2.clrDOUTbit(j,6)
            
            elif(cmd[0] == 'pinch'):
                if(cmd[1]=='on'):
                    for i in range(6):
                        DAQC2.setDOUTbit(0,i)
                        DAQC2.setDOUTbit(1,i)           
                elif(cmd[1]=='off'):
                    for i in range(6):
                        DAQC2.clrDOUTbit(0,i)
                        DAQC2.clrDOUTbit(1,i)
                else:
                    j = (int(cmd[1])-1) // 6
                    k = (int(cmd[1])-1) % 6           
                    if(cmd[2]=='on'):
                        DAQC2.setDOUTbit(j,k)
                    elif(cmd[2]=='off'):
                        DAQC2.clrDOUTbit(j,k)
                
            elif(cmd[0] == 'photo'):
                if(len(cmd)==1):
                    v=DAQC2.getADCall(0)
                    v+=DAQC2.getADCall(1)
                    print("\nvoltages:\n")
                    for i in range(6):
                        print(v[i])
                    for i in range(6):
                        print(v[8+i])
                else:
                    j = (int(cmd[1])-1) // 6
                    k = (int(cmd[1])-1) % 6           
                    v=DAQC2.getADC(j,k)
                    print("\nvoltage:\n",v)                
        
        except ValueError: 
            print("nope")

# prime and fill
def prime():
    print('\n\nEnter culture number or keep pressing enter to start/stop the pump. exit to leave. \n\n')
    
    flag=False

    condition = 1
    while(condition):
        cmd = input("\n: ")

        try: 
            
            if(cmd == 'exit'):
                condition = 0
                DAQC2.setDOUTall(0,0)
                DAQC2.setDOUTall(1,0)
                print('\n\npriming complete\n\n')
                break
                
            elif(cmd==''):
                if(flag==False):
                    DAQC2.setDOUTbit(0,6)
                    DAQC2.setDOUTbit(1,6)
                    print('pumps on')
                if(flag==True):
                    DAQC2.clrDOUTbit(0,6)
                    DAQC2.clrDOUTbit(1,6)
                    print('pumps off')
                flag=not(flag)

            elif(int(cmd)-1 in cultures):
                DAQC2.setDOUTall(0,0)
                DAQC2.setDOUTall(1,0)
                flag=False
                i=int(cmd)-1
                print('culture %s' % str(i+1))         
                DAQC2.setDOUTbit(0,i)
                DAQC2.setDOUTbit(1,i) 

        except ValueError: 
            print("nope")
    
    DAQC2.setDOUTall(0,0)
    DAQC2.setDOUTall(1,0)

def fill():    
    for i in cultures:
        DAQC2.setDOUTbit(0,i)
        time.sleep(0.2)
        DAQC2.setDOUTbit(0,6)
        time.sleep(210)
        DAQC2.clrDOUTbit(0,6)
        time.sleep(0.2)
        DAQC2.clrDOUTbit(0,i)
        
    DAQC2.setDOUTall(0,0)
    DAQC2.setDOUTall(1,0)
    print('\n\nfilling complete\n\n')

def clean():

    print('\n\nrunning alcohol\n\n')    
    for i in cultures:
        DAQC2.setDOUTbit(0,i)
        DAQC2.setDOUTbit(1,i)
        time.sleep(0.2)
        DAQC2.setDOUTbit(0,6)
        DAQC2.setDOUTbit(1,6)
        time.sleep(30)
        DAQC2.clrDOUTbit(0,6)
        DAQC2.clrDOUTbit(1,6)
        time.sleep(0.2)
        DAQC2.clrDOUTbit(0,i)
        DAQC2.clrDOUTbit(1,i)

    input("\n\npress enter to run water\n\n")
    print('\n\nrunning water\n\n')    
    for i in cultures:
        DAQC2.setDOUTbit(0,i)
        DAQC2.setDOUTbit(1,i)
        time.sleep(0.2)
        DAQC2.setDOUTbit(0,6)
        DAQC2.setDOUTbit(1,6)
        time.sleep(100)
        DAQC2.clrDOUTbit(0,6)
        DAQC2.clrDOUTbit(1,6)
        time.sleep(0.2)
        DAQC2.clrDOUTbit(0,i)
        DAQC2.clrDOUTbit(1,i)

    input("\n\npress enter to run air\n\n")
    print('\n\nrunning air\n\n')    
    for i in cultures:
        DAQC2.setDOUTbit(0,i)
        DAQC2.setDOUTbit(1,i)
        time.sleep(0.2)
        DAQC2.setDOUTbit(0,6)
        DAQC2.setDOUTbit(1,6)
        time.sleep(210)
        DAQC2.clrDOUTbit(0,6)
        DAQC2.clrDOUTbit(1,6)
        time.sleep(0.2)
        DAQC2.clrDOUTbit(0,i)
        DAQC2.clrDOUTbit(1,i)

    DAQC2.setDOUTall(0,0)
    DAQC2.setDOUTall(1,0)
    print('\n\ncleaning complete\n\n')

# measure blanks
def measure_blanks():
    DAQC2.setDOUTbit(1,7)
    time.sleep(0.5)
    for i in cultures:
        tempVolt = []
        tEnd = time.time() + 1
        while(time.time() < tEnd):
            tempVolt.append(DAQC2.getADC(0,i))
        blanks[i] = sum(tempVolt) / len(tempVolt)
    DAQC2.clrDOUTbit(1,7)
    s=""
    for item in blanks:
        s+=str("%.4f" % item)+" "
    with open("INI.txt","r+") as file:
        ini = file.read()
        ini = ini.split('\n')
        ini[7] = s[0:-1]
        file.seek(0)
        for lines in ini:
            file.write(lines)
            file.write("\n")
    print('\n\nblanks recorded\n\n')


def display_OD():
    od = read_OD()
    print("\nOD:\n")
    for i in range(6):
        print("%.4f" % od[i]) 

# measure, update past values and record
def measure_update_record(name):
    global OD, concentrations, dilutions, time_init, times, OD_past
    OD = read_OD()
    t=(time.time()-time_init)
    write_log(name,t/60,OD,concentrations,dilutions)
    concentrations = np.zeros(6)
    dilutions = np.zeros(6)
    times[0:-1]=times[1:]
    times[-1]=t
    OD_past[:,0:-1]=OD_past[:,1:]
    OD_past[:,-1]=OD
    t_last=time.time()
    return t_last

# run device
def run():

    print('\n\nStarting %s run at %s\n\n' % (mode , datetime.now().strftime('%H:%M:%S')))

    name="LOG-"+datetime.now().strftime('%Y-%m-%d_%H_%M')+".csv"
    with open(name,"a") as file:
        file.write("%s run on %s" % (mode, datetime.now().strftime('%Y/%m/%d at %H:%M:%S')))

    global OD, concentrations, dilutions, time_init, times, OD_past, checkpause, time_next, time_start, delay, flagpeak, peakdose
    time_init = time.time()
    time_limit = time.time() + experiment_duration * 3600
    
    if mode=='growthcurve':
        while time.time() < time_limit:
            t_last=measure_update_record(name)
            time.sleep(OD_interval)
    else:       
        t_last=measure_update_record(name)
        while np.mean([OD[j] for j in cultures]) < OD_target :
            time.sleep(OD_interval)
            t_last=measure_update_record(name)

        print('OD good to go')
        time_start = time.time()
        time_next = time_next + time.time() + 0 * 3600

        thread_pause = threading.Thread(target=check_pause)
        thread_pause.start()

        while time.time() < time_limit:
        
            t_cycle=time.time()
            t_last=measure_update_record(name)
            for i in cultures:
                t_pump1, t_pump2 = controller(i)
                if t_pump1+t_pump2 > 0:
                    dispense(i,t_pump1,t_pump2)
                    concentrations[i], dilutions[i] = calculate(concentrations[i],t_pump1,t_pump2)
                    t_last=measure_update_record(name)

            t_left=cycle-time.time()+t_cycle
            while t_left > OD_interval:
                time.sleep(OD_interval)
                t_last=measure_update_record(name)
                t_left=cycle-time.time()+t_cycle
            
            if t_left>0:
                time.sleep(t_left)

            if(checkpause==True):
                print('\nwaiting\n')
                direct()
                checkpause=False
                print('\nresuming\n')
    
    thread_pause.join()
    DAQC2.setDOUTall(0,0)
    DAQC2.setDOUTall(1,0)





print("\n\nAll parameters are in the INI file. All results are written to the LOG file.")

condition = 1
while(condition):

    # read INI
    with open("INI.txt","r") as file:
        ini = file.read()
    ini = ini.split('\n')
    mode = ini[1]
    ss = ini[4]
    cultures = [int(s)-1 for s in ss.split(' ')]
    ss = ini[7]
    blanks = [float(s) for s in ss.split(' ')]
    rate = float(ini[10])    
    experiment_duration = float(ini[13])
    OD_interval = float(ini[16])*60
    cycle = float(ini[19])*60
    dil_fixed = float(ini[22])/3600
    drugdose = float(ini[25])
    OD_target = float(ini[28])
    Kp = float(ini[31])
    Ki = float(ini[32])
    Kd = float(ini[33])
    ss = ini[36]
    od_gradient = [float(s) for s in ss.split(' ')]
    ss = ini[39]
    dilution_rate_gradient = [float(s) for s in ss.split(' ')]
    peakdose = float(ini[42])*np.ones(6)
    shock_cycle = float(ini[45])
    delay0 = float(ini[48])*3600
    ss = ini[51]
    conc = [float(s) for s in ss.split(' ')]
    
    time_init = time.time()
    time_start = 0
    time_next = np.zeros(6)
    delay = delay0*np.ones(6)
    flagpeak = np.ones(6)<0
    concentrations = np.zeros(6)
    dilutions = np.zeros(6)
    OD = np.zeros(6)

    times=np.zeros(180)
    OD_past=np.zeros((6,180))

    OD = read_OD()

    checkpause = False

    print("\n")
    print("0 quit")
    print("1 direct control")
    print("2 measure OD")
    print("3 prime")
    print("4 fill")
    print("5 clean")
    print("6 measure blanks")
    print("7 RUN")

    var = input(": ")
        
    try: 
        currInput = int(var)
            
        if(currInput == 0):
            print("\n\nbye")
            DAQC2.setDOUTall(0,0)
            DAQC2.setDOUTall(1,0)
            condition = 0
            break
            
        elif(currInput == 1): 
            direct()

        elif(currInput == 2):
            display_OD()
            
        elif(currInput == 3): 
            prime()

        elif(currInput == 4): 
            fill()

        elif(currInput == 5): 
            clean()
            
        elif(currInput == 6):
            measure_blanks()
                
        elif(currInput == 7):
            run()


    except ValueError:
        DAQC2.setDOUTall(0,0)
        DAQC2.setDOUTall(1,0) 
        print("\nnope")

