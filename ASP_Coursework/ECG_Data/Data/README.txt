The code necessary to segment the raw ECG signal ('RAW-ECG.mat') and convert it into RRI intervals should look similar to:

    load RAW-ECG.mat
    Trial_1 = data(770:2.5e5);
    Trial_2 = data(2.564e5:4.99e5);
    Trial_3 = data(5.228e5:7.473e5);
    [xRRI1,fsRRI1]=ECG_to_RRI(Trial_1, fs);
    [xRRI2,fsRRI2]=ECG_to_RRI(Trial_2, fs);
    [xRRI3,fsRRI3]=ECG_to_RRI(Trial_3, fs);

Alternatively, directly load the 'RRI-DATA.mat' file into MATLAB.