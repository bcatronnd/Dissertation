close all; clc; clearvars;

RunLog.TestPoint = [20210831001 20210831002 20210831003 20210831004 20210831005 20210831006 20210831007 20210831008 20210831009 20210831010 20210831011 20210831012 20210831013 20210831014 20210831015 20210831016 20210831017 20210831018 20210831019 20210831020 20210831021 20210831022 20210831023 20210901001 20210901002 20210901003 20210901004 20210901005 20210901006 20210901007 20210901008 20210901009 20210901010 20210901011 20210901012 20210901013 20210901014 20210901015 20210901016 20210901017 20210901018 20210901019 20210901020 20210901021 20210901022 20210901023];
RunLog.CameraData = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
RunLog.DaqData = [1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
RunLog.Geometry = {'Cavity - 12.75"' 'Cavity - 12.75"' 'Cavity - 12.75"' 'Cavity - 12.75"' 'Cavity - 9"' 'Cavity - 9"' 'Cavity - 9"' 'Flush' 'Flush' 'Flush' 'Flush' 'Flush' 'Flush' 'Flush' 'Cavity - 12.75"' 'Cavity - 12.75"' 'Cavity - 12.75"' 'Cavity - 9"' 'Cavity - 9"' 'Cavity - 9"' 'Flush' 'Flush' 'Flush' 'Flush' 'Flush' 'Flush' 'Flush' 'Flush' 'Flush' 'Flush' 'Cavity - 12.75"' 'Cavity - 12.75"' 'Cavity - 12.75"' 'Cavity - 9"' 'Cavity - 9"' 'Cavity - 9"' 'Cavity - 9"' 'Cavity - 9"' 'Cavity - 9"' 'Cavity - 12.75"' 'Cavity - 12.75"' 'Cavity - 12.75"' 'Flush' 'Flush' 'Flush' 'Flush'};
RunLog.Angle = [90 90 90 90 90 90 90 90 90 90 90 105 105 105 105 105 105 105 105 105 105 105 105 90 90 90 90 75 75 75 75 75 75 75 75 75 120 120 120 120 120 120 120 120 120 120];
RunLog.AngleUnit = 'Degree';
RunLog.Length = [12.75 12.75 12.75 12.75 9 9 9 NaN NaN NaN NaN NaN NaN NaN 12.75 12.75 12.75 9 9 9 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 12.75 12.75 12.75 9 9 9 9 9 9 12.75 12.75 12.75 NaN NaN NaN NaN];
RunLog.LengthUnit = 'Inch';
RunLog.Depth = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 NaN NaN NaN NaN NaN NaN NaN 0.5 0.5 0.5 0.5 0.5 0.5 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 NaN NaN NaN NaN];
RunLog.DepthUnit = 'Inch';
RunLog.Mach = [0 0.3 0.4 0.5 0.3 0.4 0.5 0 0.3 0.4 0.5 0.3 0.4 0.5 0.3 0.4 0.5 0.3 0.4 0.5 0.3 0.4 0.5 0.3 0.4 0.5 0 0.3 0.4 0.5 0.3 0.4 0.5 0.3 0.4 0.5 0.3 0.4 0.5 0.3 0.4 0.5 0.3 0.4 0.5 0];
RunLog.FanDriveHertz = [NaN 30.6 39.8 48.8 30.6 39.8 48.8 NaN 30.6 39.8 48.8 30.6 39.8 48.8 30.6 39.8 48.8 30.6 39.8 48.8 30.6 39.8 48.8 30.6 39.7 48.6 NaN 30.6 39.7 48.6 30.6 39.7 48.6 30.6 39.7 48.6 30.6 39.7 48.6 30.6 39.7 48.6 30.6 39.7 48.6 NaN];
RunLog.PressureStatic = [NaN 92755 88443 83322 92690 88400 83130 NaN 92722 88470 83222 92720 88455 83101 92666 88366 83134 92644 88344 83255 92600 88304 83162 92960 88670 83320 NaN 92960 88680 83480 92980 88622 83530 92944 88650 83570 92930 88690 83415 92966 88628 83565 92920 88644 83455 NaN];
RunLog.PressureStaticUnit = 'Pascal';
RunLog.PressureTotal = [NaN 98728 98720 98720 98636 98677 98690 NaN 98676 98740 98770 98650 98690 98700 98577 98600 98634 98532 98555 98610 98502 98540 98598 98950 98960 98960 NaN 98900 98930 98999 98900 98911 98998 98888 98910 98995 98902 98944 98952 98900 98905 98975 98886 98922 98944 NaN];
RunLog.PressureTotalUnit = 'Pascal';
RunLog.TemperatureStart = [NaN 23.9 26.5 27.4 30.1 30.1 31 NaN 25.4 25.6 26.5 25.3 25.8 27.1 25.9 26.4 27.5 28.1 28.5 29.8 29.2 29.5 30.6 24.3 24.9 26.1 NaN 26.4 27 28 27.9 28.4 29.6 29.5 29.9 30.8 28.5 28.7 29.8 29.7 30 31 29.9 30 31 NaN];
RunLog.TemperatureEnd = [NaN 25 27.4 30.2 30.1 31 32.9 NaN 25.6 26.5 28.4 25.8 27.1 28.8 26.4 27.5 29.5 28.5 29.8 31.6 29.5 30.6 32.4 24.9 26.1 28.2 NaN 27 28 29.8 28.4 29.6 31.4 29.9 30.8 32.4 28.7 29.8 31.7 30 31 32.7 30 31 32.6 NaN];
RunLog.TemperatureUnit = 'Celcius';
RunLog.Humidity = [NaN 63 52 38 43 42 37 NaN 70 64 54 69 62 47 62 53 42 48 46 38 48 42 34 64 55 44 NaN 56 50 41 50 47 33 43 40 32 50 44 38 42 37 32 43 40 33 NaN];
RunLog.HumidityUnit = 'Percent';
RunLog.PressureAmbient = [NaN 985 985 985 985 985 985 NaN 985 985 985 985 985 985 985 985 985 985 985 985 985 985 985 987 987 987 NaN 987 987 987 987 987 987 987 987 987 987 987 987 987 987 987 987 987 987 NaN];
RunLog.PressureAmbientUnit = 'Hectopascal';
RunLog.TunnelRPMHertz = 1197/60;
RunLog.TunnelBladeNumber = 32;

RunLog.CalculatedQuantities = [];
RunLog.MachCalc = real(sqrt(5*((RunLog.PressureTotal./RunLog.PressureStatic).^(2/7)-1)));
if strcmpi(RunLog.TemperatureUnit,'Celcius')
    RunLog.SpeedOfSound = sqrt(1.4*287*(mean([RunLog.TemperatureStart; RunLog.TemperatureEnd],1)+273.15));
end
RunLog.VelocityFreeStream = RunLog.MachCalc.*RunLog.SpeedOfSound;
RunLog.BladePassingFrequency = RunLog.TunnelRPMHertz*RunLog.FanDriveHertz*RunLog.TunnelBladeNumber/60;


save('RunLog.mat');