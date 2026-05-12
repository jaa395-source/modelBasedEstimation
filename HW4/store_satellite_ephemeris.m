clc; close all; clear;
%% Top Level Variables
satelliteName = "Satellite1";
hypersonicName = "HXRV_X43";
hypersonicObjectType = "Aircraft"; % Alternative is "Missile";
ephemerisFolder = pwd; % defaults to store files in folder where script is run
dataTimeStep_sec = 1;

%% Connect to STK
app = actxGetRunningServer('STK.Application');
root = app.Personality2;
scenario = root.CurrentScenario;
root.UnitPreferences.SetCurrentUnit("DateFormat", "EpSec");

%% Get handle for Satellite and Hypersonic
sat = scenario.Children.GetElements('eSatellite').Item(satelliteName);
hypersonic = scenario.Children.GetElements("e" + char(hypersonicObjectType)).Item(hypersonicName);

%% Calculate access
accessFileName = ephemerisFolder + "\" + satelliteName + "_to_" + hypersonicName + "_access_data.csv";
access = sat.GetAccessToObject(hypersonic);
aerDP = access.DataProviders.Item("AER Data").Group.Item("Default").Exec(scenario.StartTime, scenario.StopTime, dataTimeStep_sec);

times = [];
azimuths = [];
elevations = [];
ranges = [];

for intIdx = 1:aerDP.Intervals.Count
    currInt = aerDP.Intervals.Item(int32(intIdx - 1));
    timeValues = cell2mat(currInt.DataSets.GetDataSetByName("Time").GetValues);
    azimuthValues = cell2mat(currInt.DataSets.GetDataSetByName("Azimuth").GetValues);
    elevationValues = cell2mat(currInt.DataSets.GetDataSetByName("Elevation").GetValues);
    rangeValues = cell2mat(currInt.DataSets.GetDataSetByName("Range").GetValues);

    dataN = length(timeValues);
    n = length(times);

    times(n+1:n+dataN) = timeValues;
    azimuths(n+1:n+dataN) = azimuthValues;
    elevations(n+1:n+dataN) = elevationValues;
    ranges(n+1:n+dataN) = rangeValues;
end

dataTable = table(times, azimuths, elevations, ranges);
writetable(dataTable, accessFileName);
disp("Stored AER Data...");

%% Store Satellite and Hypersonic Vehicle Ephemerides
frames = ["J2000", "Fixed", "ICRF"];
platforms = {sat, hypersonic};

for platIdx = 1:length(platforms)
    currPlat = platforms{platIdx};

    for frameIdx = 1:length(frames)
        currFrame = frames(frameIdx);
        posVelFileName = ephemerisFolder + "\" + currPlat.InstanceName + "_position_velocity_" + currFrame + ".csv";
        posDP = currPlat.DataProviders.Item("Cartesian Position").Group.Item(currFrame).Exec(scenario.StartTime, scenario.StopTime, dataTimeStep_sec);
        velDP = currPlat.DataProviders.Item("Cartesian Velocity").Group.Item(currFrame).Exec(scenario.StartTime, scenario.StopTime, dataTimeStep_sec);

        posx = cell2mat(posDP.DataSets.GetDataSetByName("x").GetValues);
        posy = cell2mat(posDP.DataSets.GetDataSetByName("y").GetValues);
        posz = cell2mat(posDP.DataSets.GetDataSetByName("z").GetValues);

        velx = cell2mat(velDP.DataSets.GetDataSetByName("x").GetValues);
        vely = cell2mat(velDP.DataSets.GetDataSetByName("y").GetValues);
        velz = cell2mat(velDP.DataSets.GetDataSetByName("z").GetValues);
        
        dataTable = table(posx, posy, posz, velx, vely, velz);
        writetable(dataTable, accessFileName);
        disp("Stored ephemerides for " + currPlat.InstanceName + " in " + currFrame + " frame");
    end
end