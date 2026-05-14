clc; close all; clear;
addpath("finalProjectFunctions\");
%% Top Level Variables
satelliteName = "Satellite1";
ephemerisFolder = pwd; % defaults to store files in folder where script is run
dataTimeStep_sec = 1;

% Connect to STK
app = actxGetRunningServer('STK.Application');
root = app.Personality2;
scenario = root.CurrentScenario;
root.UnitPreferences.SetCurrentUnit("DateFormat", "EpSec");
root.UnitPreferences.SetCurrentUnit("Distance", "m");

%% Load in hypersonic platform
hypersonicEphemerisFile = pwd + "\HXRV_X43_parsed_data_ECEF.e";
hypersonicName = "HXRV_X43";
hypersonicObjectType = "Aircraft"; % Alternative is "Missile";
obj_enumeration = char("e" + hypersonicObjectType);
object_name = hypersonicName;

if ~scenario.Children.GetElements(obj_enumeration).Contains(object_name)
    hypersonic = scenario.Children.New(obj_enumeration, object_name);
else
    hypersonic = scenario.Children.GetElements(obj_enumeration).Item(object_name);
end
hypersonic.SetRouteType("ePropagatorStkExternal");
hypersonic.Route.Filename = hypersonicEphemerisFile;
hypersonic.Route.Override = 1;
hypersonic.Route.Propagate;

constraint_type = "Duration";
constraint_enumeration = strcat("eCstr", string(constraint_type));
if hypersonic.AccessConstraints.IsConstraintActive(constraint_enumeration)
    access_constraint = hypersonic.AccessConstraints.GetActiveConstraint(constraint_enumeration);
else
    access_constraint = hypersonic.AccessConstraints.AddConstraint(constraint_enumeration);
end

access_constraint.EnableMin = 1;
access_constraint.Min = 420;
disp("Defined Hypersonic Vehicle");

%% Get handle for Satellite and Hypersonic
satCount = scenario.Children.GetElements('eSatellite').Count;
allData = [];
accessCount = 1;
satNames = {};
currFrame = "Fixed";

for satIdx = 1:satCount
    sat = scenario.Children.GetElements('eSatellite').Item(int32(satIdx - 1));

    %if sat.InstanceName ~= "CPN"
    % Calculate access
    accessFileName = ephemerisFolder + "\" + satelliteName + "_to_" + hypersonicName + "_access_data.csv";
    access = sat.GetAccessToObject(hypersonic);
    aerDP = access.DataProviders.Item("AER Data").Group.Item("NorthEastDown").Exec(scenario.StartTime, scenario.StopTime, dataTimeStep_sec);

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
        speedValues = cell2mat(hypersonic.DataProviders.Item("Cartesian Velocity").Group.Item(currFrame).Exec(scenario.StartTime, scenario.StopTime, dataTimeStep_sec).DataSets.GetDataSetByName("speed").GetValues);

        dataN = length(timeValues);
        n = length(times);

        times(n+1:n+dataN) = timeValues;
        azimuths(n+1:n+dataN) = azimuthValues;
        elevations(n+1:n+dataN) = elevationValues;
        ranges(n+1:n+dataN) = rangeValues;

        posDP = sat.DataProviders.Item("Cartesian Position").Group.Item(currFrame).Exec(times(1), times(end), dataTimeStep_sec);
        velDP = sat.DataProviders.Item("Cartesian Velocity").Group.Item(currFrame).Exec(times(1), times(end), dataTimeStep_sec);
        stateDP = sat.DataProviders.Item("LLR State").Group.Item(currFrame).Exec(times(1), times(end), dataTimeStep_sec);

        posx = cell2mat(posDP.DataSets.GetDataSetByName("x").GetValues);
        posy = cell2mat(posDP.DataSets.GetDataSetByName("y").GetValues);
        posz = cell2mat(posDP.DataSets.GetDataSetByName("z").GetValues);

        rad = cell2mat(stateDP.DataSets.GetDataSetByName("Rad").GetValues);
        lon = cell2mat(stateDP.DataSets.GetDataSetByName("Lon").GetValues);
        lat = cell2mat(stateDP.DataSets.GetDataSetByName("Lat").GetValues);

        velx = cell2mat(velDP.DataSets.GetDataSetByName("x").GetValues);
        vely = cell2mat(velDP.DataSets.GetDataSetByName("y").GetValues);
        velz = cell2mat(velDP.DataSets.GetDataSetByName("z").GetValues);

        %satPositionData(:,:,accessCount) = [posx, posy, posz, velx, vely, velz];
        satPositionData(:,:,accessCount) = [rad, deg2rad(lon), deg2rad(lat), vecnorm([velx,vely,velz],2,2), zeros(length(rad),1), zeros(length(rad),1)];
        %end

        if aerDP.Intervals.Count > 0
            allData(:,:,accessCount) = [azimuths; elevations; ranges; speedValues'];
            satNames = horzcat(satNames, sat.InstanceName);
            accessCount = accessCount + 1;
        end
    end
end
disp("Stored AER Data...");

%% Store Satellite and Hypersonic Vehicle Ephemerides
frames = ["Fixed"];
platforms = {hypersonic};

for platIdx = 1:length(platforms)
    currPlat = platforms{platIdx};

    for frameIdx = 1:length(frames)
        currFrame = frames(frameIdx);
        posVelFileName = ephemerisFolder + "\" + currPlat.InstanceName + "_position_velocity_" + currFrame + ".csv";
        posDP = currPlat.DataProviders.Item("Cartesian Position").Group.Item(currFrame).Exec(scenario.StartTime, scenario.StopTime, dataTimeStep_sec);
        velDP = currPlat.DataProviders.Item("Cartesian Velocity").Group.Item(currFrame).Exec(scenario.StartTime, scenario.StopTime, dataTimeStep_sec);
        stateDP = currPlat.DataProviders.Item("LLR State").Group.Item(currFrame).Exec(scenario.StartTime, scenario.StopTime, dataTimeStep_sec);


        posx = cell2mat(posDP.DataSets.GetDataSetByName("x").GetValues);
        posy = cell2mat(posDP.DataSets.GetDataSetByName("y").GetValues);
        posz = cell2mat(posDP.DataSets.GetDataSetByName("z").GetValues);

        rad = cell2mat(stateDP.DataSets.GetDataSetByName("Rad").GetValues);
        lon = cell2mat(stateDP.DataSets.GetDataSetByName("Lon").GetValues);
        lat = cell2mat(stateDP.DataSets.GetDataSetByName("Lat").GetValues);

        velx = cell2mat(velDP.DataSets.GetDataSetByName("x").GetValues);
        vely = cell2mat(velDP.DataSets.GetDataSetByName("y").GetValues);
        velz = cell2mat(velDP.DataSets.GetDataSetByName("z").GetValues);


        for i = 1:length(posx)
            cartesianState = [posx(i), posy(i), posz(i), velx(i), vely(i), velz(i)];
            full_state(:,i) = return_relevant_state_vector(cartesianState);
        end

        hypersonicPositionData = [rad'; deg2rad(lon)'; deg2rad(lat)'; full_state(4,:); full_state(5,:); full_state(6,:)];
        disp("Stored ephemerides for " + currPlat.InstanceName + " in " + currFrame + " frame");
    end
end

save('measurement_data.mat', 'allData', 'satPositionData', 'hypersonicPositionData', 'satNames');
