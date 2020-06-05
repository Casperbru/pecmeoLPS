% clear all 
% close all
% clc 
% fig = 0;

%% Designation of maps
mainMap = 'C:\tudatBundle\tudatApplications\stepTwo\SimulationOutput\';
timeMap = '\PECMEO_30sec\';
% constellationMap = '\GalileoGLONASSBeiDouGPSISL\';
constellationMap = '\GalileoGLONASSBeiDouGPS\';
estimationMap = '\estimationMap\';
propagationMap = '\propagationResults\';

numberOfPECMEOSats = 9;

%% Loading data
PECMEOReal = cell(numberOfPECMEOSats,1);
PECMEO = cell(numberOfPECMEOSats,1);
% PECMEOError = cell(numberOfPECMEOSats,1);

SPPLow = 5; SPPLim = 5;
KODLow = 2; KODLim = 2;
SPPISLLow = 0; SPPISLLim = 0;
KODISLLow = 0; KODISLLim = 2;

periodMap = '\1d\';
% periodMap = '\4h\';
% periodMap = '\10min\';
% systemMaps = {'\NoISL\','\GNSSISL\','\GNSSISL_1C\','\LaserISL\','\LaserISL_1C\'};
errorMaps = {'\error_none\','\error_NE\','\error_EE\','\error_NE_EE\'};
systemMaps = {'\NoISL\','\GNSSISL\','\GNSSISL_1C\','\LaserISL\','\LaserISL_1C\','\KbandISL\','\KbandISL_1C\'};
systems = {' Non ISL 9C', 'GPS ISL 9C', 'GPS ISL 1C', '  LRI ISL 9C', '  LRI ISL 1C', 'KBR ISL 9C', 'KBR ISL 1C'};
errors = {'None','Noise','Ephemeris','Combined'};

% errorMaps{1,2} = '\error_none_test2\';
% errorMaps{1,2} = '\error_NE_CovRun\';
% errorMaps{1,2} = '\error_NE_rotate\';
% errorMaps{1,5} = '\error_EE_SPP\';
% errorMaps{1,3} = '\error_NE_EE_SPP\';
% errorMaps{1,5} = '\error_EE_KOD\';
% errorMaps{1,5} = '\error_NE_EE_1mEph\';
% systemMaps{1,6} = '\30KbandISL\';
% systemMaps{1,7} = '\30KbandISL_1C\';

numberOfSystems = length(systemMaps);
numberOfErrors = length(errorMaps);
PECMEOError = cell(numberOfSystems,numberOfErrors);

% for system = 1:1
for system = 1:numberOfSystems
% 	for error = 1:1
    for error = 1:numberOfErrors
        systemMap = systemMaps{system};
        errorMap = errorMaps{error};
        for sat = 1:numberOfPECMEOSats
            name = ['PECMEOSat' num2str(sat)];
            PECMEOReal(sat) = {importdata([mainMap timeMap constellationMap periodMap propagationMap name '.dat'])};
%             PECMEOReal(sat) = {importdata([mainMap timeMap constellationMap periodMap '\propagationResultsRot\' name '.dat'])};
            for PPIt = 10:10
                for KinIt = 10:10
                    for PPISLIt = 0:0
                        for KinISLIt = 0:0
%                             PECMEO{sat}{PPIt,KinIt+1}(PPISLIt+1,KinISLIt+1) = {importdata([mainMap timeMap constellationMap periodMap systemMap errorMap estimationMap name '_' num2str(PPIt) ',' num2str(KinIt) ',' num2str(PPISLIt) ',' num2str(KinISLIt) '.dat'])};
                            PECMEO{sat}{PPIt,KinIt+1}(PPISLIt+1,KinISLIt+1) = {importdata([mainMap timeMap constellationMap periodMap systemMap errorMap estimationMap name '.dat'])};
                            PECMEOError{system,error}{PPIt,KinIt+1}{PPISLIt+1,KinISLIt+1}(:,sat) = vecnorm(PECMEOReal{sat}(:,2:4)-PECMEO{sat}{PPIt,KinIt+1}{PPISLIt+1,KinISLIt+1}(:,2:4),2,2);
                            PECMEOError{system,error}{PPIt,KinIt+1}{PPISLIt+1,KinISLIt+1}(:,10) = sum(PECMEOError{system,error}{PPIt,KinIt+1}{PPISLIt+1,KinISLIt+1}(:,1:sat),2)/numberOfPECMEOSats;
                        end
                    end
                end
            end
        end
	end
end

epochs = PECMEOReal{1}(:,1);
numberOfEpochs = length(epochs);
fprintf('Loading data done!\n')

% set(groot,'defaultLineLineWidth',1.0)

sat = 1;
name = ['PECMEOSat' num2str(sat)];

%% Graph of all the cases
%{
fig=fig+1;

% for system = 1:numberOfSystems
%     for error = 1:numberOfErrors
for system = 1:1
    for error = 3:4       
        systemMap = systemMaps{system};
        errorMap = errorMaps{error};
%         for sat = 1:numberOfPECMEOSats
%             name = ['PECMEOSat' num2str(sat)];
            for PPIt = 10:10
                for KinIt = 10:10
                    for PPISLIt = 0:0
                        for KinISLIt = 0:0
                            figure(fig)
%                             errorData = PECMEOError{character}{PPIt , KinIt +1}{PPISLIt +1, KinISLIt +1}(:,1:9); % Seperate satellites
                            errorData = PECMEOError{system,error}{PPIt , KinIt +1}{PPISLIt +1, KinISLIt +1}(:,10); % Mean of satellites
%                             plot(epochs, errorData,'o','DisplayName',[  insertBefore(erase(erase(characteristicsMaps{character},"\"),"4h_LTC0_"),"_","\") ', mean: ' num2str(mean(errorData),'%5.5e') ],'LineWidth',1);
                            plot(epochs, errorData,'o','DisplayName',[ systemMap errorMap 'mean: ' num2str(mean(errorData),'%5.5e') ]);
                            hold on;
                        end
                    end
                end
            end
%         end
    end
end
figure(fig)
lgd = legend('Location','best','FontSize',8);
% title(lgd,'SPP/KOD/SPPISL/KODISL','FontSize',8)
grid on;
% ylim([10^(-5) 10^2])
% ylim([0 0.05])
for iEpoch = 1:(numberOfEpochs-1)/60
    refreshEpoch = iEpoch*1800;
    xline(refreshEpoch,'--','HandleVisibility','off');
end
xlim([0 refreshEpoch]);
xlabel('time [s]')
ylabel('error [m]')
title([ 'Estimation errors for all satellites combined for different cases' ])
% hold off;
%}

%% Graphs per error type
%
% numberOfEpochs = length(epochs);
% numberOfEpochs = 8*60*60/30;

% fig=fig+1; figure(fig);
% markerType = {'+','o','*','.','x','s','d'};
markerType = {'o','d','s','*','x','+','.'};
fontsize = 12;
for error = 1:numberOfErrors
% for error = 1:1
    errorMap = errorMaps{error};
    fig=fig+1;
    for system = 1:numberOfSystems
%     for system = 2:2
        systemMap = systemMaps{system};
      
        figure(fig)

        errorData = PECMEOError{system,error}{10 , 10 +1}{0 +1, 0 +1}(:,10); % Mean of satellites
        errorDataAll{error,system} = PECMEOError{system,error}{10 , 10 +1}{0 +1, 0 +1}(:,10); % Mean of satellites
        errorTable(system,error) = mean(errorData);
%         errorData = PECMEOError{system,error}{PPIt , KinIt +1}{PPISLIt +1, KinISLIt +1}(:,1:9); % Seperate satellites
%         errorDataAll{error,system} = PECMEOError{system,error}{10 , 10 +1}{0 +1, 0 +1}(:,1:9); % Mean of satellites
%         errorTable2{system,error}(1,:) = mean(errorData);
%         subplot(2, 2, error) ;
%         plot(epochs, errorData,'o','DisplayName',[ replace(insertAfter(insertBefore(insertBefore(insertBefore(erase(systemMap,"\"),"_","\"),"Laser"," "),"NoISL","      "),"ISL","\_9C"),"_9C\_1C","_1C") ': ' num2str(mean(errorData),'%5.5e') ]);
%         semilogy(epochs, errorData,'o','DisplayName',[ replace(insertAfter(insertBefore(insertBefore(insertBefore(erase(systemMap,"\"),"_","\"),"Laser"," "),"NoISL","      "),"ISL","\_9C"),"_9C\_1C","_1C") ': ' num2str(mean(errorData),'%5.5e') ]);
%         semilogy(epochs, errorData,markerType{system},'DisplayName',[ replace(insertAfter(insertBefore(insertBefore(insertBefore(erase(systemMap,"\"),"_","\"),"Laser"," "),"NoISL","      "),"ISL","\_9C"),"_9C\_1C","_1C") ': ' num2str(mean(errorData),'%5.5e') ]);
%         semilogy(epochs, errorData,'o','DisplayName',[ replace(insertAfter(insertBefore(insertBefore(insertBefore(erase(systemMap,"\"),"_","\"),"Laser"," "),"NoISL","      "),"ISL","\_9C"),"_9C\_1C","_1C") ': ' num2str(mean(errorData),'%5.5e') ]);
        semilogy(epochs, errorData,'o','DisplayName',[ systems{system} ': ' replace(replace(num2str(mean(errorData),'%5.5e'),"e-0","\times10^{-"),"e-10","\times10^{-10") '}' ],'LineWidth',1.5);
%         plot(epochs, errorData,'o','DisplayName',[ systems{system} ': ' replace(replace(num2str(mean(errorData),'%5.5e'),"e-0","\times10^{-"),"e-10","\times10^{-10") '}' ]);
%         yline(errorTable(system,error),'--','HandleVisibility','off');
%         colord = get(gca,'colororder'); yline(errorTable(system,error),'--','Color',colord(system,:),'LineWidth',2);
        hold on;
    end
    
    for system = 1:numberOfSystems
%         yline(errorTable(system,error),'w--','LineWidth',3,'HandleVisibility','off');
%         colord = get(gca,'colororder'); yline(errorTable(system,error),'--','Color',colord(system,:),'LineWidth',2);
%         yline(errorTable(system,error),markerType{system},'Color',[0,0,0],'map');
%         yline(errorTable(system,error),'Color',[0,0,0],'DisplayName',markerType{system});
    end
    
    figure(fig)
    lgd = legend('Location','best','FontSize',fontsize-1);
%     lgd = legend('Location','best','NumColumns',4);
    set(gcf, 'Color', 'w');
    set(gcf, 'Position',  [0, 42, 1000, 954])
%     set(gcf, 'Position',  [0, 42, 2000, 954])
    grid on;
%     title(lgd,'TypeName: average','FontSize',12)
    title(lgd,'TypeName: 1-day average')
    % ylim([10^(-5) 10^2])
    % ylim([0 0.05])
    refreshEpoch = epochs(length(epochs));
    for iEpoch = 1:(numberOfEpochs-1)/60
        refreshEpoch = iEpoch*1800;
        xline(refreshEpoch,'--','HandleVisibility','off');
    end
    xline(refreshEpoch,'--','DisplayName','   Ephemeris update')
    xlim([0 refreshEpoch]);
    xlabel('time [s]','FontSize',fontsize)
    ylabel('error [m]','FontSize',fontsize)
    set(gca,'FontSize',fontsize)
%     title([ 'Estimation errors for all satellites combined for different cases' ])
%     title([ 'Estimation errors for all satellites combined ' insertBefore(erase(errorMap,"\"),"_","\") ])
    title([ 'Satellite average estimation error per epoch for ' errors{error} ' error case' ],'FontSize',fontsize-1)
    
    figName = erase(errorMap,"\");
    export_fig(sprintf('Results/estimation/%s_log_Sq2.pdf',figName))
end

fprintf('Plotting data done!\n')
%}

%%
colors = {'red','green','blue','cyan','magenta','yellow','black','white'};
% for error = 1:numberOfErrors
for error = 4:4
    for system = 4:numberOfSystems
%         yline(errorTable(system,error),'--','HandleVisibility','off');
        yline(errorTable(system,error),'Color',colors{system},'DisplayName',[ systems{system} ]);
    end
end
legend
grid on
set(gca,'YScale','log') 

%% Box plot

errorData=[];
groups =[];
group = 1;
fig=fig+1;
% for error = 1:numberOfErrors
for error = 4:4
    for system = 1:numberOfSystems
        errorData = [ errorData; errorDataAll{error,system} ];
        groups = [ groups; group * ones(size(errorDataAll{error,system})) ];
        group = group + 1;        
    end
end
figure(fig)
% boxplot(errorDataAll{error,system})
boxplot(errorData, groups)
xlabel('All Vehicles')
ylabel('error [m]')
title('Estimation errors with all error sources')
grid on
set(gca,'XTickLabel',systems(2:7))

% A = [16 20 15 17 22 19 17]';
% B = [22 15 16 16 16 18]';
% C = [23 9 15 18 13 27 17 14 16 15 21 19 17]';
% group = [    ones(size(A));
%          2 * ones(size(B));
%          3 * ones(size(C))];
% figure
% boxplot([A; B; C],group)
% set(gca,'XTickLabel',{'A','B','C'})


%% Box plot for 9 sats

errorData=[];
groups =[];
group = 1;
fig=fig+1;
% for error = 1:numberOfErrors
for error = 4:4
%     for system = 1:numberOfSystems
    for system = 2:2
        for sat = 1:9
            currentError = PECMEOError{system,error}{10 , 10 +1}{0 +1, 0 +1}(:,sat);
            errorData = [ errorData; currentError ];
            groups = [ groups; group * ones(size(currentError)) ];
            group = group + 1;
        end
    end
end
figure(fig)
% boxplot(errorDataAll{error,system})
boxplot(errorData, groups)
xlabel('The LPS sats')
ylabel('error [m]')
title('Estimation errors with all error sources')
grid on
% set(gca,'XTickLabel',{'1','2','C'})

%% Max baseline distance
for sat = 1:numberOfPECMEOSats
    for sat2 = setdiff(1:numberOfPECMEOSats, sat)
        PECMEORealBaseline = vecnorm(PECMEOReal{sat}(:,2:4) - PECMEOReal{sat2}(:,2:4) ,2,2);
        baselineMax(sat,sat2) = max(PECMEORealBaseline);
    end
end

max(max(baselineMax))

%% Intrgral of Laser Noise
fun = @(x) sqrt(1+ (0.003./x).^2 ) .* sqrt(1+ (0.01./x).^2 );
xmin = 0.002;
xmax = 0.1;
q = integral(fun,xmin,xmax);

%%
for error = 1:4
    mean(PECMEOError{2,error}{PPIt , KinIt +1}{PPISLIt +1, KinISLIt +1}(:,10))-mean(PECMEOError{3,error}{PPIt , KinIt +1}{PPISLIt +1, KinISLIt +1}(:,10))
end
