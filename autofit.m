
% Automated fitting of aptamer traces to find best linear fit over the
% longest time range

function rates = autofit(data,fignum)

% data is a matrix with time in the first column, data in subsequent
% columns and then deviations from repeats in the following column block of
% equal size --> [time; [data columns]; [standard deviations]);

% fignum is the number of the figure where fits are plotted

warning('off','curvefit:fit:noStartPoint');

% fit type
linXtract = fittype('k*x+b');

% use to setup a color gradient; useful for titrations
start_color = [0.05 0.05 0.05];
stop_color = [0.05 0.05 0.05];

colorGRADIENTflexible = @(i,N) start_color - (start_color-stop_color)*((i-1)/(N-1));

% set data cycles
dataiter = -1; % how many data points to remove from the end of the data upon each iteration
datamin = 10; % minimum number of data points to fit

% set R-squared threshold for a satisfactory fit
good = 0.998;

% time to start the fits in time column units
startpoint = 5; 

% set last time point to consider when starting to fit
endtime = 60;

% use for biological replicates with std already determined for technical replicates
for j = 2:size(data,2)-((size(data,2)-1)/2) % loop through traces in dataset

% use for technical replicate traces
%for j = 2:size(data,2) % loop through traces in dataset

    for i = find(data(:,1)<endtime,1,'last'):dataiter:datamin+startpoint
    
       %use with error weights (standard deviations from processed technical replicates) - fit data from startpoint data point to time determined by loop
       [currentfit,gof] = fit(data(startpoint:i,1),data(startpoint:i,j),linXtract,'Weights',1./(data(startpoint:i,j-1+size(data,2)-((size(data,2)-1)/2)) +1) );

       %use when no error weights are available (i.e., processing technical replicates)
       %[currentfit,gof] = fit(data(startpoint:i,1),data(startpoint:i,j),linXtract);

        if gof.rsquare > good
            disp(['data set #' num2str(j-1)]);
            disp(['data fit from ' num2str(data(startpoint,1)) ' to ' num2str(data(i,1)) '; with ' num2str(i-startpoint) ' data points']);
            disp(['Rsquare = ' num2str(gof.rsquare)]);
            disp(' ');
            break;
        elseif i == datamin+startpoint
            disp(['data set #' num2str(j-1)]);
            disp('*** MINIMUM DATA RANGE FIT ***');
            disp(['data fit from ' num2str(data(startpoint,1)) ' to ' num2str(data(i,1)) '; with ' num2str(i-startpoint) ' data points']);
            disp(['Rsquare = ' num2str(gof.rsquare)]);
            disp(' ');
        end
    end

    figure(fignum); hold all;

    colorGRADIENTflexible(j-1,(size(data,2)-1)/2);
    plot(data(:,1),data(:,j),'color',colorGRADIENTflexible(j-1,(size(data,2)-1)/2));
    plot(data(1:10:end,1),currentfit(data(1:10:end,1)),'--','color',colorGRADIENTflexible(j-1,(size(data,2)-1)/2)); 

%  use when fitting individual technical replicates (no error estimates in array)
%     colorGRADIENTflexible(j-1,size(data,2));
%     plot(data(:,1),data(:,j),'color',colorGRADIENTflexible(j-1,size(data,2)));
%     plot(data(1:10:end,1),currentfit(data(1:10:end,1)),'--','color',colorGRADIENTflexible(j-1,size(data,2))); 

    rates(j-1) = currentfit.k;
    currentfit

end

legend('off');
