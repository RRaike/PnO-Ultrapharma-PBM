clc
close all
format shortE

%{
Samples used:
    2604: Sample11
    2904: Sample4
    0503: Sample4
%}


Data2604 = modelDataBestAllExperiments{1};
BestData2604 = Data2604{7, 10}{end};

Data2904 = modelDataBestAllExperiments{2};
BestData2904 = Data2904{2, 1}{end};

Data0305 = modelDataBestAllExperiments{3};
BestData0305 = Data0305{5, 7}{end};

bestModelData = {BestData2604, BestData2904, BestData0305};

clear BestData2604 BestData2904 BestaData0305

[spanStorageExperiment, spanSizesStorageExperiment] = DataCharacteristics(experimentalData, newBins, 0.10, 0.90);
[spanStorageModel, spanSizesStorageModel] = DataCharacteristics(bestModelData, newBins, 0.10, 0.90);

% order of birth for every best model defined above. This is the row index.
bModel = [7, 2, 5];

% order of growth for every best model defined above. This is the column
% index.
gModel = [10, 1, 7];

averageBins = AverageModelVector(newBins);

for k = 1: length(experimentalData)
    currentExperimentalData = experimentalData{k};
    currentModelData = bestModelData{k};


    b = bModel(k);
    g = gModel(k);
    Kb = KbBestAllExperiments{k}(b,g);
    Kg = KgBestAllExperiments{k}(b,g);
    RMSE = MSEBestAllExperiments{k}(b,g);
    
    currentDay = days(k);
    
    currentSpanSizesExperiment = spanSizesStorageExperiment{k};
    d10Experiment = currentSpanSizesExperiment(1);
    meanExperiment = currentSpanSizesExperiment(2);
    d90Experiment = currentSpanSizesExperiment(3);
    
    currentSpanSizesModel = spanSizesStorageModel{k};
    d10Model = currentSpanSizesModel(1);
    meanModel = currentSpanSizesModel(2);
    d90Model = currentSpanSizesModel(3);

    spanExperiment = spanStorageExperiment{k};
    spanModel = spanStorageModel{k};

    disp(currentDay)
    disp("# crystals experiment: " + string(sprintf('%.3e', sum(currentExperimentalData))))
    disp("# crystals model: " + string(sprintf('%.3e', sum(currentModelData))))
    disp(" ")
    disp("D10 of experiment: " + string(sprintf('%.3e', d10Experiment)))
    disp("Mean of experiment: " + string(sprintf('%.3e', meanExperiment)))
    disp("D90 of experiment: " + string(sprintf('%.3e', d90Experiment)))
    disp(" ")
    disp("D10 of model: " + string(sprintf('%.3e', d10Model)))
    disp("Mean of model: " + string(sprintf('%.3e', meanModel)))
    disp("D90 of model: " + string(sprintf('%.3e', d90Model)))
    disp(" ")
    disp("Span of experiment: " + string(spanExperiment))
    disp("Span of model: " + string(spanModel))
    disp(" ")
    
    %{
     h_experiment = histogram();
     h_experiment.BinCounts = currentExperimentalData;
     h_experiment.BinEdges = newBins;
     
    hold on
    
    h_model = histogram();
    %h_model = bar(averageBins, currentModelData);
    h_model.BinCounts = currentModelData;
    h_model.BinEdges = newBins;
            
    xlabel("Crystal size [m_{crystal}]")
    ylabel(["Number probability function of crystals", "[# crystals/m_{crystal}]"])
    legend(["Experimental Data", 
            "Model data"])
    
    stringInAnnotationHeaders = {"Order of birth:", 
                                 "Order of growth:",
                                 "Birth constant:",
                                 "Growth constant:",
                                 "RMSE:"};
    annotation('textbox', [.56, .5, .3, .3], 'String', stringInAnnotationHeaders, 'FitBoxToText', 'on', 'LineStyle', 'none')

    stringInAnnotationValues = {string(b),
                                string(g),
                                string(sprintf('%.3e',Kb)),
                                string(sprintf('%.3e' ,Kg)),
                                string(sprintf('%.3e', RMSE)) };
    annotation('textbox', [.76, .5, .3, .3], 'String', stringInAnnotationValues, 'FitBoxToText', 'on', 'LineStyle', 'none')
    %}

    close all
    clc
end



%{
for experimentCounter = 3: length(experimentalData)
    currentExperimentalData = experimentalData{experimentCounter};

    MSEBestCurrentExperiment = MSEBestAllExperiments{experimentCounter};
    KgBestCurrentExperiment = KgBestAllExperiments{experimentCounter};
    KbBestCurrentExperiment = KbBestAllExperiments{experimentCounter};


    for bCounter = 1 : length(bInterval)
        for gCounter = 1:length(gInterval)
            b = bInterval(bCounter);
            g = gInterval(gCounter);

            Kb = KbBestCurrentExperiment(bCounter, gCounter);
            Kg = KgBestCurrentExperiment(bCounter, gCounter);

            currentModelSet = modelDataCurrentExperiment{bCounter, gCounter};
            bestCurrentModelData = currentModelSet{end};
            
            RMSE = MSEBestCurrentExperiment(bCounter, gCounter);
            
            disp("The current Experiment data is: " + string(days(experimentCounter)));
            disp("The current b value is: " + string(b));
            disp("The current g value is: " + string(g));
            disp("The current Kb value is: " + string(Kb));
            disp("The current Kg value is: " + string(Kg));
            disp("The current RMSE value is: " + string(RMSE));
            
            h_model = histogram();
            h_model.BinCounts = abs(currentExperimentalData);
            h_model.BinEdges = newBins;
    
            hold on
        
            h_experimental = histogram();
            h_experimental.BinCounts = abs(bestCurrentModelData);
            h_experimental.BinEdges = newBins;
            
            %title(["Comparison between experimental data of case 3", ... %+ string(days(experimentCounter)), ...
             %       " and modeled data"])
            xlabel("Crystal size [m_{crystal}]")
            ylabel(["Number probability function of crystals", "[# crystals/m_{crystal}]"])

            legend(["Experimental Data of " + string(days(experimentCounter)), 
                    "Model data"])

            stringInAnnotationHeaders = {"Order of birth:", 
                                  "Order of growth:",
                                  "Birth constant:",
                                  "Growth constant:",
                                  "RMSE:"};
            annotation('textbox', [.56, .5, .3, .3], 'String', stringInAnnotationHeaders, 'FitBoxToText', 'on', 'LineStyle', 'none')

            stringInAnnotationValues = {string(b),
                                       string(g),
                                       string(sprintf('%.3e',Kb)),
                                       string(sprintf('%.3e' ,Kg)),
                                       string(sprintf('%.3e', RMSE)) };
            annotation('textbox', [.76, .5, .3, .3], 'String', stringInAnnotationValues, 'FitBoxToText', 'on', 'LineStyle', 'none')

            close all
            clc
        end
    end
end
%}