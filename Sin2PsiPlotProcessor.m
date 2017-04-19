% Results processor for BESSY Sin2Psi Plot method results from Aug 2015
% Stanley van Kempen

%Script to evaluate BESSY measurement data. The .tau files obtained from
%the measurements have been renamed so that all the metadata to identify
%the measurement is in the file name. The file names to be evaluated by the
%script are in the file Sin2plotlist.txt. Additional files can only be
%added at the end of the file and with sample numbers exceeding 27! This is
%because subsequent results processing in the BESSY Summary Excel file is
%file-order dependent due to absolute cell references in Sheet2.

clc;
close all;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%File control

%Run function allfitdist to determine the most likely PDF of the measurement data
%On => evaluMeasDistr = 1; Off => evalMeasDistr = 0;
evalMeasDistr = 0;

%Save stress plots for every single measurement point and save to path3
%location on = 1, off = 0
singleMeasEval = 1;

% Value filtering
%Remove all measurements with a mean stress greater than stressLimit from the measurement results
%Remove all measurements with a error to mean stress ratio greater than errorRatioLimit from the measurement results 
stressLimit = 500; %MPa
errorRatioLimit = 10; %ratio

%Path references for Matlab files, measurement files, and matlab result
%output
path = '\\win.iwm.rwth-aachen.de\users\User\vanKempen\Desktop\RWTH\-= Experiments =-\BESSY Results 2016\Matlab\';
path2 = '\\win.iwm.rwth-aachen.de\users\User\vanKempen\Desktop\RWTH\-= Experiments =-\BESSY Results 2016\Tau_Files_renamed\';
path3 = '\\win.iwm.rwth-aachen.de\users\User\vanKempen\Desktop\RWTH\-= Experiments =-\BESSY Results 2016\Matlab_Results\';

fileList = getFileList(strcat(path,'Sin2plotlist.txt'),1,inf);
jEnd = length(fileList);

for j = 1:1:jEnd;
    
    sampleName = char(fileList(j));
    close all;
        
    %Import data from single file
    [Tau, S11, S11_Error, S22, S22_Error] = GetFileSin2Psi(strcat(path2, sampleName));
    
    %Filter erroneous NaN values in data files and select the use of S11 or S22
    Tau(isnan(Tau)) = [];
    
    if isnan(S11(1));
        Sii = S22;
        Sii_Error = S22_Error;
        label = '\sigma_{yy}';
        Stest = S22_Error;
    else
        Sii = S11;
        Sii_Error = S11_Error;
        label = '\sigma_{xx}';
    end

    %%Get entry number of NaN value and remove Sii and Sii_Error values from vector
        [row, col] = find(isnan(Sii));
        empty = isempty(row);
        if empty == 0;
        Sii(row) = [];
        Sii_Error(row) = [];
        end
        [row, col] = find(isnan(Sii_Error));
        empty = isempty(row);
        if empty == 0; 
        Sii(row) = [];
        Sii_Error(row) = [];
        end
        
    %Filter out excessive errors and measurement values using previously
    %defined limits
        [row, col] = find(ge(abs(Sii),stressLimit));
        empty = isempty(row);
        if empty == 0;
        Sii(row) = [];
        Sii_Error(row) = [];
        end
        
        errorRatio = abs(Sii_Error./Sii);
        [row, col] = find(ge(errorRatio,errorRatioLimit));
        empty = isempty(row);
        if empty == 0;
        Sii(row) = [];
        Sii_Error(row) = [];
        end
        
     % Check if Sii has been emptied completely, if so the vector will be
     % filled with 0s to prevent the script from stopping the evaluation
        empty = isempty(Sii);
        emptySum = sum(empty);
        if emptySum > 0;
            Sii = zeros(12,1);
            Sii_Error = zeros(12,1);
        end       
        
    % Build Structure with all results for single laminates    
    % Calculate Sii mean value and error for single laminates
    Smean = sum(Sii)/length(Sii);
    Smean_Error = sqrt(sum(Sii_Error.^2));
    Sigma = std(Sii);
    standardError = 2*Sigma/sqrt(length(Sii)); % 95% probability that mean value is between plus/min this error
    
    %Decompose file name tag
    string = cell2mat(fileList(j));
    remExt = strtok(string,'.');
    splitStr = strsplit(remExt,'_');
    splitStrCell(j,:) = splitStr;
    
    %Convert letters from material tag into numbers
    splitStr(strcmp('M',splitStr))={'1'};
    splitStr(strcmp('A',splitStr))={'3'};   % Both location A and material A are changed to 3!
    splitStr(strcmp('L',splitStr))={'1'};
    splitStr(strcmp('T',splitStr))={'2'};
    splitStr(strcmp('S',splitStr))={'2'};

    %Store laminate data and measurement results for single laminates
    xresStress(j).layupNr = splitStr(1);      %see Excel file with sample layup numbers
    xresStress(j).sampleNr = splitStr(2);     %Sample nr of particular layup
    xresStress(j).measLoc = splitStr(3);      %Measurement location on laminate surface M = center of laminate; A = outer point on laminate; 2 = intermediate point (only in large laminates) 
    xresStress(j).stressDir = splitStr(4);    %Stress orientation; L = longitudinal (S11); T = transversal (S22)
    xresStress(j).measDepth = splitStr(5);    %Measurement depth; 1 = sample surface (0 - 250 µm); 2  = just above layer interface (~350µm); 3 = @ layer interface (~450µm); 4 = just below layer interface (~550µm);
    xresStress(j).layerMat = splitStr(6);     %Material at the measurement point; A = alumina; S = Spinel; In case of an interface measurement (measDepth = 3) both an A and an S evaluation have been made since both materials were in the measured volume of the beam
    xresStress(j).stress = Sii;
    xresStress(j).error = Sii_Error;
    xresStress(j).meanStress = Smean;
    xresStress(j).standardError = standardError;
    xresStress(j).tag = cell2mat(fileList(j));
    xresStress(j).label = label;
    
end

%Plot stresses and errors for single measurements
if singleMeasEval == 1;
for j = 1:1:jEnd
    SiiLength = length(xresStress(j).stress);
    Sii = xresStress(j).stress;
    Sii_Error = xresStress(j).error;
    sampleName = xresStress(j).tag;
    label = xresStress(j).label;
    xRange = linspace(1,SiiLength,SiiLength);
    
    figure (1);
    hold on
    errorbar(xRange,Sii,Sii_Error,'.','Color',[0.1,0.1,0.1]);
%     title(sampleName, 'Interpreter', 'none');
    set(gca,'fontsize',16,'YMinorTick','on');
    xlabel(strcat('Lattice Number'));
    ylabel(strcat(label,'[MPa]'));
%     axis([0,(xRange(end)+1),-stressLimit,stressLimit]);
    axis([0,12,-stressLimit,stressLimit]);
    
    %fit linear curve to data
    if SiiLength > 1
    ft = fittype( 'poly1' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    [fitresult, gof] = fit( xRange', Sii, ft, opts );
    plot(xRange, fitresult(xRange));
    plot(xRange, confint(fitresult));
    end
    
%     %Generate regression confidence interval
%     alpha = 0.05;
%     beta = 1;
%     [top_int, bot_int,X] = regression_line_ci(alpha,beta,xRange,Sii);
%      [b,bint] = regress(Sii,xRange') 
%      xval = 0:0:xRange;
%      yhat = bint(1)+bint(2)*xval;
%      ylow = bint(1,1)+bint(2,1)*xval;
%      yupp = bint(1,2)+bint(2,2)*xval;
    
    %Print and save figure to folder
    sampleName = xresStress(j).tag;
    [pathstr,name,ext] = fileparts(sampleName);
    print(strcat(path3,name),'-dpng');
    hold off
    close all
end
end


%Group measurement results for each measurement location on different
%samples with the same layup
for j = 1:1:length(xresStress)
    r = [];
    k = str2num(cell2mat(xresStress(j).layupNr));
    l = str2num(cell2mat(xresStress(j).measLoc));
    m = str2num(cell2mat(xresStress(j).stressDir));
    n = str2num(cell2mat(xresStress(j).measDepth));
    o = str2num(cell2mat(xresStress(j).sampleNr));
    p = str2num(cell2mat(xresStress(j).layerMat));
    q = xresStress(j).stress;
    
    for a = 1:1:27
        if k == a;
            for b = 1:1:3
                if l == b
                    for c = 1:1:2
                        if m == c
                            for d = 1:1:4
                                if n == d
                                   for e = 1:1:3
                                       if o == e
                                           varName = genvarname(strcat(xresStress(j).layupNr,xresStress(j).measLoc,xresStress(j).stressDir,xresStress(j).measDepth,xresStress(j).layerMat));
                                           varStruct(j).Probe = varName;
                                           r = [r; q];
                                           varStruct(j).Stress = r;
                                           varCell(j,1) = varName;
                                       end
                                   end
                                end
                            end
                        end
                    end
                end
            end
        end
    end                                                                        
end

%Identify unique measurement points and group sample measurements
C = cellfun(@char,{varStruct.Probe},'unif',0);
[~,idx] = unique(C);
idx = sort(idx);

for i = 1:1:length(idx);
    z = cell2mat(varStruct(idx(i)).Probe);
    results(i).SampleData = cell2mat(splitStrCell(idx(i),:));
    results(i).Tag = cell2mat(varStruct(idx(i)).Probe);
        validInds=[]; 
    for j=1:numel(varCell)
         ifValid=(strcmp(varCell(j),z));
         if ifValid==1
             validInds=[validInds;j];
         end
    end
    
    w = length(validInds);
    if w == 1
       w1 = varStruct(validInds(1)).Stress;
       measPointVect = [w1'];
       results(i).EvalSamples = w;
       results(i).Sample1_Mean = xresStress((validInds(1))).meanStress;
       results(i).Sample1_SE = xresStress((validInds(1))).standardError;
       results(i).Sample2_Mean = [];
       results(i).Sample2_SE = [];
       results(i).Sample3_Mean = [];
       results(i).Sample3_SE = [];
    elseif w == 2
       w1 = varStruct(validInds(1)).Stress;
       w2 = varStruct(validInds(2)).Stress;
       measPointVect = [w1' w2'];
       results(i).EvalSamples = w;
       results(i).Sample1_Mean = xresStress((validInds(1))).meanStress;
       results(i).Sample1_SE = xresStress((validInds(1))).standardError;
       results(i).Sample2_Mean = xresStress((validInds(2))).meanStress;
       results(i).Sample2_SE = xresStress((validInds(2))).standardError;
       results(i).Sample3_Mean = [];
       results(i).Sample3_SE = [];
    elseif w == 3
       w1 = varStruct(validInds(1)).Stress;
       w2 = varStruct(validInds(2)).Stress;
       w3 = varStruct(validInds(3)).Stress;
       measPointVect = [w1' w2' w3'];
       results(i).EvalSamples = w;
       results(i).Sample1_Mean = xresStress((validInds(1))).meanStress;
       results(i).Sample1_SE = xresStress((validInds(1))).standardError;
       results(i).Sample2_Mean = xresStress((validInds(2))).meanStress;
       results(i).Sample2_SE = xresStress((validInds(2))).standardError;
       results(i).Sample3_Mean = xresStress((validInds(3))).meanStress;
       results(i).Sample3_SE = xresStress((validInds(3))).standardError;
    elseif w == 4
       w1 = varStruct(validInds(1)).Stress;
       w2 = varStruct(validInds(2)).Stress;
       w3 = varStruct(validInds(3)).Stress;
       w4 = varStruct(validInds(4)).Stress;
       measPointVect = [w1' w2' w3' w4'];
       results(i).EvalSamples = w;
       results(i).Sample1_Mean = xresStress((validInds(1))).meanStress;
       results(i).Sample1_SE = xresStress((validInds(1))).standardError;
       results(i).Sample2_Mean = xresStress((validInds(2))).meanStress;
       results(i).Sample2_SE = xresStress((validInds(2))).standardError;
       results(i).Sample3_Mean = xresStress((validInds(3))).meanStress;
       results(i).Sample3_SE = xresStress((validInds(3))).standardError;
    elseif w == 5
       fprintf('W5/%n');
       w1 = varStruct(validInds(1)).Stress;
       w2 = varStruct(validInds(2)).Stress;
       w3 = varStruct(validInds(3)).Stress;
       w4 = varStruct(validInds(4)).Stress;
       w5 = varStruct(validInds(5)).Stress;
       measPointVect = [w1' w2' w3' w4' w5'];
       results(i).EvalSamples = w;
       results(i).Sample1_Mean = xresStress((validInds(1))).meanStress;
       results(i).Sample1_SE = xresStress((validInds(1))).standardError;
       results(i).Sample2_Mean = xresStress((validInds(2))).meanStress;
       results(i).Sample2_SE = xresStress((validInds(2))).standardError;
       results(i).Sample3_Mean = xresStress((validInds(3))).meanStress;
       results(i).Sample3_SE = xresStress((validInds(3))).standardError;
    elseif w == 6
       fprintf('W6/%n');  
       w1 = varStruct(validInds(1)).Stress;
       w2 = varStruct(validInds(2)).Stress;
       w3 = varStruct(validInds(3)).Stress;
       w4 = varStruct(validInds(4)).Stress;
       w5 = varStruct(validInds(5)).Stress;
       w6 = varStruct(validInds(6)).Stress;
       measPointVect = [w1' w2' w3' w4' w5' w6'];
       results(i).EvalSamples = w;
       results(i).Sample1_Mean = xresStress((validInds(1))).meanStress;
       results(i).Sample1_SE = xresStress((validInds(1))).standardError;
       results(i).Sample2_Mean = xresStress((validInds(2))).meanStress;
       results(i).Sample2_SE = xresStress((validInds(2))).standardError;
       results(i).Sample3_Mean = xresStress((validInds(3))).meanStress;
       results(i).Sample3_SE = xresStress((validInds(3))).standardError;
    elseif w >= 7
        fprintf('W7/%n');
    end
       
    lamStressVect{i,:} = measPointVect;
    
    %Fit pdf and calculate mean stress and stddev for each measurement
    %point group
    for i = 1:1:length(lamStressVect)
        stressVect = cell2mat(lamStressVect(i));
            
        %Results for normal distribution
        pd = fitdist(stressVect','Normal');
        results(i).Mu = mean(pd);
        results(i).Sigma = std(pd);
        
        %Shift stress values to positive range to allow Weibull fitting
        
        shiftVal = min(stressVect);
        if shiftVal < 1
           shiftVal = abs(shiftVal) + 1;
        else
            shiftVal = 1;
        end
        
        stressVect = stressVect + shiftVal;
        
        %Results for Weibull distribution
        pd = fitdist(stressVect','Weibull');
        results(i).MuWeibull = mean(pd) - shiftVal;
        results(i).SigmaWeibull = std(pd);


        %Evaluate measurement distribution for single sample
        if evalMeasDistr == 1
        [D PD] = allfitdist(stressVect);
        fieldVal = extractfield(D(1), 'DistName');
        results(i).fitField = fieldVal;
        end
    end
        
           
end

%Export processed measurement results to Excel file
A = struct2cell(results);
B(:,:) = A(:,1,:);

filename = 'BESSY Results Summary.xlsx';
% xlswrite(filename,fieldnames(results)');
xlswrite(filename,B','MatlabResults','A4');
