clear all;
%%%%%%%%%%%%%%%%%%%%%
% MAIN PARAMETERS
%%%%%%%%%%%%%%%%%%%%%

% INPUT FOLDER should contain .txt files with 3 columns:
% 1)-2) x and y coordinates
% 3) intensity values


%number of averaging points for smoothing. 
% bigger values correspond to more smoothing 
nAverPoints = 50; 
%pixel resolution. number of pixels in one micron.
nPixMkm = 15.8;
%threshold values of normalized intensity to determine 
%Start and End values (in range from 0 to 1)
dThreshold = 0.33;
%including tails or not in case of center moving average
bTails = true;

% END OF PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%


% open text file with intensity profile
% ask user to locate folder
filesfolder = uigetdir;
%if cancel button is not pressed
if ~isequal(filesfolder, 0)
        disp('Script started.');
        disp('Please wait......');
        
        %get the list of all files in folder
        dirListing = dir(filesfolder);
        totFiles = length(dirListing);
        
        %counting all .txt files
        totLogFiles=0;
        for nFile = 1:totFiles
                % checking is it .txt file
                filenamein = fullfile(filesfolder,dirListing(nFile).name);
                fExtension = (filenamein(length(filenamein)-3:length(filenamein)));

                if strcmpi(fExtension,'.txt')
                    totLogFiles = totLogFiles + 1;
                end
        end
        %allocating arrays
        allresults = zeros(totLogFiles, 21);
        allnormtrace = [];
        resultheader = {'Filename', 'Start from left (1), mkm', 'Start from max value(2), mkm', 'Max, mkm', 'End from right end (1), mkm','End from max value (2), mkm', 'AISlength (End1-Start1), mkm', 'AISlength (End2-Start1), mkm', 'AISlength (End1-Start2), mkm', 'AISlength (End2-Start2), mkm', '# of smooth points', 'intensity threshold', 'resolution (pixels in mkm)', 'Tails', 'Int_Left_tail1, a.u./mkm', 'Int_Left_tail2, a.u./mkm','Int_Middle_(Start1_End1)','Int_Middle_(Start2_End1)', 'Int_Middle_(Start1_End2)', 'Int_Middle_(Start2_End2)', 'Int_Right_tail1, a.u./mkm', 'Int_Right_tail2, a.u./mkm'};
        normtraceheader = cell(2,2*totLogFiles); 
        totLogFiles = 0;
        for nFile = 1:totFiles

                % checking is it .txt file
                filenamein = fullfile(filesfolder,dirListing(nFile).name);
                fExtension = (filenamein(length(filenamein)-3:length(filenamein)));

                if strcmpi(fExtension,'.txt')
                    totLogFiles = totLogFiles + 1;
                    %read intensity trace
                    inttrace = dlmread(filenamein,'\t');
                    sz = size(inttrace);
                    %total number of points
                    nTotPoints = sz(1);
                    %array containing results
                    norm_int = zeros(nTotPoints,2,'double');
                    nonnorm_int = zeros(nTotPoints,2,'double');
                    norm_int(:,1)=inttrace(:,1)*(1/nPixMkm);
                    nonnorm_int(:,1)=norm_int(:,1);
                    nonnorm_int(:,2)=inttrace(:,2);

                    if (rem(nAverPoints,2)==0)
                        nAverPoints=nAverPoints+1;
                    end
                    nHalf = (0.5*(nAverPoints-1));

                    %SMOOTHING
                    %including 'tails' in the beginning and end of the track
                    %depending on the flag
                    if(bTails)
                        %beginning
                        for i=1:nHalf
                            norm_int(i,2) = mean(inttrace(1:(i+nHalf),2));
                        end
                        %end
                        for i=(nTotPoints-nHalf+1):nTotPoints
                            norm_int(i,2) = mean(inttrace((i-nHalf):nTotPoints,2));
                        end
                    end
                    %smoothing middle part
                    for i=(nHalf+1):(nTotPoints-nHalf)
                        norm_int(i,2) = mean(inttrace((i-nHalf):(i+nHalf),2));
                    end

                    %NORMALIZATION
                    %min intensity
                    minInt = min(norm_int(:,2));
                    %max-min intensity, normalization factor
                    dSpan = max(norm_int(:,2)) - minInt;
                    %normalize intensity
                    for i=1:nTotPoints
                        norm_int(i,2) = (norm_int(i,2) - minInt) / dSpan;
                    end 

                    %SEARCH FOR REFERENCE POINTS
                    %Max
                    [~, MaxInd] = max(norm_int(:,2));
                    MaxVal = norm_int(MaxInd,1);
                    
                    %Start 1 (from left most point)
                    i=0;
                    dSpan = 0;
                    while dSpan < dThreshold
                        i = i + 1;
                        dSpan = norm_int(i, 2);
                    end
                    if i==1
                        i = 2;
                    end
                    StartVal1 = norm_int(i-1, 1);
                    StartInd1 = i-1;

                    %Start 2 (from maximum)
                    i=MaxInd;
                    dSpan = MaxVal;
                    while (dSpan > dThreshold) && (i>1)
                        i = i - 1;
                        dSpan = norm_int(i, 2);
                    end
                    StartVal2 = norm_int(i, 1);
                    StartInd2 = i;
                    
                    %End 1 (from right most point)
                    i=nTotPoints+1;
                    dSpan = 0;
                    while dSpan < dThreshold
                        i = i - 1;
                        dSpan = norm_int(i, 2);
                    end
                    EndVal1 = norm_int(i, 1);
                    End1Ind = i;
                    
                    %End 2 (from maximum)
                    i=MaxInd;
                    dSpan = 2;
                    while (dSpan > dThreshold) && (i < nTotPoints)
                        i = i + 1;
                        dSpan = norm_int(i, 2);
                    end
                    EndVal2 = norm_int(i-1, 1);
                    End2Ind = i-1;
                    
                    

                    %Calculate AIS lengths
                    LengthVal1 = EndVal1 - StartVal1;
                    LengthVal2 = EndVal2 - StartVal1;                                            
                    LengthVal3 = EndVal1 - StartVal2;
                    LengthVal4 = EndVal2 - StartVal2;
             
                    %INTENSITY QUANTIFICATION
                    bLeftTail1 = sum(nonnorm_int(1:StartInd1,2))/nonnorm_int(StartInd1,1);
                    bLeftTail2 = sum(nonnorm_int(1:StartInd2,2))/nonnorm_int(StartInd2,1);
                    bMeanInt11 = sum(nonnorm_int((StartInd1+1):End1Ind,2))/(nonnorm_int(End1Ind,1)-nonnorm_int(StartInd1,1));
                    bMeanInt21 = sum(nonnorm_int((StartInd2+1):End1Ind,2))/(nonnorm_int(End1Ind,1)-nonnorm_int(StartInd2,1));
                    bMeanInt12 = sum(nonnorm_int((StartInd1+1):End2Ind,2))/(nonnorm_int(End2Ind,1)-nonnorm_int(StartInd1,1));
                    bMeanInt22 = sum(nonnorm_int((StartInd2+1):End2Ind,2))/(nonnorm_int(End2Ind,1)-nonnorm_int(StartInd2,1));
                    bRightTail1 = sum(nonnorm_int((End1Ind+1):nTotPoints,2))/(nonnorm_int(nTotPoints,1)-nonnorm_int(End1Ind,1));
                    bRightTail2 = sum(nonnorm_int((End2Ind+1):nTotPoints,2))/(nonnorm_int(nTotPoints,1)-nonnorm_int(End2Ind,1));
                   
                    
                    %STORING RESULTS
                    allresults(totLogFiles,:) = [StartVal1, StartVal2, MaxVal, EndVal1, EndVal2, LengthVal1, LengthVal2, LengthVal3, LengthVal4, nAverPoints, dThreshold, nPixMkm, bTails, bLeftTail1, bLeftTail2, bMeanInt11, bMeanInt21, bMeanInt12, bMeanInt22, bRightTail1, bRightTail2];
                    
                    normtraceheader{1,totLogFiles*2-1} = 'File:';
                    normtraceheader{1,totLogFiles*2} = dirListing(nFile).name;  
                    normtraceheader{2,totLogFiles*2-1} = 'length,mkm';
                    normtraceheader{2,totLogFiles*2} = 'normalized intensity';  
                    resultheader{1+totLogFiles,1} = dirListing(nFile).name;
                    if (length(allnormtrace)>0)
                        allnormtrace = padadd(allnormtrace,norm_int(:,1));
                        allnormtrace = padadd(allnormtrace,norm_int(:,2));
                    else
                        allnormtrace = norm_int;
                    end
                    %informing user about progress
                    msgline = strcat(dirListing(nFile).name, ' finished');
                    disp(msgline);                    
                end
        end
        disp('Analysis done. Saving...');
        if(totLogFiles>0)
            filenamein = fullfile(filesfolder,'summary.xlsx');
             xlswrite(filenamein,resultheader,'Sheet1','A1');
             xlswrite(filenamein,allresults,'Sheet1','B2');
             xlswrite(filenamein,normtraceheader,'Sheet2','A1');
             xlswrite(filenamein,allnormtrace,'Sheet2','A3');
        end
        disp('Done.');

end