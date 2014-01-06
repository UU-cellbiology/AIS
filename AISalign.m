%initialization. reset all variables
clear all;

%parameters
%space step of final plot
dStep = 0.1;


%open Excel file with normalized intensities
[FileName,PathName] = uigetfile('*.xls','Select Excel file with data...');
FileNameFull = sprintf('%s%s',PathName,FileName);
choiceColumn = menu('Choose a column:','Start from left','Start from max value','Max' , 'End from right end','End from max value');
disp('Program started! Wait till THE END or error message.');

%check if it is an valid Excel file
[status,sheets] = xlsfinfo(FileNameFull);
if status == 'Microsoft Excel Spreadsheet'

    %analyzing sheet #2
    nTotSheets = length(sheets);
    %read sheet data from Excel file
    %distributions of fluorescence
    [num,txt,raw] = xlsread(FileNameFull, sheets{2});
    %summary of positions
    [num1,txt1,raw1] = xlsread(FileNameFull, sheets{1});
    %number of columns with numbers in sheet
    sz = size(num);
    nTotColumns = sz(2);
    %total number of pairs to align
    nTotPlots = nTotColumns/4;
    %array containing at 1st column array index of maximum intensity values of
    %reference protein
    %and at 2nd column maximum intensity values themself
    MaximumArrays = zeros(nTotPlots,2,'double');
    %array containing at 1st column array index of length of
    %reference\aligned protein
    %and at 2nd column length from maximum till right end 
    MaxLengths = zeros(nTotPlots,2,'double');
    
    %determine maximum positions and total span of final scale
    for i=1:nTotPlots
        refval = num1(i*2,choiceColumn);
        MaximumArrays(i,2)= refval;
        spaceaxis = abs(num(:,(i*4-3))-refval);
        [~, MaxInd] = min(spaceaxis);
        
        MaximumArrays(i,1)= MaxInd;
        %MaximumArrays(i,2)= num(MaxInd,(i*4-3));
        
        [~, MaxInd] = max(num(:,(i*4-3)));
        MaxLengths(i,1) = MaxInd;
        MaxLengths(i,2) = num(MaxInd,(i*4-3)) - MaximumArrays(i,2);               
    end
    [~, MaxLeft] = max(MaximumArrays(:,2)); %index of max disctance from max till left border
    [~, MaxRight] = max(MaxLengths(:,2)); %index of max disctance from max till right border
    %total span of all staining after shifting
    dLeftEnd = (-1)*MaximumArrays(MaxLeft,2);
    dRightEnd = MaxLengths(MaxRight,2);
    nTotalXSteps = floor((dRightEnd-dLeftEnd)/dStep) ;
    
   
    %nTotalXSteps = MaximumArrays(MaxLeft,1) + (MaxLengths(MaxRight,1) - MaximumArrays(MaxRight,1) );   
    %forming new x-axis, filling with numbers
    nXAxis = zeros(nTotalXSteps, 1+2*nTotPlots, 'double');
    nXAxis(:,:) = NaN;
    dLeftEnd = ceil(dLeftEnd/dStep)*dStep;
    for i=1:nTotalXSteps
        nXAxis(i,1) = dLeftEnd + i*dStep;
    end
    
    %MaxLeftIndex = MaximumArrays(MaxLeft,1);
    %for i=1:MaxLeftIndex
    %    nXAxis(i,1) = num(i,(MaxLeft*4-3))-MaximumArrays(MaxLeft,2);
    %end
    %for i=MaximumArrays(MaxRight,1):MaxLengths(MaxRight,1)
    %    nXAxis(MaxLeftIndex+i-MaximumArrays(MaxRight,1),1) = num(i,(MaxRight*4-3))-MaximumArrays(MaxRight,2);
    %end

    for i=1:nTotPlots
        %now calculating shifted values
        nNewX = num(:,i*4-3);
        nNewX = nNewX(~isnan(nNewX));
        sz = size(nNewX);
        for j=1:sz(1)
            nNewX(j)=nNewX(j)-MaximumArrays(i,2);
        end
        dCurrMin = nNewX(1);
        dCurrMax = nNewX(sz(1));
        nRef = num(:,i*4);
        nRef  = nRef(~isnan(nRef));
        nSec = num(:,i*4-2);
        nSec  = nSec(~isnan(nSec));
        nRef = horzcat(nNewX, nRef, nSec);
        for k=1:nTotalXSteps
            xCoord = nXAxis(k,1);
            if (xCoord >=dCurrMin && xCoord<=dCurrMax)
                j=2;
                while(xCoord>nRef(j,1))
                    j=j+1;
                end
                %linear approximation
                nSlope = (nRef(j,2)-nRef(j-1,2))/(nRef(j,1)-nRef(j-1,1));
                nXAxis(k,1+i) = nRef(j-1,2)+ nSlope *(xCoord-nRef(j-1,1));
                nSlope = (nRef(j,3)-nRef(j-1,3))/(nRef(j,1)-nRef(j-1,1));
                nXAxis(k,1+i+nTotPlots) = nRef(j-1,3)+ nSlope *(xCoord-nRef(j-1,1));         
            end
        end
        
    end
    
     disp('Alignment is done. Saving...');
     filenamein = fullfile(PathName,'summary_aligned.xls');
     %delete file if it already exists
     if exist(filenamein, 'file')
        delete(filenamein);
     end
     %saving headers
     sHeaders = cell(1,1+2*nTotPlots);
     sHeaders{1,1} = 'length, mkm';
     sHeaders{1,2}= 'Reference protein replicates';
     sHeaders{1,2+nTotPlots}= 'Aligned protein replicates';
     
     for i=2:nTotPlots
         sHeaders{1,1+i}= ' ';
         sHeaders{1,1+i+nTotPlots}= ' ';
     end
     xlswrite(filenamein, sHeaders,'Sheet1','A1');    
     %saving aligned data
     xlswrite(filenamein, nXAxis,'Sheet1','A2');
     disp('THE END');
end