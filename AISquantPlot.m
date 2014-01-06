%%%%%%%%%%%%%%%%%%%%%
% MAIN PARAMETERS
%%%%%%%%%%%%%%%%%%%%%

% INPUT txt file should contain 3 columns:
% 1)-2) x and y coordinates
% 3) intensity values


%number of averaging points for smoothing. 
% bigger values correspond to more smoothing 
nAverPoints = 21; 
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
[filename,filepath] = uigetfile('.log','Select txt file with x,y coordinates and intensity profile...');
%if cancel button is not pressed
if ~isequal(filename, 0)
        disp('Script started.');
        disp('Please wait......');
        filenamein = sprintf('%s%s',filepath,filename);    
        %read intensity trace
        inttrace = csvread(filenamein);
        sz = size(inttrace);
        %total number of points
        nTotPoints = sz(1);
        %array containing results
        norm_int = zeros(nTotPoints,2,'double');
        for i = 2:nTotPoints            
            %distance along track
            norm_int(i,1) =  norm_int(i-1,1) + sqrt((inttrace(i,1)-inttrace(i-1,1))^2+(inttrace(i,2)-inttrace(i-1,2))^2)/nPixMkm;
        end
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
                norm_int(i,2) = mean(inttrace(1:(i+nHalf),3));
            end
            %end
            for i=(nTotPoints-nHalf+1):nTotPoints
                norm_int(i,2) = mean(inttrace((i-nHalf):nTotPoints,3));
            end
        end
        %smoothing middle part
        for i=(nHalf+1):(nTotPoints-nHalf)
            norm_int(i,2) = mean(inttrace((i-nHalf):(i+nHalf),3));
        end
        
        %NORMALIZATION
        %min intensity
        minInt = min(norm_int(:,2));
        %max-min intensity, normalization factor
        dSpan = max(norm_int(:,2)) - minInt;
        %normalize intensity
        norm_int(:,2) = (norm_int(:,2) - minInt) / dSpan;
        
        %SEARCH FOR REFERENCE POINTS
        %Max
        [dSpan,MaxInd] = max(norm_int(:,2));
        MaxVal = norm_int(MaxInd,1);
        %Start
        i=0;
        dSpan = 0;
        while dSpan < dThreshold
            i = i + 1;
            dSpan = norm_int(i, 2);
        end
        StartVal = norm_int(i, 1);
        
        %End
        i=nTotPoints+1;
        dSpan = 0;
        while dSpan < dThreshold
            i = i - 1;
            dSpan = norm_int(i, 2);
        end
        EndVal = norm_int(i, 1);
        
        %Calculate AIS length
        LengthVal = EndVal - StartVal;
        
        %Make plot
        plot((norm_int(:,1)),(norm_int(:,2)),'Color','b');
        line([StartVal,StartVal],[0,1],'LineStyle','--', 'Color','b');
        line([EndVal,EndVal],[0,1],'Linestyle','--','Color','b');
       
        %saving results data
        filenamein =  strcat(filenamein(1:length(filenamein)-4),'_results.xls');
        headers = {'length,mkm', 'normalized intensity', 'Start, mkm','Max, mkm','End, mkm', 'AISlength, mkm', '# of smooth points', 'intensity threshold','resolution (pixels in mkm)', 'Tails'};
        xlswrite(filenamein,headers,'Sheet1','A1');
        params = [StartVal MaxVal EndVal LengthVal nAverPoints dThreshold nPixMkm bTails];
        xlswrite(filenamein,params,'Sheet1','C2');
        xlswrite(filenamein,norm_int,'Sheet1','A2');
        disp('Done.');
end