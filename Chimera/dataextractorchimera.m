%% Geometric shape parameters are calculcated for the configuration with the
%%largest mouth area.
%%Outputs [Force Tau Gamma Iter MinArea MaxArea MeanRadius StdRadius MajorRadius MinorRadius MajorAxisOrientation OpeningTime OpeningTimeNF OpeningTimeEN flag COMx COMy] as a row vector;
%%Flag outputs the following depending on the criteria:
%%If MaxArea > 1.2*(Area of outer ring), flag = -1, Hard condition
%%If MaxArea < 1.05*(Initial Mouth Area), flag = 0, Soft condition
%%If OpeningTime < dt*10 = 0.05, flag = -2, Hard condition
%%If OpeningTime > 2*T, flag = 1, Soft Condition
%%If EllipseFit reports NaN (cannot fit), flag = 3, Soft/Hard Condition
%%For all else, flag = 2: Non pathalogical cases
%%Code used for ellipse fitting: Nikolai Chernov (2021). Ellipse Fit (Direct method) (https://www.mathworks.com/matlabcentral/fileexchange/22684-ellipse-fit-direct-method), MATLAB Central File Exchange. Retrieved December 30, 2021.

function SimResults = dataextractorchimera(fullfilename,filename) 

%% Extract parameters force,tau and realization# from the filename%%
Str = filename;
Str(strfind(Str, '=')) = [];
Key = 'force';
Index = strfind(Str, Key);
Force = sscanf(Str(Index(1) + length(Key):end), '%g', 1);
Key   = 'tau';
Index = strfind(Str, Key);
Tau = sscanf(Str(Index(1) + length(Key):end), '%g', 1);
Key   = 'i';
Index = strfind(Str, Key);
Iter = sscanf(Str(Index(1) + length(Key):end), '%g', 1);
Key   = 'gam';
Index = strfind(Str, Key);
Gamma = sscanf(Str(Index(1) + length(Key):end), '%g', 1);

%% Read file data and get simulation outputs %%
load(fullfilename);

%% Obtain Areas for each timestep

    %% Get Mouth areas in each frame
    Area = polyarea(squeeze(displ(1,:,1,:)),squeeze(displ(1,:,2,:)));
    NFcoords = zeros(6,size(displ,3),size(displ,4));
    NFcoords(2:5,:,:) = squeeze(displ(1,2:5,:,:));
    NFcoords(1,:,:) = squeeze(.5.*(displ(1,1,:,:)+displ(1,2,:,:)));
    NFcoords(6,:,:) = squeeze(.5.*(displ(1,5,:,:)+displ(1,6,:,:)));
    AreaNF = polyarea(squeeze(NFcoords(:,1,:)),squeeze(NFcoords(:,2,:)));
    AreaEN = Area-AreaNF;

%% Obtain Reference areas and times
% Lower limit on Max Area of mouth:
MaxAreaLowerLimit = 1.05*Area(1); %the maximumu mouth area should atleast be 5% greater than initial mouth area
%Upper limit on Max Area of mouth (wider than outer fixed boundary ring):
x = displ(size(displ,1),:,1,1); %Coordinates of the boundary 
y = displ(size(displ,1),:,2,1);
x = x(:);
y = y(:);
MaxAreaUpperLimit = 1.2*polyarea(x,y);

% Upper and lower limits for time of opening
OpeningTimeLowerLimit = 10*(time(2)-time(1));
OpeningTimeUpperLimit = time(end)*2;

%% Obtain Min, Max area and OpeningTime for sim.
MinArea = min(Area);
MaxArea = max(Area);
IndexMaxArea = find(Area == max(Area));
rx = displ(1,:,1,IndexMaxArea);
ry = displ(1,:,2,IndexMaxArea);
rx = rx(:);
ry = ry(:);
rx = rx-mean(rx);
ry = ry-mean(ry);
r = sqrt(rx.*rx +ry.*ry);
MeanRadius = mean(r);
StdRadius = std(r);

%%Mouth opening time
NormArea = (Area-min(Area))./(MaxArea-MinArea);
[xData, yData] = prepareCurveData( time, NormArea );
if(length(yData)>50)
    
        ft = fittype( 'a+b/(1+exp(-(x-c)/d))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Lower = [-0.1 0.9 -Inf -Inf];
        [ d, ix ] = min( abs( NormArea-0.5 ) ); %%ix is index of time vector where area crosses 0.5
        [ d, ixlow ] = min( abs( NormArea-0.25 ) ); %%ix is index of time vector where area crosses 0.5
        [ d, ixhigh ] = min( abs( NormArea-0.75 ) ); %%ix is index of time vector where area crosses 0.5
        opts.StartPoint = [0 1 time(ix(1))  0.25*(time(ixhigh(1))-time(ixlow(1)))];
        opts.Upper = [0.1 1.1 Inf Inf];
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts );
        OpeningTime = fitresult.d;
        
        %%NF opening time
        NormAreaNF = (AreaNF-min(AreaNF))./(max(AreaNF)-min(AreaNF));
        [xDataNF, yDataNF] = prepareCurveData( time, NormAreaNF );
        ft = fittype( 'a+b/(1+exp(-(x-c)/d))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Lower = [-0.1 0.9 -Inf -Inf];
        [ d, ix ] = min( abs( NormAreaNF-0.5 ) ); %%ix is index of time vector where area crosses 0.5
        [ d, ixlow ] = min( abs( NormAreaNF-0.25 ) ); %%ix is index of time vector where area crosses 0.5
        [ d, ixhigh ] = min( abs( NormAreaNF-0.75 ) ); %%ix is index of time vector where area crosses 0.5
        opts.StartPoint = [0 1 time(ix(1))  0.25*(time(ixhigh(1))-time(ixlow(1)))];
        opts.Upper = [0.1 1.1 Inf Inf];
        % Fit model to data.
        [fitresult, gof] = fit( xDataNF, yDataNF, ft, opts );
        OpeningTimeNF = fitresult.d;
        
        %%EN opening time
        NormAreaEN = (AreaEN-min(AreaEN))./(max(AreaEN)-min(AreaEN));
        [xDataEN, yDataEN] = prepareCurveData( time, NormAreaEN );
        ft = fittype( 'a+b/(1+exp(-(x-c)/d))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Lower = [-0.1 0.9 -Inf -Inf];
        [ d, ix ] = min( abs( NormAreaEN-0.5 ) ); %%ix is index of time vector where area crosses 0.5
        [ d, ixlow ] = min( abs( NormAreaEN-0.25 ) ); %%ix is index of time vector where area crosses 0.5
        [ d, ixhigh ] = min( abs( NormAreaEN-0.75 ) ); %%ix is index of time vector where area crosses 0.5
        opts.StartPoint = [0 1 time(ix(1))  0.25*(time(ixhigh(1))-time(ixlow(1)))];
        opts.Upper = [0.1 1.1 Inf Inf];
        % Fit model to data.
        [fitresult, gof] = fit( xDataEN, yDataEN, ft, opts );
        OpeningTimeEN = fitresult.d;
        
        
        %% calculate COM of mouth opening in the end

        COMx = mean(squeeze(displ(1,:,1,end)));
        COMy = mean(squeeze(displ(1,:,2,end)));

        %% Check conditions on Area and OpeningTime
        %%Flag outputs the following depending on the criteria:
        %%If MaxArea > 1.2*(Area of outer ring), flag = -1, Hard condition
        %%If MaxArea < 1.05*(Initial Mouth Area), flag = 0, Soft condition
        %%If OpeningTime < dt*10 = 0.05, flag = -2, Hard condition
        %%If OpeningTime > 2*T, flag = 1, Soft Condition
        %%If EllipseFit reports NaN (cannot fit), flag = 3, Soft/Hard Condition
        %%For all else, flag = 2: Non pathalogical cases

        flag = 2; %%Assuming no pathology to begin with
        if(MaxArea > MaxAreaUpperLimit)
            flag = -1;
            SimResults = [Force Tau Gamma Iter MinArea MaxArea MeanRadius StdRadius nan nan nan nan nan nan flag COMx COMy];
            %return ;
        end

        if(OpeningTime < OpeningTimeLowerLimit)
            flag = -2;
            SimResults = [Force Tau Gamma Iter MinArea MaxArea MeanRadius StdRadius nan nan nan OpeningTime OpeningTimeNF OpeningTimeEN flag COMx COMy];
            %return ;
        end

        if(MaxArea < MaxAreaLowerLimit)
            flag = 0;   
        end

        if(OpeningTime > OpeningTimeUpperLimit)
            flag = 1;
        end

            A = EllipseDirectFit([rx ry]);
            if(~isnan(A))
                B = -det([A(1) .5*A(2) .5*A(4);.5*A(2) A(3) .5*A(5);.5*A(4) .5*A(5) A(6)]);
                J = -det([A(1) 0.5*A(2);0.5*A(2) A(3)]);
                I = A(1)+A(3);
                MajorRadius = sqrt(2*B/(J*(-I+sqrt(I^2+4*J))));
                MinorRadius = sqrt(2*B/(J*(-I-sqrt(I^2+4*J))));
                if(A(1)<A(3))
                    MajorAxisOrientation = 0.5*acot((A(1)-A(3))/A(2));
                else if(A(1)>A(3))
                    MajorAxisOrientation = 0.5*acot((A(1)-A(3))/A(2))+pi/2;
                    else
                        MajorAxisOrientation = 0;
                    end
                end
            else
                flag = 3;
                MajorRadius = nan; MinorRadius = nan; MajorAxisOrientation = nan;
            end

        SimResults = [Force Tau Gamma Iter MinArea MaxArea MeanRadius StdRadius MajorRadius MinorRadius MajorAxisOrientation OpeningTime OpeningTimeNF OpeningTimeEN flag COMx COMy];
else
       SimResults = [Force Tau Gamma Iter nan nan nan nan nan nan nan nan nan nan -20 nan nan];
end
end

