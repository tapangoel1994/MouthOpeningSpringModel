%%Function takes the name of a single .mat file containing the positions of
%%each node at each time point  and creates a plot of the area as a
%%function of time.
%%it saves the plot as a fig and as a png in a new subfolder with the same name input file
%%inputfilename is the full name of the file with full file path
function PlotSuccess = PlotTimevsArea(inputfilename)

        load(inputfilename,'time','displ');
        Area = zeros(length(time),1);
        for iter = 1:1:length(time)
                        x = displ(1,:,1,iter);
                        y = displ(1,:,2,iter);
                        x = x(:);
                        y = y(:);
                        Area(iter) = polyarea(x,y);
                end
       fig = figure('visible','off');
       
       plot(time,Area,'-o');
       xlabel('Time(s)');
       ylabel('Area(sq.um)');
       title(inputfilename)
       figfilename = [erase(inputfilename,'.mat') '.fig'];
       pngfilename = [erase(inputfilename,'.mat') '.png'];
       
       saveas(fig,figfilename);
       saveas(fig,pngfilename);

           
            
end

