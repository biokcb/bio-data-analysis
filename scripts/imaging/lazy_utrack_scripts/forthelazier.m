function forthelazier

%Trying to improve the previous script by making it more user friendly. In
%theory, you won't have to edit the script at all to get it to do the
%sequential uTrack analysis of your movies. You still have to edit the
%files "scriptDetectLazy" and "scriptTrackLazy" for the proper/desired
%parameter values. I'm also writing this with many comments like I normally
%do, in case you want to learn some basic code. 
% written 2013 by Kristen Brown

%First, let's find out if the user wants to do both detection and tracking
%for their results or just tracking. The command below sets up a dialog box
%with a question, title, and 3 options - yes, no, and cancel. The switch
%will then take the input and execute commands based on which one is
%selected.

%The default value is set to 'Yes' here, so that pressing enter will just
%go to the Detection/Tracking part of the script.
reply = questdlg('Do you want to do detection and tracking?', 'Detection/Tracking', 'Yes', 'No', 'Cancel', 'Yes');

switch reply
    case 'Cancel', %User decided to cancel the execution of the program. 
        disp('Ended program. Exciting...'); %added this dialog just for Liz :)
        quit cancel; %This will exit the run. 

%--------------------------Tracking Only----------------------------------%
%

    case 'No', %The user decided to only do tracking, now let's ask for the location of detection files and go from there.
        
%Step 1: Gather user input

    %asks the user to open the detection results, sets this =c, and the path =a.
    %c is now an array of the file names if there are >1 file open.
    [c, directory] = uigetfile('', 'Open all of the detection results you want to track', 'Multiselect', 'on'); 
    
    a = ['/' directory '/'];
    
    %This now asks the user if they want to create a new folder to save the
    %tracking results in. This would be ideal for wanting to test out
    %multiple tracking parameters for one set of detection results. If you
    %are only going to do one set of tracking results, just hit 'No' so
    %that it'll be saved in the same folder.
    folder = questdlg('Do you want to create a new folder to save the tracking results?', 'Saving Tracking Results', 'Yes', 'No', 'Cancel', 'Yes');
    
    %Determines what to do based on the button the user presses. 
    switch folder
        case 'Cancel' %once again cancel will end the run.
            disp('Ended program. Run script again to start over');
            quit cancel;
            
        case 'Yes' %Yes will open a new dialog box that asks for the user to create a new folder name
            %it would be ideal to name this folder using terms describing
            %the parameters you used for this tracking run.
            v = inputdlg('Enter name of the folder you want to create', 'Results will be saved in a new folder within detection results folder', 1);
            newfolder = (v{1});
            mkdir(directory, newfolder); %creates the new folder
            asave = [directory, '/' newfolder]; %sets the folder for the results to be saved in to the new folder
            
        case 'No' %No will just set the save folder equal to the folder that the detection results are stored in.
            disp('Results will be saved in the same folder');
            asave = a;
    end
    
%Step 2: Open the detection files and run the tracking script.

    %Figures out the number of files you selected to do tracking on by
    %looking at the size of the array the filenames are stored in. Sets the
    %value i = 1, which will be the start of the loop.
    numfiles = numel(c);
    i = 1;
    
    %The while loop will run while i is less than or equal to the size of
    %the array 'c' which is the same as the number of detection files it
    %needs to load.
    while i <= numfiles
        
        %this establishes where the files are located so we can load them
        %for tracking by letting "file" equal to the path followed by the 
        %array element of 'c' corresponding to the i we are on in the loop.
        file = ([directory '/' c{i}]);
        
        load(file); %loads the detection results based on the value of i
        
        %sets the value for the name of the saved tracking .mat file to be
        %equal to 'Tracking_firstframe-lastframe.mat'. Example:
        %Tracking_1-2000.mat, depending on what the numbers are for the
        %loaded detection file.
        d = ['Tracking_' num2str(movieParam.firstImageNum) '-' num2str(movieParam.lastImageNum) '.mat']; 
        
        %runs the tracking script.
        scriptTrackLazy;
        
        %updates i to the next value, so the loop will start over loading
        %the second detection file of the 'c' array.
        i = i + 1;
    end

    
%--------------------------Detection & Tracking --------------------------%

    case 'Yes', %The user decided to do both the tracking and detection steps for this movie. 
    
%Step 1: Gather input from the user. 

    %Ask for the folder that houses the .tif files of the movie. Set that equal
    %to 'a'. 
        directory = uigetdir('','Please select the location of the image sequence you wish to analyze: ');

    %This creates a dialog box that asks for 4 more inputs from the user
    
        x = inputdlg({'Name of folder to be created for saved results', 'File name base of image sequence', 'Total # of frames to analyze', '# of frames per run', 'First frame to analyze'}, 'User Input', 1, {'uTrack Analysis 1', 'cell_1_', '10000', '2000', '1'}); 
    
    %Name of folder to be created for saved results - default name is
    %"uTrack Analysis 1", but it is encouraged that the user use a more
    %descriptive name with parameter value descriptions. This folder will
    %then be created in the same folder as the image sequence. 
    
    %File name base of image sequence - This is the name of the files in
    %the image sequence before the numbers. The default is "cell_1_"
    %because that's one I've used before. The default can be changed in the
    %code if there is a name you use most.
    
    %Total # of frames to analyze - This is just asking for the # of frames
    %you want in the analysis. Typically we do the entire movie, but you
    %have the option of doing just a chunk of the movie. Default is 10000
    %frames
    
    %# of frames per run - This is asking for the number of frames you want
    %the movie to divided into. Default is 2000 frames, as that is what we 
    %most often do. If the total # of frames is not divisible by the number
    %you choose, that's ok because we'll establish a conditional statement
    %about this later.
    
    %Set the value of the file base name =b so detection can run.
        b = x{2};

    %Makes the new folder for BOTH detection/tracking results within the
    %folder that the images are housed.
        newfolder = (x{1});
        mkdir(directory, newfolder);
    
    %Set the value of folder for saved results equal to asave and the input
    %directory equal to the value of the input from directory.
        asave = [directory '/' newfolder];
        a = [directory '/'];

disp(asave);
%Step 2: Figure out what the value of 'h' is. Assuming that h can only be 4 or 5
%since our movies are always >1,000 frames and <100,000 frames

        if str2double(x{3}) < 10000
            h = 4;    
        else 
            h = 5;
        end

%Step 3: Run the loop until the frames are done being analyzed.

        f = str2double(x{5}); %set start values for the # of frames. First frame is set to 1 to begin with.
        g = str2double(x{4}); %Last frame is set to the number of frames the movie is divided into. Default: 2000. Detection with start with frames 1-2000


        while g < (str2double(x{3}) + 1) && f < g %loops while the value of 'g' (current max frame #) is less than or equal to the total number of frames.

            disp([f g]); %I used this for debugging the code. It seemed to not work at first due to updating g before f so f>g and it got confused.
               %I'll leave it in, so if you're checking on it, you can see which
               %frames are being analyzed. And so you know the command for basic
               %debugging - it'll display what values you ask it to while the
               %program is running (placement within the loop is important for
               %values that update). Then you can see if anything is wrong during a
               %test run like I did.

            %set the names for the files to be saved as
            %Detection_firstframe-lastframe.mat or
            %Tracking_firstframe-lastframe.mat. For example Detection_1-2000.mat
            c = ['Detection_' num2str(f) '-' num2str(g) '.mat'];
            d = ['Tracking_' num2str(f) '-' num2str(g) '.mat'];
            %Since these names aren't descriptive about the parameters, this is why
            %I made the program create new folders. All the parameters are saved in
            %the workspace regardless, but it's nice to not have to open the files
            %to find out the basic information.

            %Run the scripts...
            scriptDetectLazy;
            scriptTrackLazy;

            %Update the values of 'f' and 'g' so the loop can start again and eventually end.

            f = (g + 1); %we don't need a condition on 'f', since the loop only depends on 'g' being less than the total. 
            %This will update it to the next start value. For example, 2001 from 1.
            %Just as a note, the loop would update 'g' first, then 'f' when the
            %expression didn't have parentheses. Adding them in fixed the problem
            %and now 'f' will update first, weirdly enough. 

            %We do need a condition on 'g', in case the movie itself is not
            %divisible by the #of frames used per run. 

            %This will establish that if the total number of frames - the segment
            %is less than the value of 'g', then 'g' will update to the final value
            %Example: Say, total # of frames is 9000. 

            %If g = 2000, and the segments are 2000 frames, 9000-2000 = 7000. Which
            %is greater than g = 2000, so the condition is not met and this value
            %is skipped.

            %However, when the loop reaches g = 8000 in this case and goes to
            %update the value, 9000-2000 = 7000, which is less than 'g', so 'g'
            %becomes 9000 because the condition is met.
            if (str2double(x{3}) - str2double(x{4})) < g
                g = str2double(x{3});

            %For everything else (ie, total # of frames - the segment is greater
            %than the current value of 'g'), the value of 'g' updates to the next segment.
            %Example: from 2000 to 4000, for segments of 2000
            else
            g = g + str2double(x{4});

            end %end of if statement
        end %end of while loop
end %end of switch statement

%Just for fun. Program will display a random phrase at the end.
y = {'We did it!' 'Huzzah! Task complete!' 'Done!' 'Finished!' 'I am so happy to be done!' 'That was a tough one!' 'Man, you are being lazy today!' 'I hope you had a good day!' 'That was fun! Let''s do it again!' 'Let''s not do this again' 'I''m tired, also feed me.' 'Say hi to Phil for me!' 'Say hi to Mike for me!' 'Hi Liz!' 'Have I told you how awesome Kristen is?' 'I hope you remembered to check the parameters before starting' 'I''m running out of things to say' 'There''s someone behind you.' 'I have a good feeling about this one...' 'How much coffee have you had today?'};
ysize = numel(y); %get the number of elements in the array 'y'
idx = (1 + fix(rand(1,1)*ysize)); %pick a random index number
disp(y(idx)); %display random element of y

