function doseresponse

%generates a dose-response curve with error bars from Alamar Blue Data
%saved as an Excel Spreadsheet

%Step 1: Ask for the file with the fluoresence data.

[filename, directory] = uigetfile('', 'Open File with Alamar Blue Data');

%Step 2: Ask how the plate was split up (assuming 3 drugs/doses), read in
%the data.

plates = inputdlg({'Drug Name & Dosage', '# of Replicates', 'Drug Name & Dosage', '# of Replicates', 'Drug Name & Dosage', '# of Replicates'}, 'Plate Information', 1, {'Vincristine 10ng/ml, 1:10', '3', '', '3', '', '4'});


if str2double(plates{6}) < 4
    xlsrange = 'C3:K8';
else 
   xlsrange = 'C3:L8';
end

alb = xlsread(filename, xlsrange);

%Step 3: Average the replicates, average the control lane. 

cntrlavg = sum(alb(6,:))/sum(str2double(plates{2}), str2double(plates{4}), str2double(plates{6}));


while k < 5
    drug1(k) = sum(alb(k,1:(str2double(plates{2}))));
    
    drug2(k) = sum(alb(k,(1+str2double(plates{2})):(str2double(plates{2})+str2double(plates{4}))));
    
    drug3(k) = sum(alb(k,(1+(str2double(plates{2})+str2double(plates{4}))):(str2double(plates{2})+str2double(plates{4})+str2double(plates{6}))));
   
   
    k = k + 1; %update k

end

avg1 = drug1./str2double(plates{2});
avg2 = drug2./str2double(plates{4});
avg3 = drug3./str2double(plates{6});

%Step 4: Calculate the fraction of control (averaged drug data / averaged
%control)

fc1 = avg1./cntrlavg;
fc2 = avg2./cntrlavg;
fc3 = avg3./cntrlavg;

%Step 5: Calculate the standard deviation



%Step 6: Plot the graphs, with error bars

hold on
plot(fc1); 
errorbar();
hold off

hold on
plot(fc2);
errorbar();
hold off

hold on
plot(fc3);
errorbar();
hold off