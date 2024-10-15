%%% D3makeprh_js %%%

%% Set paths and names
tag= 'zc19_218a';
recdir='C:\Users\Jeanne\Documents\DTAG\data\zc19_218a';
prefix='zc218a';
deploy_name='zc19_218a';
settagpath('cal','C:\Users\Jeanne\Documents\DTAG\cal','prh','C:\Users\Jeanne\Documents\DTAG\prh')

%% Read in data
%X=d3readswv(recdir,prefix,df);  %often a df of 10 or 20 is good. Want final sampling rate of 10-20Hz
X=d3readswv(recdir,prefix,10);
warning('off','last') %in matlab 2016b

%% Cal

CAL=d3deployment(recdir,prefix,deploy_name);  %deploy_name = ex. gm10_234a

%% Save location and tag data
savecal(tag,'AUTHOR','jms');
savecal(tag,'TAGON',[2019 08 06 15 38 20]); %[year month day hour minute second]
savecal(tag,'GMT2LOC','4'); %number of hours from GMT time to local time
%savecal(tag,'TAGLOC',[]); %tag on position [lat, long]
%savecal(tag,'DECL', ); %local magnetic field declination in decimal degrees
savecal(tag,'TAGID','tag327');
%savecal(tag,'CAL',CAL);


%% Pressure
[p,CAL,fs,t]=d3calpressure(X,CAL,'full');

%pick all data points on the lower broken line. 
%Redo if 2nd order fitting error is more than 0-20cm

%% Improve pressure calibration (maybe?)
[p,pc]=fix_pressure(p,t,fs);

%% plot to check
plott(p,fs)

%% Acceleration
%[A,CAL,fs]=d3calacc(X,CAL,'full',min_depth); 
%mindepth usually 10 for deep divers

[A,CAL,fs]=d3calacc(X,CAL,'full',2); 
% needs to have standard deviation below 0.04

%% Magnetometer
%[M,CAL]=d3calmag(X,CAL,'full',min_depth);

[M,CAL]=d3calmag(X,CAL,'full',2);
% needs to have standard deviation below 0.5

%Final mag field intensity 43.48 uT (0.712 uT S.D.)

%% Save
d3savecal(deploy_name,'CAL',CAL)

%% PRH predictor (see prhpredictorlmml2.docx for details)
%PRH=prhpredictor(p,A,fs,[TH,METHOD,DIR]); %choose 4-10 cycles. 

PRH=prhpredictor(p,A,fs,200,2,'both');
% look at output for moves

%% For any move, plot A and P to determine if gradual or abrupt

figure
plot(p)
hold on
plot(norm2(A),'r')
axis ij
figure
plot(norm2(A),'r')
%% Orientation table
%only include moves over 10 degrees.  
%Movement table has to be in radians. 
%To convert to radians, use format shortg and then T
%moven = [t1, t2, p0, r0, h0];
%t1 and t2 are start and end times of move in seconds since tagon.
%if move is instantaneous, t1=t2
%p0, r0, h0 are new orientation angles after the move, in radians

move1=[0 0 -1.0288 0.022464 0.86546];
move2=[5377 5377 -0.95368 -0.012243 0.73318];
move3=[11536 11536 -1.1337 -0.43309 0.40839];
move4=[16166 16166 -1.0826 -0.68725 0.078242];



%% Create OTAB File
OTAB=[move1;move2;move3;move4];

%% Comments on OTAB
%comments=['js: The exact timing of the move in the second dive could not be distinguished from the sprint and was considered to be gradual.'];
%savecal(tag,'OTABCOMMENT',comments);

%% Create whale frame acceleration
[Aw, Mw]=tag2whale(A,M,OTAB,fs);


%% Plot Aw and p to check moves are now correct

figure
plot(p)
hold on
plot (norm2(Aw),'r')
axis ij

%% Run prhpredictor on Aw to test moves are correct 
Tw = prhpredictor(p,Aw,fs,200,2,'both');

%% Save OTAB
d3savecal(deploy_name,'OTAB',OTAB)

%% Make pitch, roll, and heading
[pitch, roll]=a2pr(Aw);
[head,mm,incl]=m2h(Mw,pitch,roll); %D3 tool version
%head,mm,incl]=m2h(Mw,Aw,fs); %D4 tool version

if exist('DECL')
    head=head+DECL*pi/180;
else
    disp('Warning:No declination specified')
end

%% save prh file
saveprh(tag,'p','A','M','fs','Aw','Mw','pitch','roll','head') ;

%% Make PRH file ***DOES NOT WORK***
%d3makeprhfile(recdir,prefix,tag,25) %

%d3makeprhfile(recdir,prefix,tag,10)

%% 

j=njerk(Aw,fs);
%% plots
time = (1:length(p))/fs;
figure
ax(1)=subplot(2,1,1);
plot(time,p);grid
set(gca,'Ydir','reverse')
title('Dive Profile')

ax(2)=subplot(2,1,2);
plot(time,pitch*180/pi);grid
hold on
plot(time,roll*180/pi)
plot(time,head * 180 / pi)
hline(0,'k-')
%title('PRH')
legend('Pitch','Roll')

%ime1=time(1:1003487);
%ax(3)=subplot(3,1,3);
%plot(time,j);grid
%title('Jerk')


%ax(4)=subplot(4,1,4);
%plot(e);grid
%title('ODBA')

linkaxes(ax,'x')
%% Plots
figure
ax(1)=subplot(3,1,1)
plot(pitch*180/pi);grid
title('pitch')

ax(2)=subplot(3,1,2)
plot(roll*180/pi);grid
title('roll')

ax(3)=subplot(3,1,3)
plot(head * 180 / pi); grid
title('Heading')

linkaxes(ax,'x')
