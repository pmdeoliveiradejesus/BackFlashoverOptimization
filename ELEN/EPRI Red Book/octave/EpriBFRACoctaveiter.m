clc
clear all
w=[200 225 250 263 275 300 325 350 375 400 425 450 475 500 525 550 575 600]/100;% Aislamiento en cm
 for  h=1:length(w) 
     P(1,h)=EpriBFRACoctave(20,w(h),3.6);
 end
 