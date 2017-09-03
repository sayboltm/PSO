function [ time ] = timeStamp(display)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%f = get(0, 'Format');
time = clock;
year = time(1);
month = time(2);
date = time(3);
hour = time(4);
minute = time(5);

if display == 1
    disp (['[',num2str(month),'/',num2str(date),'/',num2str(year), ...
        ' ', num2str(hour),':',num2str(minute),']'])
end

%format f
end

