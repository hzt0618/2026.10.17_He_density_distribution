function [doy,yr]=date2doy(idate)
%    ------------------------------------------------------------
%   +              18-Sep-91
%    NAME:
%      date2doy
%    PURPOSE:
%      Convert yyyymmdd into doy (day of the year).  
%      Input can be either string or integer or float.
%      The year is an Optional return parameter.
%    CALLING SEQUENCE:(µ÷ÓÃË³Ðò)
%      date2doy, idate
%    INPUT:
%      idate input format for the date: yyyymmdd.
%       Data-type can be either string or integer or float.
%    OUTPUT:
%      doy double float with the day of the year.
%      yr  double float with the year of the returned doy.
%    HISTORY:
%      written by LZW 21-Sep-2020
%   -
%    -----------------------------------------------------------------
%    ON_ERROR, 2 ;force a return to caller on error
% 
%    Check data type of input set ascII flag and convert to yyyy,mm,dd:
  
  if isa(idate,'char')    %idate string
      sdate=idate;    %scalar flag set
  elseif isa(idate,'numeric')      %idate longint
      sdate=num2str(idate);      %vector input
  end
  
  if length(sdate)==8
      yyyy=str2double(sdate(1:4));
      mm=str2double(sdate(5:6));
      dd=str2double(sdate(7:8));
      
      imonth=[0,31,28,31,30,31,30,31,31,30,31,30,31];
      
      if mod(yyyy,4)==0 && mod(yyyy,400)~=0
          imonth(3)=29;
      end
      
      %   if size(sdate,2)~=8 || mm<1 || mm>12 || dd<1 ...
      %           || dd>imonth(mm+1)
      %       disp('Error 2:force a return to caller on error');
      %       pause;
      %   end
      
      doy=sum(imonth(1:mm))+dd;
      yr=yyyy;
  end
  if length(sdate)==7
      yr=str2double(sdate(1:4));
      doy=str2double(sdate(5:7));
  end
  
end