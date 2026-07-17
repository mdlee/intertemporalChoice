function F = setFigure(F, position, dual)
%SETFIGURE sets a figure in a way that copes with dual monitor changes
%    setFigure(F, position, dual)
%    F is figure handle
%    position is normalized [left bottom width height'
%    dual is 'upDown' or 'leftRight' depending on dual monitor arrangement
[~, name] = system('hostname');
tmp = get(0, 'MonitorPositions');
%tmp = tmp./[2 1 2 1];
if strcmp(name(1:end-1), 'mdl2025XX')
 %  disp('y')
   vec(1) = (tmp(3)-tmp(1))*position(1);
   vec(2) = (tmp(4)-tmp(2))*position(2);
   vec(3) = (tmp(3)-tmp(1))*position(3);
   vec(4) = (tmp(4)-tmp(2))*position(4);
else
   
   [rows, ~] = size(tmp);

   if rows == 2
      dual = 'leave';
   elseif (tmp(3) == 1920 && tmp(4) == 1200) || (tmp(3) == 1920 && tmp(4) == 1080)
      dual = 'single';
   elseif tmp(4) == 1080 || tmp(4) == 1200 || tmp(4) == 1050 || tmp(4) == 2560 || tmp(4) == 1440
      dual = 'leftRight';
   elseif tmp(3) == 1920 || tmp(3) == 3440
      dual = 'upDown';
   else
      error('not sure about monitors');
   end
   % dual
   switch dual
      case 'single'
         vec(1) = (tmp(3)-tmp(1))*position(1);
         vec(2) = (tmp(4)-tmp(2))*position(2);
         vec(3) = (tmp(3)-tmp(1))*position(3);
         vec(4) = (tmp(4)-tmp(2))*position(4);
      case 'upDown'
         vec(1) = (tmp(3)-tmp(1))*position(1);
         vec(2) = (tmp(4)-tmp(2))*position(2)/2;
         vec(3) = (tmp(3)-tmp(1))*position(3);
         vec(4) = (tmp(4)-tmp(2))*position(4)/2;
      case 'leftRight'
         vec(1) = (tmp(3)-tmp(1))*position(1)/2;
         vec(2) = (tmp(4)-tmp(2))*position(2);
         vec(3) = (tmp(3)-tmp(1))*position(3)/2;
         vec(4) = (tmp(4)-tmp(2))*position(4);
   end
end

% if tmp(4) == 1440
% vec(1) = (tmp(3)-tmp(1))*position(1)/2;
%         vec(2) = (tmp(4)-tmp(2))*position(2);
%         vec(3) = (tmp(3)-tmp(1))*position(3)/1.5;
%         vec(4) = (tmp(4)-tmp(2))*position(4);
% end

% vec
% dual
switch dual
   case 'leave'
      set(F, ...
         'color'             , 'w'                , ...
         'units'             , 'normalized'           , ...
         'position'          ,  position               , ...
         'papersize'         , [10 6.75]          , ...
         'paperpositionmode' , 'auto'             );

   otherwise
      set(F, ...
         'color'             , 'w'                , ...
         'units'             , 'pixels'           , ...
         'position'          ,  vec               , ...
         'papersize'         , [10 6.75]          , ...
         'paperpositionmode' , 'auto'             );
end

