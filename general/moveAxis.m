function moveAxis(handle, scale, shift)
% moveAxis(handle, scale, shift)
%   scale is 4-vector of scalar multiple for left, bottom, width, height
%   shift is 4-vector of scalar addition for left, bottom, width, height
set(handle, 'position', get(handle, 'position').*scale + shift);
end

