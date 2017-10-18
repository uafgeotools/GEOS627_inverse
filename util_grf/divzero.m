% This makes dividing by matrices element-wise easier, as we don't have to
% worry about the 0-modes. If we must divide by 0, it leaves the mode
% alone. 

function a = d(b)

a = b + (b==0);

