% index =   LOOKUP(values, lookupvalues)
% index =   LOOKUP(values, lookupvalues, direction)
% This function will return the indices of lookupvalues that are closest to
% the numbers in values. The function assumes that LOOKUPVALUES are in
% order from smallest to largest.
%direction of -1 forces lookup to only look backwards in lookupvalues. If
%direction is 0 (default), however, it chooses the closest indeces. If it
%is 1 it looks foreward in lookupvalues.  If a value falls outside the
%range of lookupvalues, it always chooses the closest index, regardless of
%the value of direction.