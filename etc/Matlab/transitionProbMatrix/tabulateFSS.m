function table = tabulateFSS(x)
%TABULATE Frequency table.
%   TABLE = TABULATE(X) takes a vector of positive integers, X,  
%   and returns a matrix, TABLE. 
%   The first column of TABLE is the values of X. The second
%   is the number of instances of this value. The last column
%   contains the percentage of each value.
%   TABULATE with no output arguments returns a formatted table
%   in the command window.
%   See also PARETO.
   
%   B.A. Jones 3-5-95
%   Copyright 1993-2000 The MathWorks, Inc. 
%   $Revision: 2.9 $  $Date: 2000/05/26 18:53:44 $

if min(size(x)) > 1,
   error('Requires a vector input.');
end

y = x(find(~isnan(x)));

if any(y ~= round(y)) | any(y < 1),
   error('Requires the values of the input to be positive integers.');
end 

maxlevels = max(max(y));
[counts values] = hist(y,(1:maxlevels));

total = sum(counts);
percents = 100*counts./total;
tmp = [values; counts; percents]';

if nargout == 0
   [m,n] = size(tmp);
   s = 37;
   s = s(ones(m,1),:);
   disp('  Value    Count    Percent');
   fprintf(1,'  %5d    %5d    %6.2f%c\n',[tmp'; s']);
else
   table = tmp;
end

