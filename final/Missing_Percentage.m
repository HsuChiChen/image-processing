function percentage=Missing_Percentage(X)

[row,col,band]=size(X);
percentage=100*length(find(X==0))/(row*col*band);

end