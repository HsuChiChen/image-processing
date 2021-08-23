
function error=Error_Calculation(X1,X2)

[row,col,band]=size(X1);
L=row*col;
image=reshape(X1,band,L);
image_recover=reshape(X2,band,L);
error =(norm(image-image_recover,'fro')^2)/(L*band);

end
