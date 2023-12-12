%Coeficient matrix
function TT=tra(L,di)
            TT=[L(1)                          L(2)                             L(3);
               di(1)                        di(2)                           di(3);
               di(2)*L(3)-di(3)*L(2) -(di(1)*L(3)-di(3)*L(1)) di(1)*L(2)-di(2)*L(1)];            
end