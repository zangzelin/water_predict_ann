function [ J ] = gotJ( e,w )
%求取雅各比矩阵
%w，e输入，输出雅各比
if(length(e) == 2)
    c = 1;
else
    c = length(e);
end

for i = 1:c
    
    if i == 1
        de = e(2)-e(1);
    elseif i == length(e)
        de = e(length(e))-e(length(e)-1);
    else
        de = (e(i)-e(i-1)+e(i+1)-e(i))/2;
    end
        
    for j = 1:length(w)

            if j == 1
                dw = w(2)-w(1);
            elseif j == length(w)
                dw = w(length(w))-w(length(w)-1);
            else
                dw = (w(j)-w(j-1)+w(j+1)-w(j));
            end
            
            J(i,j) = de/dw;
        
    end
end


end

