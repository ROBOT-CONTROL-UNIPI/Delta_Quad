function v = vee(M)
    
    errID = 'vee:NotSkewSym';
    msg = 'You passed a non skew symmetric matrix.';
    baseException = MException(errID,msg);
    for i = 1:3
        for j = i:3
            if (M(i, j) ~= -M(j, i))
                throw(baseException)
            end
        end 
    end 

    v = [-M(2,3), M(1,3), -M(1,2)].';

end