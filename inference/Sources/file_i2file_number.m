function [file_number] = file_i2file_number(file_i)
    if file_i < 10
        file_number = ['00' num2str(file_i)];
    elseif file_i > 9 && file_i < 100
        file_number = ['0' num2str(file_i)];
    elseif file_i > 99 && file_i < 1000
        file_number = num2str(file_i);
    end
end