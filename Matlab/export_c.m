function export_c(x)
    fprintf('Eigen::VectorXd solution(%d); \n', length(x)); 
    for i = 1:length(x)
        fprintf('solution(%d) = %.10f; \n', i-1, x(i));
    end
end