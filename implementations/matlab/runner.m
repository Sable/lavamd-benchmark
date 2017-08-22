function runner(dim_boxes1d_arg)
% Example: runner(6);

if dim_boxes1d_arg < 0
    disp('ERROR: Wrong value to -boxes1d parameter, cannot be <=0');
    return;
end

disp('Configuration used:');
disp('cores = 1');
disp('boxes1d = '); disp(dim_boxes1d_arg);

NUMBER_PAR_PER_BOX = 100;
expected_boxes1d   = 6;

expectedAns = [4144561.0, 181665.0, -190914.0, 140373.0];
par_cpu_alpha = 0.5;
dim_number_boxes = dim_boxes1d_arg * dim_boxes1d_arg * dim_boxes1d_arg;
dim_space_elem = dim_number_boxes * NUMBER_PAR_PER_BOX;

% allocate boxes
box_cpu_x = zeros(1, dim_number_boxes);
box_cpu_y = zeros(1, dim_number_boxes);
box_cpu.z = zeros(1, dim_number_boxes);
box_cpu_nn     = ones(1, dim_number_boxes);
box_cpu_number = zeros(1, dim_number_boxes);
box_cpu_offset = zeros(1, dim_number_boxes);
box_cpu_nei_x  = zeros(dim_number_boxes, 26);
box_cpu_nei_y  = zeros(dim_number_boxes, 26);
box_cpu_nei_z  = zeros(dim_number_boxes, 26);
box_cpu_nei_number  = zeros(dim_number_boxes, 26);
box_cpu_nei_offset  = zeros(dim_number_boxes, 26);

nh = 1;
for i = 0:dim_boxes1d_arg - 1
    for j = 0:dim_boxes1d_arg - 1
        for k = 0:dim_boxes1d_arg - 1
            box_cpu_x(nh) = k;
            box_cpu_y(nh) = j;
            box_cpu_z(nh) = i;
            box_cpu_number(nh) = nh;
            box_cpu_offset(nh) = (nh - 1) * NUMBER_PAR_PER_BOX;
            
            for u = -1:1
                for m = -1:1
                    for n = -1:1
                        if and(and(i+u >= 0, j+m >= 0), k+n >= 0)
                            if and(and(i+u < dim_boxes1d_arg, j+m < dim_boxes1d_arg), k+n < dim_boxes1d_arg)
                                if not(and(and(u==0, m==0), n==0))
                                    box_cpu_nei_x(nh, box_cpu_nn(nh)) = k + n;
                                    box_cpu_nei_y(nh, box_cpu_nn(nh)) = j + m;
                                    box_cpu_nei_z(nh, box_cpu_nn(nh)) = i + u;
                                    box_cpu_nei_number(nh, box_cpu_nn(nh)) = ...
                                        box_cpu_nei_z(nh, box_cpu_nn(nh)) * dim_boxes1d_arg * dim_boxes1d_arg + ...
                                        box_cpu_nei_y(nh, box_cpu_nn(nh)) * dim_boxes1d_arg + ...
                                        box_cpu_nei_x(nh, box_cpu_nn(nh));
                                    box_cpu_nei_offset(nh, box_cpu_nn(nh)) = ...
                                        box_cpu_nei_number(nh, box_cpu_nn(nh)) * NUMBER_PAR_PER_BOX;
                                    box_cpu_nn(nh) = box_cpu_nn(nh) + 1;
                                end
                            end
                        end
                    end
                end
            end
            nh = nh + 1;
            
        end
    end
end

% input v,x,y,z
% value_range = [1 10];
% rv_cpu = randi(value_range, 4, dim_space_elem) / 10;
rv_cpu = createMatrixRandi2(10, 4, dim_space_elem) / 10;
% input (charge)
% qv_cpu = randi(value_range, 1, dim_space_elem) / 10;
qv_cpu = createMatrixRandi2(10, 1, dim_space_elem) / 10;
% output v,x,y,z
fv_cpu = zeros(4, dim_space_elem); 

disp('starting computation');
tic();
fv_cpu = kernel_cpu(par_cpu_alpha, dim_number_boxes,...
    box_cpu_offset, box_cpu_nn, box_cpu_nei_number, ...
    rv_cpu, qv_cpu, fv_cpu, NUMBER_PAR_PER_BOX);
elapsedTime = toc();

sum_cpu = round(sum(fv_cpu,2));

fprintf(1, '{');
fprintf(1, '"status": 1,');
fprintf(1, '"options": "-boxes1d 6",');
fprintf(1, '"time": %f,', elapsedTime);
fprintf(1, '"output": [ %d, %d, %d, %d ]', sum_cpu(1), sum_cpu(2), sum_cpu(3), sum_cpu(4));
fprintf(1, '}\n');

end
