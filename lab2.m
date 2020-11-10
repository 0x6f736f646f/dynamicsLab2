function lab2
    % Values to theta2 and theta4 are  provided in the question
    %------------part a-------------
    theta2  = [40; 45; 50; 55; 60]; 
    theta4  = [70; 76; 83; 91; 100];
    link_ratios = get_link_ratios(theta2, theta4);
    disp(link_ratios);
    [a,b,c,d] = get_link_lengths(link_ratios);
    sprintf("Crank: %.4f \n Coupler: %.4f \n Follower: %.4f \n  Fixed: %.4f\n",a,b,c,d);
    %----------- part b --------------
    transmission_angles  = get_transmission_angles(a,b,c,d,40,60,1);
    input_angles = 40:1:60;
    figure;
    plot(input_angles, transmission_angles);
    xlabel("Input angles");
    ylabel("Transmission angles");
   % commenting  on the quality of the transmission angles
            if (all(transmission_angles >= 40)) || (all(transmission_angles <= 140))
                sprintf("All the transmission angles guarantee a smooth rotation");
            else
                sprintf("Some transmission angles do not guarantee a smooth rotation")
            end
     %-----------part c ---------------
     
end

function link_ratios  = get_link_ratios(theta2, theta4)
   % using the freudensteins method, this method computes the link ratios
   
   % To form a complete and functional matrix the two 1D arrays  should be
   % of the same length. If not, this raises a run time error.
   if(length(theta4) ~= length(theta2))
       disp("Matrices' lengths not equal");
       quit(1);
   end
   A = [];
   b = [];
   % rows in both the matrices are added in a loop with columns in each
   for i = 1:length(theta2)
       temp1 = [(cosd(theta4(i))) (-1 * (cosd(theta2(i)))) (1)];
       temp2 = cosd(theta2(i)-theta4(i));
       A = [A; temp1];
       b = [b; temp2];
   end
   link_ratios = lsqr(A,b);
end
function [a,b,c,d] = get_link_lengths(link_ratios)
% get_link_lengths uses the link ratios and the fixed link to find the
% lengths of the other links.
% Lengths can be negatives therefore absolutes of the calculated values are
% sorted.
    d = 180;
    a = abs(d/link_ratios(1));
    c = abs(d/link_ratios(2));
    b = abs(sqrt(a^2  + c^2 + d^2 -(link_ratios(3) * 2 * a * c)));
end
function transmission_angles = get_transmission_angles(a,b,c,d,lower_limit, upper_limit, steps)
% Transmission anngles are calculated using the obtained link lengths and
% and the respective input angles 
    transmission_angles = zeros(1, ((upper_limit - lower_limit)/steps));
    j = 1;
    for i = lower_limit:steps:upper_limit
        m = acosd(((b^2 + c^2) - (a^2 + d^2) + (2 * a * d * cosd(i))) / (2 * b * c));
        transmission_angles(j) = m;
        j = j + 1;
    end
end
function structuralErrors = get_structural_errors(theta2, theta4, link_ratios)
% Structural error is basically the difference between the left side of the
% freudeinsten's equation and the right side.
    structuralErrors = zeros(1,length(theta4));% theta4 can also be used since they are of the same length
    for i = 1:length(theta2)
        er1 = link_ratios(1)*cosd(theta4(i)) - link_ratios(2)*cosd(theta2(i)) + link_ratios(3) - cosd(theta2(i) - theta4(i));
        structuralErrors(i) = er1;
    end
end