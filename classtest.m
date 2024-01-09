classdef classtest < handle
    %CLASSTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        number(1,1) double
    end
    
    methods
        function obj = classtest(number)
            %CLASSTEST Construct an instance of this class
            %   Detailed explanation goes here
            obj.number = number;
        end
        
        function [newnumber] = method1(obj,newnumber)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.number = newnumber;
        end
    end
end

