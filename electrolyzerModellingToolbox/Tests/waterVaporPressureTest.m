classdef waterVaporPressureTest < matlab.unittest.TestCase

    properties
    end

    methods (Test)

        function testWaterVaporPressure1(testCase)
            actSolution = waterVaporPressure(300);
            expSolution = 0.035358594230767;
            testCase.verifyEqual(actSolution, expSolution, 'AbsTol', 1e-10)
        end

        function testWaterVaporPressure2(testCase)
            actSolution = waterVaporPressure(300,'model',1);
            expSolution = 0.035947102222906;
            testCase.verifyEqual(actSolution, expSolution, 'AbsTol', 1e-10)
        end

        function testWaterVaporPressure3(testCase)
            actSolution = waterVaporPressure(300,'model',2);
            expSolution = 0.035358594230767;
            testCase.verifyEqual(actSolution, expSolution, 'AbsTol', 1e-10)
        end

    end

end