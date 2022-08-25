classdef funcObjectTest < matlab.unittest.TestCase

    properties
        Parameters1
        Variables1
        Constants1
        Workspace1
        funcHandle1
        func1

        Parameters2
        Variables2
        Workspace2
        funcHandle2
        func2
    end

    methods(TestMethodSetup)
        function initializeFuncs(testCase)
            testCase.Parameters1 = struct('a',3,'b',12);
            testCase.Variables1 = struct('I',[1:5]','V',[linspace(0,2,5)]');
            testCase.Constants1 = struct('twoPi',2*pi);
            testCase.Workspace1 = struct('Parameters',testCase.Parameters1,'Variables',testCase.Variables1,'Constants',testCase.Constants1);
            testCase.funcHandle1 = @(Workspace) Workspace.Constants.twoPi*Workspace.Parameters.a^2 + Workspace.Variables.I.*Workspace.Variables.V*Workspace.Parameters.b;
            testCase.func1 = func(testCase.funcHandle1,testCase.Workspace1);


            testCase.Parameters2 = struct('a',3,'b',12,'c',60);
            testCase.Variables2 = struct('T',285,'p',20);
            testCase.Workspace2 = struct('Parameters',testCase.Parameters2,'Variables',testCase.Variables2);
            testCase.funcHandle2 = @(Workspace) Workspace.Parameters.a^(1/2) + Workspace.Variables.T./Workspace.Variables.p*Workspace.Parameters.c;
            testCase.func2 = func.createEmpty;
            testCase.func2.setFuncHandle(testCase.funcHandle2);
        end
    end

    methods (Test)

        function checkFunc1(testCase)
            result = calculate(testCase.func1);
            expResult = [56.548667764616276;68.548667764616280;
                92.548667764616280;1.285486677646163e+02;
                1.765486677646163e+02];
            testCase.verifyEqual(result, expResult, 'AbsTol', 1e-10)
        end

        function checkFunc2(testCase)
            testCase.func2.setInWorkspace(testCase.Workspace2);
            result = calculate(testCase.func2);
            expResult = 8.567320508075688e+02;
            testCase.verifyEqual(result, expResult, 'AbsTol', 1e-10)
        end

        function setFuncWorkspace(testCase)
            % Make sure that the workspace is empty before starting
            bempty1 = isempty(fieldnames(testCase.func2.Workspace.Constants));
            testCase.verifyEqual(bempty1, true)
            bempty2 = isempty(fieldnames(testCase.func2.Workspace.Variables));
            testCase.verifyEqual(bempty2, true)
            bempty3 = isempty(fieldnames(testCase.func2.Workspace.Parameters));
            testCase.verifyEqual(bempty3, true)
            bempty4 = isempty(fieldnames(testCase.func2.Workspace.Dependencies));
            testCase.verifyEqual(bempty4, true)

            % Set workspace for the func object and check values
            testCase.func2.setInWorkspace(testCase.Workspace2);
            testCase.verifyEqual(testCase.Parameters2.a, testCase.func2.Workspace.Parameters.a)
            testCase.verifyEqual(testCase.Parameters2.b, testCase.func2.Workspace.Parameters.b)
            testCase.verifyEqual(testCase.Parameters2.c, testCase.func2.Workspace.Parameters.c)
            testCase.verifyEqual(testCase.Variables2.T, testCase.func2.Workspace.Variables.T)
            testCase.verifyEqual(testCase.Variables2.p, testCase.func2.Workspace.Variables.p)
        end

        function removeFromFuncWorkspace(testCase)
            % Remove value from the func workspace and check new workspace
            testCase.func2.setInWorkspace(testCase.Workspace2);
            testCase.func2.removeFromWorkspace('b');
            bisfield = isfield(testCase.func2.Workspace.Parameters, "b");
            testCase.verifyEqual(bisfield, false);
        end

        function replaceInFuncWorkspace(testCase)
            % Replace value in the func workspace and check new workspace
            testCase.func2.setInWorkspace(testCase.Workspace2);
            % Check that value is what we expect it to be
            testCase.verifyEqual(testCase.Parameters2.a, testCase.func2.Workspace.Parameters.a)
            newValue = 15;
            testCase.func2.replaceInWorkspace('a',newValue);
            % Check that value has changed
            testCase.verifyEqual(testCase.func2.Workspace.Parameters.a, newValue);
        end

        function combineFuncs(testCase)
            % Combine two func objects
            testCase.func2.setInWorkspace(testCase.Workspace2);
            
            % Check that the equations are the same
            func3 = func.add(testCase.func1,testCase.func2);
            s2 = strcat(testCase.func1.equation,'+', testCase.func2.equation);
            isSame = all(func3.equation == s2);
            testCase.verifyEqual(isSame, true)

            % Check funcHandle
            func1handleStr = func2str(testCase.func1.funcHandle);
            func2handleStr = func2str(testCase.func2.funcHandle);
            pat = "@(" + lettersPattern + ")";
            tmpstr = extract(func2handleStr, pat);
            tmpl = length(tmpstr{1});
            shortenedFunc2handleStr = func2handleStr(tmpl+1:end);
            func3handleStr = func2str(func3.funcHandle);

            s2 = strcat(func1handleStr,'+', shortenedFunc2handleStr);
            isSame = all(func3handleStr == s2);
            testCase.verifyEqual(isSame, true)

            % Check workspace (assumes that values of first added func are
            % replaced by the values in second if they overlap)
            testCase.verifyEqual(testCase.func2.Workspace.Parameters.a, func3.Workspace.Parameters.a)
            testCase.verifyEqual(testCase.func2.Workspace.Parameters.b, func3.Workspace.Parameters.b)
            testCase.verifyEqual(testCase.func2.Workspace.Parameters.c, func3.Workspace.Parameters.c)
            testCase.verifyEqual(testCase.func1.Workspace.Variables.I, func3.Workspace.Variables.I)
            testCase.verifyEqual(testCase.func1.Workspace.Variables.V, func3.Workspace.Variables.V)
            testCase.verifyEqual(testCase.func2.Workspace.Variables.T, func3.Workspace.Variables.T)
            testCase.verifyEqual(testCase.func2.Workspace.Variables.p, func3.Workspace.Variables.p)
            testCase.verifyEqual(testCase.func1.Workspace.Constants.twoPi, func3.Workspace.Constants.twoPi)

            % Check fitlims (Currently just checks that they are empty)
            testCase.verifyEqual(isempty(fieldnames(testCase.func1.Fitlims)), true)
            testCase.verifyEqual(isempty(fieldnames(testCase.func2.Fitlims)), true)
            testCase.verifyEqual(isempty(fieldnames(func3.Fitlims)), true)
        end
    end
end