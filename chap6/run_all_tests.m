% Runs all tests with testingScript.m and testingScript_times.m
% runs all tests to generate plots and table with timings


for testcase = 1:12
	% Run plot testcases
	fprintf("Running test case %d\n\n", testcase);
	testingScript(testcase);

	fprintf("\n-----------------------\n\n");

end

for testcase = 7:12
	% Run timing testcases
	fprintf("Running timing test case %d\n\n", testcase);
	testingScript_times(testcase);

	fprintf("\n-----------------------\n\n");

end

