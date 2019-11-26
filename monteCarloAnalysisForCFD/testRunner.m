function results = testRunner()
    tests       = functiontests({testMonteCarloFramework});
    results     = runtests(tests);
end