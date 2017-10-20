import os, sys
import unittest

# Import classes to test.
from test_radiationDamageMPI import radiationDamageMPITests
from test_radiationDamage    import radiationDamageTests
from test_particle           import particleTests
from test_reference          import ReferenceTest

# Setup the suite.
def suite():
    suites = [
            unittest.makeSuite(radiationDamageMPITests, 'test'),
            unittest.makeSuite(radiationDamageTests,    'test'),
            unittest.makeSuite(particleTests,           'test'),
            unittest.makeSuite(ReferenceTest,           'test'),
             ]

    return unittest.TestSuite(suites)

# If called as script, run the suite.
if __name__=="__main__":

    result = unittest.TextTestRunner(verbosity=2).run(suite())

    # If tests were successful, return with exit code 0.
    if result.wasSuccessful():
        sys.exit(0)

    # ... otherwise, exit code 1.
    sys.exit(1)
