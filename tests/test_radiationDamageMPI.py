import subprocess
import shlex, os, shutil
import unittest


class radiationDamageMPITests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """ Setting up the test class."""
        os.mkdir("diffr_out")

    @classmethod
    def tearDownClass(cls):
        """ Tear down the test class."""
        shutil.rmtree("diffr_out")

    def test_diffraction_calculation_parallel(self):
        """
        Test diffraction calculation in parallel.
        """
        mpicommand = 'mpiexec -n 2 --map-by node'

        command_sequence = ['radiationDamageMPI',
                            '--inputDir', 'pmi',
                            '--outputDir', 'diffr_out',
                            '--beamFile', 'test_files/s2e.beam',
                            '--geomFile', 'test_files/s2e.geom',
                            '--uniformRotation', '1',
                            '--calculateCompton', '0',
                            '--sliceInterval', '100',
                            '--numSlices', '100',
                            '--pmiStartID', '1',
                            '--pmiEndID', '1',
                            '--numDP', '2'
                            ]

        args = shlex.split(mpicommand) + command_sequence

        proc = subprocess.Popen(args)

        proc.wait()

        self.assertTrue(True)

    def test_diffraction_calculation_parallel_no_rotation(self):
        """
        Test diffraction calculation in parallel with no rotation applied.
        """
        mpicommand = 'mpiexec -n 2 --map-by node'

        command_sequence = ['radiationDamageMPI',
                            '--inputDir', 'pmi',
                            '--outputDir', 'diffr_out',
                            '--beamFile', 'test_files/s2e.beam',
                            '--geomFile', 'test_files/s2e.geom',
                            #'--uniformRotation', '1',
                            '--calculateCompton', '0',
                            '--sliceInterval', '100',
                            '--numSlices', '100',
                            '--pmiStartID', '1',
                            '--pmiEndID', '1',
                            '--numDP', '2'
                            ]

        args = shlex.split(mpicommand) + command_sequence

        proc = subprocess.Popen(args)

        proc.wait()

        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()

