import unittest
import shutil
from pysingfel.radiationDamage import *

class ReferenceTest(unittest.TestCase):
    """ Test class for the reference test. """


    def testDiffractionCalculation(self):
        """
        Test diffraction calculation using just one process.
        """
        diffraction_parameters = {'uniform_rotation': None,
                                  'calculate_Compton': False,
                                  'slice_interval': 100,
                                  'number_of_slices': 5,
                                  'pmi_start_ID': 1,
                                  'pmi_stop_ID': 1,
                                  'number_of_diffraction_patterns': 1,
                                  'beam_parameter_file': 'test_files/s2e.beam',
                                  'beam_geometry_file': 'test_files/s2e.geom',
                                  'number_of_MPI_processes': 1}

        uniform_rotation = diffraction_parameters['uniform_rotation']
        calculate_Compton = diffraction_parameters['calculate_Compton']
        slice_interval = diffraction_parameters['slice_interval']
        number_of_slices = diffraction_parameters['number_of_slices']
        pmi_start_ID = diffraction_parameters['pmi_start_ID']
        pmi_stop_ID = diffraction_parameters['pmi_stop_ID']
        number_of_diffraction_patterns = diffraction_parameters['number_of_diffraction_patterns']
        beam_parameter_file = diffraction_parameters['beam_parameter_file']
        beam_geometry_file = diffraction_parameters['beam_geometry_file']

        parameters = {
                      'inputDir': 'test_files/reference_input_data/Fe2O3_supercell7x7x7/pmi',
                      'outputDir': 'diffr_out',
                      'beamFile': beam_parameter_file,
                      'geomFile': beam_geometry_file,
                      'uniformRotation': uniform_rotation,
                      'calculateCompton': calculate_Compton,
                      'sliceInterval': slice_interval,
                      'numSlices': number_of_slices,
                      'pmiStartID': pmi_start_ID,
                      'pmiEndID': pmi_stop_ID,
                      'numDP': number_of_diffraction_patterns,
                      'rotationAxis': 'xyz'}

        diffract(parameters)

        # Read results and reference data.
        reference_data = h5py.File("test_files/diffr_out_Fe2O3_reference.h5","r")['data/0000001/diffr'].value
        data = h5py.File("diffr_out/diffr_out_0000001.h5","r")['data/0000001/diffr'].value

        self.assertLess( ((data - reference_data)**2).max(), 1e-29 )

        shutil.rmtree( parameters['outputDir'] )


if __name__ == '__main__':
    unittest.main()
