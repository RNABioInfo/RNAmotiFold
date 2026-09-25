import unittest

from RNAmotiFold.bgap_rna.bgap_rna import bgap_rna


class TestBgapObject(unittest.TestCase):
    def setUp(self):
        #Arrange an empty bgap_rna obj for all tests on it
        self.obj = bgap_rna()


    def test_algorithm_binary_rnamotifold_default(self):
        #Arrange
        #Do nothing to the obj
        #Act
        binary = self.obj.algorithm_binary
        #Assert
        self.assertEqual(binary,"RNAmotiFold")

    def test_algorithm_binary_rnamotifold_pfc(self):
        #Arrange
        self.obj.pfc=True
        #Act
        binary =  self.obj.algorithm_binary
        #Assert
        self.assertEqual(binary,"RNAmotiFold_pfc")

    def test_algorithm_binary_rnamotifold_subopt(self):
        #Arrange
        self.obj.subopt = True
        #Act
        binary =  self.obj.algorithm_binary
        #Assert
        self.assertEqual(binary,"RNAmotiFold_subopt")

    def test_algorithm_binary_rnamotialign_default(self):
        #Arrange
        self.obj.algorithm = "RNAmotiAlign"
        #Act
        binary = self.obj.algorithm_binary
        #Assert
        self.assertEqual(binary,"RNAmotiAlign")

    #@unittest.expectedFailure()
    #def test_algorithm_binary_rnamotialign_pfc(self):
    #    #Arrange
    #    self.obj.algorithm = "RNAmotiAlign"
    #    #Act and Assert
    #    with self.assertRaises(NotImplementedError):
    #        self.obj.pfc = True




if __name__ == "__main__":
    unittest.main()
