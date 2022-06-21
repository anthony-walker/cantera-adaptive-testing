import models
import unittest
import cantera as ct


class TestDatabaseFunctions(unittest.TestCase):

    def test_database():
        test_model = models.Hydrogen()
        test_model.append_to_db(Hydrogen.volume_problem, 1000, ct.one_atm, 1)


if __name__ == '__main__':
    unittest.main()
