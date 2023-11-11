import unittest
import files


class FileTestCases(unittest.TestCase):
    def test_load(self):
        data, loaded_files = files.load_files("patient_a")
        self.assertEqual(len(loaded_files), 7)


if __name__ == '__main__':
    unittest.main()
