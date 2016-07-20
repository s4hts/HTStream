import unittest
import os
import fnmatch
import subprocess
import filecmp


def find_fastq_files(directory, pattern):
    """Walks the directory structure, appending filenames to an array"""
    filenames = []
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                filenames.append(filename)
    filenames.sort()
    return filenames


def sub_process(command):
    """Run a shell command """
    return subprocess.check_output(
        command, stderr=subprocess.STDOUT, shell=True)


def file_compare(command, expected, returned):
    """Use this function to compare two files"""
    subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
    return filecmp.cmp(expected, returned)


def parse_fastq(filename):
    """Use this function to create a dicionary of the fastq file"""
    with open(filename) as f:
        lines = f.readlines()
    head = [item[:-1] for item in lines[0::4]]
    read = [item[:-1] for item in lines[1::4]]
    qual = [item[:-1] for item in lines[3::4]]
    return {'Headers': head, 'Sequences': read, 'QScores': qual}


class TestCase(unittest.TestCase):

    def setUp(self):
        """Setup any thing used in muliple tests"""
        myR1file = " -1 fastqFiles/testCase_1X_R1.fastq "
        myR2file = " -2 fastqFiles/testCase_1X_R2.fastq "
        additFlags = "-N -F"
        myShellCmd = "../super_deduper"
        self.myCommand = myShellCmd + myR1file + myR2file + additFlags

    def test_find_fastq_files_recursively(self):
        """Should return all fastq files from the sub directories"""
        self.assertEqual(find_fastq_files('fastqFiles', '*.fastq'),
                         ['fastqFiles/testCase_1X_R1.fastq',
                          'fastqFiles/testCase_1X_R2.fastq'],
                         "Unable to find test files")

    def test_tab_output(self):
        """Tests for 'Reads_Read' in the output"""
        self.assertIn("Reads_Read\tReads_Written\tReads_Discarded\t"
                      "Singletons\tDoubles\tThree_Plus\t"
                      "Disqualifed_Reads\tReplacements_Called\t"
                      "Reads_Per_Second\tTotal_Time",
                      sub_process(self.myCommand),
                      "Unexpected tab structure "
                      "super_deduper -1 fastqFiles/testCase_1X_R1.fastq"
                      "-2 fastqFiles/testCase_1X_R2.fastq"
                      "-F -N")

    def test_item_from_one_exists_in_two(self):
        """Tests if the first entry in the expected output is in the input"""
        sub_process(self.myCommand)
        data01 = parse_fastq("no_dup_R1.fastq")
        self.assertIn("@higher_qscore_duplicate_of_read_one",
                      data01['Headers'],
                      "The output does not contain the expected sequence")

    def test_for_reverse_complement(self):
        """Tests that the reverse complement was removed"""
        sub_process(self.myCommand)
        data01 = parse_fastq("no_dup_R1.fastq")
        self.assertNotIn("@reverse_complement_of_read_one_R1",
                         data01['Headers'],
                         "The reverse complement was not removed")


if __name__ == '__main__':
    unittest.main()
