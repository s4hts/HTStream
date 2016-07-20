#include "fileWriter.h"
#include <stdexcept>

/*Gzip output writer*/


/*make sure no one puts multiple types in*/
/*The next four functions might be a little bloated, however, it is to ensure users know when they are messing up*/
void FileWriter::setTab(bool b) {
    if (!b) {
        tab = b;
    } else if (!fastq && !interleaved) {
        tab = b;
    } else {
        fprintf(stderr, "Within fileWriter.cpp setTab() - stdout, fastq, or interleaved all ready set\n");
        fprintf(stderr, "Please, only set one type of output, either fastq, interleaved, tab, or stdout (stdout is always tab)\n");
        fprintf(stderr, "setTab: fastq=%i interleaved=%i\n", fastq, interleaved);
        exit(24);
    }

}

void FileWriter::setFastq(bool b) {
    if (!b) {
        fastq = b;
    } else if (!fastq && !interleaved && !to_stdout) {
        fastq = b;
    } else {
        fprintf(stderr, "Within fileWriter.cpp setFastq() - tab, stdout, or interleaved all ready set\n");
        fprintf(stderr, "Please, only set one type of output, either fastq, interleaved, tab, or stdout (stdout is always tab)\n");
        exit(25);
    }

}

void FileWriter::setInterleaved(bool b) {
    if (!b) {
        interleaved = b;
    } else if (!fastq && !tab && !to_stdout) {
        interleaved = b;
    } else {
        fprintf(stderr, "Within fileWriter.cpp setInterleaved() - tab, stdout, or fastq all ready set\n");
        fprintf(stderr, "Please, only set one type of output, either fastq, interleaved, tab, or stdout (stdout is always tab)\n");
        exit(26);
    }

}

void FileWriter::setStdout (bool b) {
    if (!b) {
        to_stdout = b;
    } else if (!fastq && !interleaved) {
        to_stdout = b;
    } else {
        fprintf(stderr, "Within fileWriter.cpp setStdout() - tab, stdout, or fastq all ready set\n");
        fprintf(stderr, "Please, only set one type of output, either fastq, interleaved, tab, or stdout (stdout is always tab)\n");
    }

}



FileWriter::FileWriter(bool force) {
    setForce(force);

    fOut = NULL;
    tab = false;
    fastq = false;
    interleaved = false;
    to_stdout = false;
}


void FileWriter::OpenFile(char *fName_) {

    if (to_stdout) {
        fOut = stdout;
        setFastq(false);
        setTab(true);
     } else if (fName_ != NULL) {
        /*Checks if write permissions are held for this file*/
        if (access(fName_, F_OK) != -1) {
            if (force) {
                fName = strdup(fName_);
            } else {
                fprintf(stderr, "File %s exists - please use -F option to force overwrite, -P to change the name of the file, or manually remove or move the file\n", fName_);
                exit(10);
            }
        } else {
            fName = strdup(fName_);
        }
    } else {
        fprintf(stderr, "fName within fileWriter.cpp fName in OpenFile() is NULL\n");
        exit(6);
    }


}

void FileWriter::Closer() {
    if (fOut) {
        fclose(fOut);
    }
}

void FileWriter::writeData(readInfo *R1, readInfo *R2, readInfo *R3) {
    if ((R1 && R2 && R3)) {
        fprintf(stderr, "Error in fileWriter.cpp in function writeData\n");
        fprintf(stderr, "No format takes all three reads to write\n");
        exit(7);
    }
    if ((!R1 && !R2 && !R3)) {
        fprintf(stderr, "Error in fileWriter.cpp in function writeData\n");
        fprintf(stderr, "No format takes no reads to write\n");
        exit(7);
    }

    /*this ensures opening only if need be*/
    if (fOut == NULL) {
        fOut = fopen(fName, "w");
        /*Ensures that the file actually opens*/
        if (fOut == NULL) {
            fprintf(stderr, "Within fileWriter.cpp in OpenFile() (force) fName does not have write permissions\n");
            fprintf(stderr, "File not accessable %s\n", fName);
            exit(9);
        }

    }

    if (fastq) {
        if ((R1 && R2) || (R2 && R3) || (R1 && R3)) {
            fprintf(stderr, "Error in fileWriter.cpp in function writeData\n");
            fprintf(stderr, "Fastq format only takes one read at a time to write\n");
        }
        /*Temparary for R1, R2, or R3*/
        readInfo *tmp ;

        /*Eventually we could have different header data*/
        if (R1) {
            tmp = R1;
        } else if (R2) {
            tmp = R2;
        } else if (R3) {
            tmp = R3;
        } else {
            throw std::runtime_error("all 3 reads were null");
        }

        fprintf(fOut, "%s\n%s\n+\n%s\n", tmp->getHeader(), tmp->getSeq(), tmp->getQual());
    } else if (interleaved) {
        if (!(R1 && R2) || !R3) {
            fprintf(stderr, "Error in fileWriter.cpp in function writeData\n");
            fprintf(stderr, "Interleaved format only takes two read at a time to write\n");
        }

        fprintf(fOut, "%s\n%s\n+\n%s\n%s\n%s\n+\n%s\n", R1->getHeader(), R1->getSeq(), R1->getQual(), R2->getHeader(), R2->getSeq(), R2->getQual());
    } else if (tab) {

        if (R1 && R2) {
            fprintf(fOut, "%s\t%s\t%s\t%s\t%s\n", R1->getHeader(), R1->getSeq(), R1->getQual(), R2->getSeq(), R2->getSeq());
        } else if (R1) {
            fprintf(fOut, "%s\t%s\t%s\n", R1->getHeader(), R1->getSeq(), R1->getQual());
        }
    }







}
