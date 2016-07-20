#ifndef SOURCE_FILEHELPER_H_
#define SOURCE_FILEHELPER_H_


#include <iostream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <memory>

#include "readInfo.h"

#define F_OK 0
class FileHelper {
    
    private:
        bool tab;
        bool opt_reads;
        bool interleaved;
        bool fastq;
        FILE **files;
        int fileCount;
        int currentFile;
    public:
        void ParseFileNames(char *fName);
        void FileCheck(char *fName);
        void setFastq(bool b) {fastq = b;}
        void setTab(bool b) {tab = b;}
        void setInterleaved(bool b) {interleaved = b;}
        void Strip(char *c[]);
        void Closer();

        /*fileHelper for stdin*/
        FileHelper() {
            tab = true;
            interleaved = false;
            fastq = false;
            opt_reads = false;
            files = (FILE **)malloc(1*sizeof(FILE *));
            fileCount = 1;
            currentFile = 0;
            files[0] = stdin;    
        }

        /*fileHelper for other*/
        FileHelper(char *fName) {
            tab = false;
            interleaved = false;
            fastq = false;
            opt_reads = false;
            currentFile = 0;
            fileCount = 0; 
            ParseFileNames(fName);
            
        }
        void readData(std::shared_ptr<readInfo> &r); 
        void readData(std::shared_ptr<readInfo> &r1, std::shared_ptr<readInfo> &r2); 

};

#endif 
