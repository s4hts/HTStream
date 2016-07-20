#ifndef SOURCE_ARGCOLLECTOR_H_
#define SOURCE_ARGCOLLECTOR_H_
#include <stdlib.h>
#include <getopt.h>
#include <stdio.h>
#include <unistd.h>
#include <vector>

#include "fileHelper.h"
#include "fileWriter.h"

#define F_OK 0

class argCollector {


    public:        
        FileHelper *R1_In;
        FileHelper *R2_In;
        FileHelper *SE_In;
        FileHelper *TAB_In;
        FileHelper *INTER_In;
        FileHelper *STDIN_In;

        /*R1 will either be R1 fastq, Interleaved of Tab*/
        FileWriter *R1_Out;
        /*R2 will always be R2 fastq*/
        FileWriter *R2_Out;
        /*Se will always be fastq SE */
        FileWriter *SE_Out;

        int bpWindowSize;
        int startLoc;

        FILE *log;
    
        argCollector(int argc, char *argv[]);

};

#endif
