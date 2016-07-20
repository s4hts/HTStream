#include "fileHelper.h"
/*Parses files names, and uses FileCheck to get gzipped and if they
 * are Okay to read*/

void FileHelper::ParseFileNames(char *fName) {
    /*Allocas zero memory so we can realloc later within CheckFiles*/

    files = (FILE **)malloc(0);
    if (fName != NULL) {
        /*Breaks up the commas and opens each of those files
         * and checks if there are OK to read and gzipped*/
        char *tok;
        while ((tok = strtok_r(fName, ",", &fName))) {
            FileCheck(tok);
        } 
    } else {
        fprintf(stderr, "fName within fileHelper.cpp is NULL\n");
        exit(1);
    }
}

/*Checks if file permissions and if the file is gzipped*/
void FileHelper::FileCheck(char *fName) {
        
    if (access(fName, F_OK)!= -1) {
        /*Increments file count to help with loops and memory alloc*/
        fileCount++;
        files = (FILE **)realloc(files, sizeof(FILE *) * fileCount);
        /*Must decrement by one to be in proper spot in array*/
        files[fileCount-1] = fopen(fName, "r");
        /*Checks the first two bytes within file to check it if is gzipped*/
        /*byte 1 == 0x1f and byte 2 == 0x8b*/
        unsigned char bytes[2];

        /*A decent check to see if the file is empty or not*/
        if (!fread(bytes, sizeof(unsigned char), 2, files[fileCount-1])) {
            fprintf(stderr, "Within fileHelper.cpp FileCheck() error in fread\n");
            fprintf(stderr, "It seems that the file doesn't even have two bytes in it\n");
            exit(21);
        }
            
        /*Seeks to the start of the file*/
        rewind(files[fileCount-1]);
        if (bytes[0] == 0x1f && bytes[1] == 0x8b) {
            /*So, we want to have a thread un-gzipping and piping into FILE*
             * So we are closing this file and putting a popen gzipped file*/
            
            fclose(files[fileCount-1]);
            const char *ungzip = "gunzip -c ";
            char *popenName = (char *)malloc(sizeof(char) * (strlen(ungzip) + strlen(fName) + 1));
            sprintf(popenName, "%s%s", ungzip, fName);
            files[fileCount-1] = popen(popenName, "r");
            free(popenName);
        } else {
            /*Nothing needs to happen, it is a plain text file*/
        }
    } else {
        fprintf(stderr, "Within fileHelper.cpp FileCheck() function\n");
        fprintf(stderr, "File is not accessable. File Name: '%s;\n", fName);
        exit(2);
    }
}


void FileHelper::readData(std::shared_ptr<readInfo> &r1, std::shared_ptr<readInfo> &r2) {

    if (fastq) {
        fprintf(stderr, "Error withink fileHelper.cpp in readData(readInfo**, readInfo**)\n");
        fprintf(stderr, "Fastq should never use the constructor readData(readInfo**, readInfo**)\n");
        exit(19);
    } else if (tab) {
        char c[4096];
        if (fgets(c, 4096, files[currentFile]) == NULL) {
            if (currentFile == fileCount -1) {
                return;
            } else {
                return;
            }
        }
        char *tok, *tmpPointer, data[5][4096];
        int count = 0;
        tmpPointer = &c[0];
        while ((tok = strtok_r(tmpPointer, "\t", &tmpPointer))) {
            strcpy(data[count], tok);
            count++;
        } 
        if (count == 5) {
            data[4][strlen(data[4])-1] = '\0';
            r1 =  std::make_shared<readInfo>(data[0], data[1], data[2], opt_reads );
            r2 =  std::make_shared<readInfo>(data[0], data[3], data[4], opt_reads );
        } else if (count == 3) {
            data[2][strlen(data[2])-1] = '\0';
            r1 =  std::make_shared<readInfo>(data[0], data[1], data[2], opt_reads );
            r2 =  NULL;
        } else {
            fprintf(stderr, "Error in fileHelper.cpp readData(readInfo**, readInfo**)\n");
            fprintf(stderr, "Tab delimited format is broken. Make sure tabs are tabs and not space, and make sure there is either 3 or 5 enteries per line\n");
            fprintf(stderr, "Could also be an extra line at the end of your fastq file as well\n");
            exit(23);
        }
    } else if (interleaved) {
        /*Since it is interleaved, there are 8 lines, max of 4096*/
        char c[8][4096];

        for (int f = 0; f < 8; f++) {
            if (fgets(c[f], 4096, files[currentFile]) != NULL) {
                /*Removes \n character from string*/
                c[f][strlen(c[f])-1]='\0';
            } else {
                /*If it fails at f == 0 then the value is actually NULL*/
                if (f == 0) {
                    if (currentFile == fileCount -1) {
                        r1 = NULL;
                        r2 = NULL;
                        return;
                    } else {
                        currentFile++;
                        return readData(r1, r2);
                    }
                } else if (f == 4) {
                    r1 = std::make_shared<readInfo>(c[0], c[1], c[3], opt_reads);
                    r2 = NULL;
                } else {
                    fprintf(stderr, "Within fileHelper.cpp readData(readInfo **, readInfo**) function*/");
                    fprintf(stderr, "Something wrong with fgets(), file might be uneven\n");
                    exit(22);
                }

            }
        }

        r1 =  std::make_shared<readInfo>(c[0], c[1], c[3], opt_reads );
        r2 =  std::make_shared<readInfo>(c[4], c[5], c[7], opt_reads );

    }

}


void FileHelper::Closer() {
    for (int i = 0; i < fileCount; i++) {
        if (files[i]) {
            fclose(files[i]);
        }
    }
}
/*Starts reading files*/
void FileHelper::readData(std::shared_ptr<readInfo> &r) {
    /*Loops through each FILE **/
    /*fgets(char *, size, file) == NULL bad*/
    if (tab || interleaved ) {
        fprintf(stderr, "Error withink fileHelper.cpp in readData(readInfo**, readInfo**)\n");
        fprintf(stderr, "Tab or Interleaved should never use the constructor readData(readInfo**, readInfo**)\n");
        exit(20); 
    } else if (fastq) {
        /*This memory will be free once it leaves local scope*/
        /*4 reads, max 4096 characters*/
        char c[4][4096]; 
            
        for (int f = 0; f < 4; f++) {
            if (fgets(c[f], 4096, files[currentFile]) != NULL) {
                /*Removes the newline*/
                c[f][strlen(c[f])-1] = '\0';
                /*good!*/
            } else {
                if (f == 0) {
                    /*Increments current file if there are more reads to be had
                     * the file pointer will keep track of current location so
                     * this function sort of acts like a generator.*/
                    if (currentFile == fileCount - 1) {
                        r =  NULL;
                        return;
                    } else {
                        currentFile++;
                        return readData(r);
                    }
                }

                fprintf(stderr, "Within fileHelper.cpp readData() function\n");
                fprintf(stderr, "Something went wrong with fgets (make sure there is no extra lines)\n");
                exit(3);
            }
        }
        r =  std::make_shared<readInfo>(c[0], c[1], c[3], opt_reads );
        
    }

}







