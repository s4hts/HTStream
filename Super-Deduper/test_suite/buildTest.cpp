#include "fileHelper.h"
#include "argCollector.h"
#include "fileWriter.h"
#include "binarySearch.h"

int main(int argc, char *argv[]) {
    char c[100];
    char out[100];
    sprintf(c, "R1.fastq");
    sprintf(out, "out.txt");
    FileHelper R1 = FileHelper(c);    
    R1.setFastq(true); 
    readInfo *r1 = R1.readData();
    printf("%s\n", r1->getSeq());
    printf("%s\n", r1->getQual());
    printf("%s\n", r1->getHeader());
    argCollector AC = argCollector(argc, argv);
    FileWriter *FW = new FileWriter(true);

    FW->setForce(true);
    FW->OpenFile(out);
    BinarySearchTree *bst = new BinarySearchTree();

}
