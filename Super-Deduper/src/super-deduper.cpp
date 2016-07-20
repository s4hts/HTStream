#include "argCollector.h"
#include "binarySearch.h"
#include "readInfo.h"
#include "fileHelper.h"
#include <memory>

int main(int argc, char *argv[]) {

    argCollector args(argc, argv);

	BinarySearchTree bst(args.bpWindowSize, args.startLoc);

	/*Error checking to make sure R1 and R2 exist in arrCollector*/
	if (args.R1_In && args.R2_In) {
		FileHelper *R1 = args.R1_In;
		FileHelper *R2 = args.R2_In;
		std::shared_ptr<readInfo> r1, r2;

		/*Yeah, I hate this too, I did try and do while((r1 = R1-readData) != NULL && (r2 = R2->readData....
		 * however, it would check the R1->getData, receive NULL, and then exit if statment
		 * this means r2 wouldn't also be set to null to check the size of the file R2
		 * that is why I stuck with always true with the break*/
		while(1) {
			R1->readData(r1);
			R2->readData(r2);
			if (!r1 || !r2) {
				break;
			}
			/*The only reason I'm doing this is for error checking
			 * I need to make sure the files are the same size*/
			bst.AddNode(r1, r2);
		}
		/*Means the files are different lengths*/
		if (!r1 && r2) {
			fprintf(stderr, "File R1 is shorter than File R2\n");
			exit(16);
		} else if (r1 && !r2) {
			fprintf(stderr, "File R2 is shorter than File R1\n");
			exit(17);
		}
		(args.R1_In)->Closer();
		(args.R2_In)->Closer();

	}

	if (args.SE_In) {
		FileHelper *SE = args.SE_In;
		std::shared_ptr<readInfo> se = 0;

		while(1) {
			SE->readData(se);
			if (se) {
				/*Single end reads R1 is !null R2 is null*/
				bst.AddNode(se, std::shared_ptr<readInfo>(nullptr));
			} else {
				break;
			}
		}
		(args.SE_In)->Closer();

	}

	if (args.INTER_In) {
		FileHelper *INTER = args.INTER_In;
        std::shared_ptr<readInfo> r1, r2;

		while (1) {
			INTER->readData(r1, r2);
			if (r1 && r2) {
				bst.AddNode(r1, r2);
			} else {
				break;
			}
		}
		(args.INTER_In)->Closer();

	}

	if (args.TAB_In) {
		FileHelper *TAB = args.TAB_In;
		std::shared_ptr<readInfo> r1, r2;
		while (1) {
			TAB->readData(r1, r2);
			if (r1) {
				bst.AddNode(r1, r2);
			} else {
				break;
			}
		}
		(args.TAB_In)->Closer();

	}

	if (args.STDIN_In) {
		FileHelper *STDIN = args.STDIN_In;
		std::shared_ptr<readInfo> r1, r2;
		while (1) {
			STDIN->readData(r1, r2);
			if (r1) {
				bst.AddNode(r1, r2);
			} else {
				break;
			}
		}
		(args.STDIN_In)->Closer();


	}

	bst.PrintAndDelete(args.R1_Out, args.R2_Out, args.SE_Out);

	if (args.R1_Out) {
		(args.R1_Out)->Closer();
	}
	if (args.R2_Out) {
		(args.R2_Out)->Closer();
	}
	if (args.SE_Out) {
		(args.SE_Out)->Closer();
	}

	bst.endTime();
	bst.outputStats(args.log);
}
