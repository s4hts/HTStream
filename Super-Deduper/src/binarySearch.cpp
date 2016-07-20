#include "binarySearch.h"


void BinarySearchTree::outputStats(FILE *f) {
	fprintf(f, "Reads_Read\tReads_Written\tReads_Discarded\tSingletons\tDoubles\tThree_Plus\tDisqualifed_Reads\tReplacements_Called\tReads_Per_Second\tTotal_Time\n");
	fprintf(f, "%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%f\t%ld\n",
			reads_read,
			nodesCreated,
			dup_gone,
			singletons,
			doubles,
			threeplus,
			disReads,
			replaced,
			(double)reads_read/(double)((double)time_end-(double)time_start),
			time_end-time_start);
}


/*Print and Delete*/
void BinarySearchTree::PrintAndDelete(FileWriter *R1, FileWriter *R2, FileWriter *SE) {
	PrintAndDeletePrivate(root.get(), R1, R2, SE);
}

void BinarySearchTree::PrintAndDeletePrivate(Node *n, FileWriter *R1, FileWriter *R2, FileWriter *SE) {
	if (n != NULL) {
		if (n->count == 1) {
			singletons++;
		} else if (n->count == 2) {
			doubles++;
		} else {
			threeplus++;
		}

		PrintAndDeletePrivate(n->left.get(), R1, R2, SE);
		PrintAndDeletePrivate(n->right.get(), R1, R2, SE);

		/*Only if Fastq file exists*/
		if (R2) {
			/*R1 and R2 printed Normal Fastq format*/
			if (n->R2) {
				R1->writeData(n->R1.get(), NULL, NULL);
				R2->writeData(n->R2.get(), NULL, NULL);
			} else {
				SE->writeData(n->R1.get(), NULL, NULL);
			}
		} else {
			if (n->R2 || !SE ) {
				R1->writeData(n->R1.get(), n->R2.get(), NULL);
			} else {
				SE->writeData(n->R1.get(), NULL , NULL);
			}

		}
	}
}

/*Goes throught the char * with a 33 offset of each character*/
uint32_t BinarySearchTree::qualSum(char *q) {
	uint32_t score = 0;
	uint16_t i = 0;
	/*Should end on a null character, no tab or newlines*/
	while (q[i] != '\0') {
		/*Allows error checking*/
		if (q[i] > '~' || q[i] < '!') {
			fprintf(stderr, "Quality score is not between ascii [33,126], or [\",~]\n");
			fprintf(stderr, "Bad quality string = %s\n", q);
			fprintf(stderr, "Bad character='%c'\n", q[i]);
			exit(12);
		} else {
			score += q[i]-33;
		}
		i++;
	}

	return score;
}


bool BinarySearchTree::FlipBitsCheck(char *seq, bool r2) {
	/*Another error check for the humans*/
	if (strlen(seq) < start+charLength) {
		fprintf(stderr, "Error within binarySearch.cpp in funciton FlipBitsChars() strlen of %s larger than start + length\n", r2 ? "R2" : "R1");
		fprintf(stderr, "String Len = %lu, start = %d, length = %d\n", strlen(seq), start, charLength);
		disReads++;
		return false;
	}

	return true;

}


char RC_BP(char bp) {
	if (bp == 'N') {
		return 'N';
	} else if (bp == 'A') {
		return 'T';
	} else if (bp == 'T') {
		return 'A';
	} else if (bp == 'G') {
		return 'C';
	} else if (bp == 'C') {
		return 'G';
	} else {
		return '\0';
	}
}

void RC_Read(char *&seq) {
	char *new_seq = strdup(seq);
	int len = strlen(seq);
	int i = 0;
	while (seq[i] != '\0') {
		seq[(len - 1) - i] = RC_BP(new_seq[i]);
		if (seq[(len -1) - i] == '\0') {
			fprintf(stderr, "ERROR in RC_Read in binarySearch.cpp, Illegal char in %s\n", seq);
			exit(42);
		}
		i++;
	}
	free(new_seq);
}


/*Flip bits functionality for character string*/
bool BinarySearchTree::FlipBitsChars(std::shared_ptr<readInfo> R1, std::shared_ptr<readInfo> R2, idptr &id, bool RC) {
	char *seq_1, *seq_2;

	/*Incementor for id array*/
	uint16_t idLoc = 0, bitShifts = 0;

	seq_1 = R1->getSeq();
	if (!(FlipBitsCheck(seq_1, false))) { return false; }

	if (R2) {
		seq_2 = R2->getSeq();
		if (!(FlipBitsCheck(seq_2, true))) { return false; }
	} else {
		seq_2 = NULL;
	}

    std::unique_ptr<char[]> seq;

	if (seq_2 == NULL) {
		seq = std::unique_ptr<char[]>(new char[charLength + 1]);
		sprintf(seq.get(), "%.*s", charLength, seq_1 + start);
	}
	else if (RC) {
        seq = std::unique_ptr<char[]>(new char[charLength * 2 + 1]);
		sprintf(seq.get(), "%.*s%.*s", charLength, seq_2 + start, charLength, seq_1 + start);
	}
	else {
        seq = std::unique_ptr<char[]>(new char[charLength * 2 + 1]);
		sprintf(seq.get(), "%.*s%.*s", charLength, seq_1 + start, charLength, seq_2 + start);
	}

	/*
    if (RC) {
        RC_Read(&seq);
    }
	 */

	uint16_t loc = 0;

	id[0] = 0;

	while (seq[loc] != '\0') {
		uint16_t bit_1, bit_2;

		if (seq[loc] == 'N') {
			disReads++;
			return false;
		}

        // why + 10 here? should this be the startLoc?
		bit_1 = !!((seq[loc] + 10) & (1 << 2));
		bit_2= !!((seq[loc] + 10) & (1 << 4));
		id[idLoc] ^= (bit_1 << bitShifts);
		bitShifts++;
		id[idLoc] ^= (bit_2 << bitShifts);
		bitShifts++;

		if (bitShifts == 16) {
			bitShifts = 0;
			idLoc++;
			id[idLoc] = 0;
		}
		loc++;
	}

	idLoc++;

	for (int i = idLoc; i < newsize; i++) {
		id[i] = 0;
	}

	return true;
}


/*This creates the id that the binary search tree will use to dive down the tree*/
bool BinarySearchTree::getID(std::shared_ptr<readInfo> R1, std::shared_ptr<readInfo> R2, idptr &id) {
	/*3 bit * length / 16 this was handled constructur*/
	if (!R1 and !R2) {
		fprintf(stderr, "Both reads within binarySearch.cpp function getID() are NULL\n");
		fprintf(stderr, "R1 CANNOT be NULL\n");
		exit(15);
	}

    idptr tmp_id = idptr(new uint16_t[newsize]);
    idptr tmp_id_rc = idptr(new uint16_t[newsize]);
    
	/*Not funcitonal yet*/
	if (R1->optimized) {
		fprintf(stderr, "Optimizmations no funciotnal yet\n");
		exit(-1);
	} else {
		if (FlipBitsChars(R1, R2, tmp_id, false) && FlipBitsChars(R1, R2, tmp_id_rc, true)) {
			if (GreaterThan(tmp_id.get(), tmp_id_rc.get()) > 0) {
				id = std::move(tmp_id);
			} else {
				id = std::move(tmp_id_rc);
			}
			//printStuff(tmp_id, mallocLength);
			//printStuff(tmp_id_rc, mallocLength);
			/*Add Greater Than*/
			return true;
		}
	}
	return false;
}


/*Greather than return 1 (true), 0 for equal, and -1 for less than*/
/* 1 Greater Than
 * 0 Equal To
 * -1 Less Than*/

int BinarySearchTree::GreaterThan(uint16_t *test, uint16_t *value) {
	for (int i = 0; i < newsize; i++) {
		if (test[i] > value[i]) {
			return 1;
		} else if (test[i] < value[i]) {
			return -1;
		}
	}
	return 0;
}


/*Recursive function to add Node*/
void BinarySearchTree::PrivateAddNode(std::shared_ptr<Node> &n, std::shared_ptr<readInfo> R1_, std::shared_ptr<readInfo> R2_, idptr id, uint32_t qualScore ) {
	/*Add Node condition*/
	int tmpValue = 0;
	if (n == NULL) {
		nodesCreated++;
		n = std::make_shared<Node>(R1_, R2_, std::move(id), qualScore);

		n->count = 1;
		return;
	}

	tmpValue = GreaterThan(id.get(), n->id.get());
	if (tmpValue == 1) {
		PrivateAddNode(n->left, R1_, R2_, std::move(id), qualScore);
	} else if (tmpValue == -1) {
		PrivateAddNode(n->right, R1_, R2_, std::move(id), qualScore);
		/*Nodes are equal*/
	} else {
		/*Makes sure that single ends are kept track of*/

		if ((!R2_ && (n->single)) || (R2_ && !n->single)) {
			n->count++;
			dup_gone++;
			if (qualScore > n->qualScore) {
				replaced++;
				n->Replace(R1_, R2_, qualScore);
				/*exits*/
			} else {
			}
			/*If node is single end but the read isn't*/
		} else if (R2_) {
			PrivateAddNode(n->left, R1_, R2_, std::move(id), qualScore);
		} else if (!R2_) {
			PrivateAddNode(n->right, R1_, R2_, std::move(id), qualScore);
		}
	}

}


void BinarySearchTree::AddNode(std::shared_ptr<readInfo> R1_, std::shared_ptr<readInfo> R2_) {
    idptr id = nullptr;
	uint32_t qualScore = 0;
	reads_read++;


	if (!getID(R1_, R2_, id)) {
		disReads++;
		return;
	}
	/*R1 will never be null, but checking incase something goes horriblely wrong*/
	if (R1_.get() != NULL) {
		if (qualCheck) {
			qualScore += qualSum((R1_)->getQual());
		}
	} else {
		fprintf(stderr, "Read 1 should never be NULL within BinarySearchTree.cpp PrivateAddNode()\n");
		exit(11);
	}

	if (R2_.get() != NULL) {
		if (qualCheck) {
			qualScore += qualSum((R2_)->getQual());
		}
	} else {
		/*This just means it is single end reads*/
		/*This is acceptable, no error message*/
	}
	PrivateAddNode(root, R1_, R2_, std::move(id), qualScore);

}
