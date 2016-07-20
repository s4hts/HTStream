#ifndef SOURCE_BINARYSEARCHTREE_H_
#define SOURCE_BINARYSEARCHTREE_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <memory>

#include "readInfo.h"
#include "fileWriter.h"



class BinarySearchTree {
public:
    typedef std::unique_ptr<uint16_t[]> idptr;

private:
	class Node {
	public:
		/*Search pointers for the bst*/
        std::shared_ptr<Node> left;
        std::shared_ptr<Node> right;


		/*Contains read info, it will either have R1, or both R1 and R2*/
		std::shared_ptr<readInfo> R1;
		std::shared_ptr<readInfo> R2;

		/*ID is the key for the BST*/
		idptr id;

		/*If quality checking is off, qualScore is zero*/
		uint32_t qualScore;
		uint16_t count;

		bool single;

		/*This is the only constructor needed, the only reason to create is to have reads within each node*/
		Node (std::shared_ptr<readInfo> R1_, std::shared_ptr<readInfo> R2_, idptr id_, uint32_t qualScore_) :
            left(NULL), right(NULL), R1(R1_), R2(R2_), id(std::move(id_)), qualScore(qualScore_), count(0), single(false)
        {
			if (!R2) {
				single = true;
			}
		}

		/*Replace values within the node and frees old memory*/
		void Replace(std::shared_ptr<readInfo> R1_, std::shared_ptr<readInfo> R2_, uint32_t qualScore_) {
            R1 = R1_;
            R2 = R2_;
			qualScore = qualScore_;
		}
	};
	/*A, T, G, C binary mapped to values*/
	uint8_t A;
	uint8_t T;
	uint8_t G ;
	uint8_t C;


	/*Tree Stats*/
	time_t time_start;
	time_t time_end;
	long long unsigned int disReads;
	long long unsigned int nodesCreated;
	long long unsigned int singletons;
	long long unsigned int doubles;
	long long unsigned int threeplus;
	long long unsigned int replaced;
	long long unsigned int reads_read;
	long long unsigned int dup_gone;

	bool qualCheck;


    std::shared_ptr<Node> root;


	uint32_t qualSum(char *q);
	uint16_t charLength, newsize, start;


	bool getID(std::shared_ptr<readInfo> R1, std::shared_ptr<readInfo> R2, idptr &id);
	bool FlipBitsChars(std::shared_ptr<readInfo> R1, std::shared_ptr<readInfo> R2, idptr &id, bool RC);
	int GreaterThan(uint16_t *test, uint16_t *value);
	void PrivateAddNode(std::shared_ptr<Node> &n, std::shared_ptr<readInfo> R1_, std::shared_ptr<readInfo> R2_, idptr id, uint32_t qualScore);
	void PrintAndDeletePrivate(Node *n, FileWriter *R1, FileWriter *R2, FileWriter *SE);
	bool FlipBitsCheck(char *seq, bool r2);

public:

	BinarySearchTree(uint16_t length, uint16_t startLoc):
		/*00*/
		A(0),
		/*11*/
        T(3),
		/*10*/
		G(2),
		/*01*/
		C(1),
        time_start(0),
        time_end(0),
        disReads(0),
        nodesCreated(0),
        singletons(0),
        doubles(0),
        threeplus(0),
        replaced(0),
        reads_read(0),
        dup_gone(0),
        qualCheck(true),
        root(nullptr),
        charLength(length),
        // set the length to be 3 (for 3 bit format) * length specified / 16 (num of bits) + 1 (allocate the correct amount) 
        // 2*i = # bits,
        newsize((2*length)/sizeof(uint16_t) + 1),
        /*converts human value to correct position in zero start array*/
        start(startLoc -1)   {	}
    
	void Cleaner(uint16_t *&bin);

	void AddNode(std::shared_ptr<readInfo> R1_, std::shared_ptr<readInfo> R2_);

	void setQualCheck(bool b) {qualCheck = b;}
	void PrintAndDelete(FileWriter *R1, FileWriter *R2, FileWriter *SE);
	void endTime() {time_end = time(0);}

	/*Must be called after PrintAndDelete*/
	void outputStats(FILE *f);

};


#endif

