#ifndef SOURCE_READDATA_H
#define SOURCE_READDATA_H


#include<stdint.h>
#include<stdlib.h>
#include <string.h>

class readInfo {


private:
	char *header;
	char *seq;
	char *qual;

	uint16_t *useq;
	uint16_t *uqual;

public:

	bool optimized;

	readInfo(char *head_, char *seq_, char *qual_, bool optimized_);
	~readInfo () {
		free(header);
		free(seq);
		free(qual);
		//free(useq);
		//free(uqual);
	}
	char *getSeq() {return seq;}
	char *getQual() {return qual;}
	char *getHeader() {return header;}


};



#endif
