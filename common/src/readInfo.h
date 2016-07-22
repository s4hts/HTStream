#ifndef SOURCE_READDATA_H
#define SOURCE_READDATA_H


#include <string>

class readInfo {


private:
    std::string header;
    std::string seq;
    std::string qual;

public:

	readInfo(const char *head_, const char *seq_, const char *qual_);

	const char *getSeq() {return seq.c_str();}
	const char *getQual() {return qual.c_str();}
	const char *getHeader() {return header.c_str();}


};



#endif
