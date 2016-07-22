#include "readInfo.h"

readInfo::readInfo(const std::string &head_, const std::string& seq_, const std::string &qual_) :
    header(head_),
    seq(seq_),
    qual(qual_)
{
}
