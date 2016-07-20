#include "readInfo.h"

readInfo::readInfo(char *head_, char *seq_, char *qual_, bool optimized_) {
	header = strdup(head_);
	seq = strdup(seq_);
	qual = strdup(qual_);
	useq = NULL;
	uqual = NULL;
	optimized = optimized_;
}
