#ifndef __HTS_EXCEPTION_H_
#define __HTS_EXCEPTION_H_

#include <string>

class HtsException {
public:
    HtsException(const std::string &what_): w(what_)
    { }

    const std::string& what() const {
        return w;
    }
    
private:
    std::string w;
};

class HtsIOException : public HtsException {
public:
    HtsIOException(const std::string &what_): HtsException(what_)
    {};

};

class HtsRuntimeException : public HtsException {
public:
    HtsRuntimeException(const std::string &what_): HtsException(what_)
    {};
};

#endif // __HTS_EXCEPTION_H_
