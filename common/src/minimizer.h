#ifndef __MINIMIZER_H_
#define __MINIMIZER_H_

#include "read.h"
#include "utils.h"

#include <deque>

class MinLoc {
public:
    ReferencePtr ref;
    size_t offset = 0;  // offset is in basepairs
    size_t length = 0;  // length in basepairs
    MinLoc(ReferencePtr ref_, size_t offset_, size_t length_) :
        ref(ref_), offset(offset_), length(length_) {}
};

class Minimizer {
public:
    Minimizer(size_t w_ = 10, size_t k_ = 10) : k(k_), w(w_) {
        hts_assert(k <= 32, "k must be less than or equal to 32");

        uint64_t max = std::numeric_limits<uint64_t>::max();
        mask = max >> (64-2*k);

        // reverse kmer append position
        rshift1 = (k-1)*2;
    }

    typedef std::unordered_multimap<uint64_t, MinLoc> KmerMap;

    
    // converts a string to a : A:0, T:3, C:1, G:2
    bool ascii_to_mer(uint64_t &bp, char c) {
        // todo: optimize this func?
        switch(c) {
        case 'A':
            bp = 0;
            break;
        case 'C':
            bp = 1;
            break;
        case 'G':
            bp = 2;
            break;
        case 'T':
            bp = 3;
            break;
        default:
            return false;
        }
        return true;
    }

    // convert string to kmer, used for testing mainly
    boost::optional<uint64_t> two_bit(const std::string& seq, size_t index = 0, size_t length = -1) {
        uint64_t kmer;
        uint64_t bp;
        size_t stop = length == static_cast<size_t>(-1) ? seq.size() : index+length;
        for (size_t i = index; i < stop; ++i) {
            if(ascii_to_mer(bp, seq[i])) {
                kmer = (kmer << 2 | bp) & mask;
            } else {
                return boost::none;
            }
        }
        return kmer;
    }

    bool append_two_bit(uint64_t& kmer, uint64_t& rkmer, u_char c) {
        uint64_t bp = 0;
        if(ascii_to_mer(bp, c)) {
            kmer = (kmer << 2 | bp) & mask;
            rkmer = (rkmer >> 2) | ((3u ^ bp) << rshift1);
        } else {
            kmer = 0;
            rkmer = 0;
            return false;
        }
        return true;
    }


    // thanks to https://www.biostars.org/p/113640/ <- and a slight mod!
    uint64_t reverse_complement(const uint64_t mer)
    {
        uint64_t res = ~mer;
        res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
        res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
        res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
        res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
        res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
        return (res >> (2 * (32 - k)));
    }

    void add_reference(ReferencePtr ref) {
        const std::string &seq = ref->get_seq();
        hts_assert(seq.size()/2 >= k+w-1, "sequence length must be >= k+w-1");

        uint64_t kmer = 0;
        uint64_t rkmer = 0;
        size_t ik = 0;

        std::multiset<uint64_t> min_kmer;
        std::deque<decltype(min_kmer)::iterator> kmer_queue;

        decltype(kmers)::iterator lastmin = kmers.end();

        for (size_t i = 0; i < seq.size(); ++i) {

            // todo: add template param for min function
            if (!append_two_bit(kmer, rkmer, seq[i])) {
                lastmin = kmers.end();
                kmer_queue.clear();
                min_kmer.clear();
                ik = 0;
                continue;
            }
            if (ik < k) {
                ++ik;
            }

            if (ik == k) {
                rkmer = reverse_complement(kmer);
                kmer_queue.push_back(min_kmer.insert(std::min(kmer, rkmer)));

                if(kmer_queue.size() > w) {
                    min_kmer.erase(kmer_queue.front());
                    kmer_queue.pop_front();
                }

                if (min_kmer.size() == w) {

                    if (lastmin != kmers.end() && lastmin->first == *min_kmer.begin()) {
                        ++(lastmin->second.length);
                    } else {
                        lastmin = kmers.emplace(*min_kmer.begin(), MinLoc(ref, i-(k+w-2), k+w-1));
                    }
                }
            }
        }

    }

    KmerMap& get_kmers() { return kmers; }

private:
    size_t k;
    size_t w;
    uint64_t mask;
    uint64_t rshift1;

    KmerMap kmers;

};

#endif // __MINIMIZER_H_
