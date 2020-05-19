#ifndef ADAPTER_TRIM_H
#define ADAPTER_TRIM_H
//  this is so we can implment hash function for dynamic_bitset
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
#define STARTS 4098

#include "ioHandler.h"
#include "utils.h"
#include "threadutils.h"
#include "main_template.h"

#include <map>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <algorithm>
#include <bitset>
#include <utility>

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;

class AdapterCounters : public Counters {

public:
    uint64_t SE_Adapter_Trim = 0;
    uint64_t SE_Adapter_BpTrim = 0;

    uint64_t R1_Adapter_Trim = 0;
    uint64_t R1_Adapter_BpTrim = 0;
    uint64_t R2_Adapter_Trim = 0;
    uint64_t R2_Adapter_BpTrim = 0;

    AdapterCounters(const std::string &program_name, const po::variables_map &vm) : Counters::Counters(program_name, vm) {

        se.push_back(std::forward_as_tuple("adapterTrim", SE_Adapter_Trim));
        se.push_back(std::forward_as_tuple("adapterBpTrim", SE_Adapter_BpTrim));

        r1.push_back(std::forward_as_tuple("adapterTrim", R1_Adapter_Trim));
        r1.push_back(std::forward_as_tuple("adapterBpTrim", R1_Adapter_BpTrim));
        r2.push_back(std::forward_as_tuple("adapterTrim", R2_Adapter_Trim));
        r2.push_back(std::forward_as_tuple("adapterBpTrim", R2_Adapter_BpTrim));
    }

    using Counters::output;
    void output(SingleEndRead &ser)  {
        Counters::output(ser);
        Read &one = ser.non_const_read_one();
        if (one.getLengthTrue() < one.getLength()) {
            ++SE_Adapter_Trim;
            SE_Adapter_BpTrim += (one.getLength() - one.getLengthTrue());
        }
    }

    void output(PairedEndRead &per)  {
        Counters::output(per);
        Read &one = per.non_const_read_one();
        Read &two = per.non_const_read_two();
        if (one.getLengthTrue() < one.getLength()){
            ++R1_Adapter_Trim;
            R1_Adapter_BpTrim += one.getLength() - one.getLengthTrue();
        }
        if (two.getLengthTrue() < two.getLength()) {
            ++R2_Adapter_Trim;
            R2_Adapter_BpTrim += two.getLength() - two.getLengthTrue();
        }
    }
};

class AdapterTrimmer: public MainTemplate<AdapterCounters, AdapterTrimmer> {
public:

    void add_extra_options(po::options_description &desc) {

        setThreadPoolParams(desc);
        // number-of-threads|p

        setDefaultParamsOverlapping(desc);
        // kmer|k ; kmer-offset|r ; max-mismatch-errorDensity|x
        // check-lengths|c ; min-overlap|o

        desc.add_options()
            ("no-fixbases,X", po::bool_switch()->default_value(false), "after trimming adapter, DO NOT use consensus sequence of paired reads, only trims adapter sequence");
        desc.add_options()
            ("adapter-sequence,a", po::value<std::string>()->default_value("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"),  "Primer sequence to trim in SE adapter trimming, default is truseq ht primer sequence");
    }

    AdapterTrimmer() {
        program_name = "hts_AdapterTrimmer";
        app_description =
            "Adapter Trimmer, trims off adapters by overlapping paired-end reads and\n";
        app_description += "  trimming off overhangs which by definition are adapter sequence in standard\n";
        app_description += "  libraries. SE Reads are trimmed by overlapping the adapter-sequence and trimming off the overlap.";
    }

/*If adapater trimming is turned on that means adapter trimming and do not overlap
 * so trim adapter, but don't worry about the overlap.
 * however we still need to change the overlap
 *
 * Within the overlap if they are the same bp, then add q scores
 * If they are different bp, subtract q scores and take the larger quality base*/
    unsigned int checkIfOverlap(Read &r1, Read &r2, size_t loc1, size_t loc2, const double misDensity, const size_t &mismatch, const size_t &minOverlap, bool noFixBases = false ) {
        size_t minLoc = std::min(loc1, loc2);
        int loc1_t = loc1 - minLoc;
        int loc2_t = loc2 - minLoc;
        int r1_len = r1.getLength();
        int r2_len = r2.getLength();

        size_t maxLoop = std::min(r1_len - loc1_t, r2_len - loc2_t);
        size_t maxMis = std::min(mismatch, static_cast<size_t>(maxLoop * misDensity));

        const std::string &seq1 = r1.get_seq();
        const std::string &seq2 = r2.get_seq_rc();

        auto i1 = seq1.begin();
        std::advance(i1, loc1_t);
        auto i2 = seq2.begin();
        std::advance(i2, loc2_t);
        if (maxLoop < minOverlap || !threshold_mismatches(i1, i2, maxLoop, maxMis) ) {
            // no overlap identified
            return 0;
        }
        // overlap exists thats meet maxMis criteria
        if (!noFixBases){
            size_t read1_bp;
            size_t read2_bp;

            const std::string &qual1 = r1.get_qual();
            const std::string &qual2 = r2.get_qual_rc();

            char bp;
            char qual;

            for (size_t i = 0; i < maxLoop; ++i) {
                read1_bp = loc1_t + i;
                read2_bp = loc2_t + i;
                if (seq1[read1_bp] == seq2[read2_bp] )  {
                    qual = static_cast<char>(std::min(qual1[read1_bp] + qual2[read2_bp] - 33, 40 + 33));  // MATT: I still don't agree with this, so 38 [confident] and 2 [not so confident] return a 40 [very confident]??
                    r1.changeQual(read1_bp, qual);
                    r2.changeQual( (r2_len - 1) -  read2_bp, qual);
                } else {
                    bp = qual1[read1_bp] >= qual2[read2_bp] ? seq1[read1_bp] : seq2[read2_bp];
                    qual = static_cast<char>(std::max(qual1[read1_bp] - qual2[read2_bp] + 33, 1 + 33));

                    r1.changeSeq(read1_bp, bp);
                    r1.changeQual(read1_bp, qual);
                    //Something clever with RC
                    r2.changeSeq( (r2_len - 1) - read2_bp , rc(bp) );
                    r2.changeQual( (r2_len - 1) -  read2_bp , qual);
                }
            }
        }

        if (r1_len - loc1_t > r2_len - loc2_t) {
            r1.setRCut( ( (r2_len - loc2_t) + loc1_t));  // sins??
            if (loc2_t + r2_len + loc1_t > r1_len) { //special engulf check
                r2.setRCut(r2_len - loc2_t);
            }
        }

        return 1;
    }

/*Because of the way overlapping works, you only need to check the ends of the shorter read*/
    void getOverlappedReads(Read &r1, Read &r2, const seqLookup &seq1Map,  const double misDensity, const size_t &mismatch, const size_t &minOver, const size_t &checkLengths, const size_t kmer, bool noFixBases = false ) {

        // checkLengths should at most be 1/2 read length to check all kmers in read
        std::string seq2 = r2.get_seq_rc();
        //if all we are doing is trimming adapters, why look for long inserts at all? Skip the first loop? second loop still has product > read_length
        /*Do a quick check if the shorter read kmer shows up in longer read (read 2)
         * If it does, then try the brute force approach*/
        for (size_t bp = 0; bp < (checkLengths - kmer); ++bp) {
            // check first checkLength kmers at the beginning of read2 (mins)
            auto test = seq1Map.equal_range(seq2.substr(bp, kmer));
            for (auto it = test.first; it != test.second; ++it) {
                unsigned int overlapped = checkIfOverlap(r1, r2, it->second, bp, misDensity, mismatch, minOver, noFixBases);
                if (overlapped) {
                    return;
                }
            }
        }
        for (size_t bp = seq2.length() - (checkLengths + kmer); bp <= (seq2.length() - kmer) ; ++bp) {
            // check last checkLengths kmers at the end of the read2 (sins)
            auto test = seq1Map.equal_range(seq2.substr(bp, kmer));
            for (auto it = test.first; it != test.second; ++it) {
                unsigned int overlapped = checkIfOverlap(r1, r2, it->second, bp, misDensity, mismatch, minOver, noFixBases);
                if (overlapped) {
                    return;
                }
            }
        }
    }

    ReadBasePtr check_read_pe(PairedEndReadPtr pe, const double misDensity, const size_t mismatch, const size_t minOver, const size_t checkLengths, const size_t kmer, const size_t kmerOffset, bool noFixBases = false ) {

        Read &r1 = pe->non_const_read_one();
        Read &r2 = pe->non_const_read_two();

        bool swapped = false;
        /*Read1 is always longer than Read 2)*/
        if (r1.getLength() < r2.getLength()) {
            std::swap(r1, r2);
            swapped = true;
        }
        /* checkL needs to be as long as or longer than the shortest read */
        size_t checkL = std::min(r2.getLength(), checkLengths);
        /* kmer needs to be as long as or longer than the shortest read */
        size_t kkmer = std::min(r2.getLength(), kmer);
        /* Create a map with non-overlapping kmers*/
        seqLookup mOne = readOneMap(r1.get_seq(), kkmer, kmerOffset);
        /* returns null if no much
         * r1 and r2 and passed by ref in case only adapter trimming is on */
        getOverlappedReads(r1, r2, mOne, misDensity, mismatch, minOver, checkL, kkmer, noFixBases) ;
        if (swapped) {
            std::swap(r1, r2);
        }
        return std::dynamic_pointer_cast<ReadBase>(pe);
    }

    unsigned int checkIfAdapter(Read &r1, Read &adapter, size_t loc1, size_t loc2, const double misDensity, const size_t &mismatch, const size_t &minOverlap ) {
        size_t minLoc = std::min(loc1, loc2);
        int loc1_t = loc1 - minLoc;
        int loc2_t = loc2 - minLoc;
        int r1_len = r1.getLength();
        int adapter_len = adapter.getLength();

        size_t maxLoop = std::min(r1_len - loc1_t, adapter_len - loc2_t);
        size_t maxMis = std::min(mismatch, static_cast<size_t>(maxLoop * misDensity));

        const std::string &seq1 = r1.get_seq();
        const std::string &seq_adapter = adapter.get_seq();

        auto i1 = seq1.begin();
        std::advance(i1, loc1_t);
        auto i2 = seq_adapter.begin();
        std::advance(i2, loc2_t);
        if (maxLoop < minOverlap || !threshold_mismatches(i1, i2, maxLoop, maxMis) ) {
            // no overlap identified
            return 0;
        }
        // overlap exists thats meet maxMis criteria
        r1.setRCut(loc1_t);

        return 1;
    }

    ReadBasePtr check_read_se(SingleEndReadPtr se , const double misDensity, const size_t &mismatch, const size_t &minOver, const size_t &checkLengths, const size_t kmer, const size_t kmerOffset, std::string adapter_seq = "" ) {

        Read &r1 = se->non_const_read_one();
        /* if adapter is longer than sequence, set length to same as seq */
        if (adapter_seq.length() > r1.getLength()){
            adapter_seq = adapter_seq.substr(0, r1.getLength());
        }

        Read adapter = Read(adapter_seq, "", "");

        /* checkL and kkmer cannot be larger than the adapter length */
        size_t kkmer = std::min(adapter.getLength(), kmer);
        /* Create a map with non-overlapping kmers*/
        seqLookup mOne = readOneMap(r1.get_seq(), kkmer, kmerOffset);

        /*Do a quick check if the adapter kmer shows up in the read
         * If it does, then try the brute force approach*/
        for (size_t bp = 0; bp < (checkLengths - kmer); ++bp) {
            // check first checkLength kmers at the beginning of read
            auto test = mOne.equal_range(adapter_seq.substr(bp, kmer));
            for (auto it = test.first; it != test.second; ++it) {
                unsigned int overlapped = checkIfAdapter(r1, adapter, it->second, bp, misDensity, mismatch, minOver);
                if (overlapped) {
                    return se;
                }
            }
        }
        return std::dynamic_pointer_cast<ReadBase>(se);
    }

    void writer_thread(std::shared_ptr<OutputWriter> pe,  std::shared_ptr<OutputWriter> se, AdapterCounters &counter, threadsafe_queue<std::future<ReadBasePtr>> &futures) {

        WriterHelper writer(pe, se);
        while(!futures.is_done()) {
            std::future<ReadBasePtr> fread;
            futures.wait_and_pop(fread);

            ReadBasePtr rbase = fread.get();

            // null read indicates all done
            if (!rbase.get()) {
                futures.set_done();
                return;
            }

            counter.output(*rbase);
            writer(*rbase);
        }
    }


/*This is the helper class for overlap
 * The idea is in the wet lab, they set up sequences to sequences toward each other
 * Sometimes, these squences can overlap
 * There are two cases in which they can overlap
 * They can either overlap just barely on the ends we call - these lins (long insert)
 * They can also overlap way to much to the point they have adapters
 * in the read - or a sin (short insert)
 * With a lin it is useful to have a higher confidence in the bases in the overlap and longer read
 * With a sin it is useful to have the higher confidence as well as removing the adapters*/
    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, AdapterCounters &counter, const po::variables_map &vm) {

        const double misDensity = vm["max-mismatch-errorDensity"].as<double>();
        const size_t mismatch = vm["max-mismatch"].as<size_t>();
        const size_t minOver = vm["min-overlap"].as<size_t>();
        const size_t checkLengths = vm["check-lengths"].as<size_t>();
        const size_t kmer = vm["kmer"].as<size_t>();
        const size_t kmerOffset = vm["kmer-offset"].as<size_t>();
        bool noFixBases = vm["no-fixbases"].as<bool>();
        std::string adapter = vm["adapter-sequence"].as<std::string>();
        size_t num_threads = vm["number-of-threads"].as<size_t>();

        threadsafe_queue<std::future<ReadBasePtr>> futures(50000);
        thread_pool threads(50000, num_threads);

        std::thread output_thread([=, &counter, &futures]() mutable { writer_thread(pe, se, counter, futures); });
        thread_guard tg(output_thread);

        try {
            while(reader.has_next()) {
                auto i = reader.next();
                PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());
                if (per) {
                    std::shared_ptr<PairedEndRead> sper = std::make_shared<PairedEndRead>(std::move(*per));
                    counter.input(*sper);
                    futures.push(threads.submit([=]() mutable {
                                return check_read_pe(sper, misDensity, mismatch, minOver, checkLengths, kmer, kmerOffset, noFixBases); }));

                } else {
                    SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());

                    if (ser) {
                        std::shared_ptr<SingleEndRead> sser = std::make_shared<SingleEndRead>(std::move(*ser));
                        counter.input(*sser);
                        futures.push(threads.submit([=]() mutable {
                                    return check_read_se(sser, misDensity, mismatch, minOver, checkLengths, kmer, kmerOffset, adapter); }));
                    } else {
                        throw std::runtime_error("Unknown read type");
                    }
                }
            }

            // null ptr indicates end of processing
            futures.push(threads.submit([]() { return ReadBasePtr(); }));
        } catch (...) {

            // make sure threads stop on exception
            futures.set_done();
            throw;
        }
    }
};
#endif
